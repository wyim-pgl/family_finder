"""EPA-ng phylogenetic placement adjudicator (issue #23).

The v2 tiered architecture (issue #21) nominates candidate families for a
gene by profile bits (tier 2) or structure TM-score (tier 3); when those
single-score methods are AMBIGUOUS, this module lets the tree decide —
the mechanism proven on the Ppc1 pair (#13): IQ-TREE topology unites it,
no single-score method can. Placement is an adjudication layer between
caller-supplied candidate families; it never redefines homology by itself.

Workflow (all external calls are module-level, mockable wrappers):
  1. raxml-ng --evaluate re-estimates branch lengths + model parameters on
     the family reference tree (family trees come from FastTree, which
     writes no model file; EPA-ng needs the .bestModel).
  2. hmmalign (against an hmmbuild of the reference) or mafft --add aligns
     the query into the reference alignment columns.
  3. epa-ng --split separates the combined alignment back into reference
     and query parts (EPA-ng writes only query.fasta — src/util/split.hpp).
  4. epa-ng places the query; the jplace v3 result is parsed with the
     file's own "fields" ordering (column positions are never hardcoded).

Design rules (as in steps/profile_assign.py and steps/esm.py):
  * Decisions use like_weight_ratio with a threshold AND a best-vs-second
    margin; exact ties are ALWAYS ambiguous, even at a zero margin.
    Likelihoods are stored for reporting, never compared across queries.
    The thresholds are UNCALIBRATED (no published LWR cutoff exists) —
    config defaults are starting points for calibration on the 5sp panel.
  * Abstention: placement always answers "which edge is least bad", never
    "none of the above" — config.epa_max_pendant (when set) rejects
    high-LWR placements whose pendant branch is so long the query belongs
    to no reference family.
  * Place in PROTEIN space: references are the confirmed protein
    alignments, not codon alignments — divergent queries (the whole point
    of adjudication) are exactly where codon alignment is least reliable.
  * REFERENCE TREES SHOULD NOT BE RAW FastTree TOPOLOGIES: placement
    inherits the reference topology wholesale, and FastTree was measured
    splitting the true Ppc1 pair that IQ-TREE unites (#13) — a placement
    on such a reference returns high LWR for the known failure. Prefer
    IQ-TREE (or otherwise re-estimated) references; raxml-ng --evaluate
    re-optimizes branch lengths and model on the FIXED input topology but
    cannot repair it.
  * edge_family_map is an APPROXIMATION: an edge inherits the majority
    family of the leaf clade below it, so a placement on a deep edge of a
    merged clan resolves to the family occupying that side; edges whose
    clade is tied between families are omitted (-> ambiguous).
"""

import json
import logging
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from config import Config
from utils.newick import parse_newick

logger = logging.getLogger("family_finder")

# Substitution model handed to raxml-ng --evaluate. Family references are
# protein alignments (the confirmed MAFFT alignments), so LG+G.
RAXML_EVAL_MODEL = "LG+G"

# jplace v3 per-placement fields this module needs. The file's "fields"
# array defines the column order — parse_jplace reads it, never positions.
REQUIRED_JPLACE_FIELDS = (
    "edge_num",
    "likelihood",
    "like_weight_ratio",
    "distal_length",
    "pendant_length",
)

# jplace edge markers: "{N}" after a branch length or a bare name/paren.
_EDGE_AFTER_DIST = re.compile(r"(:[^,(){};]+)\{(\d+)\}")
_EDGE_BARE = re.compile(r"\{(\d+)\}")
_NAME_MARKER = re.compile(r"@@(\d+)@@$")


@dataclass(frozen=True)
class Placement:
    """One jplace placement row for one query."""

    edge_num: int
    likelihood: float          # reporting only — never compared across queries
    like_weight_ratio: float   # the decision statistic
    distal_length: float
    pendant_length: float


@dataclass(frozen=True)
class JointReference:
    """A joint reference over a set of candidate families (merged clan).

    The caller builds it once per candidate set: `alignment` and `tree`
    span the UNION of the candidate families' confirmed members, and
    family_of_leaf maps every reference leaf back to its family.
    `alignment` is a PROTEIN alignment, and `tree` should be an IQ-TREE
    (or otherwise re-estimated) topology — see the module docstring on
    why raw FastTree references are only a knowing fallback.
    """

    alignment: Path
    tree: Path
    family_of_leaf: Dict[str, str]


# ---------------------------------------------------------------------------
# Subprocess wrappers (module-level so tests can mock them)
# ---------------------------------------------------------------------------

def _run_cmd(cmd: List[str], step: str, stdout_path: Optional[Path] = None) -> None:
    """Run one external command; optionally capture stdout to a file.

    Module-level so tests mock this one function for every wrapper.
    Raises RuntimeError on a non-zero exit.
    """
    logger.info(f"Running {step}: {' '.join(cmd)}")
    if stdout_path is not None:
        with open(stdout_path, "w") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    else:
        result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"{step} failed:\n{result.stderr}")
        raise RuntimeError(f"{step} failed with return code {result.returncode}")


def evaluate_model(
    ref_alignment: Path, ref_tree: Path, outdir: Path, config: Config
) -> Tuple[Path, Path]:
    """raxml-ng --evaluate; returns (.bestModel path, .bestTree path).

    Family trees come from FastTree, which writes no model parameters and
    whose branch lengths EPA-ng cannot trust — raxml-ng re-optimizes both
    on the FIXED topology and emits <prefix>.raxml.bestModel (the model
    file run_epa passes via --model, the recommended workflow) and
    <prefix>.raxml.bestTree (re-optimized branch lengths). PLACEMENT MUST
    RUN ON THE EVALUATED bestTree, not the input tree: branch lengths
    differ, and the evaluated tree is the one jplace edge numbering will
    refer to (edge_family_map guards against drift anyway by reading the
    tree string embedded in the jplace file itself).

    WARNING (#13 Ppc1): --evaluate cannot repair the topology. FastTree
    confirmed trees split the true Ppc1 pair that IQ-TREE unites; passing
    a raw FastTree reference is acceptable only as a fallback the caller
    knowingly chooses — prefer IQ-TREE references for adjudication.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / "eval"
    cmd = [
        config.raxml_ng_bin, "--evaluate",
        "--msa", str(ref_alignment),
        "--tree", str(ref_tree),
        "--model", RAXML_EVAL_MODEL,
        "--prefix", str(prefix),
    ]
    _run_cmd(cmd, "raxml-ng --evaluate")
    return Path(f"{prefix}.raxml.bestModel"), Path(f"{prefix}.raxml.bestTree")


def align_query(
    query_fasta: Path, ref_alignment: Path, outdir: Path, config: Config
) -> Path:
    """Align query sequences into the reference alignment columns.

    Backend chosen by config.epa_query_align:
      * "hmmalign"  — hmmbuild a profile from the reference, then hmmalign
        with --mapali so the output combines reference + query (--trim
        drops insert-state columns that would break the reference layout).
      * "mafft-add" — mafft --add query --keeplength ref (query-side
        insertions are cut, reference columns preserved).

    Returns the combined alignment path (reference + query rows).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    combined = outdir / "combined.fa"
    backend = config.epa_query_align
    if backend == "hmmalign":
        hmm_path = outdir / "ref.hmm"
        _run_cmd(
            [config.hmmbuild_bin, "--amino", str(hmm_path), str(ref_alignment)],
            "hmmbuild (EPA reference)",
        )
        _run_cmd(
            [
                config.hmmalign_bin, "--amino", "--trim",
                "--outformat", "afa",
                "-o", str(combined),
                "--mapali", str(ref_alignment),
                str(hmm_path), str(query_fasta),
            ],
            "hmmalign",
        )
    elif backend == "mafft-add":
        _run_cmd(
            [
                config.mafft_bin,
                "--add", str(query_fasta),
                "--keeplength",
                str(ref_alignment),
            ],
            "mafft --add",
            stdout_path=combined,
        )
    else:
        raise ValueError(
            f"Unknown epa_query_align backend: {backend!r} "
            f"(expected 'hmmalign' or 'mafft-add')"
        )
    return combined


def split_alignment(
    combined: Path, ref_alignment: Path, outdir: Path, config: Config
) -> Tuple[Path, Path]:
    """epa-ng --split the combined alignment; returns (ref_msa, query_msa).

    EPA-ng writes ONLY query.fasta into --outdir (src/util/split.hpp) —
    the reference half of the pair is the untouched ref_alignment, whose
    columns the query part is guaranteed to match after the split.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        config.epa_ng_bin, "--redo",
        "--split", str(ref_alignment), str(combined),
        "--outdir", str(outdir),
    ]
    _run_cmd(cmd, "epa-ng --split")
    return Path(ref_alignment), outdir / "query.fasta"


def run_epa(
    ref_tree: Path,
    ref_msa: Path,
    query_msa: Path,
    model_file: Path,
    outdir: Path,
    config: Config,
) -> Path:
    """Place queries on the reference tree; returns the epa_result.jplace path.

    ref_tree should be the raxml-ng EVALUATED bestTree from evaluate_model
    (jplace edge numbering refers to the tree actually given to epa-ng;
    downstream mapping reads the jplace's own tree string, but branch
    lengths only exist on the evaluated tree). See the evaluate_model
    warning: the topology must be trustworthy — not a raw FastTree tree —
    because placement inherits it wholesale.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        config.epa_ng_bin, "--redo",
        "--tree", str(ref_tree),
        "--ref-msa", str(ref_msa),
        "--query", str(query_msa),
        "--model", str(model_file),
        "--outdir", str(outdir),
    ]
    _run_cmd(cmd, "epa-ng")
    return outdir / "epa_result.jplace"


# ---------------------------------------------------------------------------
# jplace parsing + acceptance
# ---------------------------------------------------------------------------

def parse_jplace(path: Path) -> Dict[str, List[Placement]]:
    """Parse a jplace v3 file into query -> placements (LWR descending).

    Column positions come from the file's "fields" array — jplace writers
    are free to reorder columns, so positions are NEVER hardcoded. Both
    name conventions are handled: "n" (plain list) and "nm" ([name,
    multiplicity] pairs). Raises ValueError when a required field is
    missing from the fields array.
    """
    with open(path) as f:
        data = json.load(f)

    fields = data.get("fields", [])
    missing = [name for name in REQUIRED_JPLACE_FIELDS if name not in fields]
    if missing:
        raise ValueError(f"jplace file {path} missing required fields: {missing}")
    idx = {name: fields.index(name) for name in REQUIRED_JPLACE_FIELDS}

    result: Dict[str, List[Placement]] = {}
    for entry in data.get("placements", []):
        names = entry.get("n")
        if names is None:
            names = [pair[0] for pair in entry.get("nm", [])]
        placements = []
        for row in entry.get("p", []):
            placements.append(Placement(
                edge_num=int(row[idx["edge_num"]]),
                likelihood=float(row[idx["likelihood"]]),
                like_weight_ratio=float(row[idx["like_weight_ratio"]]),
                distal_length=float(row[idx["distal_length"]]),
                pendant_length=float(row[idx["pendant_length"]]),
            ))
        # LWR descending; edge_num as deterministic secondary key only.
        placements.sort(key=lambda p: (-p.like_weight_ratio, p.edge_num))
        for name in names:
            result[name] = placements
    return result


def accept_placement(
    placements: List[Placement], config: Config
) -> Optional[Placement]:
    """Accept the best placement iff LWR floor AND best-vs-second margin hold.

    Acceptance requires ALL of:
      * best.like_weight_ratio >= config.epa_min_lwr
      * pendant gate: config.epa_max_pendant is 0.0 (disabled) or
        best.pendant_length <= it. LWR only says which edge is LEAST BAD —
        a high-LWR placement on a huge pendant branch means the query
        belongs to no reference family, and placement never abstains on
        its own; this gate is the abstention rule.
      * single placement, or best-vs-second LWR margin >=
        config.epa_lwr_margin — exact ties are ALWAYS ambiguous, even
        with a zero margin threshold.
    The LWR thresholds are UNCALIBRATED (no published cutoff exists).

    Returns the accepted Placement, or None (ambiguous / no signal).
    """
    if not placements:
        return None
    ordered = sorted(placements, key=lambda p: (-p.like_weight_ratio, p.edge_num))
    best = ordered[0]
    if best.like_weight_ratio < config.epa_min_lwr:
        return None
    if 0.0 < config.epa_max_pendant < best.pendant_length:
        return None
    if len(ordered) == 1:
        return best
    second = ordered[1]
    has_margin = (
        best.like_weight_ratio - second.like_weight_ratio >= config.epa_lwr_margin
        and best.like_weight_ratio > second.like_weight_ratio
    )
    return best if has_margin else None


# ---------------------------------------------------------------------------
# Edge -> family mapping + adjudication
# ---------------------------------------------------------------------------

def edge_family_map(
    ref_tree_newick: str, family_of_leaf: Dict[str, str]
) -> Dict[int, str]:
    """Map jplace edge numbers to families via the reference leaf partition.

    The jplace "tree" string carries {N} edge markers after each node's
    branch length. APPROXIMATION (documented per issue #23): each edge is
    assigned the MAJORITY family among the leaves of the clade BELOW it —
    exact for edges inside a family's clade, and for the deep edges of a
    merged clan it resolves to the family occupying that side. Edges whose
    clade is tied between families (e.g. the root edge of a balanced
    two-family clan) are OMITTED from the map, so a placement there stays
    ambiguous. Leaves absent from family_of_leaf are uninformative.
    """
    # Move each {N} marker into the preceding node's name so the plain
    # utils.newick parser can read the string, then strip the markers.
    s = _EDGE_AFTER_DIST.sub(lambda m: f"@@{m.group(2)}@@{m.group(1)}", ref_tree_newick)
    s = _EDGE_BARE.sub(lambda m: f"@@{m.group(1)}@@", s)
    root = parse_newick(s)

    edge_of_node: Dict[int, int] = {}
    for node in root.postorder():
        m = _NAME_MARKER.search(node.name)
        if m:
            edge_of_node[id(node)] = int(m.group(1))
            node.name = node.name[: m.start()]

    mapping: Dict[int, str] = {}
    for node in root.postorder():
        edge_num = edge_of_node.get(id(node))
        if edge_num is None:
            continue
        counts: Dict[str, int] = {}
        for leaf in node.leaf_names():
            family = family_of_leaf.get(leaf)
            if family:
                counts[family] = counts.get(family, 0) + 1
        if not counts:
            continue
        ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
        if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
            continue  # tied clade -> ambiguous edge, never id-broken
        mapping[edge_num] = ranked[0][0]
    return mapping


def _jplace_tree(path: Path) -> str:
    """The jplace file's reference tree string (with {N} edge markers)."""
    with open(path) as f:
        return json.load(f)["tree"]


def adjudicate(
    gene_id: str,
    candidate_families: Iterable[str],
    family_refs: JointReference,
    query_fasta: Path,
    outdir: Path,
    config: Config,
) -> Optional[str]:
    """Place one ambiguous gene into the merged-clan reference; tree decides.

    For the MERGED-CLAN pattern: candidate_families come from an ambiguous
    tier-2/3 nomination (steps.profile_assign.assign / steps.esm.
    tier3_assign); the caller supplies family_refs, a JointReference whose
    alignment + tree span the union of those families' confirmed members.
    The gene is placed once on the joint tree and the accepted edge is
    mapped back to a family via edge_family_map.

    Returns the winning family id, or None when the placement is rejected
    (LWR floor / margin / tie), lands on an edge without a family majority,
    or resolves outside candidate_families (defensive: the joint reference
    should only contain candidates, but the tree is the authority).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    candidates = set(candidate_families)

    # Place on the raxml-ng EVALUATED tree (re-optimized branch lengths);
    # jplace edge numbers refer to it, and edge_family_map reads the tree
    # string back out of the jplace file so the numbering cannot drift.
    model_file, eval_tree = evaluate_model(
        family_refs.alignment, family_refs.tree, outdir, config
    )
    combined = align_query(query_fasta, family_refs.alignment, outdir, config)
    ref_msa, query_msa = split_alignment(combined, family_refs.alignment, outdir, config)
    jplace = run_epa(eval_tree, ref_msa, query_msa, model_file, outdir, config)

    placements = parse_jplace(jplace).get(gene_id)
    if not placements:
        logger.debug(f"EPA adjudication: no placements for {gene_id}")
        return None
    best = accept_placement(placements, config)
    if best is None:
        logger.debug(f"EPA adjudication: {gene_id} ambiguous (LWR floor/margin)")
        return None

    mapping = edge_family_map(_jplace_tree(jplace), family_refs.family_of_leaf)
    family = mapping.get(best.edge_num)
    if family is None:
        logger.debug(
            f"EPA adjudication: {gene_id} placed on edge {best.edge_num} "
            f"with no family majority"
        )
        return None
    if family not in candidates:
        logger.debug(
            f"EPA adjudication: {gene_id} resolved to {family}, "
            f"outside candidates {sorted(candidates)}"
        )
        return None
    return family


def write_adjudications_tsv(
    adjudications: Dict[str, Optional[str]], path: Path
) -> None:
    """Adjudication table: gene_id, winning family (empty if none), status."""
    path = Path(path)
    with open(path, "w") as f:
        f.write("gene_id\tfamily_id\tstatus\n")
        for gene_id in sorted(adjudications):
            family = adjudications[gene_id]
            status = "placed" if family is not None else "ambiguous"
            f.write(f"{gene_id}\t{family or ''}\t{status}\n")
    logger.debug(f"EPA adjudications written: {path}")
