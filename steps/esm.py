"""Tier-3 structural assignment scaffolding: ESM-2 / ESMFold / Foldseek (issue #17).

Three-tier assignment: 1) DIAMOND+MCL, 2) HMM profiles (steps/profile_assign.py),
3) this module — structure comparison for the final unplaced pool only.
Structure is conserved far longer than sequence, so ESMFold + Foldseek is the
last resort for genes with no sequence signal at all (the #8 D-class).

Local scaffolding only: every external tool call is a module-level, mockable
subprocess wrapper driven by a config command TEMPLATE (cluster tooling varies —
esm-extract vs custom scripts, GPU arrays). Heavy compute happens on the cluster;
everything here that makes a DECISION is pure Python and unit-tested.

Design rules (mirroring steps/profile_assign.py):
  * Decisions use alntmscore with a best-vs-second MARGIN; E-values are
    stored for reporting but NEVER compared between hits.
  * Equal-score ties are ambiguous, never broken by input order.
  * ESM-2 embedding cosine similarity is a tie-breaker for #13's ambiguous
    profile assignments (embedding_tiebreak) — a secondary signal, not a
    clustering criterion.
  * Low mean pLDDT genes are flagged as mispredicted/truncated model
    candidates — the 7th evidence line for pseudogene scoring.
"""

import logging
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from config import Config

logger = logging.getLogger("family_finder")

# The ONE foldseek output layout this module parses. Declared in the
# easy-search command (--format-output) and in parse_foldseek — keep in sync.
# (foldseek's own default has 12 BLAST-tab columns and no alntmscore.)
FOLDSEEK_FORMAT_OUTPUT = "query,target,evalue,bits,alntmscore"

# Maximum candidates reported per ambiguous gene (as in profile_assign).
MAX_AMBIGUOUS_CANDIDATES = 5

# Structure filename extensions stripped when mapping a foldseek target
# back to a family id (target DB entries are <family_id>.pdb representatives).
_STRUCT_SUFFIXES = (".pdb.gz", ".cif.gz", ".pdb", ".cif")


@dataclass(frozen=True)
class StructHit:
    """One gene-vs-family-representative structural hit."""

    gene_id: str
    family_id: str
    evalue: float      # significance gate / reporting only — never compared
    bits: float
    alntmscore: float


# ---------------------------------------------------------------------------
# Subprocess wrappers (module-level so tests can mock them)
# ---------------------------------------------------------------------------

def _run_template(template: str, mapping: Dict[str, str], step: str) -> None:
    """Format a config command template, shlex-split, and run it.

    Module-level so tests mock this one function for every template-driven
    step. Raises ValueError on an empty template (step not configured) and
    RuntimeError on a non-zero exit.
    """
    if not template.strip():
        raise ValueError(f"{step}: command template not configured (see config.py)")
    cmd = shlex.split(template.format(**mapping))
    logger.info(f"Running {step}: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"{step} failed:\n{result.stderr}")
        raise RuntimeError(f"{step} failed with return code {result.returncode}")


def run_esm_embed(fasta: Path, outdir: Path, config: Config) -> Path:
    """Extract ESM-2 embeddings for one proteome FASTA. Returns outdir.

    Driven by config.esm_embed_cmd, e.g.
    "esm-extract esm2_t33_650M_UR50D {fasta} {outdir} --include mean".
    Run ONCE per proteome; convert the tool's output to the TSV cache with
    save_embeddings so downstream consumers (embedding_tiebreak here,
    ESM-ECForest in steps/ecforest.py) never re-embed.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    _run_template(
        config.esm_embed_cmd,
        {"fasta": str(fasta), "outdir": str(outdir)},
        "ESM-2 embedding",
    )
    return outdir


def run_esmfold(fasta: Path, outdir: Path, config: Config) -> Path:
    """Fold sequences with ESMFold. Returns the output PDB directory.

    Fold ONLY the post-tier-2 unplaced genes plus one representative per
    family (annotate_families.py --prep-fold-list), never the full proteome.
    Driven by config.esmfold_cmd, e.g. "esm-fold -i {fasta} -o {outdir}".
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    _run_template(
        config.esmfold_cmd,
        {"fasta": str(fasta), "outdir": str(outdir)},
        "ESMFold",
    )
    return outdir


def run_foldseek(
    query_pdb_dir: Path, target_db: Path, out_tsv: Path, config: Config
) -> Path:
    """foldseek easy-search of query structures vs the family-representative DB.

    The output columns are pinned with --format-output to
    FOLDSEEK_FORMAT_OUTPUT — parse_foldseek parses exactly that layout.
    Target DB entries must be named <family_id>.pdb (one representative per
    family) so targets map back to families.
    """
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_tsv.parent / "foldseek_tmp"
    cmd = [
        config.foldseek_bin, "easy-search",
        str(query_pdb_dir), str(target_db), str(out_tsv), str(tmp_dir),
        "--format-output", FOLDSEEK_FORMAT_OUTPUT,
    ]
    logger.info(f"Running foldseek: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"foldseek failed:\n{result.stderr}")
        raise RuntimeError(f"foldseek failed with return code {result.returncode}")
    return out_tsv


# ---------------------------------------------------------------------------
# Embedding cache + cosine tie-break (the #13 ambiguous-case hook)
# ---------------------------------------------------------------------------

def save_embeddings(embeddings: Dict[str, List[float]], path: Path) -> None:
    """Persist embeddings: gene_id<TAB>comma-joined floats, with a header.

    Once-per-proteome cache — the single interchange format between the
    cluster-side extraction and every local consumer.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("gene_id\tembedding\n")
        for gene in sorted(embeddings):
            vec = ",".join(repr(x) for x in embeddings[gene])
            f.write(f"{gene}\t{vec}\n")


def load_embeddings(path: Path) -> Dict[str, List[float]]:
    embeddings: Dict[str, List[float]] = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            gene, vec = line.rstrip("\n").split("\t")
            embeddings[gene] = [float(x) for x in vec.split(",")]
    return embeddings


def cosine(a: List[float], b: List[float]) -> float:
    """Cosine similarity of two plain-list vectors (no numpy dependency).

    Zero-norm vectors get similarity 0.0 (uninformative, not an error).
    """
    if len(a) != len(b):
        raise ValueError(f"Embedding dimension mismatch: {len(a)} vs {len(b)}")
    dot = sum(x * y for x, y in zip(a, b))
    norm_a = sum(x * x for x in a) ** 0.5
    norm_b = sum(y * y for y in b) ** 0.5
    if norm_a == 0.0 or norm_b == 0.0:
        return 0.0
    return dot / (norm_a * norm_b)


def embedding_tiebreak(
    gene: str,
    candidates: Dict[str, List[str]],
    embeddings: Dict[str, List[float]],
    min_delta: float,
) -> Optional[str]:
    """Break a #13 ambiguous profile assignment with embedding similarity.

    Args:
        gene: the ambiguous gene.
        candidates: family_id -> representative gene ids (the ambiguous
            candidate families from steps.profile_assign.assign()).
        embeddings: ESM-2 embedding cache (load_embeddings).
        min_delta: required margin of the best family's MEAN cosine
            similarity over the runner-up's (config.embed_tiebreak_min_delta).

    Returns the winning family, or None when the gene/representatives lack
    embeddings, or the margin is not met (equal means stay ambiguous — this
    is a tie-BREAKER, it never invents a preference).
    """
    gene_emb = embeddings.get(gene)
    if gene_emb is None:
        return None

    scores: Dict[str, float] = {}
    for family_id, reps in candidates.items():
        sims = [
            cosine(gene_emb, embeddings[rep])
            for rep in reps
            if rep in embeddings
        ]
        if sims:
            scores[family_id] = sum(sims) / len(sims)

    if not scores:
        return None
    # Deterministic ordering: score desc, then family id (id never decides
    # the winner — an id-only difference fails the margin below).
    ranked = sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))
    if len(ranked) == 1:
        return ranked[0][0]
    best, second = ranked[0], ranked[1]
    if best[1] - second[1] >= min_delta:
        return best[0]
    return None


# ---------------------------------------------------------------------------
# pLDDT (pseudogene evidence byproduct)
# ---------------------------------------------------------------------------

def parse_plddt(pdb_path: Path) -> float:
    """Mean pLDDT of one predicted structure = mean B-factor of CA atoms.

    ESMFold (like AlphaFold) writes per-residue pLDDT into the B-factor
    column. Pure fixed-column PDB text parsing: atom name is columns 13-16,
    B-factor is columns 61-66 (1-based). Raises ValueError when the file
    has no CA atoms.
    """
    total = 0.0
    n = 0
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                total += float(line[60:66])
            except ValueError:
                continue
            n += 1
    if n == 0:
        raise ValueError(f"No CA atoms found in {pdb_path}")
    return total / n


def flag_low_plddt(
    plddt_by_gene: Dict[str, float], config: Config
) -> List[str]:
    """Genes whose mean pLDDT falls below config.plddt_flag_below.

    Sorted flag list for pseudogene evidence (steps/pseudogene.py) —
    ADVISORY input to the evidence scorer, never a hard filter.
    """
    return sorted(
        gene for gene, plddt in plddt_by_gene.items()
        if plddt < config.plddt_flag_below
    )


# ---------------------------------------------------------------------------
# Foldseek parsing + tier-3 assignment
# ---------------------------------------------------------------------------

def _target_to_family(target: str) -> str:
    """Map a foldseek target name back to its family id.

    Target DB entries are one representative structure per family named
    <family_id>.pdb; foldseek may report the name with or without the
    extension and with a trailing chain suffix (e.g. "_A" is NOT stripped —
    single-chain predicted monomers carry no chain suffix in easy-search).
    """
    for suffix in _STRUCT_SUFFIXES:
        if target.endswith(suffix):
            return target[: -len(suffix)]
    return target


def parse_foldseek(tsv: Path) -> Dict[str, List[StructHit]]:
    """Parse a foldseek easy-search TSV with FOLDSEEK_FORMAT_OUTPUT columns.

    Columns (tab-separated, no header): query, target, evalue, bits,
    alntmscore. Multiple hits to the same family keep only the best
    alntmscore. Returns gene_id -> hits sorted by alntmscore DESCENDING
    (bits as a deterministic secondary sort key only). E-values are stored
    for reporting and are never compared.
    """
    best: Dict[Tuple[str, str], StructHit] = {}
    with open(tsv) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                logger.debug(f"Skipping malformed foldseek line: {line[:80]!r}")
                continue
            try:
                hit = StructHit(
                    gene_id=fields[0],
                    family_id=_target_to_family(fields[1]),
                    evalue=float(fields[2]),
                    bits=float(fields[3]),
                    alntmscore=float(fields[4]),
                )
            except ValueError:
                logger.debug(f"Skipping malformed foldseek line: {line[:80]!r}")
                continue
            key = (hit.gene_id, hit.family_id)
            prev = best.get(key)
            if prev is None or hit.alntmscore > prev.alntmscore:
                best[key] = hit

    hits_by_gene: Dict[str, List[StructHit]] = {}
    for hit in best.values():
        hits_by_gene.setdefault(hit.gene_id, []).append(hit)
    for hits in hits_by_gene.values():
        hits.sort(key=lambda h: (-h.alntmscore, -h.bits, h.family_id))
    return hits_by_gene


def tier3_assign(
    hits_by_gene: Dict[str, List[StructHit]], config: Config
) -> Tuple[Dict[str, str], List[dict]]:
    """Decide tier-3 structural assignments from parsed foldseek hits.

    A gene is accepted into its top-alntmscore family iff BOTH:
      * best.alntmscore >= config.tier3_min_tmscore
      * no second family, or best-vs-second alntmscore margin >=
        config.tier3_margin_tm — equal-TM ties are ALWAYS ambiguous, even
        with a zero margin threshold. E-values never decide anything.

    A gene failing only the margin is reported in `ambiguous` with its top
    candidates (mirror of steps.profile_assign.assign()). Assigned genes
    still require tree validation before committing — same safeguard as
    profile assignment; this function only nominates.

    Returns (assignments: gene -> family, ambiguous: [{"gene", "candidates"}]).
    """
    assignments: Dict[str, str] = {}
    ambiguous: List[dict] = []

    for gene_id, hits in hits_by_gene.items():
        if not hits:
            continue
        best = hits[0]
        if best.alntmscore < config.tier3_min_tmscore:
            continue  # no structural signal — stays unplaced
        second = hits[1] if len(hits) > 1 else None
        has_margin = second is None or (
            best.alntmscore - second.alntmscore >= config.tier3_margin_tm
            and best.alntmscore > second.alntmscore
        )
        if has_margin:
            assignments[gene_id] = best.family_id
            continue
        ambiguous.append({
            "gene": gene_id,
            "candidates": [
                {
                    "family_id": h.family_id,
                    "alntmscore": h.alntmscore,
                    "bits": h.bits,
                }
                for h in hits[:MAX_AMBIGUOUS_CANDIDATES]
            ],
        })

    return assignments, ambiguous


def write_tier3_tsv(
    assignments: Dict[str, str],
    hits_by_gene: Dict[str, List[StructHit]],
    path: Path,
) -> None:
    """Assignment table with structure-evidence columns (issue #17 output)."""
    with open(path, "w") as f:
        f.write("gene_id\tfamily_id\talntmscore\tbits\tevalue\n")
        for gene_id in sorted(assignments):
            family_id = assignments[gene_id]
            hit = next(
                (h for h in hits_by_gene.get(gene_id, []) if h.family_id == family_id),
                None,
            )
            if hit is None:
                continue
            f.write(
                f"{gene_id}\t{family_id}\t{hit.alntmscore:.3f}\t"
                f"{hit.bits:.1f}\t{hit.evalue:.2e}\n"
            )
    logger.debug(f"Tier-3 assignments written: {path}")
