"""EC-number annotation layer via ESM-ECForest (issue #20).

ESM-ECForest (https://github.com/Xiao88866/ESM-ECForest) is a two-stage
ESM-2 + Random Forest framework: enzyme/non-enzyme classification, then
hierarchical EC prediction. It runs in its OWN conda env (Python 3.9,
PyTorch/Transformers/ESM/scikit-learn; pre-trained models from Zenodo), so
the wrapper is a config command template — and it SHARES the issue #17
ESM-2 embedding cache (steps/esm.py): embed once per proteome, never twice.

Same design position as DeepLoc (issue #16): an ANNOTATION layer over
confirmed families — annotate members and map EC onto the family tree with
the utils/newick.py Fitch machinery (exactly as steps/retargeting.py maps
localizations). NEVER a clustering criterion. A member whose predicted EC
disagrees with the family consensus is an ADVISORY review flag, never a
hard filter — EC changes within families are real biology.

Methods caveat: the model is trained on general enzyme data; plant-specific
accuracy is unvalidated. Verify the CAM spot-checks (PEPC = 4.1.1.31,
NADP-ME = 1.1.1.40, PPDK = 2.7.9.1) before trusting genome-wide.
"""

import csv
import logging
import shlex
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Set

from utils.newick import fitch, parse_newick, state_switches

logger = logging.getLogger("family_finder")

# Truthy spellings accepted for the is_enzyme column.
_TRUE_VALUES = {"1", "true", "yes", "enzyme", "y", "t"}


@dataclass(frozen=True)
class ECSwitchEvent:
    """A gene-tree branch where the Fitch-resolved EC state changes."""

    family_id: str
    parent_ec: str
    child_ec: str
    clade_leaves: FrozenSet[str]

    @property
    def clade_id(self) -> str:
        return ",".join(sorted(self.clade_leaves))


# ---------------------------------------------------------------------------
# Subprocess wrapper + parser (module-level, mockable)
# ---------------------------------------------------------------------------

def run_ecforest(input_path: Path, outdir: Path, config) -> Path:
    """Run ESM-ECForest on a FASTA or precomputed-embedding input.

    Driven by config.ecforest_cmd, e.g.
    "conda run -n enzyllm_env python predict.py --input {input} --outdir {outdir}".
    Prefer feeding the #17 embedding cache as {input} — do not embed twice.
    Returns the newest results CSV in outdir (deeploc.run_deeploc convention).
    """
    if not config.ecforest_cmd.strip():
        raise ValueError("ESM-ECForest: ecforest_cmd not configured (see config.py)")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = shlex.split(
        config.ecforest_cmd.format(input=str(input_path), outdir=str(outdir))
    )
    logger.info(f"Running ESM-ECForest: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"ESM-ECForest failed:\n{result.stderr}")
        raise RuntimeError(
            f"ESM-ECForest failed with return code {result.returncode}"
        )
    csvs = sorted(outdir.glob("*.csv"), key=lambda p: p.stat().st_mtime)
    if not csvs:
        raise FileNotFoundError(f"No ESM-ECForest results CSV found in {outdir}")
    return csvs[-1]


def parse_ecforest(csv_path: Path) -> Dict[str, dict]:
    """ESM-ECForest results CSV -> gene_id -> {is_enzyme, ec, confidence}.

    Columns are taken generically (deeploc.parse_deeploc convention): the
    first column is the gene id; a column named like "is_enzyme" gives the
    stage-1 call; a column named like "EC" gives the predicted EC string
    (empty for non-enzymes); a column named like "confidence" gives the
    prediction confidence (0.0 when absent/unparseable).
    """
    preds: Dict[str, dict] = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"Empty ESM-ECForest CSV: {csv_path}")
        id_col = reader.fieldnames[0]
        enzyme_col = _find_col(reader.fieldnames, ("is_enzyme", "enzyme"))
        ec_col = _find_col(reader.fieldnames, ("ec_number", "ec"))
        conf_col = _find_col(reader.fieldnames, ("confidence", "conf", "probability"))
        for row in reader:
            ec = (row.get(ec_col) or "").strip() if ec_col else ""
            is_enzyme = False
            if enzyme_col:
                is_enzyme = (row.get(enzyme_col) or "").strip().lower() in _TRUE_VALUES
            elif ec:
                is_enzyme = True  # no stage-1 column: a non-empty EC implies enzyme
            try:
                confidence = float(row.get(conf_col, "")) if conf_col else 0.0
            except (TypeError, ValueError):
                confidence = 0.0
            preds[row[id_col]] = {
                "is_enzyme": is_enzyme,
                "ec": ec,
                "confidence": confidence,
            }
    return preds


def _find_col(fieldnames, wanted) -> Optional[str]:
    """First column whose lowercased name matches one of `wanted`, in order."""
    lowered = {name.lower(): name for name in fieldnames}
    for w in wanted:
        if w in lowered:
            return lowered[w]
    return None


# ---------------------------------------------------------------------------
# Once-per-proteome cache (TSV with header, deeploc.save_cache convention)
# ---------------------------------------------------------------------------

def save_cache(preds: Dict[str, dict], path: Path) -> None:
    """Persist predictions: gene<TAB>is_enzyme<TAB>ec<TAB>confidence per line."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("gene_id\tis_enzyme\tec\tconfidence\n")
        for gene in sorted(preds):
            p = preds[gene]
            f.write(
                f"{gene}\t{int(p['is_enzyme'])}\t{p['ec']}\t{p['confidence']:g}\n"
            )


def load_cache(path: Path) -> Dict[str, dict]:
    preds: Dict[str, dict] = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            gene, is_enzyme, ec, confidence = line.rstrip("\n").split("\t")
            preds[gene] = {
                "is_enzyme": is_enzyme == "1",
                "ec": ec,
                "confidence": float(confidence),
            }
    return preds


# ---------------------------------------------------------------------------
# Family-level consensus + tree mapping (pure functions, unit-tested)
# ---------------------------------------------------------------------------

def _gene_ec(pred: Optional[dict], min_conf: float = 0.0) -> Optional[str]:
    """A gene's usable EC state: enzyme call with a non-empty EC clearing
    min_conf, else None (uninformative — never forced)."""
    if not pred or not pred.get("is_enzyme"):
        return None
    ec = pred.get("ec", "")
    if not ec or pred.get("confidence", 0.0) < min_conf:
        return None
    return ec


def family_ec_consensus(
    members: Set[str],
    predictions: Dict[str, dict],
    min_agree: float = 0.5,
) -> Optional[str]:
    """Majority EC among a family's EC-predicted members, else None.

    The consensus is the most common EC iff its share of VOTING members
    (those with a usable EC) is >= min_agree; a tie for the top EC is no
    consensus. Members without a usable EC abstain rather than dilute —
    non-enzyme families simply return None.
    """
    votes = Counter(
        ec for ec in (_gene_ec(predictions.get(m)) for m in members) if ec
    )
    if not votes:
        return None
    ranked = votes.most_common()
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        return None  # tied majority — never pick by dict order
    top_ec, top_n = ranked[0]
    if top_n / sum(votes.values()) < min_agree:
        return None
    return top_ec


def consensus_disagreements(
    members: Set[str],
    predictions: Dict[str, dict],
    consensus: Optional[str],
) -> List[str]:
    """Members whose predicted EC disagrees with the family consensus.

    ADVISORY review flags only (like the DeepLoc disagreement flag) — a
    profile-assigned or merged gene appearing here deserves a look, but EC
    changes within families are real biology, so this is never a filter.
    """
    if consensus is None:
        return []
    return sorted(
        m for m in members
        if _gene_ec(predictions.get(m)) not in (None, consensus)
    )


def ec_switch_events(
    family_id: str,
    newick_text: str,
    predictions: Dict[str, dict],
    exclude: Set[str] = frozenset(),
) -> List[ECSwitchEvent]:
    """EC-switch branches on one family gene tree (EC string = Fitch state).

    Exactly the steps/retargeting.py pattern on utils/newick.py: leaves get
    their predicted EC as state (None for non-enzymes, empty ECs, and the
    exclude set — truncated models fake switches just like retargeting),
    Fitch resolves ancestral ECs, and every state-changing branch becomes an
    event. An EC-switching branch + a #16 branch-site selection signal =
    enzymatic neofunctionalization candidate.
    """
    root = parse_newick(newick_text)

    states: Dict[str, Optional[str]] = {}
    for leaf in root.leaf_names():
        if leaf in exclude:
            states[leaf] = None
        else:
            states[leaf] = _gene_ec(predictions.get(leaf))

    informative = {s for s in states.values() if s}
    if len(informative) < 2:
        return []  # uniform or uninformative — no switch possible

    fitch(root, states)
    return [
        ECSwitchEvent(family_id, parent, child, leaves)
        for parent, child, leaves in state_switches(root)
    ]


def write_events_tsv(events: List[ECSwitchEvent], path) -> None:
    with open(path, "w") as f:
        f.write("family_id\tparent_ec\tchild_ec\tn_clade_genes\tclade_genes\n")
        for e in events:
            f.write(
                f"{e.family_id}\t{e.parent_ec}\t{e.child_ec}\t"
                f"{len(e.clade_leaves)}\t{e.clade_id}\n"
            )
