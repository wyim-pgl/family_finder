"""DeepLoc 2 wrapper: run once per proteome, parse, cache with probabilities.

Issue #16. DeepLoc is an annotation layer over family gene trees — NEVER a
clustering criterion (families are homology units; localization changes
within them). The probability vector is kept, not just the argmax label,
so dual targeting is not flattened away.
"""

import csv
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, Optional

from config import Config

logger = logging.getLogger("family_finder")

# DeepLoc 2 CSV columns that are not class probabilities
NON_PROB_COLS = {"Protein_ID", "Localizations", "Signals"}


def run_deeploc(fasta: Path, outdir: Path, config: Config) -> Path:
    """Run DeepLoc 2 on one proteome FASTA. Returns the results CSV path.

    Run once per proteome and cache (save_cache) — never per family.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [config.deeploc_bin, "-f", str(fasta), "-o", str(outdir),
           "-m", config.deeploc_model]
    if config.deeploc_device:
        cmd += ["-d", config.deeploc_device]

    logger.info(f"Running DeepLoc: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"DeepLoc failed:\n{result.stderr}")
        raise RuntimeError(f"DeepLoc failed with return code {result.returncode}")

    csvs = sorted(outdir.glob("*.csv"), key=lambda p: p.stat().st_mtime)
    if not csvs:
        raise FileNotFoundError(f"No DeepLoc results CSV found in {outdir}")
    return csvs[-1]


def parse_deeploc(csv_path: Path) -> Dict[str, dict]:
    """DeepLoc 2 results CSV -> gene_id -> {'label': str, 'probs': {class: p}}.

    Columns are taken generically: the first column is the ID, 'Localizations'
    is the label, every other float-parseable column is a class probability.
    """
    preds: Dict[str, dict] = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"Empty DeepLoc CSV: {csv_path}")
        id_col = reader.fieldnames[0]
        for row in reader:
            probs = {}
            for key, val in row.items():
                if key == id_col or key in NON_PROB_COLS:
                    continue
                try:
                    probs[key] = float(val)
                except (TypeError, ValueError):
                    continue
            preds[row[id_col]] = {"label": row.get("Localizations", ""), "probs": probs}
    return preds


def top_location(pred: Optional[dict], min_prob: float = 0.5) -> Optional[str]:
    """Highest-probability class if it clears min_prob, else None (ambiguous
    predictions — including split dual-targeting probabilities — are treated
    as uninformative rather than forced into a single state)."""
    if not pred:
        return None
    if not pred["probs"]:
        return pred["label"] or None
    loc, p = max(pred["probs"].items(), key=lambda kv: kv[1])
    return loc if p >= min_prob else None


def save_cache(preds: Dict[str, dict], path: Path) -> None:
    """Persist predictions: gene<TAB>label<TAB>probs-json per line."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("gene_id\tlabel\tprobs_json\n")
        for gene in sorted(preds):
            p = preds[gene]
            f.write(f"{gene}\t{p['label']}\t{json.dumps(p['probs'], sort_keys=True)}\n")


def load_cache(path: Path) -> Dict[str, dict]:
    preds: Dict[str, dict] = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            gene, label, probs_json = line.rstrip("\n").split("\t")
            preds[gene] = {"label": label, "probs": json.loads(probs_json)}
    return preds
