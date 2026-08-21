#!/usr/bin/env python3
"""Merge every annotation axis into one per-gene matrix (reusable stack).

Project goal: distant/edge members are judged by MULTI-AXIS agreement,
never by one tool. This CLI merges the gate-validated axes —

  --emapper        eggNOG-mapper .emapper.annotations   (EC/symbol, orthology)
  --clean-csv      CLEAN maxsep CSV                     (EC + confidence)
  --foldseek-tsv   structure_function_transfer.tsv      (fs_transfer.py output:
                   foldseek vs AFDB-SwissProt + UniProt join)
  --deeploc-csv    DeepLoc 2.x results CSV              (localization, signals)
  --signalp        SignalP 6 prediction_results.txt     (secretory SP)

— into annotation_matrix.tsv (one row per gene, one column block per axis)
and, when --expected-ec is given, a membership verdict per gene:

  member    >=2 axes support the expected EC
  intruder  0 support and >=1 CONFIDENT contradiction (wrong-EC foldseek hit
            with qTM >= 0.6, or wrong EC at confidence >= 0.5)
  review    everything else — abstentions (CLEAN conf ~0, foldseek no-hit,
            absent rows) are not evidence of intrusion (PEPC lesson: the
            remote clan members abstain on every axis, the one true intruder
            Ccac_g10054 is contradicted by two independent axes)

All axes are optional; whatever is given is merged. Runs anywhere (pure
Python; the heavy compute already happened on the cluster).

Example (PEPC clan):
  python annotation_matrix.py --outdir annot \\
      --emapper clan.emapper.annotations --clean-csv clan_maxsep.csv \\
      --foldseek-tsv structure_function_transfer.tsv \\
      --deeploc-csv results.csv --signalp prediction_results.txt \\
      --expected-ec 4.1.1.31
"""

import argparse
import csv
import logging
from pathlib import Path
from typing import Dict, Optional

from steps.ec_sources import parse_clean, parse_emapper

logger = logging.getLogger("family_finder")

CLEAN_CONF_MIN = 0.5   # below this CLEAN abstains (its EC string is noise)
FS_QTM_MIN = 0.6       # a structural hit this good with the wrong EC contradicts

_EMPTY_ROW = {
    "emapper_ec": "", "emapper_symbol": "", "emapper_desc": "",
    "clean_ec": "", "clean_conf": 0.0,
    "fs_uniprot": "", "fs_name": "", "fs_ec": "", "fs_qtm": 0.0, "fs_fident": 0.0,
    "deeploc_loc": "", "deeploc_signals": "",
    "signalp": "", "signalp_sp_prob": 0.0,
}


def _parse_emapper_full(path: Path) -> Dict[str, dict]:
    """emapper annotations keeping symbol/description alongside the EC."""
    rows: Dict[str, dict] = {}
    header = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##") or not line:
                continue
            if line.startswith("#query"):
                header = line.split("\t")
                continue
            if header is None:
                continue
            row = dict(zip(header, line.split("\t")))
            rows[row["#query"]] = row
    return rows


def load_axes(
    emapper: Optional[Path] = None,
    clean_csv: Optional[Path] = None,
    foldseek_tsv: Optional[Path] = None,
    deeploc_csv: Optional[Path] = None,
    signalp_txt: Optional[Path] = None,
) -> dict:
    """Parse whichever axis files are given into per-axis gene dicts."""
    axes: dict = {}
    if emapper:
        axes["emapper"] = parse_emapper(Path(emapper))
        axes["emapper_full"] = _parse_emapper_full(Path(emapper))
    if clean_csv:
        axes["clean"] = parse_clean(Path(clean_csv))
    if foldseek_tsv:
        fs: Dict[str, dict] = {}
        with open(foldseek_tsv, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                fs[row["gene_id"]] = row
        axes["foldseek"] = fs
    if deeploc_csv:
        dl: Dict[str, dict] = {}
        with open(deeploc_csv, newline="") as f:
            for row in csv.DictReader(f):
                dl[row["Protein_ID"]] = row
        axes["deeploc"] = dl
    if signalp_txt:
        sp: Dict[str, dict] = {}
        with open(signalp_txt) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                sp[p[0].split()[0]] = {
                    "prediction": p[1],
                    "sp_prob": float(p[3]) if len(p) > 3 else 0.0,
                }
        axes["signalp"] = sp
    return axes


def _canon(gene: str) -> str:
    """Canonical join key: foldseek PDB filenames (and other tools) replace
    '.' with '_' — Gene.1 and Gene_1 are the same gene."""
    return gene.replace(".", "_")


def build_matrix(axes: dict) -> Dict[str, dict]:
    """One row per gene across the union of all axes; absent values empty.

    Genes are joined on a canonical id (dots folded to underscores); the
    display id prefers the dotted (original annotation) spelling.
    """
    by_canon: Dict[str, dict] = {}   # canon -> {axis: data}
    display: Dict[str, str] = {}
    emapper_full = {_canon(g): v for g, v in axes.get("emapper_full", {}).items()}
    for key in ("emapper", "clean", "foldseek", "deeploc", "signalp"):
        for g, data in axes.get(key, {}).items():
            c = _canon(g)
            by_canon.setdefault(c, {})[key] = data
            # prefer a dotted original id for display
            if c not in display or ("." in g and "." not in display[c]):
                display[c] = g

    matrix: Dict[str, dict] = {}
    for c in sorted(by_canon):
        per_axis = by_canon[c]
        row = dict(_EMPTY_ROW)
        if "emapper" in per_axis:
            row["emapper_ec"] = per_axis["emapper"]["ec"]
            full = emapper_full.get(c, {})
            pref = (full.get("Preferred_name") or "").strip()
            row["emapper_symbol"] = "" if pref == "-" else pref
            row["emapper_desc"] = (full.get("Description") or "").strip()
        if "clean" in per_axis:
            row["clean_ec"] = per_axis["clean"]["ec"]
            row["clean_conf"] = per_axis["clean"]["confidence"]
        if "foldseek" in per_axis:
            fs = per_axis["foldseek"]
            row["fs_uniprot"] = fs.get("uniprot", "")
            row["fs_name"] = fs.get("protein_name", "")
            row["fs_ec"] = fs.get("ec", "")
            row["fs_qtm"] = float(fs.get("qtmscore") or 0.0)
            row["fs_fident"] = float(fs.get("fident") or 0.0)
        if "deeploc" in per_axis:
            dl = per_axis["deeploc"]
            row["deeploc_loc"] = dl.get("Localizations", "")
            row["deeploc_signals"] = dl.get("Signals", "")
        if "signalp" in per_axis:
            row["signalp"] = per_axis["signalp"]["prediction"]
            row["signalp_sp_prob"] = per_axis["signalp"]["sp_prob"]
        matrix[display[c]] = row
    return matrix


def membership_verdict(row: dict, expected_ec: str) -> dict:
    """Judge one gene by axis agreement against the family's expected EC."""
    support, contradict = [], []

    if row["emapper_ec"]:
        (support if row["emapper_ec"] == expected_ec else contradict).append("emapper")
    if row["clean_ec"] and row["clean_conf"] >= CLEAN_CONF_MIN:
        (support if row["clean_ec"] == expected_ec else contradict).append("clean")
    if row["fs_ec"]:
        if row["fs_ec"] == expected_ec:
            support.append("foldseek")
        elif row["fs_qtm"] >= FS_QTM_MIN:
            contradict.append("foldseek")

    if len(support) >= 2:
        verdict = "member"
    elif not support and contradict:
        verdict = "intruder"
    else:
        verdict = "review"
    return {
        "verdict": verdict,
        "n_support": len(support),
        "n_contradict": len(contradict),
        "support_axes": ",".join(support),
        "contradict_axes": ",".join(contradict),
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--emapper", default=None)
    ap.add_argument("--clean-csv", default=None)
    ap.add_argument("--foldseek-tsv", default=None)
    ap.add_argument("--deeploc-csv", default=None)
    ap.add_argument("--signalp", default=None)
    ap.add_argument("--expected-ec", default=None,
                    help="Family EC — adds membership verdict columns")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    axes = load_axes(
        emapper=args.emapper, clean_csv=args.clean_csv,
        foldseek_tsv=args.foldseek_tsv, deeploc_csv=args.deeploc_csv,
        signalp_txt=args.signalp,
    )
    if not axes:
        ap.error("No axis inputs given")
    matrix = build_matrix(axes)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "annotation_matrix.tsv"
    cols = list(_EMPTY_ROW)
    verdict_cols = ["verdict", "n_support", "n_contradict",
                    "support_axes", "contradict_axes"]
    with open(out, "w") as f:
        header = ["gene_id"] + cols + (verdict_cols if args.expected_ec else [])
        f.write("\t".join(header) + "\n")
        for g, row in matrix.items():
            fields = [g] + [str(row[c]) for c in cols]
            if args.expected_ec:
                v = membership_verdict(row, args.expected_ec)
                fields += [str(v[c]) for c in verdict_cols]
            f.write("\t".join(fields) + "\n")
    logger.info(f"{len(matrix)} genes x {len(axes)} axes -> {out}")
    if args.expected_ec:
        from collections import Counter
        counts = Counter(
            membership_verdict(r, args.expected_ec)["verdict"] for r in matrix.values()
        )
        logger.info(f"verdicts: {dict(counts)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
