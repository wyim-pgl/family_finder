#!/usr/bin/env python3
"""Is a subfamily-specific region structurally disordered — and only there?

Reads per-residue pLDDT from the B-factor column of predicted structures
(ESMFold writes it there), maps a range of ALIGNMENT columns onto each
sequence's own residues, and compares mean pLDDT inside that range against
the whole protein.

The control is the entire point. Protein termini have low pLDDT in every
structure, so "the region is low-confidence" is worthless on its own. This
script measures the same alignment columns in the NON-focal members of the
family and reports the difference between the two groups with a
Mann-Whitney U test (no scipy).

Measured on the PEPC clan, SF3 signal region (alignment columns 165-209):
SF3 mean delta -0.304 across 9/9 members, non-focal members -0.071 with a
median of +0.009, p = 1.95e-05 — the drop is subfamily-specific, not the
universal floppy N-terminus.

Usage:
  python region_disorder.py --pdb-dir pdb --alignment clan.aln \\
      --region 165-209 --members sf3_members.txt -o disorder.json

Feed the JSON to `subfamily_report.py --disorder-json` to add it to the
subfunctionalization narrative as an independent evidence axis.
"""

import argparse
import json
import math
import statistics
import sys
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
from utils.gene_ids import match_ids
import sys
from pathlib import Path
from typing import Dict, List


def canon(gene: str) -> str:
    """Structure files are named from gene ids with '.' replaced by '_'."""
    return gene.replace(".", "_")


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def residue_plddt(pdb: Path) -> Dict[int, float]:
    """Residue number -> pLDDT (B-factor of its first atom)."""
    out: Dict[int, float] = {}
    for line in open(pdb):
        if line.startswith("ATOM"):
            resi = int(line[22:26])
            if resi not in out:
                out[resi] = float(line[60:66])
    return out


def columns_to_residues(aligned_seq: str, lo: int, hi: int) -> List[int]:
    """Alignment columns [lo, hi] -> this sequence's own residue numbers."""
    resi = 0
    picked = []
    for col, ch in enumerate(aligned_seq, start=1):
        if ch == "-":
            continue
        resi += 1
        if lo <= col <= hi:
            picked.append(resi)
    return picked


def mann_whitney_p(a: List[float], b: List[float]) -> float:
    """Two-sided Mann-Whitney U p-value via the normal approximation."""
    if not a or not b:
        return float("nan")
    merged = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks: Dict[int, float] = {}
    i = 0
    while i < len(merged):
        j = i
        while j + 1 < len(merged) and merged[j + 1][0] == merged[i][0]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[k] = r
        i = j + 1
    n1, n2 = len(a), len(b)
    r1 = sum(ranks[k] for k, (_v, grp) in enumerate(merged) if grp == 0)
    u1 = r1 - n1 * (n1 + 1) / 2
    u = min(u1, n1 * n2 - u1)
    sd = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    if sd == 0:
        return 1.0
    z = (u - n1 * n2 / 2) / sd
    return math.erfc(abs(z) / math.sqrt(2))


def measure(pdb_dir: Path, alignment: Path, lo: int, hi: int,
            members: set, min_residues: int = 5) -> dict:
    raw_aln = read_fasta(alignment)
    aln = {canon(k): v for k, v in raw_aln.items()}
    focal_members = {canon(m) for m in members}

    # One shared matcher, so a '.'-vs-'_' or transcript-suffix difference is
    # reported instead of quietly reducing coverage (#42, utils/gene_ids.py).
    pdb_names = [p.stem for p in sorted(Path(pdb_dir).glob("*.pdb"))]
    matched = match_ids(pdb_names, list(raw_aln)).mapping

    focal, other, rows = [], [], []
    skipped = {"not_in_alignment": [], "region_too_short": [], "no_plddt": []}
    for pdb in sorted(Path(pdb_dir).glob("*.pdb")):
        gene = canon(pdb.stem)
        ref = matched.get(pdb.stem)
        seq = raw_aln[ref] if ref is not None else aln.get(gene)
        if seq is None:
            # An id-form mismatch between structure filenames and alignment
            # names drops genes here without a word; on the PEPC clan that once
            # left 29 of 111 matching. Record it so shrinking coverage is
            # distinguishable from an absent signal.
            skipped["not_in_alignment"].append(gene)
            continue
        region = columns_to_residues(seq, lo, hi)
        if len(region) < min_residues:
            skipped["region_too_short"].append(gene)
            continue
        plddt = residue_plddt(pdb)
        inside = [plddt[r] for r in region if r in plddt]
        if len(inside) < min_residues:
            skipped["no_plddt"].append(gene)
            continue
        whole = statistics.mean(plddt.values())
        delta = statistics.mean(inside) - whole
        (focal if gene in focal_members else other).append(delta)
        rows.append({"gene": gene,
                     "group": "focal" if gene in focal_members else "other",
                     "region_plddt": round(statistics.mean(inside), 4),
                     "whole_plddt": round(whole, 4),
                     "delta": round(delta, 4)})

    result = {
        "region": f"{lo}-{hi}",
        "n_focal": len(focal),
        "n_other": len(other),
        "delta": round(statistics.mean(focal), 4) if focal else None,
        "delta_other": round(statistics.mean(other), 4) if other else None,
        "median_other": round(statistics.median(other), 4) if other else None,
        "all_below": bool(focal) and all(d < 0 for d in focal),
        "p": mann_whitney_p(focal, other) if (focal and other) else None,
        "per_gene": rows,
        "n_skipped": sum(len(v) for v in skipped.values()),
        "skipped": skipped,
    }
    return result


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--pdb-dir", required=True,
                    help="Directory of predicted structures (pLDDT in B-factor)")
    ap.add_argument("--alignment", required=True, help="Family alignment (FASTA)")
    ap.add_argument("--region", required=True, help="Alignment columns, 'LO-HI'")
    ap.add_argument("--members", required=True,
                    help="File of focal-subfamily gene ids, one per line")
    ap.add_argument("-o", "--out", required=True, help="JSON output")
    args = ap.parse_args()

    lo, hi = (int(x) for x in args.region.split("-"))
    members = {l.strip() for l in open(args.members) if l.strip()}
    res = measure(Path(args.pdb_dir), Path(args.alignment), lo, hi, members)

    if not res["n_focal"] or not res["n_other"]:
        print("Not enough structures on one side to make the comparison — "
              "the control is what makes this evidence, so refusing to "
              "report a one-sided number.", file=sys.stderr)
    Path(args.out).write_text(json.dumps(res, indent=2))
    print(f"focal  n={res['n_focal']} mean delta {res['delta']}")
    print(f"other  n={res['n_other']} mean delta {res['delta_other']} "
          f"(median {res['median_other']})")
    print(f"Mann-Whitney p = {res['p']}")
    print(f"-> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
