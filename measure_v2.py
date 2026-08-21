#!/usr/bin/env python3
"""v2 acceptance metric: per-species unplaced rate vs the v1 baseline.

Baseline (resume.md, converged 5sp v1 run, after HMMER rescue):
  species  total   unplaced_after_R10  rate     final_placed
  Mcry     25226   3765                14.9%    88.0%
The headline number is Mcry's unplaced rate; v2 must beat 14.9%.

Usage: measure_v2.py <run_dir>
"""
import sys
from collections import Counter
from pathlib import Path

BASELINE = {  # species -> (total_genes, unplaced_v1, rate_v1)
    "Mcry": (25226, 3765, 0.149),
    "Cgig": (29163, 1698, 0.058),
    "CgigH": (27583, 3497, 0.127),
    "Ococ": (33745, 2639, 0.078),
    "Obas": (28244, 2364, 0.084),
}


def species_of(gene):
    return gene.split("_", 1)[0]


def main():
    run = Path(sys.argv[1])
    placed = Counter()
    n_fam = 0
    with open(run / "summary.tsv") as f:
        next(f)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            n_fam += 1
            for g in fields[4].split(","):
                if g:
                    placed[species_of(g)] += 1

    unplaced = Counter()
    fa = run / "unplaced_proteins.fa"
    if fa.exists():
        for line in open(fa):
            if line.startswith(">"):
                unplaced[species_of(line[1:].split()[0])] += 1

    print(f"families: {n_fam}\n")
    print(f"{'species':8s} {'total':>7s} {'placed':>7s} {'unplaced':>9s} "
          f"{'rate':>7s} {'v1 rate':>8s} {'delta':>8s}")
    for sp, (total, _u1, r1) in BASELINE.items():
        p, u = placed.get(sp, 0), unplaced.get(sp, 0)
        seen = p + u
        rate = u / total if total else 0.0
        flag = ""
        if seen and abs(seen - total) > total * 0.01:
            flag = f"  [!] accounted {seen} != {total}"
        print(f"{sp:8s} {total:7d} {p:7d} {u:9d} {rate:6.1%} {r1:7.1%} "
              f"{rate - r1:+7.1%}{flag}")

    mc = unplaced.get("Mcry", 0) / BASELINE["Mcry"][0]
    verdict = "PASS" if mc < 0.149 else "FAIL"
    print(f"\nHEADLINE Mcry unplaced {mc:.1%} vs baseline 14.9% -> {verdict}")


if __name__ == "__main__":
    main()
