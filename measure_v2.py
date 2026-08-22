#!/usr/bin/env python3
"""v2 acceptance metric: per-species unplaced rate vs the v1 baseline.

Baseline (resume.md, converged 5sp v1 run, after HMMER rescue):
  species  total   unplaced_after_R10  rate     final_placed
  Mcry     25226   3765                14.9%    88.0%
The headline number is Mcry's unplaced rate; v2 must beat 14.9%.

Issue #36: the raw per-species rate is CONFOUNDED by annotation policy.
Cgig.pep.fa has a hard floor at 151 aa (zero proteins under 100 aa) while
Mcry's shortest protein is 3 aa with 3,374 under 100 aa. Short models fail to
cluster far more often, so comparing raw rates across species compares
annotators as much as it compares the pipeline. Pass --pep-dir to get the
length-stratified rates, which are the ones that can be compared, plus the
per-species minimum protein length that makes the confound visible.

Usage: measure_v2.py <run_dir> [--pep-dir DIR] [--floor 100]
"""
import argparse
from collections import Counter
from pathlib import Path

BASELINE = {  # species -> (total_genes, unplaced_v1, rate_v1)
    "Mcry": (25226, 3765, 0.149),
    "Cgig": (29163, 1698, 0.058),
    "CgigH": (27583, 3497, 0.127),
    "Ococ": (33745, 2639, 0.078),
    "Obas": (28244, 2364, 0.084),
}

DEFAULT_FLOOR = 100
# One annotation keeping zero genes below the floor while another keeps many
# means the two strata are not the same population.
SHORT_FRACTION_RATIO_LIMIT = 3.0


def species_of(gene):
    return gene.split("_", 1)[0]


def read_lengths(pep_dir):
    """{species: {gene: protein_length}} from <pep_dir>/<species>.pep.fa."""
    lengths = {}
    for path in sorted(Path(pep_dir).glob("*.pep.fa")):
        species = path.name[: -len(".pep.fa")]
        genes = {}
        name = None
        total = 0
        with open(path) as f:
            for line in f:
                if line.startswith(">"):
                    if name:
                        genes[name] = total
                    name = line[1:].split()[0]
                    total = 0
                else:
                    total += len(line.strip().rstrip("*"))
        if name:
            genes[name] = total
        lengths[species] = genes
    return lengths


def count_placed_and_unplaced(run):
    placed, unplaced = Counter(), Counter()
    n_families = 0
    with open(run / "summary.tsv") as f:
        next(f)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            n_families += 1
            for gene in fields[4].split(","):
                if gene:
                    placed[species_of(gene)] += 1

    unplaced_genes = set()
    fasta = run / "unplaced_proteins.fa"
    if fasta.exists():
        for line in open(fasta):
            if line.startswith(">"):
                gene = line[1:].split()[0]
                unplaced_genes.add(gene)
                unplaced[species_of(gene)] += 1
    return n_families, placed, unplaced, unplaced_genes


def print_raw_table(placed, unplaced):
    print(f"{'species':8s} {'total':>7s} {'placed':>7s} {'unplaced':>9s} "
          f"{'rate':>7s} {'v1 rate':>8s} {'delta':>8s}")
    for species, (total, _unplaced_v1, rate_v1) in BASELINE.items():
        n_placed, n_unplaced = placed.get(species, 0), unplaced.get(species, 0)
        seen = n_placed + n_unplaced
        rate = n_unplaced / total if total else 0.0
        flag = ""
        if seen and abs(seen - total) > total * 0.01:
            flag = f"  [!] accounted {seen} != {total}"
        print(f"{species:8s} {total:7d} {n_placed:7d} {n_unplaced:9d} "
              f"{rate:6.1%} {rate_v1:7.1%} {rate - rate_v1:+7.1%}{flag}")


def short_fractions(lengths, floor):
    """{species: (min_length, n_short, n_total)}."""
    return {
        species: (min(genes.values()) if genes else 0,
                  sum(1 for n in genes.values() if n < floor),
                  len(genes))
        for species, genes in lengths.items()
    }


def confound_verdict(profile):
    """Reasons the raw cross-species rate cannot be compared, if any."""
    fractions = {sp: (n_short / total if total else 0.0)
                 for sp, (_m, n_short, total) in profile.items()}
    if len(fractions) < 2:
        return []
    with_short = [sp for sp, f in fractions.items() if f > 0]
    without = [sp for sp, f in fractions.items() if f == 0]
    if with_short and without:
        return [f"species {sorted(without)} have NO proteins below the floor "
                f"while {sorted(with_short)} do — their annotators applied "
                "different length policies, so the raw rates are not comparable"]
    positive = [f for f in fractions.values() if f > 0]
    if positive and max(positive) / min(positive) > SHORT_FRACTION_RATIO_LIMIT:
        return [f"short-protein fraction varies {max(positive):.1%} vs "
                f"{min(positive):.1%} across species (> "
                f"{SHORT_FRACTION_RATIO_LIMIT:g}x) — annotation policy differs"]
    return []


def print_confound_check(profile, floor):
    print(f"\n=== annotation profile (floor {floor} aa) ===")
    print(f"{'species':8s} {'min_len':>8s} {'n_short':>8s} {'total':>8s} {'pct_short':>10s}")
    for species in sorted(profile):
        min_len, n_short, total = profile[species]
        pct = n_short / total if total else 0.0
        print(f"{species:8s} {min_len:8d} {n_short:8d} {total:8d} {pct:9.1%}")
    for reason in confound_verdict(profile):
        print(f"[!] CONFOUND: {reason}")


def print_stratified(lengths, unplaced_genes, floor):
    print(f"\n=== length-stratified unplaced rate (floor {floor} aa) ===")
    print(f"{'species':8s} {'stratum':6s} {'total':>8s} {'unplaced':>9s} {'rate':>8s}")
    comparable = {}
    for species in sorted(lengths):
        genes = lengths[species]
        for stratum, keep in (("short", True), ("long", False)):
            members = [g for g, n in genes.items() if (n < floor) == keep]
            n_unplaced = sum(1 for g in members if g in unplaced_genes)
            rate = n_unplaced / len(members) if members else 0.0
            shown = f"{rate:.1%}" if members else "-"
            print(f"{species:8s} {stratum:6s} {len(members):8d} "
                  f"{n_unplaced:9d} {shown:>8s}")
            if stratum == "long" and members:
                comparable[species] = rate
    if comparable:
        summary = " | ".join(f"{sp} {r:.1%}" for sp, r in sorted(comparable.items()))
        print(f"\nCOMPARABLE (>= {floor} aa only) unplaced rate: {summary}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir")
    parser.add_argument("--pep-dir", help="proteome FASTAs, enables stratification")
    parser.add_argument("--floor", type=int, default=DEFAULT_FLOOR,
                        help="protein length separating the short/long strata")
    args = parser.parse_args(argv)

    run = Path(args.run_dir)
    n_families, placed, unplaced, unplaced_genes = count_placed_and_unplaced(run)

    print(f"families: {n_families}\n")
    print_raw_table(placed, unplaced)

    if args.pep_dir:
        lengths = read_lengths(args.pep_dir)
        print_confound_check(short_fractions(lengths, args.floor), args.floor)
        print_stratified(lengths, unplaced_genes, args.floor)
    else:
        print("\n[i] raw rates only. Pass --pep-dir to stratify by protein "
              "length — annotation length policy differs between species and "
              "confounds the raw comparison (issue #36)")

    mcry_rate = unplaced.get("Mcry", 0) / BASELINE["Mcry"][0]
    verdict = "PASS" if mcry_rate < 0.149 else "FAIL"
    print(f"\nHEADLINE Mcry unplaced {mcry_rate:.1%} vs baseline 14.9% -> {verdict}")


if __name__ == "__main__":
    main()
