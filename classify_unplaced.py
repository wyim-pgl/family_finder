#!/usr/bin/env python3
"""Classify why each unplaced gene stayed unplaced (issue #36).

The Mcry unplaced rate (13.4%) is 2-5x the cacti. resume.md already rejected
pep/CDS integrity, DIAMOND sensitivity, pruning and taxon sampling, leaving
MCL/graph fragmentation. This tool turns that verdict into per-gene counts so
the next intervention is chosen on evidence rather than on the leftover
hypothesis.

Every unplaced gene gets exactly one verdict, from the strongest evidence
found across all rounds:

  PRUNED           was in a cross-species orthogroup at or above
                   min_orthogroup_size, so it was aligned and treed and
                   pruning ejected it
  SPLINTER         was in a cross-species orthogroup BELOW the size gate --
                   it clustered with other species but the cluster was too
                   small to work on. Graph fragmentation.
  LINEAGE_ONLY     only ever clustered with its own species
  SINGLETON_HIT    never clustered, but DIAMOND found a cross-species
                   homolog. MCL cut an edge that existed.
  SINGLETON_NOHIT  never clustered and no cross-species homolog. A real
                   orphan as far as this panel can tell.

SPLINTER + SINGLETON_HIT are the graph-fragmentation share; SINGLETON_NOHIT
is the honest ceiling that no clustering change can move.
"""
import argparse
import glob
import gzip
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Set

MIN_ORTHOGROUP_SIZE = 4

VERDICTS = ("PRUNED", "SPLINTER", "LINEAGE_ONLY", "SINGLETON_HIT", "SINGLETON_NOHIT")

ROLLUP = {
    "splinter_or_graph_cut": ("SPLINTER", "SINGLETON_HIT"),
    "lineage_specific": ("LINEAGE_ONLY",),
    "true_orphan": ("SINGLETON_NOHIT",),
    "pruned": ("PRUNED",),
}


@dataclass(frozen=True)
class Observation:
    """One orthogroup a gene was seen in, summarised from that gene's view."""
    n_genes: int
    n_other_species: int


def species_of(gene_id: str) -> str:
    return gene_id.split("_", 1)[0]


# ------------------------------------------------------------ parsing ----

def parse_orthogroups(path) -> Dict[str, List[str]]:
    """OrthoFinder Orthogroups.tsv -> {orthogroup: [gene, ...]}.

    Species columns are comma-space separated and may be empty; empty cells
    must not turn into phantom members.
    """
    groups: Dict[str, List[str]] = {}
    with open(path) as f:
        f.readline()  # header: Orthogroup + one column per species
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            members: List[str] = []
            for cell in fields[1:]:
                members.extend(g.strip() for g in cell.split(",") if g.strip())
            groups[fields[0]] = members
    return groups


def observe(members: Sequence[str], species: str) -> Observation:
    return Observation(
        n_genes=len(members),
        n_other_species=sum(1 for g in members if species_of(g) != species),
    )


# ----------------------------------------------------- classification ----

def classify_gene(observations: Iterable[Observation],
                  has_cross_species_hit: bool,
                  min_orthogroup_size: int = MIN_ORTHOGROUP_SIZE) -> str:
    """Strongest evidence across rounds wins; see module docstring."""
    observations = list(observations)
    cross = [o for o in observations if o.n_other_species > 0 and o.n_genes >= 2]
    if any(o.n_genes >= min_orthogroup_size for o in cross):
        return "PRUNED"
    if cross:
        return "SPLINTER"
    if any(o.n_genes >= 2 for o in observations):
        return "LINEAGE_ONLY"
    return "SINGLETON_HIT" if has_cross_species_hit else "SINGLETON_NOHIT"


def rollup(counts: Dict[str, int]) -> Dict[str, int]:
    return {bucket: sum(counts.get(v, 0) for v in verdicts)
            for bucket, verdicts in ROLLUP.items()}


# ------------------------------------------------------------ diamond ----

def load_orthofinder_ids(working_dir):
    """-> ({internal_id: gene_id}, {species_index: species_name})."""
    genes = {}
    with open(f"{working_dir}/SequenceIDs.txt") as f:
        for line in f:
            if ": " not in line:
                continue
            iid, gid = line.rstrip("\n").split(": ", 1)
            genes[iid] = gid.split()[0]
    species = {}
    with open(f"{working_dir}/SpeciesIDs.txt") as f:
        for line in f:
            if ": " not in line:
                continue
            idx, name = line.rstrip("\n").split(": ", 1)
            species[idx] = name.split(".")[0]
    return genes, species


def cross_species_hits(working_dir, species: str, wanted: Set[str],
                       max_evalue: float) -> Set[str]:
    """Genes in `wanted` with a DIAMOND hit to a DIFFERENT species."""
    genes, species_names = load_orthofinder_ids(working_dir)
    query_idx = [i for i, name in species_names.items() if name == species]
    if not query_idx:
        raise ValueError(
            f"{species} not among clustered species {sorted(species_names.values())}"
        )
    found: Set[str] = set()
    for qi in query_idx:
        for path in sorted(glob.glob(f"{working_dir}/Blast{qi}_*.txt.gz")):
            subject_idx = os.path.basename(path).split("_")[1].split(".")[0]
            if species_names.get(subject_idx) == species:
                continue
            with gzip.open(path, "rt") as f:
                for line in f:
                    fields = line.split("\t")
                    if len(fields) < 11:
                        continue
                    gene = genes.get(fields[0])
                    if gene is None or gene in found or gene not in wanted:
                        continue
                    try:
                        if float(fields[10]) <= max_evalue:
                            found.add(gene)
                    except ValueError:
                        continue
    return found


# --------------------------------------------------------------- main ----

def gather_observations(run_dir, species, genes: Set[str]):
    """{gene: [Observation, ...]} over every round's orthogroups."""
    per_gene = defaultdict(list)
    rounds = sorted(glob.glob(f"{run_dir}/round_*/orthofinder/Results_*/Orthogroups/"
                              "Orthogroups.tsv"))
    if not rounds:
        sys.exit(f"No Orthogroups.tsv under {run_dir}/round_*/")
    for path in rounds:
        for members in parse_orthogroups(path).values():
            obs = None
            for gene in members:
                if gene in genes:
                    obs = obs or observe(members, species)
                    per_gene[gene].append(obs)
    print(f"rounds scanned: {len(rounds)}")
    return per_gene


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--species", required=True, help="species prefix to classify")
    ap.add_argument("--unplaced", help="default <run-dir>/unplaced_proteins.fa")
    ap.add_argument("--working-dir", help="round-1 OrthoFinder WorkingDirectory")
    ap.add_argument("--min-orthogroup-size", type=int, default=MIN_ORTHOGROUP_SIZE)
    ap.add_argument("--evalue", type=float, default=1e-3)
    ap.add_argument("--out", required=True, help="per-gene TSV")
    args = ap.parse_args(argv)

    unplaced_fa = args.unplaced or f"{args.run_dir}/unplaced_proteins.fa"
    genes = {line[1:].split()[0] for line in open(unplaced_fa) if line.startswith(">")}
    genes = {g for g in genes if species_of(g) == args.species}
    print(f"unplaced {args.species} genes: {len(genes)}")
    if not genes:
        sys.exit(f"No unplaced genes for {args.species} in {unplaced_fa}")

    per_gene = gather_observations(args.run_dir, args.species, genes)
    print(f"of those, seen in some orthogroup: {len(per_gene)}")

    working = args.working_dir
    if not working:
        hits = sorted(glob.glob(f"{args.run_dir}/round_01/orthofinder/Results_*/"
                                "WorkingDirectory"))
        working = hits[-1] if hits else None
    never_clustered = {g for g in genes
                       if all(o.n_genes < 2 for o in per_gene.get(g, []))}
    with_hit: Set[str] = set()
    if working and never_clustered:
        print(f"DIAMOND scan for {len(never_clustered)} never-clustered genes: {working}")
        with_hit = cross_species_hits(working, args.species, never_clustered, args.evalue)
        print(f"  with a cross-species hit (E<={args.evalue:g}): {len(with_hit)}")
    elif not working:
        print("WARNING: no WorkingDirectory found — every never-clustered gene "
              "will be reported as SINGLETON_NOHIT, which understates graph cuts")

    counts = Counter()
    with open(args.out, "w") as f:
        f.write("gene\tverdict\tmax_og_size\tmax_other_species\tn_rounds_seen\n")
        for gene in sorted(genes):
            obs = per_gene.get(gene, [])
            verdict = classify_gene(obs, gene in with_hit, args.min_orthogroup_size)
            counts[verdict] += 1
            f.write("{}\t{}\t{}\t{}\t{}\n".format(
                gene, verdict,
                max((o.n_genes for o in obs), default=0),
                max((o.n_other_species for o in obs), default=0),
                len(obs)))

    total = sum(counts.values())
    print(f"\n=== {args.species}: {total} unplaced genes ===")
    for verdict in VERDICTS:
        n = counts.get(verdict, 0)
        print(f"  {verdict:<16} {n:>7}  {100.0 * n / total:5.1f}%")
    print("\n--- rollup ---")
    for bucket, n in rollup(counts).items():
        print(f"  {bucket:<22} {n:>7}  {100.0 * n / total:5.1f}%")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
