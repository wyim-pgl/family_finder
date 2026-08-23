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

That verdict is the gene's `primary_reason`. It is not the whole story: a
60-aa model with a broken CDS and no cross-species hit is all three things at
once. Model quality therefore rides alongside as additive `flags` --
SHORT_PROTEIN (below the length floor), INVALID_CDS (the CDS is not a clean
ORF), NO_CDS (no CDS record, which is not the same as knowing it is broken).
Flags never change the primary reason; the verdict distribution recorded in
resume.md 1.5 has to stay put.
"""
import argparse
import glob
import gzip
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

MIN_ORTHOGROUP_SIZE = 4
# Protein length below which a model is too short to compare across
# annotations (issue #36: Cgig keeps nothing under 151 aa, Mcry keeps 3 aa).
DEFAULT_PROTEIN_FLOOR = 100

VERDICTS = ("PRUNED", "SPLINTER", "LINEAGE_ONLY", "SINGLETON_HIT", "SINGLETON_NOHIT")

FLAGS = ("SHORT_PROTEIN", "INVALID_CDS", "NO_CDS")

START_CODON = "ATG"
STOP_CODONS = ("TAA", "TAG", "TGA")

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


@dataclass(frozen=True)
class CdsIntegrity:
    """Whether a CDS reads as a clean ORF. `present` is false when no CDS
    record exists for the gene at all — an unknown, not a defect."""
    present: bool = False
    has_start: bool = False
    has_stop: bool = False
    in_frame: bool = False
    no_internal_stop: bool = False

    @property
    def complete(self) -> bool:
        return (self.present and self.has_start and self.has_stop
                and self.in_frame and self.no_internal_stop)

    @property
    def label(self) -> str:
        if not self.present:
            return "absent"
        return "complete" if self.complete else "incomplete"


@dataclass(frozen=True)
class Row:
    """One unplaced gene: why it stayed unplaced, and what it is."""
    gene: str
    primary_reason: str
    flags: Tuple[str, ...]
    protein_length: int
    cds_status: str
    max_og_size: int
    max_other_species: int
    n_rounds_seen: int
    comparable: bool


def species_of(gene_id: str) -> str:
    return gene_id.split("_", 1)[0]


# ------------------------------------------------------ CDS integrity ----

def assess_cds(sequence: Optional[str]) -> CdsIntegrity:
    """Start codon, terminal stop, length divisible by three, no internal stop.

    An out-of-frame CDS is still scanned for internal stops over its whole
    codons — a frameshifted model usually shows both, and reporting only the
    frame hides that.
    """
    if sequence is None:
        return CdsIntegrity()
    seq = "".join(sequence.split()).upper().replace("-", "")
    if not seq:
        return CdsIntegrity()
    codons = [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]
    body = codons[:-1] if codons and codons[-1] in STOP_CODONS else codons
    return CdsIntegrity(
        present=True,
        has_start=seq.startswith(START_CODON),
        has_stop=bool(codons) and codons[-1] in STOP_CODONS,
        in_frame=len(seq) % 3 == 0,
        no_internal_stop=not any(c in STOP_CODONS for c in body),
    )


def flags_for(protein_length: int, integrity: CdsIntegrity,
              floor: int = DEFAULT_PROTEIN_FLOOR) -> Tuple[str, ...]:
    """Additive model-quality flags. Independent of the primary reason."""
    flags = []
    if protein_length < floor:
        flags.append("SHORT_PROTEIN")
    if not integrity.present:
        flags.append("NO_CDS")
    elif not integrity.complete:
        flags.append("INVALID_CDS")
    return tuple(flags)


def is_comparable(protein_length: int, cds_complete: bool,
                  floor: int = DEFAULT_PROTEIN_FLOOR) -> bool:
    """Long enough AND a clean ORF — the only population whose unplaced rate
    can be compared across annotations (issue #36).

    Takes the CDS verdict as a bool so that callers holding one boolean per
    gene (measure_v2 over whole proteomes) and callers holding a full
    CdsIntegrity (per-gene rows here) share one definition of comparable.
    """
    return protein_length >= floor and cds_complete


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


def read_fasta(path, wanted: Optional[Set[str]] = None) -> Dict[str, str]:
    """{id: sequence}, keeping only `wanted` ids when given.

    Deliberately not Biopython: this tool runs on the cluster head node next
    to the outputs it reads, with nothing installed.
    """
    seqs: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None and (wanted is None or name in wanted):
        seqs[name] = "".join(chunks)
    return seqs


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


def classify_all(genes: Iterable[str],
                 observations: Dict[str, List[Observation]],
                 with_hit: Set[str],
                 protein_lengths: Dict[str, int],
                 cds_seqs: Dict[str, str],
                 floor: int = DEFAULT_PROTEIN_FLOOR,
                 min_orthogroup_size: int = MIN_ORTHOGROUP_SIZE) -> List[Row]:
    """One Row per gene — a primary reason for every gene, plus its flags."""
    rows = []
    for gene in sorted(genes):
        obs = observations.get(gene, [])
        integrity = assess_cds(cds_seqs.get(gene))
        length = protein_lengths.get(gene, 0)
        rows.append(Row(
            gene=gene,
            primary_reason=classify_gene(obs, gene in with_hit, min_orthogroup_size),
            flags=flags_for(length, integrity, floor),
            protein_length=length,
            cds_status=integrity.label,
            max_og_size=max((o.n_genes for o in obs), default=0),
            max_other_species=max((o.n_other_species for o in obs), default=0),
            n_rounds_seen=len(obs),
            comparable=is_comparable(length, integrity.complete, floor),
        ))
    return rows


def genes_without_a_reason(rows: Sequence[Row], genes: Iterable[str]) -> Set[str]:
    """Genes that came out with no usable primary reason.

    The acceptance criterion for issue #36 is that this is empty; it is
    counted rather than assumed.
    """
    reasoned = {r.gene for r in rows if r.primary_reason in VERDICTS}
    return set(genes) - reasoned


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
    ap.add_argument("--cds", help="CDS FASTA, default <run-dir>/unplaced_cds.fa")
    ap.add_argument("--cds-dir", help="proteome CDS dir holding <species>.cds.fa; "
                                      "used when --cds is absent")
    ap.add_argument("--floor", type=int, default=DEFAULT_PROTEIN_FLOOR,
                    help="protein length below which a model is flagged short")
    ap.add_argument("--out", required=True, help="per-gene TSV")
    args = ap.parse_args(argv)

    unplaced_fa = args.unplaced or f"{args.run_dir}/unplaced_proteins.fa"
    proteins = read_fasta(unplaced_fa)
    genes = {g for g in proteins if species_of(g) == args.species}
    lengths = {g: len(proteins[g].rstrip("*")) for g in genes}
    print(f"unplaced {args.species} genes: {len(genes)}")
    if not genes:
        sys.exit(f"No unplaced genes for {args.species} in {unplaced_fa}")

    candidates = [args.cds] if args.cds else [f"{args.run_dir}/unplaced_cds.fa"]
    if args.cds_dir:
        candidates.append(f"{args.cds_dir}/{args.species}.cds.fa")
    cds_path = next((p for p in candidates if os.path.exists(p)), None)
    cds_seqs: Dict[str, str] = {}
    if cds_path:
        cds_seqs = read_fasta(cds_path, wanted=genes)
        print(f"CDS read for {len(cds_seqs)}/{len(genes)} genes: {cds_path}")
    else:
        print("WARNING: no CDS FASTA found — every gene will be flagged NO_CDS "
              "and no comparable subset can be reported (pass --cds/--cds-dir)")

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

    rows = classify_all(genes, per_gene, with_hit, lengths, cds_seqs,
                        args.floor, args.min_orthogroup_size)
    missing = genes_without_a_reason(rows, genes)
    with open(args.out, "w") as f:
        f.write("gene\tprimary_reason\tflags\tprotein_length\tcds_status\t"
                "max_og_size\tmax_other_species\tn_rounds_seen\tcomparable\n")
        for r in rows:
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                r.gene, r.primary_reason, ",".join(r.flags) or "-",
                r.protein_length, r.cds_status, r.max_og_size,
                r.max_other_species, r.n_rounds_seen,
                "yes" if r.comparable else "no"))

    counts = Counter(r.primary_reason for r in rows)
    total = len(rows)
    print(f"\n=== {args.species}: {total} unplaced genes ===")
    for verdict in VERDICTS:
        n = counts.get(verdict, 0)
        print(f"  {verdict:<16} {n:>7}  {100.0 * n / total:5.1f}%")
    print("\n--- rollup ---")
    for bucket, n in rollup(counts).items():
        print(f"  {bucket:<22} {n:>7}  {100.0 * n / total:5.1f}%")

    print(f"\n--- flags (additive; they do NOT change the reason above) ---")
    flag_counts = Counter(flag for r in rows for flag in r.flags)
    for flag in FLAGS:
        n = flag_counts.get(flag, 0)
        print(f"  {flag:<16} {n:>7}  {100.0 * n / total:5.1f}%")

    comparable = [r for r in rows if r.comparable]
    print(f"\n--- comparable subset (>= {args.floor} aa AND complete CDS) ---")
    print(f"  {len(comparable)} of {total} unplaced genes ({100.0 * len(comparable) / total:.1f}%)")
    if comparable:
        for verdict in VERDICTS:
            n = sum(1 for r in comparable if r.primary_reason == verdict)
            print(f"  {verdict:<16} {n:>7}  {100.0 * n / len(comparable):5.1f}%")

    print(f"\nunplaced genes with no primary reason: {len(missing)} "
          f"({100.0 * len(missing) / total:.1f}%)")
    print(f"\nwrote {args.out}")
    if missing:
        sys.exit(f"ERROR: {len(missing)} genes came out with no reason, "
                 f"e.g. {sorted(missing)[:5]} — the classification is incomplete")


if __name__ == "__main__":
    main()
