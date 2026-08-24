#!/usr/bin/env python3
"""Build a codon supermatrix from a family_finder run's own single-copy markers.

Issue #14 estimated the 5-species tree this way; issue #41 needs the same for
the 15-species panel, where the committed tree still carries dummy branch
lengths (all 1.0) and would silently disable species-aware pruning.

Two ways to pick markers:

  --markers genecount   OrthoFinder round-1 orthogroups with exactly one gene
                        in every species. Independent of the species tree and
                        of pruning, so it works even when the run that produced
                        the orthogroups used a dummy tree.
  --markers summary     Confirmed families with one gene per species, read from
                        summary.tsv. This is what produced the committed 5sp
                        tree; it depends on pruning having been sane.

Per-locus codon alignments are reused from
``round_NN/orthogroups/<OG>/codon.afa`` and only realigned (MAFFT + pal2nal)
when that file is missing or unusable, which needs --pep-dir/--cds-dir.

Output: <out>/supermatrix.fa and <out>/partitions.txt, ready for
``iqtree2 -s supermatrix.fa -m GTR+F+G4``.
"""
import argparse
import os
import random
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

FAMILY_ID_RE = re.compile(r"^R(\d+)_(OG\d+)$")


# ---------------------------------------------------------------- fasta ----

def read_fasta(path) -> Dict[str, str]:
    """Minimal FASTA reader. Keeps the id up to the first whitespace."""
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def species_of(gene_id: str) -> str:
    """Gene IDs are SpeciesPrefix_GeneID everywhere in this pipeline."""
    return gene_id.split("_", 1)[0]


# ------------------------------------------------------- marker selection --

@dataclass(frozen=True)
class Marker:
    family: str
    round: int
    members: List[str]


def strict_single_copy_orthogroups(genecount_tsv, species: Sequence[str]) -> List[str]:
    """Orthogroups with exactly one gene in every requested species.

    Reads OrthoFinder's ``Orthogroups.GeneCount.tsv``. The trailing ``Total``
    column is not a species and must not be matched against.
    """
    with open(genecount_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        missing = [s for s in species if s not in header[1:]]
        if missing:
            raise ValueError(
                f"Species absent from {genecount_tsv} header: {missing} "
                f"(header has {header[1:]})"
            )
        columns = [header.index(s) for s in species]
        selected = []
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if all(fields[i] == "1" for i in columns):
                selected.append(fields[0])
    return selected


def single_copy_families(summary_tsv, species: Sequence[str]) -> List[Marker]:
    """Confirmed families holding exactly one gene from each requested species."""
    want = sorted(species)
    markers: List[Marker] = []
    with open(summary_tsv) as f:
        f.readline()  # header
        for line in f:
            if not line.strip():
                continue
            # summary_v3.tsv adds merged_from and cluster_id, so take the
            # first five and ignore the rest rather than raising on the table
            # that now ships (#50).
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            family, rnd, n_genes, n_species, gene_list = fields[:5]
            if int(n_genes) != len(want) or int(n_species) != len(want):
                continue
            members = gene_list.split(",")
            if sorted(species_of(m) for m in members) != want:
                continue
            markers.append(Marker(family, int(rnd), members))
    return markers


def orthogroup_dir(run_dir: str, family: str) -> str:
    """`R2_OG0000004` -> `<run_dir>/round_02/orthogroups/OG0000004`."""
    m = FAMILY_ID_RE.match(family)
    if not m:
        raise ValueError(f"Not a family id of the form R<round>_OG<n>: {family!r}")
    return f"{run_dir}/round_{int(m.group(1)):02d}/orthogroups/{m.group(2)}"


# ------------------------------------------------------ codon filtering ----

def drop_all_gap_codons(by_species: Dict[str, str]) -> Tuple[Dict[str, str], int]:
    """Remove codon columns that are gaps in every species.

    Returns the filtered alignment and the number of codons kept. A trailing
    partial codon is discarded. Ragged input raises rather than being
    truncated into a frame shift.
    """
    seqs = {sp: s.upper() for sp, s in by_species.items()}
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) > 1:
        raise ValueError(f"Alignment has unequal sequence lengths: {sorted(lengths)}")
    length = lengths.pop() if lengths else 0
    usable = length - length % 3
    kept = [i for i in range(0, usable, 3)
            if any(s[i:i + 3] != "---" for s in seqs.values())]
    filtered = {sp: "".join(s[i:i + 3] for i in kept) for sp, s in seqs.items()}
    return filtered, len(kept)


def concatenate(loci: Sequence[Tuple[str, Dict[str, str]]],
                species: Sequence[str]) -> Tuple[Dict[str, str], List[Tuple[str, int, int]]]:
    """Concatenate per-locus alignments, returning 1-based inclusive partitions."""
    parts: Dict[str, List[str]] = {sp: [] for sp in species}
    partitions: List[Tuple[str, int, int]] = []
    pos = 0
    for name, by_sp in loci:
        absent = [sp for sp in species if sp not in by_sp]
        if absent:
            raise ValueError(f"Locus {name} is missing species: {absent}")
        width = len(by_sp[species[0]])
        for sp in species:
            parts[sp].append(by_sp[sp])
        partitions.append((name, pos + 1, pos + width))
        pos += width
    return {sp: "".join(chunks) for sp, chunks in parts.items()}, partitions


def format_partitions(partitions: Sequence[Tuple[str, int, int]]) -> str:
    return "".join(f"DNA, {name} = {a}-{b}\n" for name, a, b in partitions)


# ------------------------------------------------------------- realign ----

class Realigner:
    """MAFFT + pal2nal fallback for loci with no usable saved codon.afa."""

    def __init__(self, pep_dir, cds_dir, species, mafft, pal2nal):
        self.pep_dir, self.cds_dir = pep_dir, cds_dir
        self.species, self.mafft, self.pal2nal = species, mafft, pal2nal
        self._pep = self._cds = None

    @property
    def available(self) -> bool:
        return bool(self.pep_dir and self.cds_dir)

    def _pools(self):
        if self._pep is None:
            self._pep, self._cds = {}, {}
            for sp in self.species:
                self._pep.update(read_fasta(f"{self.pep_dir}/{sp}.pep.fa"))
                self._cds.update(read_fasta(f"{self.cds_dir}/{sp}.cds.fa"))
        return self._pep, self._cds

    def __call__(self, members: Sequence[str]) -> Optional[Dict[str, str]]:
        if not self.available:
            return None
        pep, cds = self._pools()
        if any(m not in pep or m not in cds for m in members):
            return None
        with tempfile.TemporaryDirectory(prefix="marker_") as td:
            pf, cf, af = f"{td}/p.fa", f"{td}/c.fa", f"{td}/p.afa"
            with open(pf, "w") as f:
                f.writelines(f">{m}\n{pep[m].rstrip('*')}\n" for m in members)
            with open(cf, "w") as f:
                f.writelines(f">{m}\n{cds[m]}\n" for m in members)
            try:
                aln = subprocess.run([self.mafft, "--auto", "--quiet", pf],
                                     capture_output=True, text=True, timeout=600)
                if aln.returncode != 0:
                    return None
                with open(af, "w") as f:
                    f.write(aln.stdout)
                p2n = subprocess.run(["perl", self.pal2nal, af, cf, "-output", "fasta"],
                                     capture_output=True, text=True, timeout=600)
                if p2n.returncode != 0 or not p2n.stdout.strip():
                    return None
                out = f"{td}/codon.afa"
                with open(out, "w") as f:
                    f.write(p2n.stdout)
                return read_fasta(out)
            except subprocess.TimeoutExpired:
                return None


# ---------------------------------------------------------------- main ----

def _locus_by_species(aln: Dict[str, str], species: Sequence[str],
                      members: Optional[Sequence[str]]) -> Optional[Dict[str, str]]:
    """Subset an alignment to one sequence per species, or None if impossible."""
    if members is not None:
        if not all(g in aln for g in members):
            return None
        aln = {g: aln[g] for g in members}
    by_sp: Dict[str, str] = {}
    for gene, seq in aln.items():
        sp = species_of(gene)
        if sp in by_sp:
            return None  # more than one gene from a species: not a marker
        by_sp[sp] = seq
    if any(sp not in by_sp for sp in species):
        return None
    return {sp: by_sp[sp] for sp in species}


def collect_markers(args, species) -> List[Tuple[str, str, Optional[List[str]]]]:
    """Return (name, orthogroup_dir, members-or-None) for every candidate marker."""
    if args.markers == "genecount":
        genecount = args.genecount or (
            f"{args.run_dir}/round_01/orthofinder/Results_*/Orthogroups/"
            "Orthogroups.GeneCount.tsv")
        if "*" in genecount:
            import glob
            hits = sorted(glob.glob(genecount))
            if not hits:
                sys.exit(f"No GeneCount.tsv matched {genecount}")
            genecount = hits[-1]
        print(f"GeneCount: {genecount}")
        ogs = strict_single_copy_orthogroups(genecount, species)
        print(f"strict single-copy orthogroups: {len(ogs)}")
        return [(og, f"{args.run_dir}/round_01/orthogroups/{og}", None) for og in ogs]

    families = single_copy_families(args.summary or f"{args.run_dir}/summary.tsv", species)
    print(f"single-copy families: {len(families)}")
    return [(m.family, orthogroup_dir(args.run_dir, m.family), m.members) for m in families]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, help="pipeline output dir (output_15sp)")
    ap.add_argument("--species", required=True, help="comma-separated species prefixes")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--markers", choices=("genecount", "summary"), default="genecount")
    ap.add_argument("--genecount", help="explicit Orthogroups.GeneCount.tsv path")
    ap.add_argument("--summary", help="explicit summary.tsv path")
    ap.add_argument("--cap", type=int, default=0, help="0 = use every marker")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--pep-dir")
    ap.add_argument("--cds-dir")
    ap.add_argument("--mafft", default="mafft")
    ap.add_argument("--pal2nal", default="pal2nal.pl")
    args = ap.parse_args(argv)

    species = [s for s in args.species.split(",") if s]
    os.makedirs(args.out, exist_ok=True)
    os.environ.pop("MAFFT_BINARIES", None)  # conda conflict, see CLAUDE.md

    candidates = collect_markers(args, species)
    if args.cap and len(candidates) > args.cap:
        random.seed(args.seed)
        candidates = sorted(random.sample(candidates, args.cap))
        print(f"capped to {args.cap} (seed {args.seed})")

    realign = Realigner(args.pep_dir, args.cds_dir, species, args.mafft, args.pal2nal)
    loci: List[Tuple[str, Dict[str, str]]] = []
    reused = realigned = dropped = 0
    for name, og_dir, members in candidates:
        by_sp = None
        codon = f"{og_dir}/codon.afa"
        if os.path.exists(codon):
            by_sp = _locus_by_species(read_fasta(codon), species, members)
            if by_sp is not None:
                reused += 1
        if by_sp is None:
            fallback = realign(members) if members else None
            by_sp = _locus_by_species(fallback, species, members) if fallback else None
            if by_sp is None:
                dropped += 1
                continue
            realigned += 1
        try:
            filtered, n_codons = drop_all_gap_codons(by_sp)
        except ValueError:
            dropped += 1
            continue
        if n_codons == 0:
            dropped += 1
            continue
        loci.append((name, filtered))

    if not loci:
        sys.exit("No usable markers — refusing to write an empty supermatrix")

    concat, partitions = concatenate(loci, species)
    n_sites = len(concat[species[0]])
    print(f"loci reused from saved codon.afa: {reused}")
    print(f"loci realigned (mafft+pal2nal):   {realigned}")
    print(f"loci dropped:                     {dropped}")
    print(f"final: n_loci={len(loci)}, n_taxa={len(species)}, n_sites={n_sites}")

    with open(f"{args.out}/supermatrix.fa", "w") as f:
        for sp in species:
            f.write(f">{sp}\n{concat[sp]}\n")
    with open(f"{args.out}/partitions.txt", "w") as f:
        f.write(format_partitions(partitions))
    print(f"wrote {args.out}/supermatrix.fa and partitions.txt")


if __name__ == "__main__":
    main()
