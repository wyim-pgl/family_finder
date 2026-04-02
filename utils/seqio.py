"""FASTA I/O and CDS-protein mapping utilities."""

from pathlib import Path
from typing import Dict
from Bio import SeqIO


def read_fasta(path: str) -> Dict[str, str]:
    """Read FASTA file, return dict of id -> sequence string."""
    seqs = {}
    for record in SeqIO.parse(path, "fasta"):
        seqs[record.id] = str(record.seq)
    return seqs


def read_fasta_dir(directory: str) -> Dict[str, str]:
    """Read all FASTA files in a directory, return combined dict of id -> sequence."""
    seqs = {}
    dirpath = Path(directory)
    for ext in ("*.fa", "*.fasta", "*.faa", "*.pep"):
        for fpath in dirpath.glob(ext):
            seqs.update(read_fasta(str(fpath)))
    return seqs


def write_fasta(seqs: Dict[str, str], path: str, wrap: int = 80):
    """Write sequences dict to FASTA file."""
    outpath = Path(path)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w") as f:
        for seq_id, seq in seqs.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i : i + wrap] + "\n")


def build_seq_map(directory: str) -> Dict[str, str]:
    """Build a flat mapping of gene_id -> sequence from all FASTA files in directory."""
    return read_fasta_dir(directory)


def split_by_species(
    seqs: Dict[str, str], outdir: str, delimiter: str = "_"
) -> Path:
    """Split a pool of sequences into per-species FASTA files.

    Gene IDs are expected as SpeciesPrefix{delimiter}GeneID.
    Returns path to the output directory containing per-species files.
    """
    outpath = Path(outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    species_seqs: Dict[str, Dict[str, str]] = {}
    for gene_id, seq in seqs.items():
        species = gene_id.split(delimiter, 1)[0]
        species_seqs.setdefault(species, {})[gene_id] = seq

    for species, sseqs in species_seqs.items():
        write_fasta(sseqs, str(outpath / f"{species}.fa"))

    return outpath
