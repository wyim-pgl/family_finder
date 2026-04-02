"""Alignment wrappers: MAFFT for protein, pal2nal for codon alignment."""

import logging
import subprocess
from pathlib import Path
from typing import Dict

from config import Config
from utils.seqio import write_fasta

logger = logging.getLogger("family_finder")


def align_protein(seqs: Dict[str, str], outpath: Path, config: Config) -> Path:
    """Align protein sequences with MAFFT.

    Returns path to the aligned FASTA file.
    """
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    # Write input
    input_fa = outpath.parent / "input_proteins.fa"
    write_fasta(seqs, str(input_fa))

    # Build MAFFT command
    if config.mafft_strategy == "auto":
        cmd = [config.mafft_bin, "--auto", str(input_fa)]
    elif config.mafft_strategy == "linsi":
        cmd = [config.mafft_bin, "--localpair", "--maxiterate", "1000", str(input_fa)]
    elif config.mafft_strategy == "ginsi":
        cmd = [config.mafft_bin, "--globalpair", "--maxiterate", "1000", str(input_fa)]
    else:
        cmd = [config.mafft_bin, "--auto", str(input_fa)]

    logger.debug(f"Running MAFFT: {' '.join(cmd)}")

    with open(outpath, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logger.error(f"MAFFT failed for {input_fa}:\n{result.stderr}")
        raise RuntimeError(f"MAFFT failed with return code {result.returncode}")

    return outpath


def codon_align(
    protein_aln: Path, cds_seqs: Dict[str, str], outpath: Path, config: Config
) -> Path:
    """Generate codon alignment using pal2nal.

    Args:
        protein_aln: Path to protein alignment (FASTA).
        cds_seqs: Dict of gene_id -> CDS nucleotide sequence.
        outpath: Output path for codon alignment.
        config: Pipeline configuration.

    Returns path to codon-aligned FASTA file.
    """
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    # Write CDS sequences to temp file
    cds_fa = outpath.parent / "cds_unaligned.fa"
    write_fasta(cds_seqs, str(cds_fa))

    cmd = [
        config.pal2nal_bin,
        str(protein_aln), str(cds_fa),
        "-output", "fasta",
    ]

    logger.debug(f"Running pal2nal: {' '.join(cmd)}")

    with open(outpath, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logger.error(f"pal2nal failed:\n{result.stderr}")
        raise RuntimeError(f"pal2nal failed with return code {result.returncode}")

    return outpath
