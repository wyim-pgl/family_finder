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
        logger.debug(f"Unknown MAFFT strategy '{config.mafft_strategy}', using --auto")
        cmd = [config.mafft_bin, "--auto", str(input_fa)]

    logger.debug(f"Running MAFFT: {' '.join(cmd)}")

    with open(outpath, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logger.error(f"MAFFT failed for {input_fa}:\n{result.stderr}")
        raise RuntimeError(f"MAFFT failed with return code {result.returncode}")

    return outpath


def _filter_internal_stops(
    protein_aln: Path, cds_seqs: Dict[str, str]
) -> tuple:
    """Remove genes with internal stop codons from protein alignment and CDS.

    pal2nal fails when protein contains '*' (internal stop codons), which
    indicates CDS annotation errors (frameshifts, pseudogenes, etc.).

    Returns:
        (filtered_protein_aln_path or None, filtered_cds_seqs, removed_ids)
    """
    from utils.seqio import read_fasta

    prot_seqs = read_fasta(str(protein_aln))
    removed = set()
    for gid, seq in prot_seqs.items():
        seq_clean = seq.rstrip("*")
        if "*" in seq_clean:
            removed.add(gid)
            logger.debug(
                f"  Removing {gid}: internal stop codons in protein "
                f"({seq_clean.count('*')} stops)"
            )

    if not removed:
        return protein_aln, cds_seqs, removed

    # Filter protein alignment: write only clean sequences
    filtered_prot = {gid: seq for gid, seq in prot_seqs.items() if gid not in removed}
    filtered_cds = {gid: seq for gid, seq in cds_seqs.items() if gid not in removed}

    if len(filtered_prot) < 2:
        return None, filtered_cds, removed

    filtered_aln_path = protein_aln.parent / "proteins_nostop.afa"
    write_fasta(filtered_prot, str(filtered_aln_path))

    logger.debug(
        f"  Filtered {len(removed)} genes with internal stops, "
        f"{len(filtered_prot)} remaining"
    )
    return filtered_aln_path, filtered_cds, removed


def codon_align(
    protein_aln: Path, cds_seqs: Dict[str, str], outpath: Path, config: Config
) -> Path:
    """Generate codon alignment using pal2nal.

    Filters out genes with internal stop codons before running pal2nal,
    as these indicate CDS annotation errors that cause pal2nal to fail.

    Args:
        protein_aln: Path to protein alignment (FASTA).
        cds_seqs: Dict of gene_id -> CDS nucleotide sequence.
        outpath: Output path for codon alignment.
        config: Pipeline configuration.

    Returns path to codon-aligned FASTA file, or None on failure.
    """
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    # Filter genes with internal stop codons (annotation errors)
    filtered_aln, filtered_cds, removed = _filter_internal_stops(
        Path(protein_aln), cds_seqs
    )
    if removed:
        logger.warning(
            f"pal2nal: removed {len(removed)} genes with internal stop codons "
            f"from {protein_aln.name}: {removed}"
        )
    if filtered_aln is None:
        logger.warning(f"pal2nal: too few sequences after stop-codon filtering")
        return None

    # Write CDS sequences to temp file
    cds_fa = outpath.parent / "cds_unaligned.fa"
    write_fasta(filtered_cds, str(cds_fa))

    cmd = [
        config.pal2nal_bin,
        str(filtered_aln), str(cds_fa),
        "-output", "fasta",
    ]

    logger.debug(f"Running pal2nal: {' '.join(cmd)}")

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.warning(f"pal2nal failed:\n{result.stderr[:300]}")
        return None

    # Check output is non-empty
    if not result.stdout.strip():
        logger.warning(f"pal2nal produced empty output for {protein_aln}")
        return None

    with open(outpath, "w") as out_f:
        out_f.write(result.stdout)

    return outpath
