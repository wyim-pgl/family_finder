#!/usr/bin/env python3
"""Standalone pseudogene detection for completed family_finder runs.

Analyzes an existing pipeline output directory to identify pseudogene
candidates without re-running the full pipeline.

Usage:
    python find_pseudogenes.py \
        --protein-dir data/pep --cds-dir data/cds \
        --outdir output_5sp --species Ococ

    # Or analyze all species:
    python find_pseudogenes.py \
        --protein-dir data/pep --cds-dir data/cds \
        --outdir output_5sp
"""

import argparse
import logging
import os
import sys
from pathlib import Path

# Prevent MAFFT_BINARIES conflict
os.environ.pop("MAFFT_BINARIES", None)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect pseudogene candidates from a completed family_finder run.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--protein-dir", required=True,
        help="Directory of per-species protein FASTA files",
    )
    parser.add_argument(
        "--cds-dir", required=True,
        help="Directory of per-species CDS FASTA files",
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Pipeline output directory (must contain summary.tsv)",
    )
    parser.add_argument(
        "--species", type=str, default=None,
        help="Restrict analysis to one species prefix (e.g., Ococ)",
    )
    parser.add_argument(
        "--truncation-threshold", type=float, default=0.5,
        help="Flag genes shorter than this fraction of family median (default: 0.5)",
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output directory for pseudogene results (default: <outdir>/pseudogene_analysis)",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable debug logging",
    )
    return parser.parse_args()


def load_families_from_summary(summary_path: Path):
    """Load confirmed families from summary.tsv."""
    families = {}
    with open(summary_path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                fam_id = parts[0]
                genes = set(parts[4].split(","))
                families[fam_id] = genes
    return families


def main():
    args = parse_args()

    # Setup logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger("family_finder")

    outdir = Path(args.outdir)
    summary_path = outdir / "summary.tsv"

    if not summary_path.exists():
        logger.error(f"summary.tsv not found in {outdir}. Run the pipeline first.")
        sys.exit(1)

    # Load data
    from config import Config
    from utils.seqio import build_seq_map

    config = Config()
    config.pseudogene_truncation_threshold = args.truncation_threshold

    logger.info(f"Loading protein sequences from {args.protein_dir}")
    protein_seqs = build_seq_map(args.protein_dir)
    logger.info(f"Loaded {len(protein_seqs)} protein sequences")

    logger.info(f"Loading CDS sequences from {args.cds_dir}")
    cds_seqs = build_seq_map(args.cds_dir)
    logger.info(f"Loaded {len(cds_seqs)} CDS sequences")

    if not protein_seqs:
        logger.error(f"No protein sequences found in {args.protein_dir}")
        sys.exit(1)
    if not cds_seqs:
        logger.error(f"No CDS sequences found in {args.cds_dir}")
        sys.exit(1)

    logger.info(f"Loading families from {summary_path}")
    families = load_families_from_summary(summary_path)
    logger.info(f"Loaded {len(families)} families")

    # Run pseudogene detection
    from steps.pseudogene import (
        detect_pseudogenes,
        write_pseudogene_report,
        write_pseudogene_summary,
        write_pseudogene_fasta,
        write_family_pseudogene_report,
        write_species_comparison,
        write_pseudogene_bed,
        write_chromosomal_distribution,
    )

    species_filter = args.species
    evidence = detect_pseudogenes(
        protein_seqs=protein_seqs,
        cds_seqs=cds_seqs,
        families=families,
        outdir=outdir,
        config=config,
        species_filter=species_filter,
    )

    # Write results
    result_dir = Path(args.output) if args.output else outdir / "pseudogene_analysis"
    result_dir.mkdir(parents=True, exist_ok=True)

    suffix = f"_{species_filter}" if species_filter else ""

    # Core reports
    write_pseudogene_report(
        evidence, result_dir / f"pseudogene_candidates{suffix}.tsv", species_filter
    )

    n_analyzed = len(protein_seqs)
    if species_filter:
        n_analyzed = sum(
            1 for gid in protein_seqs
            if gid.split(config.species_delimiter, 1)[0] == species_filter
        )
    write_pseudogene_summary(
        evidence, result_dir / f"pseudogene_summary{suffix}.txt", n_analyzed, species_filter
    )
    write_pseudogene_fasta(
        evidence, protein_seqs, cds_seqs, result_dir
    )

    # Extended analyses
    write_family_pseudogene_report(
        evidence, families, result_dir / f"family_pseudogene_enrichment{suffix}.tsv"
    )

    # Cross-species comparison only makes sense when analyzing all species
    if not species_filter:
        write_species_comparison(
            evidence, protein_seqs, families,
            result_dir / "species_comparison.tsv",
            config.species_delimiter,
        )

    write_pseudogene_bed(
        evidence, result_dir / f"pseudogene_candidates{suffix}.bed"
    )
    write_chromosomal_distribution(
        evidence, protein_seqs,
        result_dir / f"chromosomal_distribution{suffix}.tsv",
        species_filter, config.species_delimiter,
    )

    # Print quick summary to stdout
    from collections import Counter
    cls_counts = Counter(ev.classification for ev in evidence.values())
    total_pseudo = sum(v for k, v in cls_counts.items() if k.startswith("pseudogene"))

    print(f"\n{'='*60}")
    print(f"Pseudogene Detection Results{' ('+species_filter+')' if species_filter else ''}")
    print(f"{'='*60}")
    print(f"Genes analyzed:     {n_analyzed}")
    print(f"Candidates found:   {total_pseudo} ({100*total_pseudo/max(n_analyzed,1):.1f}%)")
    print(f"  High confidence:  {cls_counts.get('pseudogene_high', 0)}")
    print(f"  Medium confidence:{cls_counts.get('pseudogene_medium', 0)}")
    print(f"  Low confidence:   {cls_counts.get('pseudogene_low', 0)}")
    print(f"\nOutput files in {result_dir}/:")
    for f in sorted(result_dir.iterdir()):
        if f.is_file():
            size_kb = f.stat().st_size / 1024
            print(f"  {f.name:45s} ({size_kb:.1f} KB)")


if __name__ == "__main__":
    main()
