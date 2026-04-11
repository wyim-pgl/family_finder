#!/usr/bin/env python3
"""family_finder: Iterative gene family construction pipeline.

Repeatedly runs OrthoFinder → alignment → tree building → species-aware pruning.
Outliers from each round are re-clustered in subsequent rounds, allowing
displaced sequences to find their true gene families.
"""

import argparse
import os
import sys
from pathlib import Path

# Prevent MAFFT_BINARIES conflict when running inside conda/micromamba env
os.environ.pop("MAFFT_BINARIES", None)

# Ensure orthofinder conda env bin is in PATH (for mcl, diamond, etc.)
_of_env_bin = "/data/gpfs/assoc/pgl/bin/conda/conda_envs/orthofinder/bin"
if _of_env_bin not in os.environ.get("PATH", ""):
    os.environ["PATH"] = _of_env_bin + ":" + os.environ.get("PATH", "")

from config import Config
from utils.logging_setup import setup_logging
import pipeline


def parse_args():
    parser = argparse.ArgumentParser(
        description="Iterative gene family construction with species-aware pruning.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required arguments
    parser.add_argument(
        "--protein-dir", required=True,
        help="Directory of per-species protein FASTA files",
    )
    parser.add_argument(
        "--cds-dir", required=True,
        help="Directory of per-species CDS FASTA files",
    )
    parser.add_argument(
        "--species-tree", required=True,
        help="Newick species tree file (e.g., from ASTRAL)",
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Output directory",
    )

    # Optional arguments
    parser.add_argument(
        "--config", default=None,
        help="JSON configuration file (overrides defaults)",
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume from the latest checkpoint",
    )
    parser.add_argument(
        "--max-rounds", type=int, default=None,
        help="Maximum number of iterative rounds (default: 10)",
    )
    parser.add_argument(
        "--threshold", type=float, default=None,
        help="Distance ratio threshold for pruning (default: 5.0)",
    )
    parser.add_argument(
        "--threads", type=int, default=None,
        help="Number of parallel workers (default: 8)",
    )
    parser.add_argument(
        "--tree-builder", choices=["fasttree", "iqtree"], default=None,
        help="Tree building program (default: fasttree)",
    )
    parser.add_argument(
        "--run-codeml", action="store_true",
        help="Run codeml on confirmed families after pipeline completes",
    )
    parser.add_argument(
        "--no-hmmer-rescue", action="store_true",
        help="Disable HMMER rescue step for unplaced genes",
    )
    parser.add_argument(
        "--hmmer-evalue", type=float, default=None,
        help="E-value threshold for HMMER rescue (default: 1e-5)",
    )
    parser.add_argument(
        "--no-pseudogene-detection", action="store_true",
        help="Disable post-convergence pseudogene detection",
    )
    parser.add_argument(
        "--pseudogene-species", type=str, default=None,
        help="Restrict pseudogene analysis to one species (e.g., Ococ)",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable debug logging",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Load config
    if args.config:
        config = Config.from_json(args.config)
    else:
        config = Config()

    # CLI overrides
    if args.max_rounds is not None:
        config.max_rounds = args.max_rounds
    if args.threshold is not None:
        config.distance_ratio_threshold = args.threshold
    if args.threads is not None:
        config.n_workers = args.threads
        config.orthofinder_threads = args.threads
    if args.tree_builder is not None:
        config.tree_builder = args.tree_builder
    if args.run_codeml:
        config.run_codeml = True
    if args.no_hmmer_rescue:
        config.hmmer_rescue = False
    if args.hmmer_evalue is not None:
        config.hmmer_evalue = args.hmmer_evalue
    if args.no_pseudogene_detection:
        config.pseudogene_detection = False
    if args.pseudogene_species is not None:
        config.pseudogene_species_filter = args.pseudogene_species

    # Setup logging
    import logging
    level = logging.DEBUG if args.verbose else logging.INFO
    setup_logging(Path(args.outdir), level)

    logger = logging.getLogger("family_finder")
    logger.info("family_finder pipeline starting")
    logger.info(f"Protein dir: {args.protein_dir}")
    logger.info(f"CDS dir: {args.cds_dir}")
    logger.info(f"Species tree: {args.species_tree}")
    logger.info(f"Output dir: {args.outdir}")
    logger.info(f"Max rounds: {config.max_rounds}")
    logger.info(f"Pruning threshold: {config.distance_ratio_threshold}")
    logger.info(f"Workers: {config.n_workers}")
    logger.info(f"Tree builder: {config.tree_builder}")

    # Run pipeline
    pipeline.run(
        protein_dir=args.protein_dir,
        cds_dir=args.cds_dir,
        species_tree_path=args.species_tree,
        outdir=args.outdir,
        config=config,
        resume=args.resume,
    )

    # Run codeml if requested
    if config.run_codeml:
        logger.info("Running codeml on confirmed families...")
        from steps.codeml import run_codeml_on_families
        # Reload families from summary
        families = {}
        summary = Path(args.outdir) / "summary.tsv"
        if summary.exists():
            with open(summary) as f:
                f.readline()  # skip header
                for line in f:
                    parts = line.strip().split("\t")
                    family_id = parts[0]
                    gene_list = set(parts[4].split(","))
                    families[family_id] = gene_list

            from utils.seqio import build_seq_map
            cds_pool = build_seq_map(args.cds_dir)
            run_codeml_on_families(
                families, cds_pool,
                Path(args.outdir) / "codeml", config,
            )

    logger.info("Done.")


if __name__ == "__main__":
    main()
