"""OrthoFinder wrapper and output parser."""

import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from config import Config

logger = logging.getLogger("family_finder")


def run_orthofinder(input_dir: Path, output_dir: Path, config: Config) -> Path:
    """Run OrthoFinder on per-species FASTA files.

    Args:
        input_dir: Directory containing per-species protein FASTA files.
        output_dir: Directory for OrthoFinder results.
        config: Pipeline configuration.

    Returns:
        Path to the OrthoFinder Results directory.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        config.orthofinder_bin,
        "-f", str(input_dir),
        "-t", str(config.orthofinder_threads),
        "-o", str(output_dir),
    ]
    if config.orthofinder_extra_args:
        cmd.extend(config.orthofinder_extra_args.split())

    logger.info(f"Running OrthoFinder: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"OrthoFinder failed:\n{result.stderr}")
        raise RuntimeError(f"OrthoFinder failed with return code {result.returncode}")

    # Find the Results_* directory
    results_dirs = sorted(output_dir.glob("Results_*"))
    if not results_dirs:
        # OrthoFinder sometimes nests under OrthoFinder/Results_*
        results_dirs = sorted(output_dir.glob("OrthoFinder/Results_*"))
    if not results_dirs:
        raise FileNotFoundError(f"No Results_* directory found in {output_dir}")

    results_dir = results_dirs[-1]  # latest
    logger.info(f"OrthoFinder results: {results_dir}")
    return results_dir


def parse_orthogroups(results_dir: Path) -> Dict[str, List[str]]:
    """Parse Orthogroups.tsv from OrthoFinder results.

    Returns:
        Dict mapping orthogroup ID -> list of gene IDs.
    """
    og_file = results_dir / "Orthogroups" / "Orthogroups.tsv"
    if not og_file.exists():
        raise FileNotFoundError(f"Orthogroups.tsv not found at {og_file}")

    orthogroups = {}
    with open(og_file) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            og_id = parts[0]
            genes = []
            for cell in parts[1:]:
                if cell.strip():
                    genes.extend(g.strip() for g in cell.split(",") if g.strip())
            if genes:
                orthogroups[og_id] = genes

    logger.info(f"Parsed {len(orthogroups)} orthogroups")
    return orthogroups


def get_orthofinder_species_tree(results_dir: Path) -> Optional[Path]:
    """Get the species tree inferred by OrthoFinder, if available."""
    tree_path = results_dir / "Species_Tree" / "SpeciesTree_rooted.txt"
    if tree_path.exists():
        return tree_path
    return None
