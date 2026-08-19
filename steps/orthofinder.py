"""OrthoFinder wrapper and output parser."""

import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from config import Config

logger = logging.getLogger("family_finder")


def build_orthofinder_cmd(input_dir: Path, output_dir: Path, config: Config) -> List[str]:
    """Assemble the OrthoFinder argv. Pure function — unit-testable without
    OrthoFinder installed.

    Two modes:
    - fresh run: -f <input_dir> -o <output_dir>
    - blast reuse (config.orthofinder_reuse_blast_dir set): -b <WorkingDirectory>;
      -f/-o are omitted because OrthoFinder rejects -o with -b and writes
      results next to the reused WorkingDirectory instead.
    """
    if config.orthofinder_reuse_blast_dir:
        cmd = [config.orthofinder_bin, "-b", str(config.orthofinder_reuse_blast_dir)]
    else:
        cmd = [
            config.orthofinder_bin,
            "-f", str(input_dir),
            "-o", str(output_dir),
        ]
    cmd += ["-t", str(config.orthofinder_threads)]
    cmd += ["-I", str(config.orthofinder_inflation)]
    if config.orthofinder_search_program:
        cmd += ["-S", config.orthofinder_search_program]
    if config.orthofinder_stop_after_orthogroups:
        # Accepted by v3.1.3 (run/process_args.py) though hidden from -h.
        # Skips gene trees / species tree / ortholog inference AND
        # Comparative_Genomics_Statistics/ — Orthogroups*.tsv still written.
        cmd.append("-og")
    if config.orthofinder_extra_args:
        extra = config.orthofinder_extra_args.split()
        blocked = {"-o", "--output", "-f", "--fasta", "-b", "-I", "-S", "-og", "-t"}
        for arg in extra:
            if arg in blocked:
                raise ValueError(f"orthofinder_extra_args contains blocked flag: {arg}")
        cmd.extend(extra)
    return cmd


def run_orthofinder(input_dir: Path, output_dir: Path, config: Config) -> Path:
    """Run OrthoFinder on per-species FASTA files.

    Args:
        input_dir: Directory containing per-species protein FASTA files.
        output_dir: Directory for OrthoFinder results.
        config: Pipeline configuration.

    Returns:
        Path to the OrthoFinder Results directory.
    """
    output_dir = Path(output_dir)
    if config.orthofinder_reuse_blast_dir:
        # Blast-reuse mode: never delete anything — the WorkingDirectory being
        # reused may live under a previous round's output tree.
        pass
    elif output_dir.exists():
        # OrthoFinder requires -o dir to NOT already exist
        import shutil
        shutil.rmtree(output_dir)

    cmd = build_orthofinder_cmd(input_dir, output_dir, config)

    logger.info(f"Running OrthoFinder: {' '.join(cmd)}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    if result.returncode != 0:
        logger.error(f"OrthoFinder failed:\n{result.stdout}")
        raise RuntimeError(f"OrthoFinder failed with return code {result.returncode}")

    logger.info(f"OrthoFinder stdout (last 5 lines):\n" +
                "\n".join(result.stdout.strip().splitlines()[-5:]))

    # Find the Results_* directory
    if config.orthofinder_reuse_blast_dir:
        # With -b, OrthoFinder writes results next to the reused
        # WorkingDirectory, not under our -o (which is not passed).
        blast_dir = Path(config.orthofinder_reuse_blast_dir)
        search_roots = [blast_dir, blast_dir / "OrthoFinder", blast_dir.parent]
    else:
        search_roots = [output_dir, output_dir / "OrthoFinder"]
    results_dirs = []
    for root in search_roots:
        results_dirs.extend(root.glob("Results_*"))
    if not results_dirs:
        raise FileNotFoundError(
            f"No Results_* directory found under {[str(r) for r in search_roots]}"
        )

    # Latest by modification time (name sort breaks across month boundaries
    # and when -b appends to a directory holding older Results_*)
    results_dir = max(results_dirs, key=lambda p: p.stat().st_mtime)
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
