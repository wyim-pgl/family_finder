"""Gene tree pruning using TreeShrink + species-aware distance filtering.

Two-stage pruning:
  1. TreeShrink: statistical outlier detection based on branch length distribution
  2. Species-aware ratio: compare observed vs expected distances from species tree
"""

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from statistics import median
from typing import Dict, Optional, Set, Tuple

from ete4 import Tree

from config import Config
from utils.species import get_species

logger = logging.getLogger("family_finder")


def prune_orthogroup(
    tree_path: str,
    gene_to_species: Dict[str, str],
    expected_distances: Dict[Tuple[str, str], float],
    config: Config,
) -> Tuple[Set[str], Set[str]]:
    """Prune outlier genes using TreeShrink then species-aware ratio filtering.

    Stage 1 (TreeShrink): Remove statistically long branches.
    Stage 2 (Species-aware): Remove genes whose observed/expected distance
             ratio exceeds the threshold.

    Returns:
        (confirmed_gene_ids, outlier_gene_ids)
    """
    tree = Tree(Path(tree_path).read_text().strip())
    leaves = list(tree.leaves())
    all_genes = {leaf.name for leaf in leaves}

    # Too few leaves to prune
    if len(leaves) < config.min_species_for_pruning:
        return (all_genes, set())

    # Stage 1: TreeShrink
    ts_outliers = _run_treeshrink(tree_path, config)
    stage1_confirmed = all_genes - ts_outliers

    if ts_outliers:
        logger.debug(f"  TreeShrink removed {len(ts_outliers)}: {ts_outliers}")

    # Stage 2: Species-aware ratio (on TreeShrink survivors)
    species_in_tree = {gene_to_species.get(g, "") for g in stage1_confirmed}
    species_in_tree.discard("")

    if len(species_in_tree) >= config.min_species_for_pruning and len(stage1_confirmed) >= config.min_species_for_pruning:
        ratio_outliers = _species_aware_filter(
            tree_path, stage1_confirmed, gene_to_species, expected_distances, config
        )
    else:
        ratio_outliers = set()

    # Combine outliers from both stages
    total_outliers = ts_outliers | ratio_outliers
    confirmed = all_genes - total_outliers

    logger.info(
        f"Pruning result: {len(confirmed)} confirmed, {len(total_outliers)} outliers "
        f"(treeshrink={len(ts_outliers)}, ratio={len(ratio_outliers)})"
    )

    return (confirmed, total_outliers)


def _run_treeshrink(tree_path: str, config: Config) -> Set[str]:
    """Run TreeShrink on a single gene tree.

    Returns set of outlier gene IDs identified by TreeShrink.
    """
    tmpdir = tempfile.mkdtemp(prefix="treeshrink_")
    try:
        # TreeShrink expects: indir/gene1/input.tree
        gene_dir = Path(tmpdir) / "gene1"
        gene_dir.mkdir()
        shutil.copy(tree_path, gene_dir / "input.tree")

        outdir = Path(tmpdir) / "output"

        cmd = [
            "run_treeshrink.py",
            "-i", tmpdir,
            "-o", str(outdir),
            "-q", str(config.treeshrink_quantile),
            "-c",  # centroid reroot
        ]

        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300
        )

        if result.returncode != 0:
            logger.debug(f"TreeShrink warning: {result.stderr[:200]}")
            return set()

        # Parse output: output/gene1/output.txt lists removed taxa
        removed_file = outdir / "gene1" / "output.txt"
        if not removed_file.exists():
            return set()

        outliers = set()
        with open(removed_file) as f:
            for line in f:
                # TreeShrink outputs tab-separated gene IDs per line
                for gene_id in line.strip().split("\t"):
                    gene_id = gene_id.strip()
                    if gene_id:
                        outliers.add(gene_id)
        return outliers

    except subprocess.TimeoutExpired:
        logger.debug("TreeShrink timed out, skipping")
        return set()
    except Exception as e:
        logger.debug(f"TreeShrink failed: {e}")
        return set()
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _species_aware_filter(
    tree_path: str,
    gene_ids: Set[str],
    gene_to_species: Dict[str, str],
    expected_distances: Dict[Tuple[str, str], float],
    config: Config,
) -> Set[str]:
    """Species-aware distance ratio filtering on a subset of genes.

    Returns set of outlier gene IDs.
    """
    tree = Tree(Path(tree_path).read_text().strip())
    leaves = [l for l in tree.leaves() if l.name in gene_ids]

    outlier_scores = {}
    for leaf_i in leaves:
        sp_i = gene_to_species.get(leaf_i.name)
        if sp_i is None:
            outlier_scores[leaf_i.name] = 0.0
            continue

        ratios = []
        for leaf_j in leaves:
            if leaf_i is leaf_j:
                continue
            sp_j = gene_to_species.get(leaf_j.name)
            if sp_j is None or sp_i == sp_j:
                continue

            expected = expected_distances.get((sp_i, sp_j))
            if expected is None or expected <= 0:
                continue

            observed = tree.get_distance(leaf_i, leaf_j)
            ratios.append(observed / expected)

        if ratios:
            outlier_scores[leaf_i.name] = median(ratios)
        else:
            outlier_scores[leaf_i.name] = 0.0

    outliers = set()
    for gene_id, score in outlier_scores.items():
        if score > config.distance_ratio_threshold:
            outliers.add(gene_id)
            logger.debug(
                f"  Ratio outlier: {gene_id} (species={gene_to_species.get(gene_id)}, "
                f"score={score:.2f})"
            )

    return outliers
