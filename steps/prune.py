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
from typing import Callable, Dict, List, Optional, Set, Tuple

from ete4 import Tree

from config import Config
from utils.species import get_species

logger = logging.getLogger("family_finder")

# Minimum leaves for pruning to be meaningful (independent of species count).
MIN_LEAVES_FOR_PRUNING = 3
# Floor for the internal path length after terminal-branch subtraction.
PRUNE_EPS = 1e-9


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

    # Too few leaves or too few distinct species to prune. The old code
    # compared the LEAF count against min_species_for_pruning (issue #14).
    n_species = len({gene_to_species.get(leaf.name) for leaf in leaves} - {None})
    if len(leaves) < MIN_LEAVES_FOR_PRUNING or n_species < config.min_species_for_pruning:
        return (all_genes, set())

    # Stage 1: TreeShrink
    ts_outliers = _run_treeshrink(tree_path, config)
    stage1_confirmed = all_genes - ts_outliers

    if ts_outliers:
        logger.debug(f"  TreeShrink removed {len(ts_outliers)}: {ts_outliers}")

    # Stage 2: Species-aware ratio (on TreeShrink survivors)
    species_in_tree = {gene_to_species.get(g, "") for g in stage1_confirmed}
    species_in_tree.discard("")

    if len(species_in_tree) >= config.min_species_for_pruning and len(stage1_confirmed) >= MIN_LEAVES_FOR_PRUNING:
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
            logger.warning(
                f"TreeShrink failed (rc={result.returncode}), skipping: "
                f"{result.stderr[:200]}"
            )
            return set()

        # Parse output: output/gene1/output.txt lists removed taxa
        removed_file = outdir / "gene1" / "output.txt"
        if not removed_file.exists():
            logger.warning(
                f"TreeShrink produced no output file ({removed_file}), skipping"
            )
            return set()

        outliers = set()
        with open(removed_file) as f:
            for line in f:
                # TreeShrink outputs tab-separated gene IDs per line
                for gene_id in line.strip().split("\t"):
                    gene_id = gene_id.strip()
                    if gene_id:
                        outliers.add(gene_id)
        logger.info(
            f"TreeShrink ran on {Path(tree_path).name}: "
            f"{len(outliers)} outlier(s) flagged"
        )
        return outliers

    except subprocess.TimeoutExpired:
        logger.warning("TreeShrink timed out, skipping")
        return set()
    except Exception as e:
        logger.warning(f"TreeShrink unavailable or failed, skipping: {e}")
        return set()
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def compute_pair_ratios(
    dist_fn: Callable[[str, str], float],
    leaves_meta: Dict[str, Tuple[Optional[str], float]],
    expected: Dict[Tuple[str, str], float],
    eps: float = PRUNE_EPS,
) -> Tuple[Dict[str, List[float]], List[float]]:
    """Internal-path ratios r(i,j) for all cross-species gene pairs (pure).

    r(i,j) = max(d_obs(i,j) - b_i - b_j, eps) / d_exp(sp_i, sp_j)

    Subtracting both TERMINAL branch lengths removes the focal gene's own
    evolutionary rate from its score: a fast-evolving but correctly placed
    gene no longer looks misplaced (issue #14).

    Args:
        dist_fn: (gene_i, gene_j) -> observed tree distance.
        leaves_meta: gene_id -> (species or None, terminal branch length).
        expected: (sp_i, sp_j) -> expected species-tree distance (both orders).

    Returns:
        (ratios_by_gene, all_ratios). Every gene in leaves_meta has an entry
        in ratios_by_gene (possibly empty). Pairs with unknown species,
        same-species, or missing/non-positive expected distance are skipped —
        the same skips as the legacy criterion. all_ratios holds each
        unordered pair once.
    """
    genes = list(leaves_meta)
    ratios_by_gene: Dict[str, List[float]] = {g: [] for g in genes}
    all_ratios: List[float] = []

    for idx, gene_i in enumerate(genes):
        sp_i, b_i = leaves_meta[gene_i]
        if sp_i is None:
            continue
        for gene_j in genes[idx + 1:]:
            sp_j, b_j = leaves_meta[gene_j]
            if sp_j is None or sp_i == sp_j:
                continue
            d_exp = expected.get((sp_i, sp_j))
            if d_exp is None or d_exp <= 0:
                continue
            d_obs = dist_fn(gene_i, gene_j)
            r = max(d_obs - (b_i or 0.0) - (b_j or 0.0), eps) / d_exp
            ratios_by_gene[gene_i].append(r)
            ratios_by_gene[gene_j].append(r)
            all_ratios.append(r)

    return ratios_by_gene, all_ratios


def relative_scores(
    ratios_by_gene: Dict[str, List[float]],
    all_ratios: List[float],
) -> Dict[str, Tuple[float, float]]:
    """Per-gene (S_raw, S_norm) from pairwise internal-path ratios (pure).

    S_raw(i)  = median_j r(i,j) over the gene's cross-species pairs.
    M_F       = median over ALL computed r(i,j) pairs in the family.
    S_norm(i) = S_raw(i) / M_F.

    Dividing by the family median cancels the species-tree/codon-tree unit
    mismatch by construction. Genes with no ratios, or families where
    M_F <= 0 / no pairs, get (0.0, 0.0) = not prunable.
    """
    family_median = median(all_ratios) if all_ratios else 0.0
    scores: Dict[str, Tuple[float, float]] = {}
    for gene, ratios in ratios_by_gene.items():
        if not ratios or family_median <= 0:
            scores[gene] = (0.0, 0.0)
        else:
            s_raw = median(ratios)
            scores[gene] = (s_raw, s_raw / family_median)
    return scores


def decide_relative_outliers(
    scores: Dict[str, Tuple[float, float]],
    relative_threshold: float,
    score_floor: float,
) -> Set[str]:
    """Prune iff S_norm > relative_threshold AND S_raw > score_floor (pure).

    Both must hold: the raw-score floor prevents pruning the worst gene of a
    perfectly fine (tight) family just because it is relatively worst.
    """
    return {
        gene
        for gene, (s_raw, s_norm) in scores.items()
        if s_norm > relative_threshold and s_raw > score_floor
    }


def compute_absolute_scores(
    dist_fn: Callable[[str, str], float],
    leaves_meta: Dict[str, Tuple[Optional[str], float]],
    expected: Dict[Tuple[str, str], float],
) -> Dict[str, float]:
    """Legacy scoring (pure): S(i) = median_j d_obs(i,j) / d_exp(sp_i, sp_j).

    Exact pre-#14 behavior: no terminal-branch subtraction, no family
    normalization; unknown species or no usable pairs -> score 0.0.
    """
    genes = list(leaves_meta)
    scores: Dict[str, float] = {}
    for gene_i in genes:
        sp_i, _ = leaves_meta[gene_i]
        if sp_i is None:
            scores[gene_i] = 0.0
            continue
        ratios = []
        for gene_j in genes:
            if gene_j == gene_i:
                continue
            sp_j, _ = leaves_meta[gene_j]
            if sp_j is None or sp_i == sp_j:
                continue
            d_exp = expected.get((sp_i, sp_j))
            if d_exp is None or d_exp <= 0:
                continue
            ratios.append(dist_fn(gene_i, gene_j) / d_exp)
        scores[gene_i] = median(ratios) if ratios else 0.0
    return scores


def _species_aware_filter(
    tree_path: str,
    gene_ids: Set[str],
    gene_to_species: Dict[str, str],
    expected_distances: Dict[Tuple[str, str], float],
    config: Config,
) -> Set[str]:
    """Species-aware distance ratio filtering on a subset of genes.

    config.prune_criterion selects the scoring:
      "relative" (default) — internal-path ratio with per-family
        normalization (compute_pair_ratios / relative_scores).
      "absolute" — legacy median(d_obs/d_exp) vs distance_ratio_threshold.

    Returns set of outlier gene IDs.
    """
    tree = Tree(Path(tree_path).read_text().strip())
    # Prune tree to surviving genes so distances aren't distorted by removed taxa
    surviving = [l for l in tree.leaves() if l.name in gene_ids]
    if len(surviving) < 2:
        return set()
    tree.prune(surviving, preserve_branch_length=True)
    leaves = list(tree.leaves())

    name_to_leaf = {leaf.name: leaf for leaf in leaves}
    leaves_meta = {
        leaf.name: (gene_to_species.get(leaf.name), leaf.dist or 0.0)
        for leaf in leaves
    }

    def dist_fn(name_a: str, name_b: str) -> float:
        return tree.get_distance(name_to_leaf[name_a], name_to_leaf[name_b])

    if config.prune_criterion == "relative":
        ratios_by_gene, all_ratios = compute_pair_ratios(
            dist_fn, leaves_meta, expected_distances
        )
        scores = relative_scores(ratios_by_gene, all_ratios)
        outliers = decide_relative_outliers(
            scores, config.prune_relative_threshold, config.prune_score_floor
        )
        for gene_id in outliers:
            s_raw, s_norm = scores[gene_id]
            logger.debug(
                f"  Ratio outlier: {gene_id} (species={gene_to_species.get(gene_id)}, "
                f"S_raw={s_raw:.2f}, S_norm={s_norm:.2f})"
            )
        return outliers

    # Legacy "absolute" criterion (pre-#14 behavior, kept for reproducibility)
    outlier_scores = compute_absolute_scores(dist_fn, leaves_meta, expected_distances)

    outliers = set()
    for gene_id, score in outlier_scores.items():
        if score > config.distance_ratio_threshold:
            outliers.add(gene_id)
            logger.debug(
                f"  Ratio outlier: {gene_id} (species={gene_to_species.get(gene_id)}, "
                f"score={score:.2f})"
            )

    return outliers
