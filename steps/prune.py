"""Species-aware branch length pruning.

For each gene in a gene tree, computes an outlier score based on the ratio
of observed pairwise distances (from the gene tree) to expected distances
(from the species tree). Genes with a median ratio exceeding the threshold
are classified as outliers.
"""

import logging
from statistics import median
from typing import Dict, Set, Tuple

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
    """Prune outlier genes from a gene tree based on species-aware distance ratios.

    For each leaf gene_i, computes:
        outlier_score = median over all inter-species gene_j of:
            observed_distance(gene_i, gene_j) / expected_distance(species_i, species_j)

    Genes with outlier_score > config.distance_ratio_threshold are outliers.

    Args:
        tree_path: Path to Newick gene tree.
        gene_to_species: Mapping of gene_id -> species name.
        expected_distances: Pairwise species distances from species tree.
        config: Pipeline configuration.

    Returns:
        (confirmed_gene_ids, outlier_gene_ids)
    """
    tree = Tree(open(tree_path).read().strip())
    leaves = tree.leaves()

    # Too few leaves to prune meaningfully
    if len(leaves) < config.min_species_for_pruning:
        all_genes = {leaf.name for leaf in leaves}
        return (all_genes, set())

    # Check species diversity
    species_in_tree = {gene_to_species.get(leaf.name, "") for leaf in leaves}
    species_in_tree.discard("")
    if len(species_in_tree) < config.min_species_for_pruning:
        all_genes = {leaf.name for leaf in leaves}
        return (all_genes, set())

    # Compute outlier score for each leaf
    outlier_scores = {}
    for leaf_i in leaves:
        sp_i = gene_to_species.get(leaf_i.name)
        if sp_i is None:
            logger.warning(f"No species mapping for {leaf_i.name}, keeping it")
            outlier_scores[leaf_i.name] = 0.0
            continue

        ratios = []
        for leaf_j in leaves:
            if leaf_i is leaf_j:
                continue
            sp_j = gene_to_species.get(leaf_j.name)
            if sp_j is None or sp_i == sp_j:
                continue  # skip within-species (paralog) comparisons

            expected = expected_distances.get((sp_i, sp_j))
            if expected is None or expected <= 0:
                continue

            observed = leaf_i.get_distance(leaf_j)
            ratios.append(observed / expected)

        if ratios:
            outlier_scores[leaf_i.name] = median(ratios)
        else:
            # Can't evaluate (e.g., all other genes are same species)
            outlier_scores[leaf_i.name] = 0.0

    # Classify
    confirmed = set()
    outliers = set()
    for gene_id, score in outlier_scores.items():
        if score > config.distance_ratio_threshold:
            outliers.add(gene_id)
            logger.debug(
                f"  Outlier: {gene_id} (species={gene_to_species.get(gene_id)}, "
                f"score={score:.2f})"
            )
        else:
            confirmed.add(gene_id)

    logger.info(
        f"Pruning result: {len(confirmed)} confirmed, {len(outliers)} outliers "
        f"(threshold={config.distance_ratio_threshold})"
    )

    return (confirmed, outliers)
