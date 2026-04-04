"""Species tree utilities: loading and pairwise distance computation."""

from pathlib import Path
from typing import Dict, Tuple
from ete4 import Tree


def load_species_tree(path: str) -> Tree:
    """Load a species tree from a Newick file."""
    return Tree(Path(path).read_text().strip())


def compute_pairwise_distances(species_tree: Tree) -> Dict[Tuple[str, str], float]:
    """Compute pairwise distances between all leaf pairs in the species tree.

    Returns dict mapping (sp_a, sp_b) -> distance for all ordered pairs.
    """
    leaves = list(species_tree.leaves())
    distances = {}
    for i, leaf_a in enumerate(leaves):
        for leaf_b in leaves[i + 1 :]:
            dist = species_tree.get_distance(leaf_a, leaf_b)
            distances[(leaf_a.name, leaf_b.name)] = dist
            distances[(leaf_b.name, leaf_a.name)] = dist
    return distances


def get_species(gene_id: str, delimiter: str = "_") -> str:
    """Extract species prefix from a gene ID."""
    return gene_id.split(delimiter, 1)[0]
