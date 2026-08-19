"""Species tree utilities: loading and pairwise distance computation."""

from pathlib import Path
from typing import Dict, List, Set, Tuple
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


# Heuristic bounds for substitution-rate branch lengths (issue #14).
# Below the minimum the tree likely carries coalescent units (e.g. ASTRAL);
# above the maximum the units are likely wrong (e.g. mya, gene counts).
MIN_PLAUSIBLE_MAX_DISTANCE = 0.05
MAX_PLAUSIBLE_MAX_DISTANCE = 10.0


def validate_species_tree(species_tree, expected_species: Set[str]) -> List[str]:
    """Validate a loaded species tree against the species present in the data.

    Checks (issue #14):
      - leaf-name set vs expected_species (both directions reported)
      - non-positive or missing branch lengths
      - max pairwise distance implausibly small (coalescent/ASTRAL-like
        units) or implausibly large (wrong units)

    Returns a list of human-readable problem strings; empty list = OK.
    """
    problems: List[str] = []

    leaves = list(species_tree.leaves())
    leaf_names = {leaf.name for leaf in leaves}

    missing = expected_species - leaf_names
    extra = leaf_names - expected_species
    if missing:
        problems.append(
            f"Species in data but missing from species tree: {sorted(missing)} "
            "— pruning is silently disabled for genes of these species"
        )
    if extra:
        problems.append(
            f"Species tree leaves absent from the data: {sorted(extra)}"
        )

    bad_branches = []
    for node in species_tree.traverse():
        if node.up is None:  # root has no branch above it
            continue
        if node.dist is None or node.dist <= 0:
            bad_branches.append(node.name or "<internal>")
    if bad_branches:
        problems.append(
            f"Non-positive or missing branch lengths on {len(bad_branches)} "
            f"node(s): {bad_branches} — distance ratios will be unreliable"
        )

    if len(leaves) >= 2:
        distances = compute_pairwise_distances(species_tree)
        max_dist = max(distances.values())
        if max_dist < MIN_PLAUSIBLE_MAX_DISTANCE:
            problems.append(
                f"Max pairwise species distance is {max_dist:.4g} "
                f"(< {MIN_PLAUSIBLE_MAX_DISTANCE}) — branch lengths look like "
                "coalescent units (e.g. ASTRAL). Use substitution-rate branch "
                "lengths (e.g. IQ-TREE concatenation)"
            )
        elif max_dist > MAX_PLAUSIBLE_MAX_DISTANCE:
            problems.append(
                f"Max pairwise species distance is {max_dist:.4g} "
                f"(> {MAX_PLAUSIBLE_MAX_DISTANCE}) — branch lengths are likely "
                "in the wrong units for substitution-rate distance ratios"
            )

    return problems
