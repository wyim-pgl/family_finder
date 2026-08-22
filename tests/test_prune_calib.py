"""Tests for the calibrated relative pruning criterion and species tree
validation (issue #14).

The pruning math lives in pure functions (steps.prune.compute_pair_ratios,
relative_scores, decide_relative_outliers, compute_absolute_scores) that
operate on plain dicts, so no ete4 is needed. ete4 is stubbed before import,
same pattern as tests/test_size_gate.py.
"""

import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# steps.prune and utils.species import ete4 at module level; stub it so
# these tests run without the conda environment.
if "ete4" not in sys.modules:
    ete4_stub = types.ModuleType("ete4")
    ete4_stub.Tree = object
    sys.modules["ete4"] = ete4_stub

from config import Config  # noqa: E402
from steps.prune import (  # noqa: E402
    compute_absolute_scores,
    compute_pair_ratios,
    decide_relative_outliers,
    relative_scores,
)
from utils.species import validate_species_tree  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_dist_fn(pairwise):
    """Symmetric lookup over a dict keyed by unordered gene-name pairs."""
    def dist_fn(a, b):
        return pairwise[(a, b)] if (a, b) in pairwise else pairwise[(b, a)]
    return dist_fn


class FakeNode:
    def __init__(self, name, dist, up=None):
        self.name = name
        self.dist = dist
        self.up = up


class FakeTree:
    """Duck-typed species tree: leaves(), traverse(), get_distance()."""

    def __init__(self, leaf_dists, pairwise, internal_dists=()):
        self.root = FakeNode(None, 0.0, up=None)
        self._internal = [FakeNode(None, d, up=self.root) for d in internal_dists]
        self._leaves = [FakeNode(n, d, up=self.root) for n, d in leaf_dists.items()]
        self._pairwise = pairwise

    def leaves(self):
        return list(self._leaves)

    def traverse(self):
        return [self.root] + self._internal + self._leaves

    def get_distance(self, a, b):
        key = (a.name, b.name)
        if key in self._pairwise:
            return self._pairwise[key]
        return self._pairwise[(b.name, a.name)]


# Family fixture: 3 species (S1..S3), expected distance 1.0 for every pair.
# g1/g2/g3 are well behaved (internal path 0.8, terminal branches 0.1).
# gfast (S1) is fast-evolving but CORRECTLY placed: terminal branch 2.0,
#   same 0.8 internal path to everyone.
# gbad (S3) is MISPLACED: normal terminal branch 0.1, inflated internal
#   path 4.0 to everyone.
EXPECTED = {
    (a, b): 1.0
    for a in ("S1", "S2", "S3")
    for b in ("S1", "S2", "S3")
    if a != b
}

LEAVES_META = {
    "S1_g1": ("S1", 0.1),
    "S2_g2": ("S2", 0.1),
    "S3_g3": ("S3", 0.1),
    "S1_gfast": ("S1", 2.0),
    "S3_gbad": ("S3", 0.1),
}

PAIRWISE = {
    # well-behaved cross-species pairs: 0.8 internal + both terminals
    ("S1_g1", "S2_g2"): 1.0,
    ("S1_g1", "S3_g3"): 1.0,
    ("S2_g2", "S3_g3"): 1.0,
    # gfast: correct placement, long terminal branch
    ("S1_gfast", "S2_g2"): 2.9,   # 0.8 + 2.0 + 0.1
    ("S1_gfast", "S3_g3"): 2.9,
    ("S1_gfast", "S3_gbad"): 4.0 + 2.0 + 0.1,  # bad internal path dominates
    # gbad: misplaced, inflated internal path 4.0
    ("S3_gbad", "S1_g1"): 4.2,    # 4.0 + 0.1 + 0.1
    ("S3_gbad", "S2_g2"): 4.2,
    # same-species pairs present but must be skipped
    ("S1_g1", "S1_gfast"): 99.0,
    ("S3_g3", "S3_gbad"): 99.0,
}


def _scores(pairwise=PAIRWISE, leaves_meta=LEAVES_META):
    ratios_by_gene, all_ratios = compute_pair_ratios(
        make_dist_fn(pairwise), leaves_meta, EXPECTED
    )
    return relative_scores(ratios_by_gene, all_ratios)


# ---------------------------------------------------------------------------
# Relative criterion: pure math
# ---------------------------------------------------------------------------

def test_terminal_branch_subtraction_rescues_fast_gene_flags_misplaced_gene():
    # Act
    scores = _scores()

    # Assert: fast-but-correct gene normalizes to ~1; misplaced gene is high
    s_raw_fast, s_norm_fast = scores["S1_gfast"]
    s_raw_bad, s_norm_bad = scores["S3_gbad"]
    assert s_norm_fast == pytest.approx(1.0)
    assert s_raw_fast == pytest.approx(0.8)
    assert s_norm_bad == pytest.approx(5.0)   # 4.0 / M_F(0.8)
    assert s_raw_bad == pytest.approx(4.0)

    # Decision with default thresholds: only the misplaced gene is pruned
    cfg = Config()
    outliers = decide_relative_outliers(
        scores, cfg.prune_relative_threshold, cfg.prune_score_floor
    )
    assert outliers == {"S3_gbad"}


def test_same_species_pairs_are_skipped():
    # Arrange / Act: the 99.0 same-species distances must never be used
    ratios_by_gene, _ = compute_pair_ratios(
        make_dist_fn(PAIRWISE), LEAVES_META, EXPECTED
    )
    # Assert: g1 has pairs to g2, g3, gbad only (gfast is same-species)
    assert len(ratios_by_gene["S1_g1"]) == 3


def test_per_family_normalization_cancels_global_scale():
    # Arrange: multiply every observed distance AND terminal branch by 5
    # (a uniform unit change of the gene tree)
    scaled_pairwise = {k: v * 5.0 for k, v in PAIRWISE.items()}
    scaled_meta = {g: (sp, b * 5.0) for g, (sp, b) in LEAVES_META.items()}

    # Act
    base = _scores()
    scaled = _scores(scaled_pairwise, scaled_meta)

    # Assert: S_norm identical for every gene; S_raw scales by 5
    for gene in LEAVES_META:
        assert scaled[gene][1] == pytest.approx(base[gene][1])
        assert scaled[gene][0] == pytest.approx(base[gene][0] * 5.0)


def test_score_floor_protects_tight_family():
    # Arrange: a tight family — relative worst gene has S_norm 5 but tiny S_raw
    scores = {
        "A_g1": (0.1, 1.0),
        "B_g2": (0.1, 1.0),
        "C_g3": (0.5, 5.0),  # S_norm > 3.0 but S_raw < 2.0
    }

    # Act
    outliers = decide_relative_outliers(scores, 3.0, 2.0)

    # Assert: floor prevents pruning the worst of a perfectly fine family
    assert outliers == set()


def test_no_pairs_or_zero_family_median_means_not_prunable():
    # Arrange: gene with no cross-species pairs at all
    scores = relative_scores({"A_g1": []}, [])
    # Assert
    assert scores["A_g1"] == (0.0, 0.0)


# ---------------------------------------------------------------------------
# Legacy absolute criterion: unchanged behavior
# ---------------------------------------------------------------------------

def test_absolute_scores_match_legacy_math():
    # Arrange: legacy S(i) = median_j(d_obs/d_exp), no terminal subtraction
    leaves_meta = {
        "S1_a": ("S1", 0.1),
        "S2_b": ("S2", 0.1),
        "S3_c": ("S3", 0.1),
        "Xx_unknown": (None, 0.1),  # unknown species -> score 0.0
    }
    pairwise = {
        ("S1_a", "S2_b"): 2.0,
        ("S1_a", "S3_c"): 4.0,
        ("S2_b", "S3_c"): 6.0,
    }

    # Act
    scores = compute_absolute_scores(make_dist_fn(pairwise), leaves_meta, EXPECTED)

    # Assert: exact legacy values (d_exp = 1.0 everywhere)
    assert scores["S1_a"] == pytest.approx(3.0)   # median(2.0, 4.0)
    assert scores["S2_b"] == pytest.approx(4.0)   # median(2.0, 6.0)
    assert scores["S3_c"] == pytest.approx(5.0)   # median(4.0, 6.0)
    assert scores["Xx_unknown"] == 0.0


def test_config_defaults():
    cfg = Config()
    assert cfg.prune_criterion == "relative"
    assert cfg.prune_relative_threshold == 3.0
    assert cfg.prune_score_floor == 2.0
    assert cfg.distance_ratio_threshold == 5.0  # legacy knob unchanged


# ---------------------------------------------------------------------------
# validate_species_tree
# ---------------------------------------------------------------------------

def _good_tree():
    return FakeTree(
        leaf_dists={"S1": 0.5, "S2": 0.5, "S3": 0.2},
        pairwise={
            ("S1", "S2"): 1.0,
            ("S1", "S3"): 1.0,
            ("S2", "S3"): 0.7,
        },
        internal_dists=(0.3,),
    )


def test_validate_passes_good_tree():
    # Act
    problems = validate_species_tree(_good_tree(), {"S1", "S2", "S3"})
    # Assert
    assert problems == []


def test_validate_catches_name_mismatch_both_directions():
    # Arrange: data has S4 (not in tree); tree has S3 (not in data)
    problems = validate_species_tree(_good_tree(), {"S1", "S2", "S4"})

    # Assert: both differences reported
    assert any("S4" in p and "missing from species tree" in p for p in problems)
    assert any("S3" in p and "absent from the data" in p for p in problems)


def test_validate_flags_astral_like_tiny_branch_lengths():
    # Arrange: coalescent-unit-scale tree (max pairwise distance < 0.05)
    tiny = FakeTree(
        leaf_dists={"S1": 0.005, "S2": 0.005, "S3": 0.002},
        pairwise={
            ("S1", "S2"): 0.01,
            ("S1", "S3"): 0.01,
            ("S2", "S3"): 0.007,
        },
    )

    # Act
    problems = validate_species_tree(tiny, {"S1", "S2", "S3"})

    # Assert
    assert any("coalescent" in p for p in problems)


def test_validate_flags_nonpositive_branch_lengths():
    # Arrange: one leaf with zero branch length
    tree = FakeTree(
        leaf_dists={"S1": 0.5, "S2": 0.0, "S3": 0.2},
        pairwise={
            ("S1", "S2"): 0.5,
            ("S1", "S3"): 0.7,
            ("S2", "S3"): 0.2,
        },
    )

    # Act
    problems = validate_species_tree(tree, {"S1", "S2", "S3"})

    # Assert
    assert any("Non-positive or missing branch lengths" in p for p in problems)


def test_validate_flags_a_topology_only_tree_whose_branches_are_all_equal():
    # Arrange: data_15sp/species_tree.nwk shipped with every branch = 1.0
    # (issue #41). The old bounds let it through — its max pairwise distance
    # lands at exactly 10.0 and the guard was `> 10.0`.
    dummy = FakeTree(
        leaf_dists={"S1": 1.0, "S2": 1.0, "S3": 1.0},
        pairwise={("S1", "S2"): 2.0, ("S1", "S3"): 4.0, ("S2", "S3"): 4.0},
        internal_dists=(1.0, 1.0),
    )

    # Act
    problems = validate_species_tree(dummy, {"S1", "S2", "S3"})

    # Assert
    assert any("identical" in p for p in problems)


def test_uniform_branch_length_check_is_independent_of_scale():
    # Arrange: all 0.1 — max pairwise 0.4 sits inside the plausible band and
    # would never trip the distance bounds, yet carries no rate information.
    dummy = FakeTree(
        leaf_dists={"S1": 0.1, "S2": 0.1, "S3": 0.1},
        pairwise={("S1", "S2"): 0.2, ("S1", "S3"): 0.4, ("S2", "S3"): 0.4},
        internal_dists=(0.1,),
    )

    # Act
    problems = validate_species_tree(dummy, {"S1", "S2", "S3"})

    # Assert
    assert any("identical" in p for p in problems)


def test_an_estimated_tree_is_not_mistaken_for_a_topology_only_tree():
    # Act
    problems = validate_species_tree(_good_tree(), {"S1", "S2", "S3"})

    # Assert
    assert not any("identical" in p for p in problems)


def test_two_branch_tree_is_too_small_to_call_uniformity_suspicious():
    # Arrange: (S1:0.1,S2:0.1) is a legitimate two-taxon estimate
    pair = FakeTree(
        leaf_dists={"S1": 0.1, "S2": 0.1},
        pairwise={("S1", "S2"): 0.2},
    )

    # Act
    problems = validate_species_tree(pair, {"S1", "S2"})

    # Assert
    assert not any("identical" in p for p in problems)
