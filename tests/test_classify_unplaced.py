"""Tests for classify_unplaced.py (issue #36).

Pure classification logic only — no OrthoFinder, no DIAMOND. The point of the
module is to say WHY a gene stayed unplaced, so every verdict is pinned here.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import classify_unplaced as cu


# ------------------------------------------------------------ parsing ----

ORTHOGROUPS = "\n".join([
    "Orthogroup\tCgig\tMcry\tObas",
    "OG0000001\tCgig_a, Cgig_b\tMcry_x\tObas_p",
    "OG0000002\t\tMcry_y, Mcry_z\t",
    "OG0000003\tCgig_c\t\t",
]) + "\n"


def test_orthogroup_members_are_collected_across_species_columns(tmp_path):
    p = tmp_path / "og.tsv"
    p.write_text(ORTHOGROUPS)

    groups = cu.parse_orthogroups(p)

    assert groups["OG0000001"] == ["Cgig_a", "Cgig_b", "Mcry_x", "Obas_p"]
    assert groups["OG0000002"] == ["Mcry_y", "Mcry_z"]
    assert groups["OG0000003"] == ["Cgig_c"]


def test_empty_species_cells_do_not_become_phantom_members(tmp_path):
    p = tmp_path / "og.tsv"
    p.write_text(ORTHOGROUPS)

    groups = cu.parse_orthogroups(p)

    assert "" not in groups["OG0000002"]
    assert len(groups["OG0000002"]) == 2


# -------------------------------------------------------- observations ----

def test_observation_counts_own_and_other_species_members():
    obs = cu.observe(["Mcry_x", "Mcry_y", "Cgig_a"], "Mcry")

    assert obs == cu.Observation(n_genes=3, n_other_species=1)


def test_observation_of_a_species_specific_cluster_has_no_other_species():
    obs = cu.observe(["Mcry_x", "Mcry_y"], "Mcry")

    assert obs == cu.Observation(n_genes=2, n_other_species=0)


# ----------------------------------------------------- classification ----

def test_gene_in_a_cross_species_orthogroup_above_the_size_gate_was_pruned():
    obs = [cu.Observation(n_genes=6, n_other_species=4)]

    assert cu.classify_gene(obs, has_cross_species_hit=True) == "PRUNED"


def test_cross_species_orthogroup_below_the_size_gate_is_a_splinter():
    obs = [cu.Observation(n_genes=3, n_other_species=2)]

    assert cu.classify_gene(obs, has_cross_species_hit=True) == "SPLINTER"


def test_size_gate_boundary_is_inclusive_at_min_orthogroup_size():
    assert cu.classify_gene([cu.Observation(4, 2)], True) == "PRUNED"
    assert cu.classify_gene([cu.Observation(3, 2)], True) == "SPLINTER"


def test_clusters_holding_only_the_genes_own_species_are_lineage_specific():
    obs = [cu.Observation(n_genes=7, n_other_species=0)]

    assert cu.classify_gene(obs, has_cross_species_hit=False) == "LINEAGE_ONLY"


def test_never_clustered_but_diamond_found_a_cross_species_homolog():
    assert cu.classify_gene([], has_cross_species_hit=True) == "SINGLETON_HIT"
    assert cu.classify_gene([cu.Observation(1, 0)], True) == "SINGLETON_HIT"


def test_never_clustered_and_no_homolog_anywhere_is_a_true_orphan():
    assert cu.classify_gene([], has_cross_species_hit=False) == "SINGLETON_NOHIT"


def test_the_strongest_evidence_across_rounds_wins():
    # unassigned in round 1, splinter in round 2, properly clustered in round 3
    obs = [cu.Observation(1, 0), cu.Observation(3, 1), cu.Observation(9, 5)]

    assert cu.classify_gene(obs, has_cross_species_hit=True) == "PRUNED"


def test_a_cross_species_splinter_outranks_a_larger_lineage_only_cluster():
    # being cut away from other species is the more informative failure
    obs = [cu.Observation(20, 0), cu.Observation(2, 1)]

    assert cu.classify_gene(obs, has_cross_species_hit=True) == "SPLINTER"


# ------------------------------------------------------------ rollup ----

def test_rollup_maps_five_verdicts_onto_the_three_buckets_the_issue_asked_for():
    counts = {"PRUNED": 10, "SPLINTER": 20, "SINGLETON_HIT": 30,
              "LINEAGE_ONLY": 40, "SINGLETON_NOHIT": 50}

    rolled = cu.rollup(counts)

    assert rolled == {"splinter_or_graph_cut": 50,   # SPLINTER + SINGLETON_HIT
                      "lineage_specific": 40,
                      "true_orphan": 50,
                      "pruned": 10}


def test_rollup_reports_zero_for_verdicts_that_never_occurred():
    rolled = cu.rollup({"SINGLETON_NOHIT": 3})

    assert rolled == {"splinter_or_graph_cut": 0, "lineage_specific": 0,
                      "true_orphan": 3, "pruned": 0}
