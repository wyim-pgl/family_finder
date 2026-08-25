"""Tests for budget-capped nomination clustering.

Pure-dict/list fixtures only. No external tools, no ete4, no pipeline wiring.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from config import Config  # noqa: E402
from steps.nominate_clusters import (  # noqa: E402
    advance_deferred_edges,
    collapse_vote_edges,
    nominate_clusters,
)


def _edge(a, b, votes, size, frac):
    return (a, b, votes, size, frac)


def _edge_keys(records):
    return {tuple(rec["edge_key"]) for rec in records}


def test_under_cap_component_passes_through_untouched():
    edges = [
        _edge("A", "B", 9, 10, 0.9),
        _edge("B", "A", 4, 8, 0.5),
        _edge("B", "C", 4, 8, 0.5),
        _edge("C", "B", 2, 6, 0.333),
    ]
    family_sizes = {"A": 10, "B": 8, "C": 6}

    out = nominate_clusters(edges, family_sizes, Config(max_cluster_genes=50))

    assert len(out["clusters"]) == 1
    cluster = out["clusters"][0]
    assert cluster["families"] == ("A", "B", "C")
    assert cluster["total_genes"] == 24
    assert cluster["min_admitted_weight"] is None
    assert cluster["decomposed"] is False
    assert _edge_keys(out["pairwise_nominations"]) == set()
    assert _edge_keys(out["deferred_edges"]) == set()
    assert set(cluster["edge_keys"]) == {("A", "B"), ("B", "C")}


def test_hairball_decomposes_and_every_cut_edge_becomes_pairwise_or_deferred():
    edges = [
        _edge("A", "B", 90, 100, 0.9),
        _edge("B", "A", 80, 100, 0.8),
        _edge("B", "C", 80, 100, 0.8),
        _edge("C", "B", 70, 100, 0.7),
        _edge("C", "D", 70, 100, 0.7),
        _edge("D", "C", 60, 100, 0.6),
        _edge("D", "E", 60, 100, 0.6),
        _edge("E", "D", 50, 100, 0.5),
        _edge("A", "E", 40, 100, 0.4),
    ]
    family_sizes = {fam: 100 for fam in "ABCDE"}

    out = nominate_clusters(edges, family_sizes, Config(max_cluster_genes=250))

    assert [c["families"] for c in out["clusters"]] == [("A", "B"), ("C", "D")]
    assert [c["min_admitted_weight"] for c in out["clusters"]] == [0.9, 0.7]
    assert _edge_keys(out["pairwise_nominations"]) == {("B", "C"), ("D", "E"), ("A", "E")}
    assert _edge_keys(out["deferred_edges"]) == set()
    # All undirected links are accounted exactly once.
    accounted = set()
    for cluster in out["clusters"]:
        for edge_key in cluster["edge_keys"]:
            assert edge_key not in accounted
            accounted.add(edge_key)
    for group in (out["pairwise_nominations"], out["deferred_edges"]):
        for rec in group:
            assert tuple(rec["edge_key"]) not in accounted
            accounted.add(tuple(rec["edge_key"]))
    assert accounted == set(collapse_vote_edges(edges))


def test_mcry_flagship_shape_is_clustered_under_cap_and_survives_as_pairwise_when_target_component_is_oversize():
    edges = [
        _edge("R1_OG0000440", "R1_OG0013449", 59, 63, 0.937),
        _edge("R1_OG0008467", "R1_OG0013449", 17, 18, 0.944),
        _edge("R1_OG0009826", "R1_OG0000440", 15, 17, 0.882),
        _edge("R1_OG0013449", "R1_OG0019668", 15, 15, 1.000),
        _edge("R2_OG0000359", "R1_OG0019668", 2, 9, 0.222),
    ]
    family_sizes = {
        "R1_OG0000440": 63,
        "R1_OG0008467": 18,
        "R1_OG0009826": 17,
        "R1_OG0013449": 15,
        "R1_OG0019668": 6,
        "R2_OG0000359": 9,
    }

    under_cap = nominate_clusters(edges, family_sizes, Config(max_cluster_genes=200))
    assert [c["families"] for c in under_cap["clusters"]] == [tuple(sorted(family_sizes))]
    assert _edge_keys(under_cap["pairwise_nominations"]) == set()

    oversize = nominate_clusters(edges, family_sizes, Config(max_cluster_genes=25))
    assert ("R1_OG0019668", "R2_OG0000359") in _edge_keys(oversize["pairwise_nominations"])
    pair = next(
        rec for rec in oversize["pairwise_nominations"]
        if tuple(rec["edge_key"]) == ("R1_OG0019668", "R2_OG0000359")
    )
    assert pair["total_genes"] == 15
    assert pair["reciprocal"] is False


def test_star_hub_refusals_become_pairwise_nominations():
    edges = [
        _edge("A", "Hub", 60, 60, 1.0),
        _edge("B", "Hub", 60, 60, 1.0),
        _edge("C", "Hub", 60, 60, 1.0),
        _edge("D", "Hub", 60, 60, 1.0),
    ]
    family_sizes = {"A": 60, "B": 60, "C": 60, "D": 60, "Hub": 120}

    out = nominate_clusters(edges, family_sizes, Config(max_cluster_genes=180))

    assert [c["families"] for c in out["clusters"]] == [("A", "Hub")]
    assert _edge_keys(out["pairwise_nominations"]) == {
        ("B", "Hub"), ("C", "Hub"), ("D", "Hub"),
    }


def test_exact_frac_ties_are_deterministic_and_reported():
    edges = [
        _edge("B", "C", 30, 100, 0.3),
        _edge("A", "C", 30, 100, 0.3),
        _edge("A", "B", 30, 100, 0.3),
    ]
    family_sizes = {"A": 100, "B": 100, "C": 100}

    out = nominate_clusters(edges, family_sizes, Config(max_cluster_genes=150))

    assert out["tie_groups"] == [
        {"weight": 0.3, "edges": (("A", "B"), ("A", "C"), ("B", "C"))}
    ]
    assert _edge_keys(out["pairwise_nominations"]) == set()
    assert _edge_keys(out["deferred_edges"]) == {
        ("A", "B"), ("A", "C"), ("B", "C"),
    }


def test_pair_sum_exceeding_cap_is_deferred_even_when_neither_endpoint_exceeds_it():
    edges = [_edge("Mid1", "Mid2", 120, 300, 0.4)]
    family_sizes = {"Mid1": 300, "Mid2": 300}

    out = nominate_clusters(edges, family_sizes, Config(max_cluster_genes=500))

    assert out["clusters"] == []
    assert out["pairwise_nominations"] == []
    assert out["deferred_edges"][0]["reason"] == "pair_exceeds_cap"
    assert out["deferred_edges"][0]["retry_count"] == 0
    assert out["deferred_edges"][0]["terminal"] is False


def _deferred_fixture():
    return [{
        "edge_key": ("Big1", "Big2"),
        "families": ("Big1", "Big2"),
        "weight": 0.4,
        "reason": "pair_exceeds_cap",
        "retry_count": 0,
        "terminal": False,
    }]


def test_deferred_edge_goes_terminal_when_its_own_families_are_unchanged():
    # An unrelated family changing must NOT keep the edge alive: the check
    # is per-edge, or a busy run would retry the same deferral forever.
    advanced = advance_deferred_edges(_deferred_fixture(),
                                      changed_families={"Unrelated"})

    assert advanced[0]["retry_count"] == 1
    assert advanced[0]["terminal"] is True
    assert advanced[0]["reason"] == "unchanged_after_retry"


def test_deferred_edge_is_retried_while_one_of_its_endpoints_changes():
    advanced = advance_deferred_edges(_deferred_fixture(),
                                      changed_families={"Big2"})

    assert advanced[0]["retry_count"] == 0
    assert advanced[0]["terminal"] is False
    assert advanced[0]["reason"] == "pair_exceeds_cap"


def test_duplicate_directed_edge_raises():
    edges = [
        _edge("A", "B", 8, 10, 0.8),
        _edge("A", "B", 3, 10, 0.3),
    ]

    with pytest.raises(ValueError, match="duplicate directed edge"):
        collapse_vote_edges(edges)


def test_weight_and_reciprocity_survive_the_3_decimal_tsv_round_trip():
    # frac serialized as %.3f: 1/2500 lands on disk as 0.000 and a strong
    # reverse edge would read as non-reciprocal. The exact fraction is
    # recomputed from votes/from_size, so neither happens.
    edges = [
        _edge("Tiny", "Huge", 1, 2500, 0.000),
        _edge("Huge", "Tiny", 9, 10, 0.9),
    ]
    links = collapse_vote_edges(edges)

    link = links[("Huge", "Tiny")]
    assert link["weight"] == 0.9
    assert link["reciprocal"] is True
    assert min(link["frac_ab"], link["frac_ba"]) == pytest.approx(1 / 2500)


def test_edge_from_size_disagreeing_with_family_table_raises():
    edges = [_edge("A", "B", 5, 10, 0.5)]
    family_sizes = {"A": 99, "B": 10}

    with pytest.raises(ValueError, match="different runs"):
        nominate_clusters(edges, family_sizes, Config(max_cluster_genes=500))
