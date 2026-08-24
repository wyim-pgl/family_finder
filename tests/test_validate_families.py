"""The family-boundary campaign, in the repository instead of on one cluster.

`fragment_verdict` had the same disease `detect_merge_candidates` had: it
existed, it was tested, and nothing outside its own tests ever called it. The
5,618-cluster campaign that produced summary_v3.tsv ran from ad-hoc scripts on
pronghorn (`judge_one.py`, `array.sbatch`), so the published family table could
not be reproduced, re-run at a different threshold, or applied to any other
dataset - which is the whole point of a pipeline.

This module holds the pure parts, so they test without MAFFT, FastTree,
pal2nal, or a cluster: read the clusters, apply the verdicts to the family
table, and prove the gene set survives. The subprocess parts stay in thin
wrappers that get monkeypatched, per the repository's test conventions.
"""
import sys
import types
from pathlib import Path

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from validate_families import (apply_merges, load_clusters, load_families,
                               merge_groups_from_rows)


def _write(path, rows):
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))


# ---------------------------------------------------------------------------
# reading the run's own outputs
# ---------------------------------------------------------------------------

def test_load_families_reads_members_and_sizes(tmp_path):
    # Arrange
    _write(tmp_path / "summary.tsv", [
        ("family_id", "round", "n_genes", "n_species", "gene_list"),
        ("R1_OG0000001", 1, 2, 2, "Sp1_a,Sp2_b"),
        ("R2_OG0000002", 2, 1, 1, "Sp1_c"),
    ])

    # Act
    fams = load_families(tmp_path / "summary.tsv")

    # Assert
    assert fams == {"R1_OG0000001": ["Sp1_a", "Sp2_b"],
                    "R2_OG0000002": ["Sp1_c"]}


def test_load_clusters_reads_membership(tmp_path):
    # Arrange
    _write(tmp_path / "fragmentation_clusters.tsv", [
        ("cluster_id", "n_families", "n_genes", "families"),
        ("C0001", 2, 5, "R1_OG0000001,R2_OG0000002"),
    ])

    # Act
    clusters = load_clusters(tmp_path / "fragmentation_clusters.tsv")

    # Assert
    assert clusters == {"C0001": ["R1_OG0000001", "R2_OG0000002"]}


# ---------------------------------------------------------------------------
# verdict rows -> merge groups
# ---------------------------------------------------------------------------

def test_merge_groups_come_only_from_interleaved_rows():
    # Arrange: MONOPHYLETIC is explicitly undecided, never merged
    rows = [
        ("C1", "A", "INTERLEAVED", 5, 0, "A+B", "members mix with: B"),
        ("C1", "B", "INTERLEAVED", 4, 0, "A+B", "members mix with: A"),
        ("C1", "M", "MONOPHYLETIC", 3, 0, "", "topology alone cannot..."),
    ]

    # Act
    groups = merge_groups_from_rows(rows)

    # Assert
    assert groups == [["A", "B"]]


def test_merge_groups_ignore_singleton_groups():
    # Arrange: a group naming only itself merges nothing
    rows = [("C1", "A", "INTERLEAVED", 5, 0, "A", "members mix with: X")]

    # Act / Assert
    assert merge_groups_from_rows(rows) == []


def test_merge_groups_are_deduplicated_across_rows():
    # Arrange: every member of a group carries the same group string
    rows = [
        ("C1", "A", "INTERLEAVED", 5, 0, "A+B+C", ""),
        ("C1", "B", "INTERLEAVED", 4, 0, "A+B+C", ""),
        ("C1", "C", "INTERLEAVED", 3, 0, "A+B+C", ""),
    ]

    # Act
    groups = merge_groups_from_rows(rows)

    # Assert
    assert groups == [["A", "B", "C"]]


# ---------------------------------------------------------------------------
# applying merges to the family table
# ---------------------------------------------------------------------------

def test_apply_merges_folds_members_into_the_first_family():
    # Arrange
    families = {"A": ["Sp1_a"], "B": ["Sp2_b"], "C": ["Sp3_c"]}

    # Act
    merged, prov = apply_merges(families, [["A", "B"]])

    # Assert
    assert merged == {"A": ["Sp1_a", "Sp2_b"], "C": ["Sp3_c"]}
    assert prov["A"] == "A+B"


def test_apply_merges_preserves_every_gene_exactly_once():
    # Arrange: the invariant the campaign asserted - no gene lost, none twice
    families = {"A": ["g1", "g2"], "B": ["g3"], "C": ["g4"], "D": ["g5"]}

    # Act
    merged, _ = apply_merges(families, [["A", "B"], ["C", "D"]])

    # Assert
    before = sorted(g for m in families.values() for g in m)
    after = sorted(g for m in merged.values() for g in m)
    assert after == before
    assert len(after) == len(set(after))


def test_apply_merges_reduces_the_count_by_group_size_minus_one():
    # Arrange
    families = {c: [f"g{c}"] for c in "ABCDEF"}

    # Act
    merged, _ = apply_merges(families, [["A", "B", "C"], ["D", "E"]])

    # Assert
    assert len(merged) == len(families) - (3 - 1) - (2 - 1)


def test_apply_merges_is_a_no_op_without_groups():
    # Arrange
    families = {"A": ["g1"], "B": ["g2"]}

    # Act
    merged, prov = apply_merges(families, [])

    # Assert
    assert merged == families
    assert prov == {}


def test_apply_merges_refuses_a_group_naming_an_unknown_family():
    # Arrange: a verdict file from a different run must not silently drop
    # genes or invent a family
    families = {"A": ["g1"]}

    # Act / Assert
    try:
        apply_merges(families, [["A", "GHOST"]])
    except KeyError as e:
        assert "GHOST" in str(e)
    else:
        raise AssertionError("expected a KeyError naming the unknown family")


def test_apply_merges_handles_a_family_in_two_groups():
    # Arrange: overlapping groups are one component, not two merges - taking
    # them at face value would copy genes into both.
    families = {"A": ["g1"], "B": ["g2"], "C": ["g3"]}

    # Act
    merged, _ = apply_merges(families, [["A", "B"], ["B", "C"]])

    # Assert
    after = sorted(g for m in merged.values() for g in m)
    assert after == ["g1", "g2", "g3"]
    assert len(merged) == 1
