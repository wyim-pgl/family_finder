"""HOG levels as the SUBFAMILY axis - measured, not assumed (issue #47).

OrthoFinder already wrote a hierarchy: one HOG table per internal node of the
species tree, N0 (root) through N13 in the 15-species round-1 run. It was
tempting to read that as the family-rank hierarchy the pipeline lacks, so that
"family" could be a declared cut instead of a per-cluster judgement.

It is not, and the run's own numbers say so:

    OG 28,150  ->  N0 HOG 28,161      (more, not fewer)
    OGs split into >1 HOG at N0:  3,153
    HOGs holding genes from >1 OG:    0   <- decisive

A HOG never crosses an orthogroup boundary. The hierarchy SUBDIVIDES; it
cannot merge, and merging is what fragmentation needs. Worse for that purpose,
the tables are per ROUND: the Ococ PEPC flagship is in round 1's N0 while the
Mcry flagship is a round-2 gene absent from that table entirely, so no common
level can hold both.

What it is exactly right for is the level BELOW family. Once a family is
settled - by an outgroup, as PTPC was - these levels give nested subfamily
partitions at every species-tree node, already computed, with the depth as a
free parameter. That is the axis this module exposes.
"""
import sys
import types
from pathlib import Path

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from steps.hog import (hog_level_files, nested_partitions, parse_hog_table,
                       refines, subfamily_partition)


def _hog_table(path, rows, species=("Sp1", "Sp2")):
    header = ["HOG", "OG", "Gene Tree Parent Clade", *species]
    lines = ["\t".join(header)]
    for hog, og, clade, *cols in rows:
        lines.append("\t".join([hog, og, clade, *cols]))
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------

def test_parse_hog_table_maps_every_gene_to_its_hog(tmp_path):
    # Arrange: species columns hold comma-space separated gene lists
    p = tmp_path / "N0.tsv"
    _hog_table(p, [
        ("N0.HOG0000001", "OG0000001", "n1", "Sp1_a, Sp1_b", "Sp2_c"),
        ("N0.HOG0000002", "OG0000002", "n1", "Sp1_d", ""),
    ])

    # Act
    got = parse_hog_table(p)

    # Assert
    assert got == {"Sp1_a": "N0.HOG0000001", "Sp1_b": "N0.HOG0000001",
                   "Sp2_c": "N0.HOG0000001", "Sp1_d": "N0.HOG0000002"}


def test_parse_hog_table_ignores_empty_cells(tmp_path):
    # Arrange
    p = tmp_path / "N0.tsv"
    _hog_table(p, [("N0.HOG0000001", "OG0000001", "n1", "", "  ")])

    # Act / Assert
    assert parse_hog_table(p) == {}


def test_parse_hog_table_refuses_a_gene_in_two_hogs(tmp_path):
    # Arrange: a gene in two HOGs at the same level would make the partition
    # ill-defined and silently double-count it downstream
    p = tmp_path / "N0.tsv"
    _hog_table(p, [
        ("N0.HOG0000001", "OG0000001", "n1", "Sp1_a", ""),
        ("N0.HOG0000002", "OG0000002", "n1", "Sp1_a", ""),
    ])

    # Act / Assert
    try:
        parse_hog_table(p)
    except ValueError as e:
        assert "Sp1_a" in str(e)
    else:
        raise AssertionError("expected a duplicated gene to be refused")


def test_hog_level_files_are_ordered_by_node_number(tmp_path):
    # Arrange: N10 must not sort before N2
    for name in ("N0.tsv", "N2.tsv", "N10.tsv", "N0.ids.tsv", "notes.txt"):
        (tmp_path / name).write_text("HOG\tOG\tGene Tree Parent Clade\n")

    # Act
    got = [p.name for p in hog_level_files(tmp_path)]

    # Assert: the .ids variant uses internal ids and is not a level
    assert got == ["N0.tsv", "N2.tsv", "N10.tsv"]


# ---------------------------------------------------------------------------
# partitions
# ---------------------------------------------------------------------------

def test_subfamily_partition_groups_a_family_by_hog():
    # Arrange
    genes = ["Sp1_a", "Sp1_b", "Sp2_c"]
    hog = {"Sp1_a": "H1", "Sp1_b": "H1", "Sp2_c": "H2"}

    # Act
    got = subfamily_partition(genes, hog)

    # Assert
    assert got == {"H1": ["Sp1_a", "Sp1_b"], "H2": ["Sp2_c"]}


def test_subfamily_partition_reports_genes_with_no_hog():
    # Arrange: HOG tables are per ROUND, so a family merged across rounds has
    # members no single table covers. Counting them is the point - a silent
    # drop would understate the family.
    genes = ["Sp1_a", "Sp2_missing"]
    hog = {"Sp1_a": "H1"}

    # Act
    got = subfamily_partition(genes, hog)

    # Assert
    assert got == {"H1": ["Sp1_a"], None: ["Sp2_missing"]}


def test_refines_accepts_a_deeper_level_that_only_splits():
    # Arrange: nesting is the property that makes a declared depth meaningful
    shallow = {"a": "H1", "b": "H1", "c": "H2"}
    deeper = {"a": "D1", "b": "D2", "c": "D3"}

    # Act / Assert
    assert refines(deeper, shallow) is True


def test_refines_rejects_a_level_that_merges_across_the_shallow_one():
    # Arrange: if a deeper level joined two shallow HOGs it would not be a
    # refinement, and "subfamily depth" would stop being a single number
    shallow = {"a": "H1", "b": "H2"}
    deeper = {"a": "D1", "b": "D1"}

    # Act / Assert
    assert refines(deeper, shallow) is False


def test_nested_partitions_returns_one_partition_per_level():
    # Arrange
    genes = ["a", "b", "c"]
    levels = [("N0", {"a": "H1", "b": "H1", "c": "H1"}),
              ("N1", {"a": "D1", "b": "D1", "c": "D2"})]

    # Act
    got = nested_partitions(genes, levels)

    # Assert
    assert got == [("N0", {"H1": ["a", "b", "c"]}),
                   ("N1", {"D1": ["a", "b"], "D2": ["c"]})]


def test_nested_partitions_is_empty_for_an_empty_family():
    # Arrange / Act / Assert
    assert nested_partitions([], [("N0", {"a": "H1"})]) == [("N0", {})]


def test_truncated_row_is_refused_not_silently_shortened(tmp_path):
    import pytest

    from steps.hog import parse_hog_table

    p = tmp_path / "N1.tsv"
    p.write_text(
        "HOG\tOG\tGene Tree Parent Clade\tSpA\tSpB\n"
        "N1.HOG0000001\tOG0000001\tn0\tSpA_g1\n"   # SpB column lost
    )
    with pytest.raises(ValueError, match="columns"):
        parse_hog_table(p)
