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

import pytest

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from validate_families import (apply_merges, cluster_membership_fingerprint,
                               judge_cluster, load_clusters, load_families,
                               load_sequence_pool, merge_groups_from_rows,
                               validate_verdict_coverage,
                               validate_verdict_file_stamps,
                               verdict_file_is_current, write_summary_v3,
                               write_verdict_rows)


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


def test_load_families_refuses_duplicate_rows_instead_of_losing_the_first(tmp_path):
    # Arrange: dict assignment used to replace the first row silently, making
    # Sp1_a disappear before the later conservation comparison was computed.
    _write(tmp_path / "summary.tsv", [
        ("family_id", "round", "n_genes", "n_species", "gene_list"),
        ("R1_OG1", 1, 1, 1, "Sp1_a"),
        ("R1_OG1", 1, 1, 1, "Sp2_b"),
    ])

    # Act / Assert
    with pytest.raises(ValueError, match="duplicate family id"):
        load_families(tmp_path / "summary.tsv")


def test_load_clusters_refuses_a_family_in_two_components(tmp_path):
    # Arrange: cluster_of would otherwise keep only the last component while
    # the judge/apply paths consumed both.
    _write(tmp_path / "fragmentation_clusters.tsv", [
        ("cluster_id", "n_families", "n_genes", "families"),
        ("C1", 2, 2, "A,B"),
        ("C2", 2, 2, "B,C"),
    ])

    # Act / Assert
    with pytest.raises(ValueError, match="both C1 and C2"):
        load_clusters(tmp_path / "fragmentation_clusters.tsv")


def test_sequence_pool_refuses_duplicate_ids_across_fastas(tmp_path):
    # Arrange: dict.update previously made the selected peptide depend on file
    # name order while every expected id still appeared present.
    (tmp_path / "a.fa").write_text(">Sp1_g\nMAAA\n")
    (tmp_path / "b.fasta").write_text(">Sp1_g\nMBBB\n")

    # Act / Assert
    with pytest.raises(ValueError, match="duplicate peptide id 'Sp1_g'"):
        load_sequence_pool(tmp_path, "peptide")


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


def test_validate_verdict_coverage_refuses_missing_fragment_rows():
    # Arrange: one verdict file exists for the cluster, but one family row is
    # missing. Counting files alone would accept this as a finished campaign.
    clusters = {"C1": ["A", "B"]}
    rows = [("C1", "A", "INTERLEAVED", 2, 0, "A+B", "members mix with: B")]

    # Act / Assert
    try:
        validate_verdict_coverage(rows, clusters)
    except ValueError as e:
        assert "missing fragment rows" in str(e)
        assert "C1: B" in str(e)
    else:
        raise AssertionError("expected truncated verdict coverage to be refused")


def test_validate_verdict_coverage_refuses_a_cross_cluster_merge_group():
    rows = [
        ("C1", "A", "INTERLEAVED", "1", "0", "A+C", ""),
        ("C1", "B", "MONOPHYLETIC", "1", "0", "", ""),
        ("C2", "C", "INTERLEAVED", "1", "0", "A+C", ""),
        ("C2", "D", "MONOPHYLETIC", "1", "0", "", ""),
    ]

    with pytest.raises(ValueError, match="crosses the cluster boundary"):
        validate_verdict_coverage(rows, {"C1": ["A", "B"], "C2": ["C", "D"]})


def test_judge_cluster_uses_protein_tree_when_cds_is_partial(tmp_path, monkeypatch):
    # Arrange: passing a partial CDS set to pal2nal would let the codon tree
    # omit the missing members, silently changing the verdict's leaf set.
    captured = {}

    def fake_build_cluster_tree(pep, cds, workdir, cfg):
        captured["cds"] = cds
        tree = workdir / "tree.nwk"
        tree.write_text("((Sp1_a,Sp1_b),(Sp2_c,Sp2_d));")
        return tree

    monkeypatch.setattr("validate_families.build_cluster_tree", fake_build_cluster_tree)

    families = {"A": ["Sp1_a", "Sp1_b"], "B": ["Sp2_c", "Sp2_d"]}
    pep_pool = {g: "MPEP" for genes in families.values() for g in genes}
    cds_pool = {"Sp1_a": "ATG", "Sp1_b": "ATG", "Sp2_c": "ATG"}  # Sp2_d missing
    cfg = type("Cfg", (), {"mafft_bin": "mafft", "fasttree_bin": "FastTree",
                           "pal2nal": "pal2nal.pl", "threads": 1})()

    # Act
    rows = judge_cluster("C1", ["A", "B"], families, pep_pool, cds_pool,
                         tmp_path, cfg)

    # Assert
    assert captured["cds"] is None
    assert [r[1] for r in rows] == ["A", "B"]


def test_judge_cluster_refuses_a_family_member_absent_from_pep_pool(
    tmp_path, monkeypatch,
):
    # Arrange: the old path wrote a subset FASTA and let fragment_verdict call
    # the omitted member merely "missing", producing a plausible judgement.
    monkeypatch.setattr(
        "validate_families.build_cluster_tree",
        lambda *a, **k: pytest.fail("tree building must not start"),
    )
    families = {"A": ["Sp1_a", "Sp2_missing"], "B": ["Sp3_b"]}
    pep_pool = {"Sp1_a": "MPEP", "Sp3_b": "MPEP"}
    cfg = type("Cfg", (), {})()

    # Act / Assert
    with pytest.raises(ValueError, match="Sp2_missing"):
        judge_cluster("C1", ["A", "B"], families, pep_pool, {}, tmp_path, cfg)


def test_judge_cluster_refuses_a_tree_that_drops_an_input_leaf(tmp_path, monkeypatch):
    # Arrange: even with a complete peptide pool, a partial pal2nal/FastTree
    # result must not decide a family boundary from the surviving subset.
    def partial_tree(*args):
        tree = tmp_path / "tree.nwk"
        tree.write_text("(Sp1_a,Sp3_b);")
        return tree

    monkeypatch.setattr("validate_families.build_cluster_tree", partial_tree)
    families = {"A": ["Sp1_a", "Sp2_a"], "B": ["Sp3_b"]}
    pep_pool = {g: "MPEP" for members in families.values() for g in members}

    # Act / Assert
    with pytest.raises(RuntimeError, match="tree leaf set"):
        judge_cluster("C1", ["A", "B"], families, pep_pool, {}, tmp_path,
                      type("Cfg", (), {})())


def test_stamped_resume_file_must_be_complete(tmp_path):
    # Arrange: a matching fingerprint alone is insufficient if a killed writer
    # left only one fragment row behind.
    rows = [
        ("C1", "A", "MONOPHYLETIC", "1", "0", "", "topology alone"),
        ("C1", "B", "MONOPHYLETIC", "1", "0", "", "topology alone"),
    ]
    families = {"A": ["Sp1_a"], "B": ["Sp2_b"]}
    path = tmp_path / "C1.rows.tsv"
    membership = cluster_membership_fingerprint("C1", ["A", "B"], families)
    write_verdict_rows(path, rows, "input-v1", membership)
    assert verdict_file_is_current(path, "input-v1", "C1", ["A", "B"], families)

    # A file from the short-lived pre-membership stamp implementation must be
    # refreshed by --judge, otherwise a later --apply can only reject it.
    without_membership = "\n".join(
        line for line in path.read_text().splitlines()
        if not line.startswith("# judge_membership\t")
    ) + "\n"
    path.write_text(without_membership)
    assert not verdict_file_is_current(
        path, "input-v1", "C1", ["A", "B"], families,
    )
    write_verdict_rows(path, rows, "input-v1", membership)

    # Act: retain the stamps but truncate one data row.
    path.write_text("\n".join(path.read_text().splitlines()[:-1]) + "\n")

    # Assert
    assert not verdict_file_is_current(
        path, "input-v1", "C1", ["A", "B"], families,
    )


def test_apply_stamp_refuses_same_family_counts_with_changed_membership(tmp_path):
    # Arrange: the old and new runs use the same cluster/family ids and each
    # family still has one gene, so row coverage and n_members both pass. The
    # tree verdict belongs to Sp1_old, however, not its replacement Sp1_new.
    clusters = {"C1": ["A", "B"]}
    old_families = {"A": ["Sp1_old"], "B": ["Sp2_b"]}
    new_families = {"A": ["Sp1_new"], "B": ["Sp2_b"]}
    rows = [
        ("C1", "A", "INTERLEAVED", "1", "0", "A+B", ""),
        ("C1", "B", "INTERLEAVED", "1", "0", "A+B", ""),
    ]
    path = tmp_path / "C1.rows.tsv"
    old_membership = cluster_membership_fingerprint(
        "C1", clusters["C1"], old_families,
    )
    write_verdict_rows(path, rows, "old-full-input", old_membership)
    validate_verdict_file_stamps(path, rows, clusters, old_families)

    # Act / Assert
    with pytest.raises(ValueError, match="stale verdict membership"):
        validate_verdict_file_stamps(
            path, rows, clusters, new_families,
        )


def test_write_summary_uses_the_source_runs_species_delimiter(tmp_path):
    original = tmp_path / "summary.tsv"
    _write(original, [
        ("family_id", "round", "n_genes", "n_species", "gene_list"),
        ("R1_OG1", 1, 2, 1, "Sp1|a,Sp1|b"),
    ])
    out = tmp_path / "summary_v3.tsv"

    # Act
    write_summary_v3(
        {"R1_OG1": ["Sp1|a", "Sp1|b"]}, {}, original, out, {}, "|",
    )

    # Assert: hardcoding '_' counted each complete gene id as a species.
    assert out.read_text().splitlines()[1].split("\t")[3] == "1"


def test_write_summary_refuses_missing_round_metadata(tmp_path):
    original = tmp_path / "summary.tsv"
    _write(original, [
        ("family_id", "n_genes", "n_species", "gene_list"),
        ("R1_OG1", 1, 1, "Sp1_a"),
    ])

    with pytest.raises(ValueError, match="round metadata"):
        write_summary_v3(
            {"R1_OG1": ["Sp1_a"]}, {}, original,
            tmp_path / "summary_v3.tsv", {},
        )


def test_write_summary_refuses_ids_without_a_species_delimiter(tmp_path):
    original = tmp_path / "summary.tsv"
    _write(original, [
        ("family_id", "round", "n_genes", "n_species", "gene_list"),
        ("R1_OG1", 1, 1, 1, "no_species_prefix"),
    ])

    with pytest.raises(ValueError, match="lacks delimiter"):
        write_summary_v3(
            {"R1_OG1": ["missingdelimiter"]}, {}, original,
            tmp_path / "summary_v3.tsv", {},
        )


# ---------------------------------------------------------------------------
# outgroup selection
# ---------------------------------------------------------------------------
# The tree can only say "one family" when it has something outside the family
# to compare against, and vote_edges already ranks every family's neighbours -
# so the nearest families OUTSIDE the cluster are the outgroup, for free.

def test_outgroup_families_are_the_nearest_non_members():
    # Arrange: the cluster is {A,B}; C and D are ranked neighbours outside it
    from validate_families import pick_outgroup_families
    edges = [("A", "C", 8, 10, 0.8), ("B", "D", 5, 10, 0.5),
             ("A", "B", 9, 10, 0.9)]

    # Act
    out = pick_outgroup_families(["A", "B"], edges, n=2)

    # Assert: highest-frac non-members first, cluster members never chosen
    assert out == ["C", "D"]


def test_outgroup_selection_never_returns_a_cluster_member():
    # Arrange: every edge points back inside the cluster
    from validate_families import pick_outgroup_families
    edges = [("A", "B", 9, 10, 0.9), ("B", "A", 8, 10, 0.8)]

    # Act
    out = pick_outgroup_families(["A", "B"], edges, n=2)

    # Assert: an outgroup drawn from inside the family would be circular
    assert out == []


def test_outgroup_selection_is_capped_and_ordered():
    # Arrange
    from validate_families import pick_outgroup_families
    edges = [("A", "X", 1, 10, 0.1), ("A", "Y", 5, 10, 0.5),
             ("A", "Z", 3, 10, 0.3)]

    # Act
    out = pick_outgroup_families(["A"], edges, n=2)

    # Assert
    assert out == ["Y", "Z"]


def test_outgroup_selection_deduplicates():
    # Arrange: two cluster members can name the same neighbour
    from validate_families import pick_outgroup_families
    edges = [("A", "X", 9, 10, 0.9), ("B", "X", 8, 10, 0.8)]

    # Act
    out = pick_outgroup_families(["A", "B"], edges, n=3)

    # Assert
    assert out == ["X"]


def _cfg():
    return type("Cfg", (), {"mafft_bin": "mafft", "fasttree_bin": "FastTree",
                            "pal2nal": "pal2nal.pl", "threads": 1})


def test_judge_cluster_puts_the_outgroup_in_the_tree_and_the_verdict(tmp_path, monkeypatch):
    # Arrange: with an outgroup present the verdict must be the outgroup one,
    # and the outgroup sequences must actually reach the alignment - a verdict
    # naming an outgroup that never entered the tree would be a silent lie.
    import validate_families as vf
    seen = {}

    def fake_build(pep, cds, workdir, cfg):
        from utils.seqio import read_fasta
        seen["pep_ids"] = sorted(read_fasta(str(pep)))
        tree = workdir / "tree.nwk"
        tree.write_text("(((Sp1_a,Sp2_b),og1),(og2,Sp3_c));")
        return tree

    monkeypatch.setattr(vf, "build_cluster_tree", fake_build)
    families = {"A": ["Sp1_a"], "B": ["Sp2_b"], "C": ["Sp3_c"],
                "OG": ["og1", "og2"]}
    pool = {g: "MPEP" for m in families.values() for g in m}

    # Act
    rows = vf.judge_cluster("C1", ["A", "B", "C"], families, pool, {},
                            tmp_path, _cfg(), outgroup_families=["OG"])

    # Assert
    assert seen["pep_ids"] == ["Sp1_a", "Sp2_b", "Sp3_c", "og1", "og2"]
    status = {r[1]: r[2] for r in rows}
    assert status == {"A": "SAME_FAMILY", "B": "SAME_FAMILY",
                      "C": "SEPARATE"}
    assert not any(r[1] == "OG" for r in rows)


def test_judge_cluster_without_an_outgroup_is_unchanged(tmp_path, monkeypatch):
    # Arrange
    import validate_families as vf

    def fake_build(pep, cds, workdir, cfg):
        tree = workdir / "tree.nwk"
        tree.write_text("((Sp1_a,Sp2_b),(Sp3_c,Sp4_d));")
        return tree

    monkeypatch.setattr(vf, "build_cluster_tree", fake_build)
    families = {"A": ["Sp1_a", "Sp2_b"], "B": ["Sp3_c", "Sp4_d"]}
    pool = {g: "MPEP" for m in families.values() for g in m}

    # Act
    rows = vf.judge_cluster("C1", ["A", "B"], families, pool, {}, tmp_path,
                            _cfg())

    # Assert
    assert {r[2] for r in rows} == {"MONOPHYLETIC"}


def test_the_edges_that_built_a_cluster_cannot_supply_its_outgroup():
    # Arrange: a cluster IS a connected component of the vote edges, so by
    # construction no edge leaves it. Measured on the shipped 15sp file: all
    # four edges out of C0297's five PEPC fragments point back inside, and the
    # file's minimum frac is exactly 0.600 - the cut. An outgroup therefore
    # has to come from edges BELOW the component threshold, which only the
    # uncut dump (vote_edges at min_frac=0.0) contains.
    from validate_families import pick_outgroup_families
    component_edges = [("A", "B", 9, 10, 0.9), ("B", "C", 8, 10, 0.8)]
    with_subthreshold = component_edges + [("A", "X", 2, 10, 0.2)]

    # Act
    from_component = pick_outgroup_families(["A", "B", "C"], component_edges, 3)
    from_full = pick_outgroup_families(["A", "B", "C"], with_subthreshold, 3)

    # Assert
    assert from_component == []
    assert from_full == ["X"]


def test_explicit_outgroup_families_are_used_verbatim():
    # Arrange: the defensible outgroup is an externally justified one, not a
    # below-cut graph neighbour. For PEPC that is the BTPC group: scored
    # against the four Arabidopsis anchors on 15-species sequences alone, the
    # five PTPC fragments sit 700-1000 bits nearer PPC1/2/3 while
    # R1_OG0009826 sits 812 bits nearer PPC4. One sign flip, not a gradient.
    from validate_families import resolve_outgroup
    edges = [("A", "X", 2, 10, 0.2)]

    # Act
    explicit = resolve_outgroup(["A", "B"], edges, n_from_edges=3,
                                explicit=["OG_BTPC"])
    from_edges = resolve_outgroup(["A", "B"], edges, n_from_edges=3,
                                  explicit=None)

    # Assert: an explicit panel wins over anything the graph would pick
    assert explicit == ["OG_BTPC"]
    assert from_edges == ["X"]


def test_explicit_outgroup_may_not_name_a_cluster_member():
    # Arrange
    from validate_families import resolve_outgroup

    # Act / Assert
    try:
        resolve_outgroup(["A", "B"], [], n_from_edges=0, explicit=["B"])
    except ValueError as e:
        assert "B" in str(e)
    else:
        raise AssertionError("expected a cluster member to be refused")


# ---------------------------------------------------------------------------
# the outgroup verdicts must actually reach the family table
# ---------------------------------------------------------------------------
# Judging emitted SAME_FAMILY/SEPARATE/OUTGROUP_SPLIT while apply accepted only
# INTERLEAVED/MONOPHYLETIC, so an outgroup campaign would have been rejected or
# silently merged nothing. The tests below are the seam nobody was testing.

def test_same_family_groups_are_applied():
    # Arrange
    rows = [("C1", "A", "SAME_FAMILY", 5, 0, "A+B", "groups with B"),
            ("C1", "B", "SAME_FAMILY", 4, 0, "A+B", "groups with A")]

    # Act
    groups = merge_groups_from_rows(rows)

    # Assert
    assert groups == [["A", "B"]]


def test_separate_and_outgroup_split_never_merge():
    # Arrange: SEPARATE is a decision NOT to merge; OUTGROUP_SPLIT is a
    # conflict report. Neither may put a fragment into a group.
    rows = [("C1", "A", "SEPARATE", 5, 0, "", "outgroup lies between"),
            ("C1", "B", "OUTGROUP_SPLIT", 4, 0, "", "both sides")]

    # Act / Assert
    assert merge_groups_from_rows(rows) == []


def test_apply_accepts_the_outgroup_statuses():
    # Arrange
    from validate_families import validate_verdict_coverage
    clusters = {"C1": ["A", "B", "C"]}
    rows = [("C1", "A", "SAME_FAMILY", 2, 0, "A+B", ""),
            ("C1", "B", "SAME_FAMILY", 2, 0, "A+B", ""),
            ("C1", "C", "SEPARATE", 1, 0, "", "")]

    # Act / Assert: must not raise on a status it has never seen before
    validate_verdict_coverage(rows, clusters)
