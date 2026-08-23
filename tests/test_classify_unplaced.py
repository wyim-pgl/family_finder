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


# ----------------------------------------------------- CDS integrity ----
#
# Issue #36: a short model with a broken CDS and no cross-species hit is all
# three things at once. Integrity is measured per gene so that "true orphan"
# can be told apart from "the model itself is not a gene".

def _cds(*codons):
    return "".join(codons)


def test_a_cds_with_start_stop_frame_and_no_internal_stop_is_complete():
    integrity = cu.assess_cds(_cds("ATG", "GCT", "TAA"))

    assert integrity.has_start and integrity.has_stop
    assert integrity.in_frame and integrity.no_internal_stop
    assert integrity.complete


def test_a_cds_without_a_start_codon_is_incomplete():
    integrity = cu.assess_cds(_cds("GCT", "GCT", "TAA"))

    assert not integrity.has_start
    assert not integrity.complete


def test_a_cds_without_a_terminal_stop_codon_is_incomplete():
    integrity = cu.assess_cds(_cds("ATG", "GCT", "GCT"))

    assert not integrity.has_stop
    assert not integrity.complete


def test_a_cds_whose_length_is_not_divisible_by_three_is_out_of_frame():
    integrity = cu.assess_cds(_cds("ATG", "GCT", "TAA") + "G")

    assert not integrity.in_frame
    assert not integrity.complete


def test_an_internal_stop_codon_makes_the_cds_incomplete():
    integrity = cu.assess_cds(_cds("ATG", "TGA", "GCT", "TAA"))

    assert not integrity.no_internal_stop
    assert not integrity.complete


def test_the_terminal_stop_is_not_counted_as_an_internal_stop():
    integrity = cu.assess_cds(_cds("ATG", "GCT", "TGA"))

    assert integrity.no_internal_stop
    assert integrity.complete


def test_all_three_stop_codons_are_recognised():
    for stop in ("TAA", "TAG", "TGA"):
        assert cu.assess_cds(_cds("ATG", "GCT", stop)).complete, stop


def test_case_and_whitespace_do_not_change_the_verdict():
    assert cu.assess_cds("atg gct taa\n").complete


def test_an_absent_cds_is_neither_complete_nor_reported_as_broken():
    integrity = cu.assess_cds(None)

    assert not integrity.present
    assert not integrity.complete
    assert integrity.label == "absent"


def test_integrity_labels_separate_absent_from_incomplete():
    assert cu.assess_cds(_cds("ATG", "GCT", "TAA")).label == "complete"
    assert cu.assess_cds(_cds("GCT", "GCT", "TAA")).label == "incomplete"


# ------------------------------------------------------------- flags ----

COMPLETE = _cds("ATG", "GCT", "TAA")
BROKEN = _cds("ATG", "TGA", "GCT", "TAA")


def test_a_protein_below_the_floor_is_flagged_short():
    flags = cu.flags_for(60, cu.assess_cds(COMPLETE), floor=100)

    assert flags == ("SHORT_PROTEIN",)


def test_the_floor_is_inclusive_so_a_gene_at_the_floor_is_not_short():
    assert cu.flags_for(100, cu.assess_cds(COMPLETE), floor=100) == ()


def test_a_broken_cds_is_flagged_invalid():
    assert cu.flags_for(300, cu.assess_cds(BROKEN), floor=100) == ("INVALID_CDS",)


def test_a_missing_cds_is_flagged_separately_from_a_broken_one():
    # not knowing is not the same as knowing it is broken
    assert cu.flags_for(300, cu.assess_cds(None), floor=100) == ("NO_CDS",)


def test_one_gene_can_carry_several_flags_at_once():
    flags = cu.flags_for(60, cu.assess_cds(BROKEN), floor=100)

    assert set(flags) == {"SHORT_PROTEIN", "INVALID_CDS"}


def test_a_long_gene_with_a_complete_cds_carries_no_flags():
    assert cu.flags_for(300, cu.assess_cds(COMPLETE), floor=100) == ()


def test_comparable_requires_both_the_length_floor_and_a_complete_cds():
    assert cu.is_comparable(300, cu.assess_cds(COMPLETE).complete, floor=100)
    assert not cu.is_comparable(60, cu.assess_cds(COMPLETE).complete, floor=100)
    assert not cu.is_comparable(300, cu.assess_cds(BROKEN).complete, floor=100)
    assert not cu.is_comparable(300, cu.assess_cds(None).complete, floor=100)


# ------------------------------------- primary reason plus flags rows ----

def _rows(genes, observations=None, with_hit=(), lengths=None, cds=None):
    return cu.classify_all(
        genes,
        observations or {},
        set(with_hit),
        lengths or {g: 300 for g in genes},
        cds if cds is not None else {g: COMPLETE for g in genes},
    )


def test_every_unplaced_gene_comes_out_with_a_primary_reason():
    genes = {"Mcry_a", "Mcry_b", "Mcry_c", "Mcry_d"}
    rows = _rows(genes,
                 observations={"Mcry_a": [cu.Observation(6, 3)],
                               "Mcry_b": [cu.Observation(3, 1)],
                               "Mcry_c": [cu.Observation(5, 0)]},
                 with_hit={"Mcry_d"})

    assert cu.genes_without_a_reason(rows, genes) == set()
    assert {r.gene for r in rows} == genes
    assert all(r.primary_reason in cu.VERDICTS for r in rows)


def test_a_gene_missing_from_the_rows_is_reported_not_silently_dropped():
    rows = _rows({"Mcry_a"})

    assert cu.genes_without_a_reason(rows, {"Mcry_a", "Mcry_b"}) == {"Mcry_b"}


def test_flags_are_additive_and_never_change_the_primary_reason():
    """The known answer in resume.md 1.5 is a verdict distribution. Flags must
    not move it, so the same observations must classify identically whether
    the models are pristine or wrecked."""
    genes = {"Mcry_a", "Mcry_b", "Mcry_c", "Mcry_d"}
    observations = {"Mcry_a": [cu.Observation(6, 3)],
                    "Mcry_b": [cu.Observation(3, 1)],
                    "Mcry_c": [cu.Observation(5, 0)]}

    pristine = _rows(genes, observations, with_hit={"Mcry_d"})
    wrecked = _rows(genes, observations, with_hit={"Mcry_d"},
                    lengths={g: 30 for g in genes},
                    cds={g: BROKEN for g in genes})

    assert ({r.gene: r.primary_reason for r in pristine}
            == {r.gene: r.primary_reason for r in wrecked})
    assert all(r.flags == () for r in pristine)
    assert all(set(r.flags) == {"SHORT_PROTEIN", "INVALID_CDS"} for r in wrecked)


def test_a_row_carries_the_evidence_behind_its_reason_and_its_flags():
    rows = _rows({"Mcry_a"}, {"Mcry_a": [cu.Observation(6, 3), cu.Observation(2, 1)]},
                 lengths={"Mcry_a": 250})

    row = rows[0]
    assert row.primary_reason == "PRUNED"
    assert (row.max_og_size, row.max_other_species, row.n_rounds_seen) == (6, 3, 2)
    assert row.protein_length == 250
    assert row.cds_status == "complete"
    assert row.comparable


def test_a_gene_with_no_cds_record_is_not_counted_as_comparable():
    rows = _rows({"Mcry_a"}, cds={})

    assert rows[0].flags == ("NO_CDS",)
    assert not rows[0].comparable


# --------------------------------------------------------- CLI output ----

def _make_run(tmp_path, proteins, cds):
    """proteins/cds: {gene: sequence}. One round, no WorkingDirectory."""
    (tmp_path / "unplaced_proteins.fa").write_text(
        "".join(f">{g}\n{s}\n" for g, s in proteins.items()))
    if cds is not None:
        (tmp_path / "unplaced_cds.fa").write_text(
            "".join(f">{g}\n{s}\n" for g, s in cds.items()))
    og = tmp_path / "round_01" / "orthofinder" / "Results_x" / "Orthogroups"
    og.mkdir(parents=True)
    (og / "Orthogroups.tsv").write_text(
        "Orthogroup\tCgig\tMcry\n"
        "OG0000001\tCgig_a, Cgig_b, Cgig_c\tMcry_pruned\n"
        "OG0000002\tCgig_d\tMcry_splinter\n"
    )
    return tmp_path


def _run_cli(tmp_path, *extra):
    run = _make_run(tmp_path, {
        "Mcry_pruned": "M" * 300,
        "Mcry_splinter": "M" * 300,
        "Mcry_short": "M" * 60,
        "Mcry_nocds": "M" * 300,
    }, {
        "Mcry_pruned": COMPLETE,
        "Mcry_splinter": BROKEN,
        "Mcry_short": COMPLETE,
    })
    out = run / "genes.tsv"
    cu.main(["--run-dir", str(run), "--species", "Mcry", "--out", str(out), *extra])
    rows = [line.split("\t") for line in
            out.read_text().rstrip("\n").split("\n")[1:]]
    return {r[0]: r for r in rows}


def test_cli_writes_a_primary_reason_and_flags_per_gene(tmp_path, capsys):
    rows = _run_cli(tmp_path)
    header_free = capsys.readouterr().out

    assert rows["Mcry_pruned"][1] == "PRUNED" and rows["Mcry_pruned"][2] == "-"
    assert rows["Mcry_splinter"][1] == "SPLINTER"
    assert rows["Mcry_splinter"][2] == "INVALID_CDS"
    assert rows["Mcry_short"][2] == "SHORT_PROTEIN"
    assert rows["Mcry_nocds"][2] == "NO_CDS"
    assert "no primary reason: 0 (0.0%)" in header_free


def test_cli_reports_the_comparable_subset_and_the_flag_counts(tmp_path, capsys):
    _run_cli(tmp_path)
    out = capsys.readouterr().out

    assert "comparable subset (>= 100 aa AND complete CDS)" in out
    assert "1 of 4" in out           # only Mcry_pruned is long and clean
    assert "SHORT_PROTEIN" in out and "INVALID_CDS" in out


def test_cli_floor_is_configurable(tmp_path):
    rows = _run_cli(tmp_path, "--floor", "50")

    assert rows["Mcry_short"][2] == "-"      # 60 aa clears a 50 aa floor
