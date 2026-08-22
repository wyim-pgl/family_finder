"""Tests for steps/subfamily.py — generalized subfamily diagnostics.

Three pure layers, each testable with synthetic dicts (house convention):
  * sdp_scan          — subfamily-diagnostic residues on any alignment
  * taxonomic_composition — is a split taxonomy-driven or a paralog split?
  * structure_coherence   — foldseek all-vs-all within/between coherence
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from steps.subfamily import (  # noqa: E402
    clade_rank_label,
    load_taxonomy,
    parse_pairwise_scores,
    sdp_scan,
    species_monophyly,
    structure_coherence,
    taxonomic_composition,
)
from utils.newick import parse_newick  # noqa: E402


# ---------------------------------------------------------------------------
# sdp_scan
# ---------------------------------------------------------------------------

# 8 seqs x 6 cols; groups A (4) and B (4). Col 1 (0-based 0) separates the
# groups perfectly (A=K, B=S); col 2 is globally conserved; col 3 is noisy.
ALN = {
    "Sp1_a1": "KAAWTR", "Sp2_a2": "KACWTR", "Sp3_a3": "KAAWTR", "Sp4_a4": "KADWTR",
    "Sp1_b1": "SAAWTR", "Sp2_b2": "SAEWTR", "Sp3_b3": "SAAWTR", "Sp4_b4": "SAFWTR",
}
GROUPS = {
    "SF_A": ["Sp1_a1", "Sp2_a2", "Sp3_a3", "Sp4_a4"],
    "SF_B": ["Sp1_b1", "Sp2_b2", "Sp3_b3", "Sp4_b4"],
}


def test_sdp_scan_finds_the_separating_column():
    rows = sdp_scan(ALN, GROUPS, min_group=4)
    cols = {(r["subfamily"], r["aln_col"]): r for r in rows}
    a = cols[("SF_A", 1)]
    assert a["sf_residue"] == "K" and a["rest_freq_of_sf_residue"] == 0.0
    b = cols[("SF_B", 1)]
    assert b["sf_residue"] == "S"


def test_sdp_scan_skips_globally_conserved_and_noisy_columns():
    rows = sdp_scan(ALN, GROUPS, min_group=4)
    cols = {(r["subfamily"], r["aln_col"]) for r in rows}
    assert ("SF_A", 2) not in cols  # A everywhere: conserved clan-wide
    assert ("SF_A", 3) not in cols  # too variable inside the group


def test_sdp_scan_reference_numbering_skips_ref_gaps():
    aln = {k: v for k, v in ALN.items()}
    aln["REF_x"] = "K-AWTR"  # gap at col 2 -> ref position shifts
    rows = sdp_scan(aln, GROUPS, min_group=4, ref_seq_id="REF_x")
    a = next(r for r in rows if r["subfamily"] == "SF_B" and r["aln_col"] == 1)
    assert a["ref_pos"] == 1  # first ungapped ref residue


def test_sdp_scan_small_groups_skipped():
    rows = sdp_scan(ALN, {"tiny": ["Sp1_a1"], **GROUPS}, min_group=4)
    assert all(r["subfamily"] != "tiny" for r in rows)


def test_sdp_scan_ragged_alignment_raises():
    bad = dict(ALN, Sp1_a1="KAAWT")  # one seq shorter
    with pytest.raises(ValueError, match="ragged"):
        sdp_scan(bad, GROUPS, min_group=4)


# ---------------------------------------------------------------------------
# taxonomic_composition
# ---------------------------------------------------------------------------

TAXONOMY = {
    "Sp1": {"genus": "Opuntia", "family": "Cactaceae", "order": "Caryophyllales"},
    "Sp2": {"genus": "Carnegiea", "family": "Cactaceae", "order": "Caryophyllales"},
    "Sp3": {"genus": "Talinum", "family": "Talinaceae", "order": "Caryophyllales"},
    "Sp4": {"genus": "Arabidopsis", "family": "Brassicaceae", "order": "Brassicales"},
}


def test_cross_order_subfamily_is_a_paralog_split():
    # SF_A has members from 2 orders -> the split cannot be taxonomy-driven.
    rows = taxonomic_composition(GROUPS, TAXONOMY)
    a = next(r for r in rows if r["subfamily"] == "SF_A")
    assert a["n_species"] == 4
    assert a["order_purity"] < 1.0
    assert a["verdict"] == "paralog-split (crosses order)"


def test_single_species_subfamily_is_lineage_specific():
    groups = {"SF_C": ["Sp1_x1", "Sp1_x2", "Sp1_x3"]}
    rows = taxonomic_composition(groups, TAXONOMY)
    c = rows[0]
    assert c["n_species"] == 1 and c["species_purity"] == 1.0
    assert c["verdict"] == "lineage-specific (species)"


def test_genus_pure_subfamily_reports_rank():
    # Two species, same genus entry: fabricate Sp5 sharing Opuntia.
    taxonomy = dict(TAXONOMY, Sp5={"genus": "Opuntia", "family": "Cactaceae",
                                   "order": "Caryophyllales"})
    groups = {"SF_D": ["Sp1_y1", "Sp5_y2"]}
    d = taxonomic_composition(groups, taxonomy)[0]
    assert d["verdict"] == "lineage-specific (genus)"


def test_family_pure_but_multi_genus_reports_family_rank():
    groups = {"SF_E": ["Sp1_z1", "Sp2_z2"]}  # both Cactaceae, 2 genera
    e = taxonomic_composition(groups, TAXONOMY)[0]
    assert e["verdict"] == "lineage-specific (family)"


def test_unknown_species_counts_but_never_crashes():
    groups = {"SF_F": ["Zzz_q1", "Sp1_q2"]}
    f = taxonomic_composition(groups, TAXONOMY)[0]
    assert f["n_species"] == 2
    assert "unknown" in f["notes"]


def test_load_taxonomy_roundtrip(tmp_path):
    p = tmp_path / "tax.tsv"
    p.write_text("species\tgenus\tfamily\torder\n"
                 "Sp1\tOpuntia\tCactaceae\tCaryophyllales\n")
    tax = load_taxonomy(p)
    assert tax["Sp1"]["order"] == "Caryophyllales"


# ---------------------------------------------------------------------------
# species-tree attribution (issue #27 — default verdict path)
# ---------------------------------------------------------------------------

# ((Sp1,Sp2),(Sp3,Sp4)) — Sp1+Sp2 are sisters (both Cactaceae in TAXONOMY).
SPECIES_TREE = "((Sp1:1,Sp2:1)n12:1,(Sp3:1,Sp4:1)n34:1);"


def _tree():
    return parse_newick(SPECIES_TREE)


def test_monophyletic_subset_is_lineage_specific_with_mrca():
    groups = {
        "SF_CACTUS": ["Sp1_a", "Sp2_b"],
        "SF_ALL": ["Sp1_c", "Sp2_d", "Sp3_e", "Sp4_f"],
    }
    rows = taxonomic_composition(groups, species_tree=_tree())
    cactus = next(r for r in rows if r["subfamily"] == "SF_CACTUS")
    assert cactus["monophyletic"] is True
    assert cactus["mrca_name"] == "n12" and cactus["mrca_depth"] == 1
    assert cactus["verdict"] == "lineage-specific (clade)"


def test_family_spanning_subfamilies_are_paralog_splits():
    # Two ancient paralogs each kept in every species (the PEPC Ppc1/Ppc2
    # shape): both species sets are trivially the root clade, so the
    # divergence predates every speciation -> paralog split, not lineage.
    groups = {
        "SF_P1": ["Sp1_a", "Sp2_b", "Sp3_c", "Sp4_d"],
        "SF_P2": ["Sp1_e", "Sp2_f", "Sp3_g", "Sp4_h"],
    }
    rows = taxonomic_composition(groups, species_tree=_tree())
    assert all(r["verdict"] == "paralog-split (spans family root)"
               for r in rows)


def test_non_monophyletic_set_is_paralog_split_with_interleaving_note():
    groups = {"SF_A": ["Sp1_a", "Sp3_b"], "SF_B": ["Sp2_x", "Sp4_y"]}
    rows = taxonomic_composition(groups, species_tree=_tree())
    a = next(r for r in rows if r["subfamily"] == "SF_A")
    assert a["monophyletic"] is False
    assert a["verdict"] == "paralog-split (non-monophyletic)"
    assert "interleaved species: Sp2,Sp4" in a["notes"]


def test_species_missing_from_tree_judged_on_intersection_and_noted():
    groups = {"SF_M": ["Sp1_a", "Sp2_b", "Zzz_q"], "SF_O": ["Sp3_x"]}
    rows = taxonomic_composition(groups, species_tree=_tree())
    m = next(r for r in rows if r["subfamily"] == "SF_M")
    assert m["n_in_tree"] == 2 and m["monophyletic"] is True
    assert m["verdict"] == "lineage-specific (clade)"
    assert "not in species tree: Zzz" in m["notes"]


def test_no_species_in_tree_gives_unknown_verdict():
    groups = {"SF_Z": ["Zzz_q1", "Yyy_q2"]}
    rows = taxonomic_composition(groups, species_tree=_tree())
    z = rows[0]
    assert z["monophyletic"] is None and z["n_in_tree"] == 0
    assert z["verdict"] == "unknown (no species in species tree)"


def test_taxonomy_is_a_label_layer_on_the_tree_verdict():
    groups = {
        "SF_CACTUS": ["Sp1_a", "Sp2_b"],
        "SF_ALL": ["Sp1_c", "Sp2_d", "Sp3_e", "Sp4_f"],
    }
    # absent -> works, label empty
    bare = taxonomic_composition(groups, species_tree=_tree())
    cactus = next(r for r in bare if r["subfamily"] == "SF_CACTUS")
    assert cactus["clade_label"] == ""
    assert "species_purity" not in cactus
    # present -> same verdict, rank label + purity columns added
    labeled = taxonomic_composition(groups, TAXONOMY, species_tree=_tree())
    cactus_l = next(r for r in labeled if r["subfamily"] == "SF_CACTUS")
    assert cactus_l["verdict"] == cactus["verdict"]
    assert cactus_l["clade_label"] == "Cactaceae (family)"
    assert cactus_l["family_purity"] == 1.0


def test_species_monophyly_pure_helper():
    mono = species_monophyly({"Sp3", "Sp4"}, _tree())
    assert mono["monophyletic"] is True and mono["mrca_name"] == "n34"
    mixed = species_monophyly({"Sp1", "Sp3"}, _tree())
    assert mixed["monophyletic"] is False
    assert mixed["intruders"] == ["Sp2", "Sp4"]


def test_clade_rank_label_single_species_and_no_shared_rank():
    assert clade_rank_label({"Sp1"}, {}) == "Sp1 (species)"
    assert clade_rank_label({"Sp1", "Sp4"}, TAXONOMY) == ""  # crosses order
    assert clade_rank_label({"Sp1", "Zzz"}, TAXONOMY) == ""  # incomplete


# ---------------------------------------------------------------------------
# structure_coherence (foldseek all-vs-all)
# ---------------------------------------------------------------------------

def _pairs(*entries):
    scores = {}
    for a, b, s in entries:
        scores[frozenset((a, b))] = s
    return scores


def test_structure_coherence_within_vs_between():
    scores = _pairs(
        ("Sp1_a1", "Sp2_a2", 900.0), ("Sp1_a1", "Sp3_a3", 880.0),
        ("Sp2_a2", "Sp3_a3", 910.0),                       # within A
        ("Sp1_b1", "Sp2_b2", 950.0),                       # within B
        ("Sp1_a1", "Sp1_b1", 400.0), ("Sp2_a2", "Sp2_b2", 420.0),  # A-B
    )
    groups = {"SF_A": ["Sp1_a1", "Sp2_a2", "Sp3_a3"],
              "SF_B": ["Sp1_b1", "Sp2_b2"]}
    rows = structure_coherence(groups, scores)
    a = next(r for r in rows if r["subfamily"] == "SF_A")
    assert a["mean_within"] == pytest.approx((900 + 880 + 910) / 3)
    assert a["mean_between"] == pytest.approx(410.0)
    assert a["coherent"] is True


def test_structure_coherence_missing_pairs_reported_not_invented():
    scores = _pairs(("Sp1_a1", "Sp2_a2", 900.0))
    groups = {"SF_A": ["Sp1_a1", "Sp2_a2", "Sp3_a3"], "SF_B": ["Sp1_b1"]}
    a = structure_coherence(groups, scores)[0]
    assert a["n_within_pairs"] == 1      # only the observed pair counts
    assert a["mean_between"] is None     # no observed A-B pair -> None, not 0
    assert a["coherent"] is None


def test_parse_pairwise_scores_keeps_best_and_symmetric(tmp_path):
    tsv = tmp_path / "aln.tsv"
    tsv.write_text(
        "Sp1_a1\tSp2_a2\t1e-50\t800\t0.9\n"
        "Sp2_a2\tSp1_a1\t1e-50\t850\t0.92\n"   # reverse direction, better
        "Sp1_a1\tSp1_a1\t0\t2000\t1.0\n"       # self hit ignored
    )
    scores = parse_pairwise_scores(tsv)
    assert scores == {frozenset(("Sp1_a1", "Sp2_a2")): 850.0}


# ---------------------------------------------------------------------------
# coverage suppression reporting (issue #42, found via #40)
# ---------------------------------------------------------------------------

def test_reports_columns_where_a_group_is_suppressed_by_the_coverage_filter():
    # `big` covers column 3 in only 1 of 4 members (0.25), below min_cover, so
    # sdp_scan skips that column for `big` and says nothing about it. On the
    # PEPC clan this made a whole N-terminal region look signal-free when it
    # was simply filtered out.
    from steps.subfamily import coverage_suppressed
    aln = {
        "A_1": "MKW", "A_2": "MK-", "A_3": "MK-", "A_4": "MK-",
        "B_1": "MQW", "B_2": "MQW",
    }
    groups = {"big": ["A_1", "A_2", "A_3", "A_4"], "small": ["B_1", "B_2"]}

    report = coverage_suppressed(aln, groups, min_cover=0.7)

    assert report["big"]["n_suppressed"] == 1
    assert report["big"]["columns"] == [3]
    assert report["small"]["n_suppressed"] == 0


def test_suppression_report_records_the_coverage_that_caused_it():
    from steps.subfamily import coverage_suppressed
    aln = {"A_1": "MW", "A_2": "M-", "B_1": "MW", "B_2": "MW"}
    groups = {"a": ["A_1", "A_2"], "b": ["B_1", "B_2"]}

    report = coverage_suppressed(aln, groups, min_cover=0.7)

    assert report["a"]["coverage"][3 - 1 - 1] == pytest.approx(0.5)
    assert report["a"]["min_coverage"] == pytest.approx(0.5)


def test_a_group_covering_everything_reports_no_suppression():
    from steps.subfamily import coverage_suppressed
    aln = {"A_1": "MKW", "A_2": "MKW", "B_1": "MQW", "B_2": "MQW"}
    groups = {"a": ["A_1", "A_2"], "b": ["B_1", "B_2"]}

    report = coverage_suppressed(aln, groups, min_cover=0.7)

    assert all(r["n_suppressed"] == 0 for r in report.values())


def test_groups_below_min_group_are_still_reported_as_skipped_entirely():
    from steps.subfamily import coverage_suppressed
    aln = {"A_1": "MK", "A_2": "MK", "B_1": "MQ"}
    groups = {"a": ["A_1", "A_2"], "tiny": ["B_1"]}

    report = coverage_suppressed(aln, groups, min_cover=0.7, min_group=2)

    assert report["tiny"]["skipped_too_small"] is True
    assert report["a"]["skipped_too_small"] is False


# ---------------------------------------------------------------------------
# SDP vs the reference-free core (issue #42 A/C)
# ---------------------------------------------------------------------------

def _core_fixture():
    """10 universally invariant columns, then 10 group-diagnostic ones."""
    core = "M" * 10
    a = {f"a{i}": core + "K" * 10 for i in range(5)}
    b = {f"b{i}": core + "Q" * 10 for i in range(5)}
    aln = {**a, **b}
    groups = {"a": sorted(a), "b": sorted(b)}
    return aln, groups


def test_diagnostic_columns_can_never_be_invariant_columns():
    from steps.subfamily import sdp_core_relationship
    aln, groups = _core_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=5)

    assert rel["overlap"] == 0


def test_zero_overlap_is_flagged_as_forced_so_it_is_not_read_as_evidence():
    # An SDP means one group differs; an invariant column means none do. The
    # two sets are disjoint by definition, so a zero overlap says nothing
    # about biology and must not be reported as though it did.
    from steps.subfamily import sdp_core_relationship
    aln, groups = _core_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=5)

    assert rel["overlap_is_by_construction"] is True


def test_diagnostic_columns_clustered_away_from_the_core_are_called_avoiding():
    from steps.subfamily import sdp_core_relationship
    aln, groups = _core_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=5, seed=0)

    assert rel["verdict"] == "avoids_core"
    assert rel["observed_median_distance"] > rel["null_median_distance"]


def test_no_invariant_column_means_no_interpretation_rather_than_a_verdict():
    from steps.subfamily import sdp_core_relationship
    # every column is group-diagnostic; there is no agreed core at all
    aln = {**{f"a{i}": "KKKK" for i in range(5)},
           **{f"b{i}": "QQQQ" for i in range(5)}}
    groups = {"a": [f"a{i}" for i in range(5)], "b": [f"b{i}" for i in range(5)]}

    rel = sdp_core_relationship(aln, groups, min_group=5)

    assert rel["verdict"] == "no_interpretation_available"
    assert "invariant" in rel["reason"]


def test_no_diagnostic_column_means_no_interpretation_rather_than_a_verdict():
    from steps.subfamily import sdp_core_relationship
    aln = {f"s{i}": "MMMM" for i in range(10)}
    groups = {"a": [f"s{i}" for i in range(5)], "b": [f"s{i}" for i in range(5, 10)]}

    rel = sdp_core_relationship(aln, groups, min_group=5)

    assert rel["verdict"] == "no_interpretation_available"
    assert "diagnostic" in rel["reason"]


def test_groups_too_small_to_scan_are_reported_not_silently_dropped():
    from steps.subfamily import sdp_core_relationship
    aln, groups = _core_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=99)

    assert rel["verdict"] == "no_interpretation_available"
    assert rel["skipped_groups"] == ["a", "b"]


def test_the_null_is_seeded_so_the_verdict_is_reproducible():
    from steps.subfamily import sdp_core_relationship
    aln, groups = _core_fixture()

    first = sdp_core_relationship(aln, groups, min_group=5, seed=7)
    second = sdp_core_relationship(aln, groups, min_group=5, seed=7)

    assert first["null_median_distance"] == second["null_median_distance"]
    assert first["p_value"] == second["p_value"]
