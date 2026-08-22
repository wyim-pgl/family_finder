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


# ---------------------------------------------------------------------------
# group id resolution (issue #42)
#
# sdp_scan used to keep only the members already spelled exactly as the
# alignment spells them (`[m for m in members if m in alignment]`). Everything
# else vanished without a word, which is how `Ococ_OcoChr10G09070.t1` failed to
# match groups.json (#34) and how a subfamily can be reported as having no
# diagnostic residues when in fact none of its members were ever looked at.
# ---------------------------------------------------------------------------

def test_group_members_with_a_transcript_suffix_still_reach_the_alignment():
    from steps.subfamily import resolve_groups

    aln = {"Ococ_OcoChr10G09070": "KAAWTR", "Ococ_OcoChr10G09071": "SAAWTR"}
    groups = {"SF": ["Ococ_OcoChr10G09070.t1", "Ococ_OcoChr10G09071.t1"]}

    resolved, report = resolve_groups(aln, groups)

    assert resolved["SF"] == ["Ococ_OcoChr10G09070", "Ococ_OcoChr10G09071"]
    assert report["groups"]["SF"]["unmatched"] == []
    assert report["level"] == "isoform"


def test_members_missing_from_the_alignment_are_named_not_dropped():
    from steps.subfamily import resolve_groups

    resolved, report = resolve_groups(ALN, {"SF_A": GROUPS["SF_A"] + ["ghost_g1"]})

    assert "ghost_g1" not in resolved["SF_A"]
    assert report["groups"]["SF_A"]["unmatched"] == ["ghost_g1"]
    assert report["n_unmatched"] == 1
    assert report["unmatched_fraction"] > 0


def test_a_cap_turns_a_silent_coverage_loss_into_a_failure():
    from steps.subfamily import resolve_groups

    with pytest.raises(ValueError, match="unmatched"):
        resolve_groups(ALN, {"SF_A": ["ghost1", "ghost2"]}, max_unmatched=0.1)


def test_an_underscore_gene_number_is_never_read_as_an_isoform_here_either():
    # Cgig_..._1338_000001 / _000002 are distinct loci. Loosening far enough to
    # merge them would be worse than leaving one unmatched, so the collision
    # guard in match_ids must still fire through this path.
    from steps.subfamily import resolve_groups

    aln = {"Cgig_1338_000001": "KAAWTR", "Cgig_1338_000002": "SAAWTR"}
    resolved, report = resolve_groups(aln, {"SF": ["Cgig_1338_000001.pdb",
                                                   "Cgig_1338_000002.pdb"]})

    assert resolved["SF"] == ["Cgig_1338_000001", "Cgig_1338_000002"]


def test_sdp_scan_matches_group_ids_through_the_shared_normaliser():
    aln = {k.replace("_", "_") : v for k, v in ALN.items()}
    # the groups arrive with transcript suffixes the alignment does not carry
    groups = {name: [m + ".t1" for m in members]
              for name, members in GROUPS.items()}

    rows = sdp_scan(aln, groups, min_group=4)

    assert {r["subfamily"] for r in rows} == {"SF_A", "SF_B"}


def test_sdp_scan_can_be_told_to_refuse_rather_than_scan_a_shrunken_group():
    groups = {"SF_A": ["ghost1", "ghost2", "ghost3", "ghost4"]}

    with pytest.raises(ValueError, match="unmatched"):
        sdp_scan(ALN, groups, min_group=4, max_unmatched=0.1)


def test_coverage_suppressed_uses_the_same_matcher_as_the_scan():
    from steps.subfamily import coverage_suppressed

    groups = {name: [m + ".t1" for m in members]
              for name, members in GROUPS.items()}

    report = coverage_suppressed(ALN, groups, min_cover=0.7, min_group=4)

    assert report["SF_A"]["n_members"] == 4
    assert report["SF_A"]["skipped_too_small"] is False


# ---------------------------------------------------------------------------
# sequence-identity control for structural coherence (issue #39 item 3)
#
# foldseek bits scales with sequence similarity, so "within > between" may be
# nothing but "closer relatives look more alike". The measured PEPC result (5
# of 6 subfamilies coherent, ratio median 1.24) is not citable until identity
# is either regressed out or the number carries a warning saying it was not.
# ---------------------------------------------------------------------------

def _ident(*entries):
    return {frozenset((a, b)): f for a, b, f in entries}


def test_an_uncontrolled_coherence_row_says_so():
    scores = _pairs(("Sp1_a1", "Sp2_a2", 900.0), ("Sp1_a1", "Sp1_b1", 400.0))
    groups = {"SF_A": ["Sp1_a1", "Sp2_a2"], "SF_B": ["Sp1_b1"]}

    row = structure_coherence(groups, scores)[0]

    assert row["sequence_controlled"] is False
    assert "sequence" in row["warning"].lower()
    assert row["coherent_controlled"] is None


def test_the_metric_the_number_came_from_is_recorded():
    scores = _pairs(("Sp1_a1", "Sp2_a2", 900.0), ("Sp1_a1", "Sp1_b1", 400.0))
    groups = {"SF_A": ["Sp1_a1", "Sp2_a2"], "SF_B": ["Sp1_b1"]}

    row = structure_coherence(groups, scores, metric="alntmscore")[0]

    assert row["metric"] == "alntmscore"


def _overlapping(within_lift=0.0, n=30, slope=1000.0):
    """within- and between-pairs spanning the SAME identity range.

    The comparison only means anything where both kinds of pair exist, so
    every fixture that expects a verdict has to build that overlap
    deliberately — the real PEPC data does not have it (see below).
    """
    ident, scores = {}, {}
    groups = {"SF_A": [], "SF_B": []}
    for i in range(n):
        f = 0.30 + 0.02 * i
        wp = frozenset((f"a{i}", f"a{i + 100}"))
        bp = frozenset((f"a{i}", f"b{i}"))
        ident[wp] = ident[bp] = f
        scores[wp] = slope * f + within_lift
        scores[bp] = slope * f
        groups["SF_A"] += [f"a{i}", f"a{i + 100}"]
        groups["SF_B"].append(f"b{i}")
    return groups, scores, ident


def test_identity_driven_coherence_disappears_once_identity_is_regressed_out():
    # every pair sits exactly on score = 1000 * fident, so the structural
    # signal is PURELY sequence identity and no residual separates the groups
    groups, scores, ident = _overlapping(within_lift=0.0)

    row = next(r for r in structure_coherence(groups, scores, identities=ident)
               if r["subfamily"] == "SF_A")

    assert row["sequence_controlled"] is True
    assert row["linear_fit_adequate"] is True
    assert row["mean_within_residual"] == pytest.approx(0.0, abs=1e-6)
    assert row["mean_between_residual"] == pytest.approx(0.0, abs=1e-6)
    assert row["verdict"] == "not_coherent"
    assert row["coherent_controlled"] is False  # nothing beyond sequence


def test_coherence_beyond_sequence_identity_survives_the_control():
    # same identities, but the within-A pairs score ABOVE the identity line
    groups, scores, ident = _overlapping(within_lift=300.0)

    row = next(r for r in structure_coherence(groups, scores, identities=ident)
               if r["subfamily"] == "SF_A")

    assert row["sequence_controlled"] is True
    assert row["mean_within_residual"] > row["mean_between_residual"]
    assert row["verdict"] == "coherent"
    assert row["coherent_controlled"] is True


def test_identities_too_few_to_regress_fall_back_to_a_warning():
    ident = _ident(("Sp1_a1", "Sp2_a2", 0.9), ("Sp1_a1", "Sp1_b1", 0.4))
    scores = {pair: 1000.0 * f for pair, f in ident.items()}
    groups = {"SF_A": ["Sp1_a1", "Sp2_a2"], "SF_B": ["Sp1_b1"]}

    row = structure_coherence(groups, scores, identities=ident)[0]

    assert row["sequence_controlled"] is False
    assert "regress" in row["warning"].lower() or "identit" in row["warning"].lower()


def test_constant_identity_cannot_be_regressed_out_and_says_so():
    ident = _ident(("Sp1_a1", "Sp2_a2", 0.5), ("Sp1_a1", "Sp3_a3", 0.5),
                   ("Sp2_a2", "Sp3_a3", 0.5), ("Sp1_a1", "Sp1_b1", 0.5))
    scores = _pairs(("Sp1_a1", "Sp2_a2", 900.0), ("Sp1_a1", "Sp3_a3", 880.0),
                    ("Sp2_a2", "Sp3_a3", 910.0), ("Sp1_a1", "Sp1_b1", 400.0))
    groups = {"SF_A": ["Sp1_a1", "Sp2_a2", "Sp3_a3"], "SF_B": ["Sp1_b1"]}

    row = structure_coherence(groups, scores, identities=ident)[0]

    assert row["sequence_controlled"] is False


def test_sequence_confounding_reports_how_strong_the_coupling_is():
    from steps.subfamily import sequence_confounding

    ident = _ident(("a", "b", 0.9), ("a", "c", 0.5), ("b", "c", 0.3),
                   ("a", "d", 0.7))
    scores = {pair: 1000.0 * f for pair, f in ident.items()}

    report = sequence_confounding(scores, ident)

    assert report["n_pairs"] == 4
    assert report["r"] == pytest.approx(1.0)
    assert report["slope"] == pytest.approx(1000.0)


def test_sequence_confounding_on_unrelated_numbers_is_near_zero():
    from steps.subfamily import sequence_confounding

    ident = _ident(("a", "b", 0.1), ("a", "c", 0.9), ("b", "c", 0.1),
                   ("a", "d", 0.9))
    scores = _pairs(("a", "b", 500.0), ("a", "c", 500.0),
                    ("b", "c", 700.0), ("a", "d", 700.0))

    report = sequence_confounding(scores, ident)

    assert report["r"] == pytest.approx(0.0, abs=1e-9)


# --- reading fident out of the foldseek table --------------------------------

_WIDE = (
    "query,target,fident,alnlen,evalue,bits,qtmscore,alntmscore"
)


def test_pairwise_scores_can_be_read_on_a_named_metric(tmp_path):
    p = tmp_path / "pairs.tsv"
    p.write_text("g1\tg2\t0.90\t900\t0.0\t800\t0.95\t0.93\n")

    by_bits = parse_pairwise_scores(p, metric="bits", columns=_WIDE)
    by_tm = parse_pairwise_scores(p, metric="alntmscore", columns=_WIDE)

    assert by_bits[frozenset(("g1", "g2"))] == 800.0
    assert by_tm[frozenset(("g1", "g2"))] == pytest.approx(0.93)


def test_identities_come_from_the_same_alignment_as_the_score(tmp_path):
    from steps.subfamily import parse_pair_identities

    p = tmp_path / "pairs.tsv"
    p.write_text("g1\tg2\t0.90\t900\t0.0\t800\t0.95\t0.93\n"
                 "g2\tg1\t0.40\t100\t0.0\t100\t0.50\t0.40\n")

    ident = parse_pair_identities(p, columns=_WIDE)

    # the better-scoring direction wins, and its OWN fident is the one kept
    assert ident[frozenset(("g1", "g2"))] == pytest.approx(0.90)


def test_a_table_that_does_not_match_the_declared_columns_is_refused(tmp_path):
    p = tmp_path / "pairs.tsv"
    p.write_text("g1\tg2\t0.0\t800\t0.93\n")   # the 5-column layout

    with pytest.raises(ValueError, match="column"):
        parse_pairwise_scores(p, columns=_WIDE)


def test_asking_for_a_column_the_table_does_not_have_is_refused(tmp_path):
    from steps.subfamily import parse_pair_identities

    p = tmp_path / "pairs.tsv"
    p.write_text("g1\tg2\t0.0\t800\t0.93\n")

    with pytest.raises(ValueError, match="fident"):
        parse_pair_identities(p)


def test_coherence_matches_group_ids_against_the_pair_table_ids():
    # groups arrive in alignment spelling, the pair table in PDB-filename
    # spelling. Comparing them literally finds no pairs at all and reports
    # every subfamily as having none observed — silently.
    scores = _pairs(("Obas_X_1_9", "Obas_Y_1_9", 900.0),
                    ("Obas_X_1_9", "Obas_Z_1_9", 400.0))
    groups = {"SF_A": ["Obas_X.1_9", "Obas_Y.1_9"], "SF_B": ["Obas_Z.1_9"]}

    row = structure_coherence(groups, scores)[0]

    assert row["n_within_pairs"] == 1
    assert row["n_between_pairs"] == 1


def test_group_members_with_no_structure_are_counted_per_subfamily():
    scores = _pairs(("Sp1_a1", "Sp2_a2", 900.0), ("Sp1_a1", "Sp1_b1", 400.0))
    groups = {"SF_A": ["Sp1_a1", "Sp2_a2", "ghost1"], "SF_B": ["Sp1_b1"]}

    row = structure_coherence(groups, scores)[0]

    assert row["n_members_without_structure"] == 1


# ---------------------------------------------------------------------------
# taxonomic_composition, branch by branch (issue #39 item 1)
#
# The function was refactored from 85 lines / nesting 5 into three helpers.
# The refactor was proved output-identical first (sha256 576796d8... over 14
# synthetic cases covering all nine verdicts); these tests keep each branch
# pinned from here on, so a later change to one path cannot quietly move
# another.
# ---------------------------------------------------------------------------

# ((A,B),(C,D)),(E,F) — lets a species set be a clade, span the root, or
# interleave with another subfamily's species.
TC_TREE = "(((A,B)ab,(C,D)cd)abcd,(E,F)ef)root;"

TC_TAX = {
    "A": {"genus": "Ga", "family": "Fa", "order": "Oa"},
    "B": {"genus": "Ga", "family": "Fa", "order": "Oa"},
    "C": {"genus": "Gc", "family": "Fa", "order": "Oa"},
    "D": {"genus": "Gd", "family": "Fd", "order": "Oa"},
    "E": {"genus": "Ge", "family": "Fe", "order": "Oe"},
    "F": {"genus": "Ge", "family": "Fe", "order": "Oe"},
}


def _tc(groups, **kw):
    if "species_tree" in kw:
        kw["species_tree"] = parse_newick(kw["species_tree"])
    return {r["subfamily"]: r for r in taxonomic_composition(groups, **kw)}


def test_a_clade_strictly_inside_the_family_span_is_lineage_specific():
    rows = _tc({"sf_ab": ["A_g1", "B_g1"], "sf_ef": ["E_g1", "F_g1"],
                "sf_cd": ["C_g1", "D_g1"]}, species_tree=TC_TREE)

    assert rows["sf_ab"]["verdict"] == "lineage-specific (clade)"
    assert rows["sf_ab"]["monophyletic"] is True


def test_a_set_whose_mrca_is_the_family_span_is_a_paralog_split():
    # ancient paralogs kept in every sampled species are trivially
    # monophyletic — the root-span rule is what stops that reading as
    # lineage-specific. Here the family spans only A+B, and sf_all covers
    # exactly that span while sf_a sits strictly inside it.
    rows = _tc({"sf_all": ["A_g1", "B_g1"], "sf_a": ["A_g2"]},
               species_tree=TC_TREE)

    assert rows["sf_all"]["monophyletic"] is True
    assert rows["sf_all"]["verdict"] == "paralog-split (spans family root)"
    assert rows["sf_a"]["verdict"] == "lineage-specific (clade)"


def test_an_interleaved_species_set_is_a_paralog_split_and_names_the_intruders():
    rows = _tc({"sf_ac": ["A_g1", "C_g1"], "sf_ef": ["E_g1", "F_g1"]},
               species_tree=TC_TREE)

    assert rows["sf_ac"]["verdict"] == "paralog-split (non-monophyletic)"
    assert "interleaved species: B,D" in rows["sf_ac"]["notes"]


def test_species_absent_from_the_tree_are_named_not_ignored():
    rows = _tc({"sf_mixed": ["A_g1", "Z_g1"], "sf_offtree": ["Y_g1", "Z_g2"]},
               species_tree=TC_TREE)

    assert "not in species tree: Z" in rows["sf_mixed"]["notes"]
    assert rows["sf_offtree"]["verdict"] == "unknown (no species in species tree)"


def test_taxonomy_only_labels_the_tree_verdict_it_does_not_decide_it():
    # sf_ef widens the family span, so sf_ab sits strictly inside it and the
    # root-span rule does not fire — the verdict under test is the clade one
    rows = _tc({"sf_ab": ["A_g1", "B_g1"], "sf_odd": ["A_g2", "Z_g1"],
                "sf_ef": ["E_g1"]},
               taxonomy=TC_TAX, species_tree=TC_TREE)

    assert rows["sf_ab"]["clade_label"] == "Ga (genus)"
    assert rows["sf_ab"]["genus_purity"] == 1.0        # label layer present
    assert rows["sf_ab"]["verdict"] == "lineage-specific (clade)"
    assert "unknown taxonomy for: Z" in rows["sf_odd"]["notes"]


def test_without_a_tree_the_lowest_pure_rank_decides():
    assert _tc({"s": ["A_g1", "A_g2"]}, taxonomy=TC_TAX)["s"]["verdict"] == \
        "lineage-specific (species)"
    assert _tc({"s": ["A_g1", "B_g1"]}, taxonomy=TC_TAX)["s"]["verdict"] == \
        "lineage-specific (genus)"
    assert _tc({"s": ["A_g1", "B_g1", "C_g1"]}, taxonomy=TC_TAX)["s"]["verdict"] == \
        "lineage-specific (family)"
    assert _tc({"s": ["A_g1", "B_g1", "C_g1", "D_g1"]},
               taxonomy=TC_TAX)["s"]["verdict"] == "lineage-specific (order)"


def test_without_a_tree_a_set_crossing_order_is_a_paralog_split():
    rows = _tc({"sf_ae": ["A_g1", "E_g1"]}, taxonomy=TC_TAX)

    assert rows["sf_ae"]["verdict"] == "paralog-split (crosses order)"


def test_a_species_missing_from_the_taxonomy_still_counts_as_itself():
    # identity comes from the gene-id prefix, so the row must not vanish
    rows = _tc({"sf_z": ["Z_g1", "A_g1"]}, taxonomy=TC_TAX)

    assert rows["sf_z"]["n_species"] == 2
    assert "unknown taxonomy for: Z" in rows["sf_z"]["notes"]


def test_the_species_prefix_delimiter_is_configurable():
    rows = _tc({"sf_a": ["A|g1", "A|g2"]}, taxonomy=TC_TAX, delimiter="|")

    assert rows["sf_a"]["species_dominant"] == "A"


def test_no_groups_is_an_empty_table_not_a_crash():
    assert taxonomic_composition({}, species_tree=parse_newick(TC_TREE)) == []


def test_the_tree_path_and_the_legacy_path_emit_different_columns():
    tree_row = _tc({"s": ["A_g1", "B_g1"]}, species_tree=TC_TREE)["s"]
    legacy_row = _tc({"s": ["A_g1", "B_g1"]}, taxonomy=TC_TAX)["s"]

    assert list(tree_row) == ["subfamily", "n_members", "n_species",
                             "n_in_tree", "monophyletic", "mrca_name",
                             "mrca_depth", "clade_label", "verdict", "notes"]
    assert "monophyletic" not in legacy_row
    assert legacy_row["order_dominant"] == "Oa"


# ---------------------------------------------------------------------------
# the identity control has to prove it is SUPPORTED before it may conclude
#
# Measured on the real PEPC pair table (5,343 pairs): the linear fit reaches
# only r2 = 0.430 (bits) / 0.497 (alntmscore), and the binned residual means
# run +0.006 / -0.061 / -0.068 / +0.074 / ... — systematically non-zero.
# alntmscore saturates at 1.0 while fident keeps climbing, so high-identity
# pairs sit BELOW the line whatever their structure. Within-pairs are almost
# all high-identity, so mean_within_residual < mean_between_residual comes out
# of the misspecification, not out of the structures. Reporting that as
# "not coherent" is the artifact-as-conclusion failure this issue is about.
# ---------------------------------------------------------------------------

def _saturating(n=60):
    """score = min(1.0, 1.6 * fident) — the real alntmscore/fident shape."""
    scores, ident = {}, {}
    for i in range(n):
        f = 0.30 + 0.011 * i
        pair = frozenset((f"q{i}", f"t{i}"))
        ident[pair] = f
        scores[pair] = min(1.0, 1.6 * f)
    return scores, ident


def test_a_saturating_relationship_is_reported_as_uninterpretable_not_negative():
    from steps.subfamily import residual_diagnostics

    scores, ident = _saturating()
    diag = residual_diagnostics(scores, ident)

    assert diag["linear_fit_adequate"] is False
    assert diag["max_abs_z"] > 3.0
    assert "non-monotonic" in diag["reason"] or "systematic" in diag["reason"]


def test_a_genuinely_linear_relationship_passes_the_diagnostic():
    from steps.subfamily import residual_diagnostics

    ident = {frozenset((f"q{i}", f"t{i}")): 0.3 + 0.011 * i for i in range(60)}
    # exactly linear plus an alternating wobble that carries no trend
    scores = {p: 2.0 * f + (0.001 if i % 2 else -0.001)
              for i, (p, f) in enumerate(sorted(ident.items(),
                                                key=lambda kv: kv[1]))}
    diag = residual_diagnostics(scores, ident)

    assert diag["linear_fit_adequate"] is True
    assert diag["r2"] > 0.99


def _saturating_groups():
    """The real PEPC layout: within-pairs bunched at high identity, a
    saturating score, and between-pairs spread lower down."""
    ident, scores = {}, {}
    groups = {"SF_A": [f"a{i}" for i in range(12)],
              "SF_B": [f"b{i}" for i in range(12)]}
    for i in range(11):                                  # within SF_A, high id
        f = 0.80 + 0.018 * i
        ident[frozenset((f"a{i}", f"a{i + 1}"))] = f
    for i in range(11):                                  # within SF_B, high id
        f = 0.82 + 0.016 * i
        ident[frozenset((f"b{i}", f"b{i + 1}"))] = f
    for i in range(12):                                  # between, low id
        for j in (0, 1, 2):
            ident[frozenset((f"a{i}", f"b{(i + j) % 12}"))] = 0.30 + 0.011 * (3 * i + j)
    for pair, f in ident.items():
        scores[pair] = min(1.0, 1.6 * f)                 # saturates at 1.0
    return groups, scores, ident


def test_a_misspecified_fit_cannot_be_used_to_call_a_subfamily_incoherent():
    groups, scores, ident = _saturating_groups()

    row = next(r for r in structure_coherence(groups, scores, identities=ident)
               if r["subfamily"] == "SF_A")

    assert row["linear_fit_adequate"] is False
    assert row["verdict"] == "no_interpretation_available"
    assert row["coherent_controlled"] is None
    assert "r2" in row["reason"]


def test_within_and_between_that_never_share_an_identity_range_cannot_be_compared():
    # OG8 on the real data: zero fident bins hold both a within- and a
    # between-pair, so there is nothing the residual comparison is comparing
    ident, scores = {}, {}
    for i in range(20):                      # A-vs-B pairs: identity 0.30-0.49
        p = frozenset((f"a{i % 7}", f"b{i}"))
        ident[p] = 0.30 + 0.01 * i
    for i in range(6):                       # within A: identity 0.90-0.95
        p = frozenset((f"a{i}", f"a{i + 1}"))
        ident[p] = 0.90 + 0.01 * i
    for p, f in ident.items():               # strictly linear: the fit is fine
        scores[p] = 2.0 * f
    groups = {"SF_A": [f"a{i}" for i in range(7)],
              "SF_B": [f"b{i}" for i in range(20)]}

    row = next(r for r in structure_coherence(groups, scores, identities=ident)
               if r["subfamily"] == "SF_A")

    assert row["linear_fit_adequate"] is True   # the fit is not the problem

    assert row["n_shared_identity_bins"] == 0
    assert row["verdict"] == "no_interpretation_available"
    assert row["coherent_controlled"] is None
    assert "identity range" in row["reason"]


def test_a_supported_control_still_reaches_a_verdict():
    # within and between overlap across the identity range and the fit is
    # linear, so the residual comparison means something
    ident, scores = {}, {}
    members = {"SF_A": [], "SF_B": []}
    for i in range(30):
        f = 0.30 + 0.02 * i
        wp = frozenset((f"a{i}", f"a{i + 100}"))
        bp = frozenset((f"a{i}", f"b{i}"))
        ident[wp], ident[bp] = f, f
        scores[wp], scores[bp] = 2.0 * f + 0.5, 2.0 * f      # within sits above
        members["SF_A"] += [f"a{i}", f"a{i + 100}"]
        members["SF_B"].append(f"b{i}")

    row = next(r for r in structure_coherence(members, scores, identities=ident)
               if r["subfamily"] == "SF_A")

    assert row["linear_fit_adequate"] is True
    assert row["n_shared_identity_bins"] >= 3
    assert row["verdict"] == "coherent"
    assert row["coherent_controlled"] is True


def test_an_uncontrolled_row_reports_no_interpretation_rather_than_a_verdict():
    scores = _pairs(("Sp1_a1", "Sp2_a2", 900.0), ("Sp1_a1", "Sp1_b1", 400.0))
    groups = {"SF_A": ["Sp1_a1", "Sp2_a2"], "SF_B": ["Sp1_b1"]}

    row = structure_coherence(groups, scores)[0]

    assert row["verdict"] == "no_interpretation_available"
    assert row["coherent_controlled"] is None
    assert "UNCONTROLLED" in row["reason"]


def test_the_fit_quality_travels_with_every_row():
    groups, scores, ident = _saturating_groups()

    row = structure_coherence(groups, scores, identities=ident)[0]

    assert row["fit_r2"] is not None
    assert 0.0 <= row["fit_r2"] <= 1.0


def test_sequence_controlled_means_the_control_held_not_that_it_was_attempted():
    # Passing identities= is not the same as the control working. A row whose
    # reason says UNCONTROLLED must not be counted as a controlled result by
    # anyone filtering on this column.
    groups, scores, ident = _saturating_groups()

    rows = structure_coherence(groups, scores, identities=ident)

    for row in rows:
        assert row["verdict"] == "no_interpretation_available"
        assert row["sequence_controlled"] is False
        assert "UNCONTROLLED" in row["reason"]


def test_the_flag_and_the_verdict_never_disagree():
    cases = [_saturating_groups()]
    for lift in (300.0, 0.0):
        g, s, i = _overlapping(within_lift=lift)
        cases.append((g, s, i))
    for groups, scores, ident in cases:
        for row in structure_coherence(groups, scores, identities=ident):
            decided = row["verdict"] in ("coherent", "not_coherent")
            assert row["sequence_controlled"] is decided
            assert (row["coherent_controlled"] is not None) is decided
