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
    load_taxonomy,
    parse_pairwise_scores,
    sdp_scan,
    structure_coherence,
    taxonomic_composition,
)


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
