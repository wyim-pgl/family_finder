"""Tests for steps/family_naming.py — annotate families from ONE species.

Given one annotated species (gene_id + description, e.g. the iceplant/Mcry
functional annotation), map each annotated gene to its pipeline family
(and subfamily when a groups table is supplied) and name each family by
its members' descriptions. Pure dict-in/dict-out (house convention).
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from steps.family_naming import (  # noqa: E402
    index_families,
    map_annotations,
    match_gene,
    name_groups,
)

FAMILIES = {
    "R1_OG0000001": {"Mcry_g1", "Mcry_g2", "Ococ_x1"},
    "R1_OG0000002": {"Mcry_g3.t1", "Ococ_x2"},
}
GROUPS = {"Mcry_g1": "SF_A", "Mcry_g2": "SF_B", "Ococ_x1": "SF_A"}


def test_index_and_exact_match():
    index = index_families(FAMILIES)
    assert match_gene("Mcry_g1", index) == ("Mcry_g1", "R1_OG0000001")


def test_match_adds_species_prefix_when_missing():
    # Annotation tables often carry bare ids (g1) without the pipeline's
    # SpeciesPrefix_ convention.
    index = index_families(FAMILIES)
    assert match_gene("g1", index, species="Mcry") == ("Mcry_g1", "R1_OG0000001")


def test_match_tries_isoform_suffix_variants():
    index = index_families(FAMILIES)
    # Annotation says g3, pipeline stored g3.t1
    assert match_gene("g3", index, species="Mcry") == ("Mcry_g3.t1", "R1_OG0000002")
    # And the reverse: annotation has the isoform, pipeline the locus.
    assert match_gene("Mcry_g1.t1", index) == ("Mcry_g1", "R1_OG0000001")


def test_match_returns_none_for_unknown():
    assert match_gene("nope", index_families(FAMILIES), species="Mcry") is None


def test_map_annotations_rows_and_unmatched():
    annotations = {"g1": "PEPC 1", "g2": "PEPC 2", "zzz": "mystery"}
    rows, unmatched = map_annotations(
        annotations, FAMILIES, groups=GROUPS, species="Mcry"
    )
    by_gene = {r["annotation_gene"]: r for r in rows}
    assert by_gene["g1"]["family_id"] == "R1_OG0000001"
    assert by_gene["g1"]["subfamily"] == "SF_A"
    assert by_gene["g2"]["subfamily"] == "SF_B"
    assert unmatched == ["zzz"]


def test_name_groups_majority_description():
    rows = [
        {"family_id": "F1", "subfamily": "", "description": "PEPC 1"},
        {"family_id": "F1", "subfamily": "", "description": "PEPC 1"},
        {"family_id": "F1", "subfamily": "", "description": "kinase"},
        {"family_id": "F2", "subfamily": "", "description": "kinase"},
    ]
    named = {r["group_id"]: r for r in name_groups(rows, key="family_id")}
    assert named["F1"]["name"] == "PEPC 1"
    assert named["F1"]["support"] == pytest.approx(2 / 3)
    assert named["F1"]["n_annotated"] == 3
    assert named["F2"]["name"] == "kinase"


def test_name_groups_by_subfamily_skips_blank():
    rows = [
        {"family_id": "F1", "subfamily": "SF_A", "description": "PEPC 1"},
        {"family_id": "F1", "subfamily": "", "description": "PEPC 2"},
    ]
    named = name_groups(rows, key="subfamily")
    assert len(named) == 1 and named[0]["group_id"] == "SF_A"


def test_ambiguous_suffix_match_is_unmatched_not_first_hit(caplog):
    # Both Mcry_g9.t1 and Mcry_g9.1 exist, in DIFFERENT families; a curated
    # table saying just g9 must not be resolved by suffix trial order.
    families = {
        "R1_OG0000010": {"Mcry_g9.t1"},
        "R1_OG0000011": {"Mcry_g9.1"},
    }
    index = index_families(families)

    import logging
    with caplog.at_level(logging.WARNING, logger="family_finder"):
        assert match_gene("g9", index, species="Mcry") is None
    assert any("ambiguous" in rec.message for rec in caplog.records)


def test_exact_match_beats_suffix_ambiguity():
    families = {
        "R1_OG0000010": {"Mcry_g9"},
        "R1_OG0000011": {"Mcry_g9.1"},
    }
    index = index_families(families)
    assert match_gene("g9", index, species="Mcry") == ("Mcry_g9", "R1_OG0000010")


def test_combine_sources_surfaces_contradicting_orthology():
    from steps.family_naming import combine_sources

    direct = [{"pipeline_gene": "Mcry_g1", "description": "PEPC",
               "source": "direct"}]
    ortho_agree = {"pipeline_gene": "Mcry_g1", "description": "PEPC",
                   "source": "orthology"}
    ortho_conflict = {"pipeline_gene": "Mcry_g1", "description": "PTPC",
                      "source": "orthology"}
    ortho_other = {"pipeline_gene": "Mcry_g2", "description": "MYB",
                   "source": "orthology"}

    rows, conflicts = combine_sources(
        direct, [ortho_agree, ortho_conflict, ortho_other])

    assert ortho_other in rows                       # unsuppressed passes through
    assert ortho_agree not in rows                   # echo is dropped quietly
    assert len(conflicts) == 1                       # contradiction is surfaced
    assert conflicts[0]["description"] == "PTPC"
    assert conflicts[0]["direct_description"] == "PEPC"
