"""One canonical form for gene identifiers (issue #42).

Three separate silent coverage losses in one day traced to the same cause:
structure filenames replace '.' with '_', transcript suffixes appear on one
side and not the other, and every module carries its own normalisation. Each
mismatch drops genes without an error — historically 111 structures matched 29
alignment names, and a region_disorder rerun lost 16 of 111 the same way.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.gene_ids import canon_gene_id, match_ids, strip_isoform


# ------------------------------------------------------------- canonical ----

def test_a_pdb_extension_is_not_part_of_the_identifier():
    assert canon_gene_id("Mcry_Mcr8G11630.pdb") == "Mcry_Mcr8G11630"


def test_dots_become_underscores_because_structure_filenames_do_that():
    assert canon_gene_id("Obas__JBFLFP010000003.1_000519") == \
        "Obas__JBFLFP010000003_1_000519"


def test_canonical_form_is_idempotent():
    once = canon_gene_id("Ococ_OcoChr01G12850.t1")

    assert canon_gene_id(once) == once


def test_surrounding_whitespace_is_not_an_identifier():
    assert canon_gene_id("  Mcry_Mcr8G11630 ") == "Mcry_Mcr8G11630"


# --------------------------------------------------------------- isoform ----

def test_transcript_suffixes_are_stripped_only_when_asked():
    assert strip_isoform("Ococ_OcoChr01G12850_t1") == "Ococ_OcoChr01G12850"
    assert strip_isoform("Ococ_OcoChr01G12850.t1") == "Ococ_OcoChr01G12850"


def test_a_numeric_suffix_that_is_part_of_the_locus_survives():
    # G09070 is the locus, not an isoform — only _t<N> / .<N> style suffixes go
    assert strip_isoform("Ococ_OcoChr10G09070") == "Ococ_OcoChr10G09070"


def test_stripping_a_trailing_dot_number_isoform():
    assert strip_isoform("Ahyp_AH008008.v2.1") == "Ahyp_AH008008.v2"


# ----------------------------------------------------------------- match ----

def test_exact_matches_are_preferred_and_reported_as_such():
    result = match_ids(["a", "b"], ["a", "b", "c"])

    assert result.mapping == {"a": "a", "b": "b"}
    assert result.level == "exact"
    assert result.unmatched == []


def test_falls_back_to_canonical_form_when_exact_fails():
    result = match_ids(["Obas_X.1_9.pdb"], ["Obas_X_1_9"])

    assert result.mapping == {"Obas_X.1_9.pdb": "Obas_X_1_9"}
    assert result.level == "canonical"


def test_falls_back_to_isoform_stripping_last():
    result = match_ids(["Ococ_G1.t1"], ["Ococ_G1"])

    assert result.mapping == {"Ococ_G1.t1": "Ococ_G1"}
    assert result.level == "isoform"


def test_unmatched_queries_are_named_not_counted_away():
    result = match_ids(["a", "ghost"], ["a"])

    assert result.unmatched == ["ghost"]
    assert result.unmatched_fraction == pytest.approx(0.5)


def test_matching_raises_once_too_much_of_the_input_is_lost():
    # 16 of 111 vanished silently in the real run; a ceiling makes that loud
    with pytest.raises(ValueError, match="unmatched"):
        match_ids(["a", "g1", "g2"], ["a"], max_unmatched=0.1)


def test_the_ceiling_is_not_applied_unless_requested():
    result = match_ids(["a", "g1", "g2"], ["a"])

    assert result.unmatched_fraction == pytest.approx(2 / 3)


def test_an_empty_query_set_is_not_a_division_by_zero():
    result = match_ids([], ["a"])

    assert result.unmatched_fraction == 0.0
    assert result.mapping == {}


def test_a_canonical_collision_is_an_error_not_a_silent_overwrite():
    # Two distinct references collapsing to one canonical form make the mapping
    # ambiguous. The query must actually need that level: an exact hit is used
    # as-is, so a collision further down is irrelevant and not raised.
    with pytest.raises(ValueError, match="collide"):
        match_ids(["X.1.pdb"], ["X.1", "X_1"])


def test_a_collision_at_an_unused_level_does_not_block_an_exact_match():
    result = match_ids(["X.1"], ["X.1", "X_1"])

    assert result.mapping == {"X.1": "X.1"}
    assert result.level == "exact"


def test_an_underscore_gene_number_is_not_an_isoform_suffix():
    # Cgig_..._SGP5p_1338_000001 and _000002 are distinct loci; stripping the
    # trailing _NNNNNN merged them and the collision guard caught it on real
    # structure filenames.
    assert strip_isoform("Cgig_v2_SGP5p_1338_000001") == "Cgig_v2_SGP5p_1338_000001"
    assert strip_isoform("Cgig_v2_SGP5p_1338_000002") == "Cgig_v2_SGP5p_1338_000002"


def test_neighbouring_gene_numbers_do_not_collide_when_matching():
    result = match_ids(["Cgig_1338_000001.pdb"],
                       ["Cgig_1338_000001", "Cgig_1338_000002"])

    assert result.mapping == {"Cgig_1338_000001.pdb": "Cgig_1338_000001"}
