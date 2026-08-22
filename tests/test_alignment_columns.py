"""Column selection and occupancy reporting (issue #42).

A global occupancy threshold cannot preserve a region unique to a minority
subfamily: such a region's global occupancy is capped at that subfamily's
share of the alignment. These tests pin the group-aware rule that removes
that ceiling, and the column map that makes trimmed coordinates recoverable.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.alignment import (
    alignment_id,
    apply_columns,
    column_map,
    column_occupancy,
    group_occupancy,
    select_columns,
    translate_columns,
)

# 5 sequences. Columns 1-2 are universal, column 3 is covered ONLY by the
# 2-member minority group, column 4 is empty everywhere.
ALN = {
    "A_1": "MK-–".replace("–", "-"),
    "A_2": "MK--",
    "A_3": "MK--",
    "B_1": "MKW-",
    "B_2": "MKW-",
}
GROUPS = {"big": ["A_1", "A_2", "A_3"], "small": ["B_1", "B_2"]}


def test_global_occupancy_is_the_fraction_of_non_gap_rows():
    assert column_occupancy(ALN) == [1.0, 1.0, 0.4, 0.0]


def test_group_occupancy_is_reported_per_group():
    occ = group_occupancy(ALN, GROUPS)

    assert occ["big"] == [1.0, 1.0, 0.0, 0.0]
    assert occ["small"] == [1.0, 1.0, 1.0, 0.0]


def test_global_selection_deletes_a_region_unique_to_the_minority_group():
    # column 3 is fully occupied within `small` but only 0.4 overall
    assert select_columns(ALN, threshold=0.5) == [1, 2]


def test_group_aware_selection_keeps_it():
    assert select_columns(ALN, threshold=0.5, groups=GROUPS) == [1, 2, 3]


def test_group_aware_selection_still_drops_columns_empty_in_every_group():
    kept = select_columns(ALN, threshold=0.5, groups=GROUPS)

    assert 4 not in kept


def test_all_gap_columns_are_dropped_even_at_threshold_zero():
    assert select_columns(ALN, threshold=0.0) == [1, 2, 3]


def test_groups_naming_sequences_outside_the_alignment_are_tolerated():
    groups = {"small": ["B_1", "B_2", "GHOST"]}

    assert select_columns(ALN, threshold=0.5, groups=groups) == [1, 2, 3]


def test_groups_matching_nothing_at_all_raise_instead_of_emptying_the_matrix():
    # Every group empty would silently select zero columns. A caller who passed
    # the wrong ids must hear about it, not receive an empty alignment.
    groups = {"phantom": ["NOPE"], "ghost": ["ALSO_NOPE"]}

    with pytest.raises(ValueError, match="no members"):
        select_columns(ALN, threshold=0.5, groups=groups)


def test_apply_columns_subsets_every_sequence_to_the_kept_columns():
    out = apply_columns(ALN, [1, 3])

    assert out == {"A_1": "M-", "A_2": "M-", "A_3": "M-",
                   "B_1": "MW", "B_2": "MW"}


def test_apply_columns_rejects_a_column_outside_the_alignment():
    with pytest.raises(ValueError):
        apply_columns(ALN, [1, 99])


def test_column_map_ties_kept_position_to_original_column_and_reference_residue():
    # reference B_1 = "MKW-": residues M=1, K=2, W=3; column 4 is a gap
    mapping = column_map([1, 3, 4], ALN["B_1"])

    assert mapping == [(1, 1, 1), (2, 3, 3), (3, 4, None)]


def test_column_map_reports_none_where_the_reference_is_gapped():
    mapping = column_map([3], ALN["A_1"])  # A_1 = "MK--", column 3 is a gap

    assert mapping == [(1, 3, None)]


def test_ragged_alignments_are_rejected_rather_than_silently_truncated():
    with pytest.raises(ValueError, match="unequal"):
        column_occupancy({"a": "MK", "b": "MKW"})


def test_empty_alignment_yields_no_columns():
    assert column_occupancy({}) == []
    assert select_columns({}, threshold=0.5) == []


# ---------------------------------------------------------------------------
# family-wide invariant columns (issue #42 A)
# ---------------------------------------------------------------------------

from utils.alignment import choose_reference, invariant_columns  # noqa: E402

# col1 universal M; col2 every group fixed but on DIFFERENT residues;
# col3 group "a" fixed, group "b" split; col4 universal but poorly covered in b
INV = {
    "a1": "MKWA", "a2": "MKWA", "a3": "MKWA",
    "b1": "MQWA", "b2": "MQY-", "b3": "MQW-",
}
INV_GROUPS = {"a": ["a1", "a2", "a3"], "b": ["b1", "b2", "b3"]}


def test_invariant_column_needs_every_group_fixed_on_the_same_residue():
    inv = invariant_columns(INV, INV_GROUPS, threshold=1.0, min_cover=0.5)

    assert 1 in inv and inv[1] == "M"


def test_a_column_where_groups_are_each_fixed_but_disagree_is_not_invariant():
    # col2: group a is all K, group b is all Q. Both conserved, not shared.
    inv = invariant_columns(INV, INV_GROUPS, threshold=1.0, min_cover=0.5)

    assert 2 not in inv


def test_a_column_one_group_fails_to_conserve_is_not_invariant():
    # col3: a is WWW, b is WYW -> b at 2/3
    inv = invariant_columns(INV, INV_GROUPS, threshold=1.0, min_cover=0.5)

    assert 3 not in inv


def test_threshold_is_applied_within_each_group_not_across_the_alignment():
    # at 0.66 group b's 2/3 W passes, so col3 becomes invariant
    inv = invariant_columns(INV, INV_GROUPS, threshold=0.66, min_cover=0.5)

    assert inv[3] == "W"


def test_a_group_below_min_cover_disqualifies_the_column():
    # col4: group b covers 1/3 -> cannot be called invariant on b's behalf
    inv = invariant_columns(INV, INV_GROUPS, threshold=1.0, min_cover=0.5)

    assert 4 not in inv


def test_invariant_columns_ignores_sequences_outside_the_groups():
    aln = dict(INV, outsider="XXXX")

    inv = invariant_columns(aln, INV_GROUPS, threshold=1.0, min_cover=0.5)

    assert inv[1] == "M"


def test_invariant_columns_rejects_groups_that_match_nothing():
    with pytest.raises(ValueError, match="no members"):
        invariant_columns(INV, {"ghost": ["NOPE"]}, threshold=1.0)


# ---------------------------------------------------------------------------
# reference auto-selection (issue #42 B)
# ---------------------------------------------------------------------------

def test_reference_is_the_most_complete_sequence():
    aln = {"short": "M--A", "long": "MKWA", "mid": "MK-A"}

    assert choose_reference(aln) == "long"


def test_reference_tie_is_broken_by_identifier_so_the_choice_is_deterministic():
    aln = {"zeta": "MKWA", "alpha": "MKWA", "mid": "MK-A"}

    assert choose_reference(aln) == "alpha"


def test_reference_selection_can_be_restricted_to_named_sequences():
    aln = {"short": "M--A", "long": "MKWA", "mid": "MK-A"}

    assert choose_reference(aln, candidates=["short", "mid"]) == "mid"


def test_reference_selection_on_an_empty_alignment_raises():
    with pytest.raises(ValueError):
        choose_reference({})


# --------------------------------------------------------------------------
# coordinate stamping and translation (issue #42)
#
# Residue-level artifacts carry column numbers, and nothing currently records
# WHICH alignment those columns belong to. beb_cross crosses codeml BEB sites
# (family codon alignment, untrimmed) with signal windows (clan protein
# alignment, trimAl-trimmed) by plain interval overlap, so a mismatch yields 0
# overlaps -- which classify() then reports as evidence against
# neofunctionalization. Detection is not enough: the two alignments differ
# legitimately, so they must be translatable.
# --------------------------------------------------------------------------

ALN_A = {"g1": "MK-W", "g2": "MKYW"}
ALN_B = {"g1": "--MKW", "g2": "XXMKW"}   # g1 ungapped MKW in both


def test_alignment_id_is_stable_for_identical_content():
    assert alignment_id(ALN_A) == alignment_id(dict(ALN_A))


def test_alignment_id_ignores_dictionary_order():
    reordered = {"g2": ALN_A["g2"], "g1": ALN_A["g1"]}

    assert alignment_id(reordered) == alignment_id(ALN_A)


def test_alignment_id_changes_when_a_column_is_trimmed():
    trimmed = {k: v[:3] for k, v in ALN_A.items()}

    assert alignment_id(trimmed) != alignment_id(ALN_A)


def test_alignment_id_reports_shape_alongside_the_digest():
    stamp = alignment_id(ALN_A)

    assert stamp.n_columns == 4
    assert stamp.n_sequences == 2
    assert len(stamp.digest) >= 8


# ------------------------------------------------------------ translation --

def test_translating_into_the_same_alignment_is_the_identity():
    assert translate_columns([1, 2, 4], ALN_A, ALN_A, via="g1") == [1, 2, 4]


def test_columns_translate_through_a_shared_sequence():
    # g1: ALN_A "MK-W" cols 1,2,4 hold M,K,W; ALN_B "--MKW" holds them at 3,4,5
    assert translate_columns([1, 2, 4], ALN_A, ALN_B, via="g1") == [3, 4, 5]


def test_a_column_where_the_bridge_sequence_is_gapped_cannot_be_translated():
    # ALN_A col 3 is a gap in g1
    assert translate_columns([3], ALN_A, ALN_B, via="g1") == [None]


def test_translation_refuses_a_bridge_absent_from_either_alignment():
    with pytest.raises(ValueError, match="g9"):
        translate_columns([1], ALN_A, ALN_B, via="g9")


def test_translation_refuses_a_bridge_whose_residues_disagree():
    # same id, different underlying protein -> translating would be nonsense
    other = {"g1": "MKQW-"}

    with pytest.raises(ValueError, match="differ"):
        translate_columns([1], ALN_A, other, via="g1")


def test_translation_rejects_a_column_outside_the_source_alignment():
    with pytest.raises(ValueError, match="outside"):
        translate_columns([99], ALN_A, ALN_B, via="g1")
