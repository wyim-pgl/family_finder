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
    apply_columns,
    column_map,
    column_occupancy,
    group_occupancy,
    select_columns,
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
