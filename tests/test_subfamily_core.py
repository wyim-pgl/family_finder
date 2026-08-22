"""The invariant core, the coordinate it is measured in, and who it belongs to.

Issue #42, three defects that all read the same way from outside — a number
that is present, plausible and wrong about what it measured:

  * `sdp_core_relationship` computed the invariant core but never said which
    columns `min_cover` kept it from examining. A region gappy in one subfamily
    came back as "not invariant" when it was never looked at — the same silence
    `coverage_suppressed` already removed from `sdp_scan`.
  * the distance to that core was measured in ALIGNMENT COLUMNS. Columns are
    not residues: a gapped stretch pushes two columns apart without moving a
    single residue, so the observed distance and the null were both inflated by
    indels that no protein contains.
  * the coordinates were reported without naming the sequence they are in, and
    the reference was whatever the caller happened to pass.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from steps.subfamily import (  # noqa: E402
    resolve_reference,
    sdp_core_relationship,
    sdp_scan,
)
from utils.alignment import invariant_suppressed  # noqa: E402


# ---------------------------------------------------------------------------
# min_cover suppression of the invariant core
# ---------------------------------------------------------------------------

def _gappy_core_fixture():
    """Columns 5-6 are the SAME residue in both groups but group b is gappy.

    Four of b's five members have no residue there, so coverage is 0.2 and
    `invariant_columns` drops the pair. Agreement was never in question — the
    columns were never examined.
    """
    aln = {f"a{i}": "MMMM" + "WW" + "KKKK" for i in range(5)}
    aln["b0"] = "MMMM" + "WW" + "QQQQ"
    for i in range(1, 5):
        aln[f"b{i}"] = "MMMM" + "--" + "QQQQ"
    groups = {"a": [f"a{i}" for i in range(5)], "b": [f"b{i}" for i in range(5)]}
    return aln, groups


def test_columns_min_cover_kept_out_of_the_core_are_listed():
    aln, groups = _gappy_core_fixture()

    report = invariant_suppressed(aln, groups, min_cover=0.7)

    assert report["columns"] == [5, 6]
    assert report["n_suppressed"] == 2


def test_the_suppressing_group_is_named_not_just_counted():
    aln, groups = _gappy_core_fixture()

    report = invariant_suppressed(aln, groups, min_cover=0.7)

    assert report["by_group"]["b"] == [5, 6]
    assert report["by_group"]["a"] == []


def test_core_relationship_reports_the_columns_it_never_examined():
    aln, groups = _gappy_core_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=5)

    assert rel["n_core_suppressed_columns"] == 2
    assert rel["core_suppressed_columns"] == [5, 6]
    assert rel["core_suppressed_by_group"] == {"a": 0, "b": 2}


def test_an_empty_core_says_how_much_of_it_was_never_looked_at():
    # Every column outside the diagnostic block is gappy in one group, so
    # there is no invariant core AND no evidence that the groups disagree.
    # "No core" and "never examined" must not print the same sentence.
    aln = {f"a{i}": "MMMM" + "KKKK" for i in range(5)}
    aln["b0"] = "MMMM" + "QQQQ"
    for i in range(1, 5):
        aln[f"b{i}"] = "----" + "QQQQ"
    groups = {"a": [f"a{i}" for i in range(5)], "b": [f"b{i}" for i in range(5)]}

    rel = sdp_core_relationship(aln, groups, min_group=5)

    assert rel["verdict"] == "no_interpretation_available"
    assert rel["n_core_suppressed_columns"] == 4
    assert "never examined" in rel["reason"]


# ---------------------------------------------------------------------------
# residue-space distance instead of column-space distance
# ---------------------------------------------------------------------------

def _indel_fixture():
    """Core, then a six-column insertion no group member carries, then SDPs.

    In column units the diagnostic block sits 7-11 columns from the core; in
    the residue coordinates of any member it sits 1-5 residues away, because
    the intervening columns hold no residue of theirs at all.
    """
    aln = {f"a{i}": "MMMMM" + "------" + "KKKKK" for i in range(5)}
    aln.update({f"b{i}": "MMMMM" + "------" + "QQQQQ" for i in range(5)})
    # Present only in a sequence that is shorter overall, so it cannot become
    # the automatic reference and cannot rescue the columns for any group.
    aln["insert_only"] = "-----" + "IIIIII" + "-----"
    groups = {"a": [f"a{i}" for i in range(5)], "b": [f"b{i}" for i in range(5)]}
    return aln, groups


def test_distance_is_measured_in_residues_not_columns():
    aln, groups = _indel_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=5, seed=0)

    # residues 6-10 against core residues 1-5: 1,2,3,4,5 -> median 3
    assert rel["observed_median_distance"] == 3


def test_the_column_number_would_have_been_the_wrong_answer():
    # Same alignment, old convention: the six insertion columns add six to
    # every distance although no member of either group has a residue in them.
    aln, groups = _indel_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=5, seed=0,
                                distance_space="column")

    assert rel["observed_median_distance"] == 9


def test_the_convention_and_its_reference_sequence_are_reported():
    aln, groups = _indel_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=5, seed=0)

    assert rel["distance_convention"] == "residue:a0"
    assert rel["distance_reference"] == "a0"


def test_the_column_convention_says_so_and_names_no_reference():
    aln, groups = _indel_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=5, seed=0,
                                distance_space="column")

    assert rel["distance_convention"] == "column"
    assert rel["distance_reference"] is None


def test_the_verdict_vocabulary_is_unchanged_by_the_new_convention():
    aln, groups = _indel_fixture()

    rel = sdp_core_relationship(aln, groups, min_group=5, seed=0)

    assert rel["verdict"] in {"avoids_core", "at_core", "indistinguishable",
                             "no_interpretation_available"}


def test_positions_the_reference_cannot_place_are_counted():
    aln, groups = _indel_fixture()

    # Forcing the insertion-only sequence as the reference leaves every core
    # and diagnostic column without a residue to be numbered against.
    rel = sdp_core_relationship(aln, groups, min_group=5,
                                ref_seq_id="insert_only")

    assert rel["n_sdp_untranslatable"] == 5
    assert rel["n_invariant_untranslatable"] == 5
    assert rel["verdict"] == "no_interpretation_available"
    assert "insert_only" in rel["reason"]


def test_an_unknown_distance_space_is_refused_rather_than_guessed():
    aln, groups = _indel_fixture()

    with pytest.raises(ValueError, match="distance_space"):
        sdp_core_relationship(aln, groups, min_group=5, distance_space="codon")


# ---------------------------------------------------------------------------
# reference auto-selection
# ---------------------------------------------------------------------------

REF_ALN = {
    "Sp1_a1": "KAAWTR", "Sp2_a2": "KACWTR", "Sp3_a3": "KAAWTR", "Sp4_a4": "KADWTR",
    "Sp1_b1": "SAAWTR", "Sp2_b2": "SAEWTR", "Sp3_b3": "SAAWTR", "Sp4_b4": "SAFWTR",
    "ATH_PPC1": "SA-WTR",  # characterised, but one residue short of the rest
}
REF_GROUPS = {
    "SF_A": ["Sp1_a1", "Sp2_a2", "Sp3_a3", "Sp4_a4"],
    "SF_B": ["Sp1_b1", "Sp2_b2", "Sp3_b3", "Sp4_b4"],
}


def test_a_characterised_member_is_preferred_over_the_longest_sequence():
    choice = resolve_reference(REF_ALN, characterised=["ATH_PPC1"])

    assert choice.seq_id == "ATH_PPC1"
    assert choice.source == "characterised"


def test_without_a_characterised_id_the_automatic_representative_is_used():
    choice = resolve_reference(REF_ALN)

    assert choice.seq_id == "Sp1_a1"  # most complete, ties broken on the id
    assert choice.source == "automatic"


def test_a_characterised_id_absent_from_the_alignment_falls_back_and_says_so():
    choice = resolve_reference(REF_ALN, characterised=["ZEA_PPC7"])

    assert choice.source == "automatic"
    assert choice.unmatched_characterised == ["ZEA_PPC7"]


def test_a_characterised_id_spelled_with_a_transcript_suffix_still_matches():
    choice = resolve_reference(REF_ALN, characterised=["ATH_PPC1.t1"])

    assert choice.seq_id == "ATH_PPC1"
    assert choice.source == "characterised"


def test_an_explicit_reference_outranks_a_characterised_one():
    choice = resolve_reference(REF_ALN, ref_seq_id="Sp4_b4",
                               characterised=["ATH_PPC1"])

    assert choice.seq_id == "Sp4_b4"
    assert choice.source == "explicit"


def test_an_explicit_reference_that_is_not_there_is_refused():
    with pytest.raises(ValueError, match="ZEA_PPC7"):
        resolve_reference(REF_ALN, ref_seq_id="ZEA_PPC7")


def test_every_scanned_row_names_the_sequence_its_positions_are_in():
    rows = sdp_scan(REF_ALN, REF_GROUPS, min_group=4,
                    characterised=["ATH_PPC1"])

    assert rows
    assert {r["ref_seq"] for r in rows} == {"ATH_PPC1"}


def test_positions_follow_the_characterised_reference_own_numbering():
    # ATH_PPC1 is gapped at column 3, so column 4 is its third residue while
    # it is the fourth of every other sequence.
    rows = sdp_scan(REF_ALN, REF_GROUPS, min_group=4,
                    characterised=["ATH_PPC1"])
    by_col = {r["aln_col"]: r for r in rows}

    assert by_col[1]["ref_pos"] == 1


def test_a_scan_with_no_reference_at_all_still_names_one():
    rows = sdp_scan(REF_ALN, REF_GROUPS, min_group=4)

    assert {r["ref_seq"] for r in rows} == {"Sp1_a1"}
    assert all(r["ref_pos"] is not None for r in rows)
