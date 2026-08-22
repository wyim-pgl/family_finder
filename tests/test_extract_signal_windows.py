"""extract_signal_windows.py — DeepLoc attention windows mapped to alignment
columns (issue #24 step 2, brought into the repo for #42).

The script existed only on the GPU box, so `aln_col_start` had a provenance
nobody could check: `beb_cross.py` asserted in its docstring that those columns
were the same coordinate system as codeml's BEB sites, and nothing could
confirm or refute it. Every window row now carries the `alignment_id()` stamp
of the alignment its columns came from.

Synthetic CSV/FASTA fixtures only — DeepLoc is never invoked (house rule).
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from extract_signal_windows import (  # noqa: E402
    extract_windows,
    map_windows_to_columns,
    read_attention,
    write_windows,
)
from utils.alignment import alignment_id  # noqa: E402


# ---------------------------------------------------------------------------
# window extraction: runs above mean + 2SD, gaps <= 2 merged, >= 5 aa kept
# ---------------------------------------------------------------------------

# Keep the high fraction away from 0.2: for a two-valued profile with exactly
# 20% peaks, mean + 2SD lands ON the peak value and `v > threshold` is a
# coin toss decided by floating point.
def _profile(peak_at, n=60, base=0.001, peak=0.5, width=8):
    vals = [base] * n
    for i in range(peak_at - 1, peak_at - 1 + width):
        vals[i] = peak
    return vals


def test_a_flat_profile_has_no_windows():
    assert extract_windows([0.01] * 30) == []


def test_a_sustained_attention_run_becomes_one_window():
    windows = extract_windows(_profile(11))

    assert windows == [(11, 18)]


def test_a_run_shorter_than_the_minimum_is_dropped():
    vals = [0.001] * 60
    vals[10:13] = [0.5, 0.5, 0.5]   # only 3 aa above threshold

    assert extract_windows(vals) == []


def test_two_runs_separated_by_a_short_gap_are_merged():
    vals = [0.001] * 60
    vals[10:14] = [0.5] * 4
    vals[16:20] = [0.5] * 4         # gap of 2 residues

    assert extract_windows(vals) == [(11, 20)]


def test_two_runs_separated_by_a_long_gap_stay_apart():
    vals = [0.001] * 80
    vals[10:16] = [0.5] * 6
    vals[30:36] = [0.5] * 6

    assert extract_windows(vals) == [(11, 16), (31, 36)]


def test_positions_are_one_based():
    vals = [0.5] * 6 + [0.001] * 60

    assert extract_windows(vals)[0][0] == 1


# ---------------------------------------------------------------------------
# residue -> alignment column
# ---------------------------------------------------------------------------

ALN = {"g1": "--MKWQRSTV", "g2": "XXMKWQRSTV"}


def test_windows_are_reported_in_alignment_columns():
    rows = map_windows_to_columns({"g1": [(1, 5)]}, ALN)

    # residue 1 of g1 is column 3, residue 5 is column 7
    assert rows[0]["aln_col_start"] == 3
    assert rows[0]["aln_col_end"] == 7


def test_every_row_carries_the_stamp_of_the_alignment_it_used():
    rows = map_windows_to_columns({"g1": [(1, 5)]}, ALN)

    assert rows[0]["aln_stamp"] == str(alignment_id(ALN))


def test_a_sequence_absent_from_the_alignment_is_reported_not_blanked():
    # the original script wrote an empty aln_col_start and moved on, which is
    # the silent coverage loss this issue is about
    rows = map_windows_to_columns({"ghost": [(1, 5)]}, ALN)

    assert rows == []


def test_the_unmapped_sequences_are_returned_alongside_the_rows():
    rows, skipped = map_windows_to_columns({"ghost": [(1, 5)], "g1": [(1, 5)]},
                                           ALN, report_skipped=True)

    assert [r["seq_id"] for r in rows] == ["g1"]
    assert skipped == {"ghost": "not in the alignment"}


def test_a_window_running_past_the_sequence_is_reported_not_truncated():
    rows, skipped = map_windows_to_columns({"g1": [(1, 99)]}, ALN,
                                           report_skipped=True)

    assert rows == []
    assert "g1" in skipped


# ---------------------------------------------------------------------------
# attention CSV
# ---------------------------------------------------------------------------

def test_read_attention_parses_the_deeploc_two_column_csv(tmp_path):
    p = tmp_path / "alpha_g1.csv"
    p.write_text("AA,Alpha\nM,0.01\nK,0.02\nW,0.03\n")

    assert read_attention(p) == [0.01, 0.02, 0.03]


def test_read_attention_refuses_an_unexpected_layout(tmp_path):
    p = tmp_path / "alpha_g1.csv"
    p.write_text("AA,Alpha,Extra\nM,0.01,x\n")

    with pytest.raises(ValueError, match="two columns"):
        read_attention(p)


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

def test_written_rows_keep_the_stamp_column(tmp_path):
    rows = map_windows_to_columns({"g1": [(1, 5)]}, ALN)
    out = tmp_path / "signal_windows.tsv"

    write_windows(rows, out)

    header, row = out.read_text().splitlines()[:2]
    assert "aln_stamp" in header.split("\t")
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert fields["aln_stamp"] == str(alignment_id(ALN))


def test_an_empty_result_still_writes_a_header(tmp_path):
    out = tmp_path / "signal_windows.tsv"

    write_windows([], out)

    assert out.read_text().strip().split("\t")[0] == "seq_id"
