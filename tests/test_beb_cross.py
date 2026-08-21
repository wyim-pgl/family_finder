"""beb_cross.py — codeml branch-site BEB sites x DeepLoc signal windows
(issue #24 final step). Text fixtures only, no codeml invoked.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from beb_cross import cross_windows, parse_beb

MLC_TEXT = """\
some codeml preamble
lnL(ntime: 20  np: 25):  -12345.678
Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)
Positive sites for foreground lineages Prob(w>1):
   170 K 0.712
   200 R 0.981*
   500 E 0.994**

The grid ...
"""

WINDOWS_TSV = """\
seq_id\tstart\tend\tpeak_pos\tpeak_attention\tmean_attention\taln_col_start\taln_col_end\tdeeploc_signals
GeneA\t13\t30\t27\t0.019\t0.013\t165\t202\tNuclear localization signal
GeneB\t400\t420\t410\t0.02\t0.011\t480\t510\t
"""


def test_parse_beb_extracts_sites():
    sites = parse_beb(MLC_TEXT)
    assert sites == [(170, "K", 0.712), (200, "R", 0.981), (500, "E", 0.994)]


def test_parse_beb_empty_when_no_section():
    assert parse_beb("no beb here\n") == []


def test_parse_beb_stops_at_blank_line():
    # nothing after the section's terminating blank line is a site
    text = MLC_TEXT + "   999 Q 0.9\n"
    sites = parse_beb(text)
    assert (999, "Q", 0.9) not in sites


def test_cross_windows_flags_overlap(tmp_path):
    p = tmp_path / "w.tsv"
    p.write_text(WINDOWS_TSV)
    rows = cross_windows(parse_beb(MLC_TEXT), p)
    by_site = {r["site"]: r for r in rows}
    # 170 and 200 fall inside GeneA's aln_col 165-202 window
    assert "GeneA" in by_site[170]["window_seqs"]
    assert "GeneA" in by_site[200]["window_seqs"]
    # 500 overlaps GeneB's 480-510 window
    assert by_site[500]["window_seqs"] == "GeneB"


def test_cross_windows_site_outside_all_windows(tmp_path):
    p = tmp_path / "w.tsv"
    p.write_text(WINDOWS_TSV)
    rows = cross_windows([(300, "A", 0.9)], p)
    assert rows[0]["window_seqs"] == ""
    assert rows[0]["n_windows"] == 0


def test_cross_windows_signal_types_aggregated(tmp_path):
    p = tmp_path / "w.tsv"
    p.write_text(WINDOWS_TSV)
    rows = cross_windows(parse_beb(MLC_TEXT), p)
    by_site = {r["site"]: r for r in rows}
    assert "Nuclear localization signal" in by_site[170]["signals"]
