"""beb_cross.py — codeml branch-site BEB sites x DeepLoc signal windows
(issue #24 final step). Text fixtures only, no codeml invoked.
"""

import sys
from pathlib import Path

import pytest

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
    rows = cross_windows(parse_beb(MLC_TEXT), p, allow_unverified=True)
    by_site = {r["site"]: r for r in rows}
    # 170 and 200 fall inside GeneA's aln_col 165-202 window
    assert "GeneA" in by_site[170]["window_seqs"]
    assert "GeneA" in by_site[200]["window_seqs"]
    # 500 overlaps GeneB's 480-510 window
    assert by_site[500]["window_seqs"] == "GeneB"


def test_cross_windows_site_outside_all_windows(tmp_path):
    p = tmp_path / "w.tsv"
    p.write_text(WINDOWS_TSV)
    rows = cross_windows([(300, "A", 0.9)], p, allow_unverified=True)
    assert rows[0]["window_seqs"] == ""
    assert rows[0]["n_windows"] == 0


def test_cross_windows_signal_types_aggregated(tmp_path):
    p = tmp_path / "w.tsv"
    p.write_text(WINDOWS_TSV)
    rows = cross_windows(parse_beb(MLC_TEXT), p, allow_unverified=True)
    by_site = {r["site"]: r for r in rows}
    assert "Nuclear localization signal" in by_site[170]["signals"]


# --------------------------------------------------------------------------
# coordinate verification (issue #42)
#
# BEB site numbers are columns of the family's untrimmed codon alignment;
# signal windows were written in clan protein-alignment columns, which for the
# PEPC work was trimAl-trimmed from 1,468 to 876. Overlapping them by interval
# alone always yields a number, and zero overlaps is what
# subfunctionalization.classify() reports as evidence against
# neofunctionalization.
# --------------------------------------------------------------------------

WINDOWS_STAMPED = (
    "seq_id\tstart\tend\tpeak_pos\tpeak_attention\tmean_attention\t"
    "aln_col_start\taln_col_end\tdeeploc_signals\taln_stamp\n"
    "g1\t1\t10\t5\t0.9\t0.5\t100\t200\tNucleus\tabc123:5x876\n"
)


def test_cross_refuses_when_the_two_coordinate_systems_disagree(tmp_path):
    w = tmp_path / "w.tsv"
    w.write_text(WINDOWS_STAMPED)

    with pytest.raises(ValueError, match="coordinate"):
        cross_windows([(150, "A", 0.9)], w, site_stamp="different:5x1468")


def test_cross_proceeds_when_the_stamps_match(tmp_path):
    w = tmp_path / "w.tsv"
    w.write_text(WINDOWS_STAMPED)

    rows = cross_windows([(150, "A", 0.9)], w, site_stamp="abc123:5x876")

    assert rows[0]["n_windows"] == 1
    assert rows[0]["coordinates_verified"] is True


def test_cross_refuses_an_unstamped_window_file_unless_told_to(tmp_path):
    w = tmp_path / "w.tsv"
    w.write_text(WINDOWS_STAMPED.replace("\taln_stamp", "\tunused")
                                .replace("\tabc123:5x876", "\t-"))

    with pytest.raises(ValueError, match="coordinate"):
        cross_windows([(150, "A", 0.9)], w, site_stamp="abc123:5x876")


def test_an_unverified_cross_is_marked_as_such_in_every_row(tmp_path):
    w = tmp_path / "w.tsv"
    w.write_text(WINDOWS_STAMPED.replace("\taln_stamp", "\tunused")
                                .replace("\tabc123:5x876", "\t-"))

    rows = cross_windows([(150, "A", 0.9)], w, allow_unverified=True)

    assert rows[0]["coordinates_verified"] is False


def test_cross_translates_between_alignments_when_both_are_supplied(tmp_path):
    # BEB sites live in the family codon alignment; windows in the clan
    # protein alignment. Given both and a bridge sequence, the cross should
    # translate rather than refuse (#42).
    w = tmp_path / "w.tsv"
    w.write_text(
        "seq_id\tstart\tend\tpeak_pos\tpeak_attention\tmean_attention\t"
        "aln_col_start\taln_col_end\tdeeploc_signals\n"
        "g1\t1\t3\t2\t0.9\t0.5\t3\t5\tNucleus\n"
    )
    site_aln = {"g1": "MK-W", "g2": "MKYW"}
    window_aln = {"g1": "--MKW", "g2": "XXMKW"}

    # site 4 of the family alignment is residue 3 of g1 -> column 5 of the clan
    rows = cross_windows([(4, "W", 0.9)], w,
                         site_alignment=site_aln, window_alignment=window_aln,
                         bridge="g1")

    assert rows[0]["n_windows"] == 1
    assert rows[0]["coordinates_verified"] is True
    assert rows[0]["translated_site"] == 5


def test_a_site_the_bridge_cannot_carry_is_reported_not_dropped(tmp_path):
    w = tmp_path / "w.tsv"
    w.write_text(
        "seq_id\tstart\tend\tpeak_pos\tpeak_attention\tmean_attention\t"
        "aln_col_start\taln_col_end\tdeeploc_signals\n"
        "g1\t1\t3\t2\t0.9\t0.5\t3\t5\tNucleus\n"
    )
    site_aln = {"g1": "MK-W", "g2": "MKYW"}
    window_aln = {"g1": "--MKW", "g2": "XXMKW"}

    # column 3 is a gap in the bridge, so it has no clan-alignment position
    rows = cross_windows([(3, "Y", 0.9)], w,
                         site_alignment=site_aln, window_alignment=window_aln,
                         bridge="g1")

    assert rows[0]["translated_site"] is None
    assert rows[0]["n_windows"] == 0
    assert rows[0]["untranslatable"] is True


# --------------------------------------------------------------------------
# CLI wiring (issue #42)
#
# cross_windows() gained the stamp check and the translation path, but main()
# passed neither, so the command line could only ever raise. The translation
# has to be reachable from the CLI or the pipeline path has no way to use it.
# --------------------------------------------------------------------------

def _mlc(tmp_path):
    p = tmp_path / "mlc"
    p.write_text(
        "Bayes Empirical Bayes (BEB) analysis\n"
        "Positive sites for foreground lineages Prob(w>1):\n"
        "     4 W 0.994**\n"
        "\n"
    )
    return p


def _windows(tmp_path):
    p = tmp_path / "w.tsv"
    p.write_text(
        "seq_id\tstart\tend\tpeak_pos\tpeak_attention\tmean_attention\t"
        "aln_col_start\taln_col_end\tdeeploc_signals\n"
        "g1\t1\t3\t2\t0.9\t0.5\t3\t5\tNucleus\n"
    )
    return p


def _fasta(tmp_path, name, seqs):
    p = tmp_path / name
    p.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    return p


def test_cli_refuses_an_unverifiable_cross(tmp_path, monkeypatch):
    from beb_cross import main

    monkeypatch.setattr("sys.argv", [
        "beb_cross.py", "--mlc", str(_mlc(tmp_path)),
        "--windows", str(_windows(tmp_path)), "-o", str(tmp_path / "o.tsv"),
    ])
    with pytest.raises(ValueError, match="coordinate"):
        main()


def test_cli_translates_when_given_both_alignments_and_a_bridge(tmp_path, monkeypatch):
    from beb_cross import main

    site_aln = _fasta(tmp_path, "s.fa", {"g1": "MK-W", "g2": "MKYW"})
    window_aln = _fasta(tmp_path, "w.fa", {"g1": "--MKW", "g2": "XXMKW"})
    out = tmp_path / "o.tsv"
    monkeypatch.setattr("sys.argv", [
        "beb_cross.py", "--mlc", str(_mlc(tmp_path)),
        "--windows", str(_windows(tmp_path)), "-o", str(out),
        "--site-alignment", str(site_aln),
        "--window-alignment", str(window_aln), "--bridge", "g1",
    ])

    assert main() == 0
    header, row = out.read_text().splitlines()[:2]
    fields = dict(zip(header.split("\t"), row.split("\t")))
    # site 4 is residue 3 of g1 -> window-alignment column 5, inside 3-5
    assert fields["translated_site"] == "5"
    assert fields["n_windows"] == "1"
    assert fields["coordinates_verified"] == "True"


def test_cli_can_be_told_to_proceed_unverified_and_marks_every_row(tmp_path, monkeypatch):
    from beb_cross import main

    out = tmp_path / "o.tsv"
    monkeypatch.setattr("sys.argv", [
        "beb_cross.py", "--mlc", str(_mlc(tmp_path)),
        "--windows", str(_windows(tmp_path)), "-o", str(out),
        "--allow-unverified",
    ])

    assert main() == 0
    header, row = out.read_text().splitlines()[:2]
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert fields["coordinates_verified"] == "False"


def test_a_codon_alignment_is_collapsed_to_the_columns_beb_numbers_use():
    # steps/codeml.py runs on confirmed_codon.afa -- nucleotides -- but BEB
    # site numbers are amino-acid positions. Without this step the family
    # codeml result cannot be translated against a protein alignment at all.
    from beb_cross import codon_alignment_as_protein

    protein = codon_alignment_as_protein({"g1": "ATGAAA---TGG",
                                          "g2": "ATGAAATATTGG"})

    assert protein == {"g1": "MK-W", "g2": "MKYW"}


def test_a_codon_alignment_whose_width_is_not_a_multiple_of_three_is_refused():
    from beb_cross import codon_alignment_as_protein

    with pytest.raises(ValueError, match="codon"):
        codon_alignment_as_protein({"g1": "ATGAA"})


def test_a_partial_codon_is_a_gap_not_a_guessed_residue():
    from beb_cross import codon_alignment_as_protein

    protein = codon_alignment_as_protein({"g1": "ATGA--", "g2": "ATGAAA"})

    assert protein["g1"] == "MX"


def test_cli_accepts_the_family_codon_alignment_directly(tmp_path, monkeypatch):
    from beb_cross import main

    site_aln = _fasta(tmp_path, "confirmed_codon.afa",
                      {"g1": "ATGAAA---TGG", "g2": "ATGAAATATTGG"})
    window_aln = _fasta(tmp_path, "clan.aln", {"g1": "--MKW", "g2": "XXMKW"})
    out = tmp_path / "o.tsv"
    monkeypatch.setattr("sys.argv", [
        "beb_cross.py", "--mlc", str(_mlc(tmp_path)),
        "--windows", str(_windows(tmp_path)), "-o", str(out),
        "--site-alignment", str(site_aln), "--site-alignment-is-codon",
        "--window-alignment", str(window_aln), "--bridge", "g1",
    ])

    assert main() == 0
    header, row = out.read_text().splitlines()[:2]
    fields = dict(zip(header.split("\t"), row.split("\t")))
    assert fields["translated_site"] == "5"
    assert fields["coordinates_verified"] == "True"
