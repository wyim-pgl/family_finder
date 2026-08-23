"""The rescue tier must decide by bit score, not by E-value.

Issue #45. `_parse_hmmsearch_tblout` kept the hit with the lowest E-value and
broke ties with `<`, so the family seen FIRST won. First is tblout order, which
is profile-database order, which is `sorted(glob("*.hmm"))` - and '0' < '_', so
`R10_*` sorts before `R1_*`. Measured on output_5sp_v2: of 29,943 rescued genes
1,036 (3.5%) had a tie at the best E-value, 1,030 of those at E=0 (underflow),
and bit scores would have separated 1,035 of them. Twenty-four genes sit in a
family that the bit scores do not support.

CLAUDE.md already forbids this for the per-round tier. The rule is the same
here.
"""
import sys
import types
from pathlib import Path

import pytest

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import Config
from steps.hmmer_rescue import _parse_hmmsearch_tblout


def _tbl(tmp_path, rows):
    """rows: (gene, family, evalue, bits) in the order hmmsearch emitted them."""
    path = tmp_path / "hits.tblout"
    with open(path, "w") as fh:
        fh.write("#  target  acc  query  acc  E-value  score\n")
        for gene, fam, ev, bits in rows:
            fh.write(f"{gene} - {fam} - {ev} {bits} 0.0 1.0 1.0 1.0 - -\n")
    return path


CFG = Config()


def test_the_higher_bit_score_wins_when_evalues_are_equal(tmp_path):
    """Both underflow to 0; only the bits distinguish them."""
    tbl = _tbl(tmp_path, [
        ("gene1", "R10_OG0000001", "0.0", "120.0"),   # emitted first
        ("gene1", "R1_OG0000002", "0.0", "980.0"),    # far better fit
    ])
    hits = _parse_hmmsearch_tblout(tbl, CFG.hmmer_evalue, CFG)
    assert hits["gene1"].family == "R1_OG0000002"


def test_emission_order_cannot_decide(tmp_path):
    """The same two hits in the opposite order must give the same answer."""
    rows = [("g", "R1_A", "0.0", "980.0"), ("g", "R10_B", "0.0", "120.0")]
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir(); b.mkdir()
    first = _parse_hmmsearch_tblout(_tbl(a, rows), CFG.hmmer_evalue, CFG)
    second = _parse_hmmsearch_tblout(
        _tbl(b, list(reversed(rows))), CFG.hmmer_evalue, CFG)
    assert first["g"].family == second["g"].family == "R1_A"


def test_a_true_tie_is_reported_not_broken(tmp_path):
    """Equal bits too: there is no evidence for either, and saying so is the
    only honest option. Picking one by sort order is what this issue is about."""
    tbl = _tbl(tmp_path, [
        ("g", "R10_A", "0.0", "500.0"),
        ("g", "R1_B", "0.0", "500.0"),
    ])
    hit = _parse_hmmsearch_tblout(tbl, CFG.hmmer_evalue, CFG)["g"]
    assert hit.grade == "UNRESOLVED"
    assert "tie" in hit.reason.lower()


def test_a_clear_margin_is_graded_high(tmp_path):
    tbl = _tbl(tmp_path, [("g", "A", "1e-40", "800.0"), ("g", "B", "1e-20", "300.0")])
    hit = _parse_hmmsearch_tblout(tbl, CFG.hmmer_evalue, CFG)["g"]
    assert hit.grade == "HIGH"
    assert hit.margin == pytest.approx(500.0)


def test_a_thin_margin_is_graded_provisional_not_silently_equal(tmp_path):
    """The 1,747 genes (5.8%) the per-round tier would have refused must not
    leave the rescue looking identical to the confident ones."""
    tbl = _tbl(tmp_path, [("g", "A", "1e-40", "500.0"), ("g", "B", "1e-39", "498.0")])
    hit = _parse_hmmsearch_tblout(tbl, CFG.hmmer_evalue, CFG)["g"]
    assert hit.grade == "PROVISIONAL"


def test_a_single_hit_has_no_margin_and_is_still_high(tmp_path):
    tbl = _tbl(tmp_path, [("g", "A", "1e-40", "800.0")])
    hit = _parse_hmmsearch_tblout(tbl, CFG.hmmer_evalue, CFG)["g"]
    assert hit.grade == "HIGH" and hit.margin is None


def test_the_evalue_cutoff_still_filters(tmp_path):
    tbl = _tbl(tmp_path, [("g", "A", "1.0", "800.0")])
    assert _parse_hmmsearch_tblout(tbl, 1e-5, CFG) == {}


def test_the_margin_threshold_scales_with_the_best_score(tmp_path):
    """max(margin_bits, margin_frac * best): at 10,000 bits, 10 is nothing."""
    tbl = _tbl(tmp_path, [("g", "A", "0.0", "10000.0"), ("g", "B", "0.0", "9950.0")])
    hit = _parse_hmmsearch_tblout(tbl, CFG.hmmer_evalue, CFG)["g"]
    assert hit.grade == "PROVISIONAL"   # 50 bits < 0.05 * 10000


def test_the_evalue_is_still_reported(tmp_path):
    """It is no longer the decision, but it stays in the output."""
    tbl = _tbl(tmp_path, [("g", "A", "3.2e-40", "800.0")])
    assert _parse_hmmsearch_tblout(tbl, CFG.hmmer_evalue, CFG)["g"].evalue == 3.2e-40
