"""Chunked ESMFold structure prediction for tier-3 (issue #34).

ESMFold was chosen over ProstT5 because tier3_assign gates on alntmscore, which
a ProstT5 database cannot produce. These tests pin the two properties measured
on 2026-08-22 that decide how the screen is run: cost rises steeply with length
(t = 1.68 s at 210 aa, exponent 2.55), so a handful of very long sequences
dominates the bill; and a truncated or empty PDB must fail loudly rather than
silently reducing coverage.

No GPU, no model, no ESMFold — the runner is monkeypatched and the parsers are
fed synthetic PDB text.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from steps import esmfold_chunks as ef


def atom(serial, resseq, plddt, name="CA", res="ALA"):
    return (f"ATOM  {serial:>5}  {name:<3} {res} A{resseq:>4}    "
            f"   1.000   2.000   3.000  1.00{plddt:>6.2f}           C")


def pdb(residues, plddt=85.0):
    lines = [atom(i, i, plddt) for i in range(1, residues + 1)]
    return "\n".join(lines) + "\nTER\nEND\n"


# ------------------------------------------------------ length partition ----

def test_sequences_over_the_cap_are_separated_rather_than_dropped():
    seqs = {"a": "M" * 200, "big": "M" * 1500, "b": "M" * 900}

    foldable, oversize = ef.partition_by_length(seqs, max_len=1100)

    assert sorted(foldable) == ["a", "b"]
    assert oversize == {"big": 1500}


def test_a_sequence_exactly_at_the_cap_is_foldable():
    foldable, oversize = ef.partition_by_length({"a": "M" * 1100}, max_len=1100)

    assert list(foldable) == ["a"]
    assert oversize == {}


# ----------------------------------------------------------- cost model ----

def test_runtime_estimate_reproduces_the_measured_anchor():
    # 20 sequences of 210 aa measured at 1.68 s each
    assert ef.estimate_runtime([210] * 20) == pytest.approx(33.6, rel=0.02)


def test_runtime_estimate_is_superlinear_in_length():
    single = ef.estimate_runtime([210])
    double = ef.estimate_runtime([420])

    assert double > 4 * single  # exponent 2.55 => ~5.9x


def test_oversize_sequences_dominate_the_bill_and_the_cap_shows_it():
    # 100 short plus one very long: the cap must change the estimate materially
    lengths = [210] * 100 + [3956]

    uncapped = ef.estimate_runtime(lengths, max_len=None)
    capped = ef.estimate_runtime(lengths, max_len=1100)

    assert uncapped > 5 * capped


# --------------------------------------------------------- pdb checking ----

def test_verify_pdb_returns_mean_plddt_from_the_b_factor_column(tmp_path):
    p = tmp_path / "g.pdb"
    p.write_text(pdb(3, plddt=80.0).replace("  1.00 80.00", "  1.00 90.00", 1))

    mean = ef.verify_pdb(p)

    assert mean == pytest.approx((90.0 + 80.0 + 80.0) / 3)


def test_verify_pdb_rejects_a_missing_file(tmp_path):
    with pytest.raises(ValueError, match="missing"):
        ef.verify_pdb(tmp_path / "nope.pdb")


def test_verify_pdb_rejects_a_file_with_no_atom_records(tmp_path):
    p = tmp_path / "g.pdb"
    p.write_text("HEADER something\nEND\n")

    with pytest.raises(ValueError, match="no ATOM"):
        ef.verify_pdb(p)


def test_verify_pdb_rejects_an_empty_file(tmp_path):
    p = tmp_path / "g.pdb"
    p.write_text("")

    with pytest.raises(ValueError, match="empty"):
        ef.verify_pdb(p)


def test_verify_pdb_rejects_out_of_range_plddt(tmp_path):
    # a B-factor column holding something other than pLDDT is a silent
    # wrong-format failure, not a crash
    p = tmp_path / "g.pdb"
    p.write_text(pdb(2, plddt=999.0))

    with pytest.raises(ValueError, match="pLDDT"):
        ef.verify_pdb(p)


def test_verify_pdb_rejects_a_truncated_structure(tmp_path):
    p = tmp_path / "g.pdb"
    p.write_text(pdb(5))

    with pytest.raises(ValueError, match="residues"):
        ef.verify_pdb(p, expected_residues=10)


def test_verify_pdb_accepts_the_expected_residue_count(tmp_path):
    p = tmp_path / "g.pdb"
    p.write_text(pdb(10))

    assert ef.verify_pdb(p, expected_residues=10) == pytest.approx(85.0)


# ------------------------------------------------------------ folding ------

def test_fold_chunk_verifies_every_structure_and_writes_a_sentinel(tmp_path, monkeypatch):
    seqs = {"a": "M" * 5, "b": "M" * 7}
    out = tmp_path / "pdb"

    def fake_run(sequences, out_dir, **kw):
        out_dir.mkdir(parents=True, exist_ok=True)
        for name, s in sequences.items():
            (out_dir / f"{name}.pdb").write_text(pdb(len(s)))

    monkeypatch.setattr(ef, "run_esmfold", fake_run)

    report = ef.fold_chunk(seqs, out)

    assert report["n_folded"] == 2
    assert (out / "chunk.done").exists()
    assert report["plddt"]["a"] == pytest.approx(85.0)


def test_fold_chunk_fails_hard_when_a_structure_is_missing(tmp_path, monkeypatch):
    seqs = {"a": "M" * 5, "b": "M" * 7}
    out = tmp_path / "pdb"

    def fake_partial(sequences, out_dir, **kw):
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "a.pdb").write_text(pdb(5))  # b never written

    monkeypatch.setattr(ef, "run_esmfold", fake_partial)

    with pytest.raises(ValueError, match="missing"):
        ef.fold_chunk(seqs, out)
    assert not (out / "chunk.done").exists()


def test_fold_chunk_skips_structures_that_already_exist(tmp_path, monkeypatch):
    seqs = {"a": "M" * 5}
    out = tmp_path / "pdb"
    out.mkdir()
    (out / "a.pdb").write_text(pdb(5))
    called = []

    monkeypatch.setattr(ef, "run_esmfold",
                        lambda s, d, **k: called.append(sorted(s)))

    report = ef.fold_chunk(seqs, out)

    assert called == [[]]           # nothing left to fold
    assert report["n_skipped"] == 1


def test_completed_chunks_are_recognised_by_their_sentinel(tmp_path):
    a, b = tmp_path / "c0", tmp_path / "c1"
    for d in (a, b):
        d.mkdir()
    (a / "chunk.done").write_text("")

    assert ef.completed_chunks([a, b]) == [a]
