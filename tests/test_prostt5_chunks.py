"""steps/prostt5_chunks.py — chunked ProstT5 3Di DB construction (issue #34).

ProstT5 conversion runs at ~4.3 s/sequence on CPU and 0.034 s/seq with
`--gpu 1` on a -DENABLE_CUDA=1 build (2026-08-25 measurement; the flag decides
the path there), and createdb has no internal checkpoint, so a genome-scale
screen must be chunked.

The correctness rule is the trap this module exists for: **passing a `.gguf`
FILE to `--prostt5-model` produces no 3Di at all.** createdb exits in 0.00 s
with no error and no warning, and the resulting database holds amino acids.
A screen run that way is a sequence search wearing a structure search's name.
Checking `ls <db>_ss` by hand does not survive dozens of chunks, so the check
belongs in code and must fail hard.

No foldseek is ever invoked here.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from steps.prostt5_chunks import (
    build_chunk_db,
    merge_search_results,
    split_fasta,
    verify_3di_db,
)

FASTA = (">a\nMKTAY\nIAKQR\n>b\nMKV\n>c\nMKTAYIAKQRQISFVK\n"
         ">d\nMM\n>e\nMKTA\n")


def _db(tmp_path, name="db", ss="DDPLVVVLLVVVDDD", main="MDAWETDLLIDYGGS",
        ss_index="0\t0\t16\n1\t16\t16\n", index="0\t0\t16\n1\t16\t16\n"):
    """Write a foldseek-shaped DB pair. ss=None omits the 3Di database."""
    base = tmp_path / name
    base.write_text(main)
    (tmp_path / f"{name}.index").write_text(index)
    if ss is not None:
        (tmp_path / f"{name}_ss").write_text(ss)
        (tmp_path / f"{name}_ss.index").write_text(ss_index)
    return base


# ---------------------------------------------------------------- the guard --

def test_missing_3di_database_is_the_gguf_trap_and_must_raise(tmp_path):
    # exactly what `--prostt5-model <file>.gguf` leaves behind
    db = _db(tmp_path, ss=None)

    with pytest.raises(RuntimeError, match="_ss"):
        verify_3di_db(db)


def test_empty_3di_database_raises(tmp_path):
    db = _db(tmp_path, ss="")

    with pytest.raises(RuntimeError, match="empty"):
        verify_3di_db(db)


def test_3di_database_identical_to_the_amino_acid_database_raises(tmp_path):
    # a build that copied residues into _ss is still not a structure search
    db = _db(tmp_path, ss="MDAWETDLLIDYGGS", main="MDAWETDLLIDYGGS")

    with pytest.raises(RuntimeError, match="amino acid"):
        verify_3di_db(db)


def test_3di_entry_count_must_match_the_sequence_database(tmp_path):
    db = _db(tmp_path, ss_index="0\t0\t16\n")  # one entry vs two

    with pytest.raises(RuntimeError, match="entries"):
        verify_3di_db(db)


def test_a_well_formed_pair_passes_and_returns_the_entry_count(tmp_path):
    db = _db(tmp_path)

    assert verify_3di_db(db) == 2


def test_missing_3di_index_raises(tmp_path):
    db = _db(tmp_path)
    (tmp_path / "db_ss.index").unlink()

    with pytest.raises(RuntimeError, match="_ss.index"):
        verify_3di_db(db)


# ------------------------------------------------------------ fasta split --

def test_split_fasta_puts_chunk_size_records_in_each_chunk(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(FASTA)

    chunks = split_fasta(src, tmp_path / "chunks", chunk_size=2)

    assert [p.name for p in chunks] == ["chunk_0001.fa", "chunk_0002.fa",
                                        "chunk_0003.fa"]
    assert chunks[0].read_text().count(">") == 2
    assert chunks[2].read_text().count(">") == 1


def test_split_fasta_never_cuts_a_record_in_half(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(FASTA)

    chunks = split_fasta(src, tmp_path / "chunks", chunk_size=2)

    rejoined = "".join(p.read_text() for p in chunks)
    assert rejoined.count(">") == FASTA.count(">")
    # record "a" is wrapped across two lines; both must land in the same chunk
    assert ">a\nMKTAY\nIAKQR\n" in chunks[0].read_text()
    assert not any(">" in p.read_text().split("\n")[-2] for p in chunks)


def test_split_fasta_rejects_a_nonpositive_chunk_size(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(FASTA)

    with pytest.raises(ValueError):
        split_fasta(src, tmp_path / "chunks", chunk_size=0)


def test_split_fasta_rejects_an_input_with_no_records(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text("\n")

    with pytest.raises(ValueError, match="no sequences"):
        split_fasta(src, tmp_path / "chunks", chunk_size=2)


# ------------------------------------------------------------- build path --

def test_build_chunk_db_verifies_3di_after_running_createdb(tmp_path, monkeypatch):
    src = tmp_path / "c.fa"
    src.write_text(">a\nMK\n")
    calls = []

    def fake_createdb(fasta, db, model_dir, foldseek_bin, threads, gpu=False):
        calls.append((Path(fasta).name, Path(db).name, str(model_dir), threads))
        _db(tmp_path, name=Path(db).name)

    monkeypatch.setattr("steps.prostt5_chunks.run_createdb", fake_createdb)

    out = build_chunk_db(src, tmp_path / "db", "/w/prostt5_weights",
                         foldseek_bin="foldseek", threads=16)

    assert calls == [("c.fa", "db", "/w/prostt5_weights", 16)]
    assert out == tmp_path / "db"


def test_build_chunk_db_fails_loudly_when_the_chunk_produced_no_3di(tmp_path, monkeypatch):
    src = tmp_path / "c.fa"
    src.write_text(">a\nMK\n")

    def fake_createdb(fasta, db, model_dir, foldseek_bin, threads, gpu=False):
        _db(tmp_path, name=Path(db).name, ss=None)  # the trap

    monkeypatch.setattr("steps.prostt5_chunks.run_createdb", fake_createdb)

    with pytest.raises(RuntimeError, match="db"):
        build_chunk_db(src, tmp_path / "db", "/w/prostt5_weights",
                       foldseek_bin="foldseek", threads=16)


def test_build_chunk_db_rejects_a_model_path_that_is_a_file(tmp_path, monkeypatch):
    # the trap at its source: --prostt5-model must be the weights DIRECTORY
    src = tmp_path / "c.fa"
    src.write_text(">a\nMK\n")
    gguf = tmp_path / "prostt5-f16.gguf"
    gguf.write_text("weights")
    monkeypatch.setattr("steps.prostt5_chunks.run_createdb",
                        lambda *a, **k: pytest.fail("createdb must not be reached"))

    with pytest.raises(ValueError, match="directory"):
        build_chunk_db(src, tmp_path / "db", gguf,
                       foldseek_bin="foldseek", threads=16)


# ----------------------------------------------------------------- merge --

def _result(tmp_path, name, text="q\tt\t1e-9\t100\n", done=True):
    p = tmp_path / name
    p.write_text(text)
    if done:
        (tmp_path / f"{name}.done").write_text("ok\n")
    return p


def test_merge_concatenates_completed_chunk_results(tmp_path):
    _result(tmp_path, "chunk_0001.m8", "a\tb\t1e-9\t100\n")
    _result(tmp_path, "chunk_0002.m8", "c\td\t1e-8\t90\n")
    out = tmp_path / "all.m8"

    assert merge_search_results(tmp_path, out) == 2
    assert out.read_text() == "a\tb\t1e-9\t100\nc\td\t1e-8\t90\n"


def test_merge_refuses_a_chunk_with_no_completion_sentinel(tmp_path):
    # foldseek writes no terminator of its own, so a killed task leaves a
    # valid-looking partial file; the sentinel is what distinguishes them
    _result(tmp_path, "chunk_0001.m8")
    _result(tmp_path, "chunk_0002.m8", done=False)

    with pytest.raises(RuntimeError, match="incomplete"):
        merge_search_results(tmp_path, tmp_path / "all.m8")


def test_merge_refuses_when_fewer_chunks_than_expected_are_present(tmp_path):
    _result(tmp_path, "chunk_0001.m8")

    with pytest.raises(RuntimeError, match="expected 2"):
        merge_search_results(tmp_path, tmp_path / "all.m8", expected=2)


def test_merge_refuses_an_empty_directory(tmp_path):
    with pytest.raises(RuntimeError, match="no chunk"):
        merge_search_results(tmp_path, tmp_path / "all.m8")


def test_a_chunk_that_legitimately_found_nothing_still_merges(tmp_path):
    # an empty result set is not the same as a truncated one
    _result(tmp_path, "chunk_0001.m8", text="")
    _result(tmp_path, "chunk_0002.m8", "c\td\t1e-8\t90\n")
    out = tmp_path / "all.m8"

    assert merge_search_results(tmp_path, out) == 2
    assert out.read_text() == "c\td\t1e-8\t90\n"


def test_run_createdb_gpu_flag_reaches_the_command(monkeypatch):
    # Without --gpu 1 a CUDA build silently takes the CPU ProstT5ForkRunner
    # (which can hang forever in msgsnd); the flag must reach the argv.
    import steps.prostt5_chunks as pc

    captured = {}

    def fake_run(cmd, capture_output, text):
        captured["cmd"] = cmd

        class R:
            returncode = 0
            stderr = ""
        return R()

    monkeypatch.setattr(pc.subprocess, "run", fake_run)

    pc.run_createdb("in.fa", "db", "/w/weights", "foldseek", 16, gpu=True)
    gpu_cmd = captured["cmd"]
    assert gpu_cmd[gpu_cmd.index("--gpu") + 1] == "1"

    pc.run_createdb("in.fa", "db", "/w/weights", "foldseek", 16)
    assert "--gpu" not in captured["cmd"]


def test_config_entry_point_forwards_the_device_choice(tmp_path, monkeypatch):
    # The reviewer's point: a dead `gpu` parameter changes nothing. The
    # config entry point must forward config.prostt5_gpu into the argv.
    import steps.prostt5_chunks as pc
    from config import Config

    src = tmp_path / "c.fa"
    src.write_text(">a\nMK\n")
    seen = {}

    def fake_createdb(fasta, db, model_dir, foldseek_bin, threads, gpu=False):
        seen["gpu"] = gpu
        seen["threads"] = threads
        seen["bin"] = foldseek_bin
        _db(tmp_path, name=Path(db).name)

    monkeypatch.setattr("steps.prostt5_chunks.run_createdb", fake_createdb)

    pc.build_chunk_db_from_config(src, tmp_path / "db", "/w/prostt5_weights",
                                  Config(n_workers=7))
    assert seen == {"gpu": True, "threads": 7, "bin": "foldseek"}

    pc.build_chunk_db_from_config(src, tmp_path / "db2", "/w/prostt5_weights",
                                  Config(prostt5_gpu=False))
    assert seen["gpu"] is False
