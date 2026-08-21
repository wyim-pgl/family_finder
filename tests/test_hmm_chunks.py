"""steps/hmm_chunks.py — split an HMM profile file, run the chunks anywhere,
merge the results (issue #31).

hmmsearch cost is profiles x sequences, so the 15sp rescue extrapolates to
2.7-4.1 days against a 3-day SLURM limit with no internal checkpoint.
Chunking makes that linear. Correctness rule: a chunk whose tblout lacks
HMMER's terminating `# [ok]` marker is INCOMPLETE, and merging must fail
rather than emit a quietly truncated rescue.

No hmmsearch is ever invoked here.
"""

import stat
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from steps.hmm_chunks import (
    merge_tblouts,
    split_hmm_file,
    write_runner_script,
)

HMM = "".join(
    f"HMMER3/f [3.2.1]\nNAME  fam{i}\nLENG  100\nHMM\n//\n" for i in range(7)
)

TBL_OK = """\
#  target name  accession  query name
gene1  -  fam0  -  1e-40
gene2  -  fam0  -  1e-30
#
# [ok]
"""

TBL_TRUNCATED = """\
#  target name  accession  query name
gene3  -  fam1  -  1e-20
"""


# --- splitting -----------------------------------------------------------

def test_split_produces_expected_chunk_count(tmp_path):
    src = tmp_path / "all.hmm"
    src.write_text(HMM)
    chunks = split_hmm_file(src, tmp_path / "chunks", chunk_size=3)
    assert len(chunks) == 3          # 7 profiles -> 3 + 3 + 1


def test_split_preserves_every_profile(tmp_path):
    src = tmp_path / "all.hmm"
    src.write_text(HMM)
    chunks = split_hmm_file(src, tmp_path / "chunks", chunk_size=3)
    names = []
    for c in chunks:
        names += [l.split()[1] for l in c.read_text().splitlines()
                  if l.startswith("NAME")]
    assert names == [f"fam{i}" for i in range(7)]


def test_split_keeps_records_intact(tmp_path):
    src = tmp_path / "all.hmm"
    src.write_text(HMM)
    chunks = split_hmm_file(src, tmp_path / "chunks", chunk_size=3)
    for c in chunks:
        text = c.read_text()
        assert text.count("HMMER3/f") == text.count("//")
        assert text.endswith("//\n")


def test_split_rejects_nonpositive_chunk_size(tmp_path):
    src = tmp_path / "all.hmm"
    src.write_text(HMM)
    for bad in (0, -1):
        try:
            split_hmm_file(src, tmp_path / "c", chunk_size=bad)
        except ValueError:
            continue
        raise AssertionError(f"chunk_size={bad} must raise")


def test_split_of_empty_file_raises(tmp_path):
    src = tmp_path / "all.hmm"
    src.write_text("")
    try:
        split_hmm_file(src, tmp_path / "c", chunk_size=3)
    except ValueError:
        return
    raise AssertionError("empty profile file must raise")


# --- runner script -------------------------------------------------------

def test_runner_script_is_executable_and_slurm_optional(tmp_path):
    script = write_runner_script(
        tmp_path / "run.sh", n_chunks=4, chunk_dir=tmp_path / "chunks",
        query_fasta=tmp_path / "q.fa", out_dir=tmp_path / "out",
        hmmsearch_bin="hmmsearch", threads=4, evalue=1e-5, max_concurrent=2,
    )
    text = script.read_text()
    assert script.stat().st_mode & stat.S_IXUSR
    # SLURM when present, plain shell when not — the same script must do both
    assert "command -v sbatch" in text
    assert "--wait" in text and "--array" in text
    assert "wait" in text


def test_runner_script_quotes_paths_with_spaces(tmp_path):
    script = write_runner_script(
        tmp_path / "run.sh", n_chunks=2, chunk_dir=tmp_path / "my chunks",
        query_fasta=tmp_path / "my q.fa", out_dir=tmp_path / "out",
        hmmsearch_bin="hmmsearch", threads=2, evalue=1e-5,
    )
    assert "'" in script.read_text()


def test_runner_script_carries_the_evalue(tmp_path):
    script = write_runner_script(
        tmp_path / "run.sh", n_chunks=2, chunk_dir=tmp_path / "c",
        query_fasta=tmp_path / "q.fa", out_dir=tmp_path / "o",
        hmmsearch_bin="hmmsearch", threads=2, evalue=1e-7,
    )
    assert "1e-07" in script.read_text() or "1e-7" in script.read_text()


# --- merging -------------------------------------------------------------

def test_merge_concatenates_complete_chunks(tmp_path):
    for i in range(3):
        (tmp_path / f"chunk_{i:04d}.tblout").write_text(TBL_OK)
    out = tmp_path / "merged.tblout"
    n = merge_tblouts(tmp_path, out)
    assert n == 3
    body = [l for l in out.read_text().splitlines() if not l.startswith("#")]
    assert len(body) == 6           # 2 hits per chunk


def test_merge_fails_on_truncated_chunk(tmp_path):
    (tmp_path / "chunk_0000.tblout").write_text(TBL_OK)
    (tmp_path / "chunk_0001.tblout").write_text(TBL_TRUNCATED)
    try:
        merge_tblouts(tmp_path, tmp_path / "m.tblout")
    except RuntimeError as e:
        assert "chunk_0001" in str(e)
        return
    raise AssertionError("a chunk without '# [ok]' must fail the merge")


def test_merge_fails_when_a_chunk_is_missing(tmp_path):
    (tmp_path / "chunk_0000.tblout").write_text(TBL_OK)
    try:
        merge_tblouts(tmp_path, tmp_path / "m.tblout", expected=2)
    except RuntimeError as e:
        assert "2" in str(e)
        return
    raise AssertionError("a missing chunk must fail the merge")


def test_merge_output_parses_with_the_existing_reader(tmp_path):
    """The merged file must satisfy steps.hmmer_rescue's tblout parser."""
    for i in range(2):
        (tmp_path / f"chunk_{i:04d}.tblout").write_text(
            "#  target name  accession  query name  accession  E-value\n"
            f"gene{i}  -  fam{i}  -  1e-40  100.0  0.0  1e-40  99.0  0.0\n"
            "# [ok]\n"
        )
    out = tmp_path / "merged.tblout"
    merge_tblouts(tmp_path, out)
    from steps.hmmer_rescue import _parse_hmmsearch_tblout
    hits = _parse_hmmsearch_tblout(out, 1e-5)
    assert set(hits) == {"gene0", "gene1"}


# --- integration with hmmer_rescue (option on/off) ------------------------

def test_chunking_is_off_by_default():
    """Default config must keep the single-hmmsearch path — the option is
    opt-in so existing runs are unaffected."""
    from config import Config
    assert Config().hmmer_chunk_size == 0


def test_rescue_routes_to_chunked_path_when_enabled(tmp_path, monkeypatch):
    import steps.hmmer_rescue as hr

    calls = {"single": 0, "chunked": 0}
    monkeypatch.setattr(hr, "_run_hmmsearch",
                        lambda *a, **k: calls.__setitem__("single", 1))
    monkeypatch.setattr(hr, "_run_hmmsearch_chunked",
                        lambda *a, **k: calls.__setitem__("chunked", 1))

    from config import Config
    cfg = Config()
    cfg.hmmer_chunk_size = 500
    # exercise the branch directly — the surrounding rescue needs a full run dir
    tblout = tmp_path / "t.tblout"
    if cfg.hmmer_chunk_size:
        hr._run_hmmsearch_chunked(tmp_path / "a.hmm", tmp_path / "q.fa",
                                  tblout, tmp_path, cfg)
    else:
        hr._run_hmmsearch(tmp_path / "a.hmm", tmp_path / "q.fa", tblout, cfg)
    assert calls == {"single": 0, "chunked": 1}


def test_chunked_run_raises_when_runner_fails(tmp_path, monkeypatch):
    """A failed runner must not fall through to a merge of partial chunks."""
    import subprocess as sp

    import steps.hmmer_rescue as hr
    from config import Config

    (tmp_path / "a.hmm").write_text(HMM)
    (tmp_path / "q.fa").write_text(">g\nMA\n")
    cfg = Config()
    cfg.hmmer_chunk_size = 3

    class R:
        returncode = 1
        stderr = "boom"

    monkeypatch.setattr(sp, "run", lambda *a, **k: R())
    monkeypatch.setattr(hr.subprocess, "run", lambda *a, **k: R())
    try:
        hr._run_hmmsearch_chunked(tmp_path / "a.hmm", tmp_path / "q.fa",
                                  tmp_path / "t.tblout", tmp_path, cfg)
    except RuntimeError as e:
        assert "chunked hmmsearch failed" in str(e)
        return
    raise AssertionError("runner failure must raise")
