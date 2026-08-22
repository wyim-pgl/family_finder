"""A missing tool must fail by name, not as a bare OSError from subprocess.

Issue #44. `codeml_bin` defaulted to the bare name "codeml", nothing on the
pipeline's PATH provided one, and no configuration pinned a path - so the
selection analysis could not run at all and nothing said so until a subprocess
died. Worse, a bare name means "whichever binary a shell happens to expose",
which puts irreproducibility in the configuration itself: seven PAML
installations spanning 4.9i to 4.10.10 exist across the two hosts and codeml
prints no version banner, so a completed run cannot be attributed afterwards.
"""
import os
import stat
import sys
import types

import pytest

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import Config, ToolNotFoundError, resolve_tool


def _executable(tmp_path, name):
    path = tmp_path / name
    path.write_text("#!/bin/sh\nexit 0\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def test_an_absolute_path_that_exists_resolves_to_itself(tmp_path):
    tool = _executable(tmp_path, "codeml")
    assert resolve_tool(str(tool), "codeml_bin") == str(tool)


def test_a_missing_absolute_path_names_the_config_key(tmp_path):
    with pytest.raises(ToolNotFoundError) as excinfo:
        resolve_tool(str(tmp_path / "nope"), "codeml_bin")
    assert "codeml_bin" in str(excinfo.value)


def test_a_bare_name_resolves_through_path(tmp_path, monkeypatch):
    _executable(tmp_path, "codeml")
    monkeypatch.setenv("PATH", str(tmp_path))
    assert resolve_tool("codeml", "codeml_bin") == str(tmp_path / "codeml")


def test_a_bare_name_that_is_not_on_path_fails_by_name(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(ToolNotFoundError) as excinfo:
        resolve_tool("codeml", "codeml_bin")
    message = str(excinfo.value)
    assert "codeml_bin" in message and "PATH" in message


def test_the_error_says_a_bare_name_is_not_reproducible(tmp_path, monkeypatch):
    """The message has to explain WHY pinning matters, not just that it failed.

    Someone who hits this will otherwise 'fix' it by putting any codeml on
    PATH, which reintroduces exactly the ambiguity issue #44 records.
    """
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(ToolNotFoundError) as excinfo:
        resolve_tool("codeml", "codeml_bin")
    assert "version" in str(excinfo.value).lower()


def test_config_still_carries_the_bare_default():
    """The default stays a bare name so nothing breaks for callers that do
    have the tool on PATH; resolution is what turns it into a failure."""
    assert Config().codeml_bin == "codeml"


# ---------------------------------------------------------------------------
# run_codeml: judge by output, not by exit status
# ---------------------------------------------------------------------------
# codeml writes its results and THEN exits 1 with "error: end of tree file"
# (reproduced three times; recorded in resume.md). Raising on the return code
# therefore discards completed analyses, and a caller that catches the raise
# reads a finished run as a failed one.

import textwrap

from steps.codeml import run_codeml

COMPLETE = textwrap.dedent("""\
    seed used = 1
    lnL(ntime: 20  np: 24):  -49044.147942      +0.000000
    Bayes Empirical Bayes (BEB) analysis
        11 S    0.987*
    """)


def _fake_codeml(tmp_path, exit_code, results_text):
    script = tmp_path / "codeml"
    script.write_text(
        "#!/bin/sh\n"
        f"cat <<'EOT' > results.txt\n{results_text}EOT\n"
        f"echo 'error: end of tree file' >&2\n"
        f"exit {exit_code}\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return script


def test_a_complete_run_that_exits_nonzero_is_accepted(tmp_path):
    work = tmp_path / "w"
    work.mkdir()
    (work / "codeml.ctl").write_text("seqfile = x\n")
    cfg = Config()
    cfg.codeml_bin = str(_fake_codeml(tmp_path, 1, COMPLETE))
    assert run_codeml(work / "codeml.ctl", work, cfg).exists()


def test_a_truncated_run_that_exits_nonzero_is_rejected(tmp_path):
    work = tmp_path / "w"
    work.mkdir()
    (work / "codeml.ctl").write_text("seqfile = x\n")
    cfg = Config()
    cfg.codeml_bin = str(_fake_codeml(tmp_path, 1, "seed used = 1\n"))
    with pytest.raises(RuntimeError) as excinfo:
        run_codeml(work / "codeml.ctl", work, cfg)
    assert "lnL" in str(excinfo.value)


def test_a_missing_codeml_fails_by_config_key(tmp_path):
    work = tmp_path / "w"
    work.mkdir()
    (work / "codeml.ctl").write_text("seqfile = x\n")
    cfg = Config()
    cfg.codeml_bin = str(tmp_path / "absent")
    with pytest.raises(ToolNotFoundError) as excinfo:
        run_codeml(work / "codeml.ctl", work, cfg)
    assert "codeml_bin" in str(excinfo.value)
