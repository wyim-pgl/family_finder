"""Unit tests for build_orthofinder_cmd (issue #10) — pure argv assembly,
no OrthoFinder needed."""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from config import Config  # noqa: E402
from steps.orthofinder import build_orthofinder_cmd  # noqa: E402

IN = Path("/in")
OUT = Path("/out")


def test_default_cmd_has_inflation_and_og_no_search_program():
    # Arrange
    cfg = Config()
    # Act
    cmd = build_orthofinder_cmd(IN, OUT, cfg)
    # Assert
    assert cmd[1:5] == ["-f", "/in", "-o", "/out"]
    assert cmd[cmd.index("-I") + 1] == "1.2"  # OrthoFinder v3 default, explicit
    assert "-og" in cmd
    assert "-S" not in cmd
    assert cmd[cmd.index("-t") + 1] == "8"


def test_search_program_passed_when_set():
    cfg = Config(orthofinder_search_program="diamond_ultra_sens")
    cmd = build_orthofinder_cmd(IN, OUT, cfg)
    assert cmd[cmd.index("-S") + 1] == "diamond_ultra_sens"


def test_stop_after_orthogroups_off_drops_og():
    cfg = Config(orthofinder_stop_after_orthogroups=False)
    assert "-og" not in build_orthofinder_cmd(IN, OUT, cfg)


def test_blast_reuse_uses_b_and_drops_f_o():
    # -o is rejected by OrthoFinder when combined with -b
    cfg = Config(orthofinder_reuse_blast_dir="/wd/WorkingDirectory")
    cmd = build_orthofinder_cmd(IN, OUT, cfg)
    assert cmd[1:3] == ["-b", "/wd/WorkingDirectory"]
    assert "-f" not in cmd and "-o" not in cmd
    assert "-I" in cmd  # inflation still sweepable in reuse mode


def test_extra_args_blocklist_rejects_typed_flags():
    for bad in ["-I 1.5", "-S blast", "-b /x", "-og", "-t 4", "-o /x", "-f /x"]:
        cfg = Config(orthofinder_extra_args=bad)
        with pytest.raises(ValueError, match="blocked flag"):
            build_orthofinder_cmd(IN, OUT, cfg)


def test_extra_args_appended_last():
    cfg = Config(orthofinder_extra_args="-X 5")
    cmd = build_orthofinder_cmd(IN, OUT, cfg)
    assert cmd[-2:] == ["-X", "5"]


def test_inflation_sweep_values_render():
    for i in [1.05, 1.1, 1.15, 1.3]:
        cfg = Config(orthofinder_inflation=i)
        cmd = build_orthofinder_cmd(IN, OUT, cfg)
        assert cmd[cmd.index("-I") + 1] == str(i)
