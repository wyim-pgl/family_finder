"""parse_orthogroups must not let one gene ID reach two orthogroups.

User-reported failure mode: after OrthoFinder clustering, a gene ID can
appear more than once (across orthogroups or twice inside one cell). A
duplicated ID gets processed in two OGs and can end up confirmed in two
families — the per-round conservation audit checks placed-or-recycled,
not uniqueness, so it never catches this. The parser is the choke point:
keep the first occurrence, drop the rest, and log loudly.
"""

import logging
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from steps.orthofinder import parse_orthogroups  # noqa: E402


def _write_tsv(tmp_path, rows):
    og_dir = tmp_path / "Orthogroups"
    og_dir.mkdir()
    lines = ["Orthogroup\tSp1\tSp2"] + rows
    (og_dir / "Orthogroups.tsv").write_text("\n".join(lines) + "\n")
    return tmp_path


def test_clean_file_parses_unchanged(tmp_path):
    results = _write_tsv(tmp_path, [
        "OG0000001\tSp1_a, Sp1_b\tSp2_c",
        "OG0000002\tSp1_d\t",
    ])
    ogs = parse_orthogroups(results)
    assert ogs == {
        "OG0000001": ["Sp1_a", "Sp1_b", "Sp2_c"],
        "OG0000002": ["Sp1_d"],
    }


def test_gene_in_two_orthogroups_kept_only_in_first(tmp_path, caplog):
    results = _write_tsv(tmp_path, [
        "OG0000001\tSp1_a, Sp1_b\t",
        "OG0000002\tSp1_b\tSp2_c",   # Sp1_b again
    ])
    with caplog.at_level(logging.ERROR, logger="family_finder"):
        ogs = parse_orthogroups(results)
    assert ogs["OG0000001"] == ["Sp1_a", "Sp1_b"]
    assert ogs["OG0000002"] == ["Sp2_c"]  # duplicate dropped, OG survives
    assert any("Sp1_b" in r.message for r in caplog.records)


def test_duplicate_within_one_orthogroup_collapsed(tmp_path, caplog):
    results = _write_tsv(tmp_path, [
        "OG0000001\tSp1_a, Sp1_a\tSp2_c",
    ])
    with caplog.at_level(logging.ERROR, logger="family_finder"):
        ogs = parse_orthogroups(results)
    assert ogs["OG0000001"] == ["Sp1_a", "Sp2_c"]


def test_orthogroup_left_empty_by_dedup_is_dropped(tmp_path, caplog):
    results = _write_tsv(tmp_path, [
        "OG0000001\tSp1_a\t",
        "OG0000002\tSp1_a\t",   # only member is a duplicate
    ])
    with caplog.at_level(logging.ERROR, logger="family_finder"):
        ogs = parse_orthogroups(results)
    assert "OG0000002" not in ogs
