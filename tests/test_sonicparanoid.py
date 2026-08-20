"""Tests for steps/sonicparanoid.py (issue #22 adapter): groups parsing
(incl. "*" empty cells and metadata columns), OrthoFinder-format conversion
round-trip through steps.orthofinder.parse_orthogroups, the template
runner, and config defaults.

No sonicparanoid required: the template runner is monkeypatched and the
parser is fed a synthetic ortholog_groups.tsv fixture (house convention).
"""

import os
import sys
import time
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from config import Config  # noqa: E402
from steps import sonicparanoid  # noqa: E402
from steps.orthofinder import parse_orthogroups  # noqa: E402
from steps.sonicparanoid import (  # noqa: E402
    find_groups_file,
    parse_groups,
    run_sonicparanoid,
    write_orthogroups_tsv,
)


# ---------------------------------------------------------------------------
# Fixtures (documented SonicParanoid2 ortholog_groups.tsv layout)
# ---------------------------------------------------------------------------

GROUPS_FILE_EMBEDDED_HEADERS = (
    "group_id\tgroup_size\tMcry.faa\tOcoc.faa\tCgig.faa\n"
    "1\t4\tMcry_g1,Mcry_g2\tOcoc_g1\tCgig_g1\n"
    "2\t2\t*\tOcoc_g2,Ococ_g3\t*\n"
    "3\t1\t*\t*\tCgig_g2\n"
)

GROUPS_FILE_PLAIN_HEADERS = (
    "group_id\tMcry\tOcoc\n"
    "10\tMcry_g9\t*\n"
)


def _write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


# ---------------------------------------------------------------------------
# parse_groups
# ---------------------------------------------------------------------------

def test_parse_groups_embedded_faa_headers_and_star_cells(tmp_path):
    path = _write(tmp_path / "ortholog_groups.tsv", GROUPS_FILE_EMBEDDED_HEADERS)
    groups = parse_groups(path)
    assert groups == {
        "1": ["Mcry_g1", "Mcry_g2", "Ococ_g1", "Cgig_g1"],
        "2": ["Ococ_g2", "Ococ_g3"],
        "3": ["Cgig_g2"],
    }


def test_parse_groups_skips_metadata_columns(tmp_path):
    # group_size is a known metadata column: its integer values must never
    # be ingested as gene ids.
    path = _write(tmp_path / "ortholog_groups.tsv", GROUPS_FILE_EMBEDDED_HEADERS)
    all_genes = [g for genes in parse_groups(path).values() for g in genes]
    assert "4" not in all_genes and "2" not in all_genes


def test_parse_groups_plain_headers_and_verbatim_ids(tmp_path):
    path = _write(tmp_path / "ortholog_groups.tsv", GROUPS_FILE_PLAIN_HEADERS)
    assert parse_groups(path) == {"10": ["Mcry_g9"]}  # nothing stripped


def test_parse_groups_all_star_row_dropped(tmp_path):
    text = "group_id\tMcry\tOcoc\n7\t*\t*\n"
    path = _write(tmp_path / "ortholog_groups.tsv", text)
    assert parse_groups(path) == {}


def test_parse_groups_empty_header_raises(tmp_path):
    path = _write(tmp_path / "ortholog_groups.tsv", "\n")
    with pytest.raises(ValueError, match="Empty header"):
        parse_groups(path)


# ---------------------------------------------------------------------------
# write_orthogroups_tsv: OrthoFinder-format conversion round-trip
# ---------------------------------------------------------------------------

SPECIES_ORDER = ["Cgig", "Mcry", "Ococ"]
GROUPS = {
    "1": ["Mcry_g1", "Mcry_g2", "Ococ_g1", "Cgig_g1"],
    "2": ["Ococ_g2", "Ococ_g3"],
}
ALL_GENES = ["Mcry_g1", "Mcry_g2", "Mcry_g3", "Ococ_g1", "Ococ_g2",
             "Ococ_g3", "Cgig_g1", "Cgig_g9"]


def test_roundtrip_readable_by_pipeline_parser(tmp_path):
    # Write into a Results-style dir and read back with the SAME parser
    # pipeline.run() uses (steps.orthofinder.parse_orthogroups) to prove
    # drop-in compatibility.
    results_dir = tmp_path / "Results_Test"
    out_tsv = results_dir / "Orthogroups" / "Orthogroups.tsv"
    write_orthogroups_tsv(GROUPS, SPECIES_ORDER, out_tsv,
                          tmp_path / "Orthogroups_UnassignedGenes.tsv",
                          ALL_GENES)
    parsed = parse_orthogroups(results_dir)
    assert {og: sorted(genes) for og, genes in parsed.items()} == {
        og: sorted(genes) for og, genes in GROUPS.items()
    }


def test_orthogroups_tsv_layout(tmp_path):
    out_tsv = tmp_path / "Orthogroups.tsv"
    write_orthogroups_tsv(GROUPS, SPECIES_ORDER, out_tsv,
                          tmp_path / "unassigned.tsv", ALL_GENES)
    lines = out_tsv.read_text().splitlines()
    assert lines[0] == "Orthogroup\tCgig\tMcry\tOcoc"
    assert lines[1] == "1\tCgig_g1\tMcry_g1, Mcry_g2\tOcoc_g1"
    assert lines[2] == "2\t\t\tOcoc_g2, Ococ_g3"


def test_unassigned_is_all_genes_minus_grouped(tmp_path):
    unassigned_out = tmp_path / "unassigned.tsv"
    write_orthogroups_tsv(GROUPS, SPECIES_ORDER, tmp_path / "Orthogroups.tsv",
                          unassigned_out, ALL_GENES)
    lines = unassigned_out.read_text().splitlines()
    assert lines[0] == "Orthogroup\tCgig\tMcry\tOcoc"
    rows = [line.split("\t") for line in lines[1:]]
    genes = {cell for row in rows for cell in row[1:] if cell}
    assert genes == {"Mcry_g3", "Cgig_g9"}
    # One singleton per row, gene in its own species column.
    cgig_row = next(row for row in rows if row[1])
    assert cgig_row[1] == "Cgig_g9" and cgig_row[2] == "" and cgig_row[3] == ""


def test_unknown_species_prefix_warns_but_still_grouped(tmp_path, caplog):
    groups = {"1": ["Mcry_g1", "Zzzz_g1"]}
    out_tsv = tmp_path / "Orthogroups.tsv"
    unassigned_out = tmp_path / "unassigned.tsv"
    with caplog.at_level("WARNING", logger="family_finder"):
        write_orthogroups_tsv(groups, ["Mcry"], out_tsv, unassigned_out,
                              ["Mcry_g1", "Zzzz_g1"])
    assert "Zzzz_g1" in caplog.text
    # Grouped (even if unwritable into a species column) => NOT unassigned.
    assert "Zzzz_g1" not in unassigned_out.read_text()


# ---------------------------------------------------------------------------
# Template runner + groups-file discovery
# ---------------------------------------------------------------------------

def test_run_sonicparanoid_formats_template(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(
        sonicparanoid, "_run_template",
        lambda template, mapping, step: calls.append((template, mapping, step)),
    )
    config = Config(sonicparanoid_cmd="sonicparanoid -i {pep_dir} -o {outdir} -t 16")
    out = run_sonicparanoid(tmp_path / "pep", tmp_path / "sp_out", config)
    assert out == tmp_path / "sp_out" and out.is_dir()
    template, mapping, step = calls[0]
    assert mapping == {"pep_dir": str(tmp_path / "pep"),
                       "outdir": str(tmp_path / "sp_out")}
    assert step == "SonicParanoid2"


def test_run_sonicparanoid_unconfigured_template_raises(tmp_path):
    with pytest.raises(ValueError, match="not configured"):
        run_sonicparanoid(tmp_path / "pep", tmp_path / "sp_out", Config())


def test_find_groups_file_newest_wins(tmp_path):
    old = _write(tmp_path / "runs" / "run1" / "ortholog_groups" /
                 "ortholog_groups.tsv", "group_id\tMcry\n")
    new = _write(tmp_path / "runs" / "run2" / "ortholog_groups" /
                 "ortholog_groups.tsv", "group_id\tMcry\n")
    past = time.time() - 100
    os.utime(old, (past, past))
    assert find_groups_file(tmp_path) == new


def test_find_groups_file_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        find_groups_file(tmp_path)


# ---------------------------------------------------------------------------
# Config defaults
# ---------------------------------------------------------------------------

def test_clustering_config_defaults():
    config = Config()
    assert config.clustering_method == "orthofinder"
    assert config.sonicparanoid_cmd == ""
