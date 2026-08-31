"""End-to-end tests of the adjudicate_labels.py CLI on tiny fixtures."""

import json
from pathlib import Path

import pytest

import adjudicate_labels as cli

OVERRIDES = Path(__file__).resolve().parents[1] / "data" / "symbol_stem_overrides_v1.tsv"


def write(path, lines):
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture()
def tiny_run(tmp_path):
    fams = write(tmp_path / "families.tsv", [
        "family_id\tgene_list",
        "F1\tMcry_a,Mcry_b",
        "F2\tMcry_c,Ococ_d",
        "F3\tMcry_e",
    ])
    curated = write(tmp_path / "curated.tsv", [
        "family_id\tgene\tsymbol",
        "F1\tMcry_a\tBam1",
        "F1\tMcry_b\tBam3",
        "F2\tMcry_c\tFbp1",
    ])
    ath = write(tmp_path / "ath.tsv", [
        "family_id\tgene\tsymbol\tgate_passed\tgate_reason",
        "F2\tOcoc_d\tSbp1\tTrue\t",
        "F3\tMcry_e\tTpt1\tFalse\tlow margin",
    ])
    out = tmp_path / "out"
    argv = ["--families", str(fams),
            "--calls", "mcry_curated=%s" % curated,
            "--calls", "ath_diamond=%s" % ath,
            "--overrides", str(OVERRIDES),
            "--min-label-coverage", "0.0",
            "--out", str(out)]
    return argv, out


def read_labels(out):
    lines = (out / "labels.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def test_cli_end_to_end(tiny_run, capsys):
    argv, out = tiny_run
    assert cli.main(argv) == 0
    rows = read_labels(out)
    by_family = {}
    for r in rows:
        by_family.setdefault(r["target_id"], []).append(r)
    f1 = [r for r in by_family["F1"] if r["transfer_tier"] == "stem_auto"][0]
    assert (f1["stem"], f1["verdict"]) == ("Bam", "SINGLE_STEM_MULTI_SUFFIX")
    suffixes = sorted(r["suffix"] for r in by_family["F1"]
                      if r["transfer_tier"] == "suffix_needs_tree_gate")
    assert suffixes == ["1", "3"]  # suffix column is present in the artifact
    f2 = [r for r in by_family["F2"] if r["level"] == "family"][0]
    assert f2["verdict"] == "NEEDS_REVIEW_MULTI_STEM"
    assert "F3" not in by_family  # NO_EVIDENCE emits no label rows...
    ev = (out / "evidence.tsv").read_text()
    assert '"verdict": "NO_EVIDENCE"' in ev  # ...but leaves a positive record
    assert "abstention" in ev and "low margin" in ev
    vc = (out / "verdicts.tsv").read_text()
    assert "NO_EVIDENCE\t1" in vc


def test_cli_rejects_unknown_source(tiny_run, tmp_path):
    argv, out = tiny_run
    bad = write(tmp_path / "typo.tsv", [
        "family_id\tgene\tsymbol",
        "F1\tMcry_a\tBam1",
    ])
    argv2 = argv + ["--calls", "afdb_swissport=%s" % bad]  # transposed letter
    with pytest.raises(SystemExit, match="unknown call source"):
        cli.main(argv2)


def test_cli_rejects_garbled_boolean(tiny_run, tmp_path):
    argv, out = tiny_run
    bad = write(tmp_path / "garbled.tsv", [
        "family_id\tgene\tsymbol\tgate_passed",
        "F1\tMcry_a\tBam1\tmaybe",
    ])
    argv2 = list(argv)
    argv2[argv2.index("--calls") + 1] = "mcry_curated=%s" % bad
    with pytest.raises(SystemExit, match="boolean"):
        cli.main(argv2)


def test_cli_refuses_resume_after_input_change(tiny_run, tmp_path):
    argv, out = tiny_run
    assert cli.main(argv) == 0
    fams = Path(argv[1])
    fams.write_text(fams.read_text() + "F9\tMcry_z\n")
    with pytest.raises(SystemExit, match="changed"):
        cli.main(argv)


def test_cli_rejects_calls_for_unknown_family(tiny_run, tmp_path):
    argv, out = tiny_run
    bad = write(tmp_path / "bad.tsv", [
        "family_id\tgene\tsymbol",
        "F_UNKNOWN\tMcry_a\tBam1",
    ])
    argv2 = argv + ["--calls", "ath_diamond=%s" % bad]
    with pytest.raises(SystemExit, match="unknown families"):
        cli.main(argv2)


def test_manifest_records_inputs_and_policy(tiny_run):
    argv, out = tiny_run
    cli.main(argv)
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["parser_version"] == cli.PARSER_VERSION
    assert manifest["min_label_coverage"] == 0.0
    assert len(manifest["inputs"]) == 4  # families, overrides, 2 calls files
