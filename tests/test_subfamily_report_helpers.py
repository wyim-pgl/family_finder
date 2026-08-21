"""subfamily_report.py evidence-collection helpers.

These were extracted from a 143-line main() so they could be tested at all
(the refactor was verified byte-identical on the real PEPC report first).
"""
import json
import sys
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from subfamily_report import (
    _region_bounds,
    collect_selection_evidence,
    read_branch_name_map,
    write_subfunctionalization,
)


def _args(**kw):
    base = dict(expression_share=None, signal_partition=None,
                retargeting_events=0, relax_json=None, meme_json=None,
                meme_region=None, absrel_json=None, branch_name_map=None,
                branchsite_mlc=None, branchsite_lnl=None, disorder_json=None,
                family_name="fam", focal_subfamily="SF")
    base.update(kw)
    return types.SimpleNamespace(**base)


def test_region_bounds():
    assert _region_bounds("165-209") == (165, 209)


def test_read_branch_name_map_skips_blank_lines(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("g101\tTfru_real\n\ng002\tOcoc_real\n")
    assert read_branch_name_map(p) == {"g101": "Tfru_real", "g002": "Ococ_real"}


def test_collect_evidence_with_no_axes_is_empty_but_valid():
    ev = collect_selection_evidence(_args())
    assert ev["signal_partition"] == ""
    assert "relax" not in ev and "branchsite" not in ev


def test_collect_evidence_counts_meme_sites_in_region(tmp_path):
    j = tmp_path / "meme.json"
    j.write_text(json.dumps({"MLE": {
        "headers": [["p-value", ""]],
        "content": {"0": [[0.01], [1.0], [0.02]]},   # sites 1 and 3 significant
    }}))
    ev = collect_selection_evidence(_args(meme_json=str(j), meme_region="3-5"))
    assert ev["meme_sites_total"] == 2
    assert ev["meme_sites_in_region"] == 1          # only site 3


def test_collect_evidence_splits_absrel_stem_from_terminal(tmp_path):
    j = tmp_path / "absrel.json"
    j.write_text(json.dumps({"branch attributes": {"0": {
        "Node7": {"Corrected P-value": 1e-4},
        "g101": {"Corrected P-value": 1e-9},
        "Node9": {"Corrected P-value": 0.9},
    }}}))
    ev = collect_selection_evidence(_args(absrel_json=str(j)))
    assert [b for b, _ in ev["absrel_stem"]] == ["Node7"]
    assert [b for b, _ in ev["absrel_terminal"]] == ["g101"]


def test_collect_evidence_applies_branch_name_map(tmp_path):
    j = tmp_path / "absrel.json"
    j.write_text(json.dumps({"branch attributes": {"0": {
        "g101": {"Corrected P-value": 1e-9}}}}))
    m = tmp_path / "m.tsv"
    m.write_text("g101\tTfru_contig_062_000131\n")
    ev = collect_selection_evidence(
        _args(absrel_json=str(j), branch_name_map=str(m)))
    assert ev["absrel_terminal"][0][0] == "Tfru_contig_062_000131"


def test_collect_evidence_reads_disorder_json(tmp_path):
    j = tmp_path / "d.json"
    j.write_text(json.dumps({"delta": -0.3, "n_focal": 9, "n_other": 76,
                             "p": 1e-5, "all_below": True, "per_gene": []}))
    ev = collect_selection_evidence(_args(disorder_json=str(j)))
    assert ev["disorder"]["n_other"] == 76
    assert "per_gene" not in ev["disorder"]      # bulk payload not carried


def test_write_subfunctionalization_emits_both_files(tmp_path):
    verdict = {"verdict": "subfunctionalization",
               "evidence_for": ["a", "b"], "evidence_against": ["c"]}
    write_subfunctionalization(_args(), verdict, {}, tmp_path)
    md = (tmp_path / "subfunctionalization.md").read_text()
    tsv = (tmp_path / "subfunctionalization.tsv").read_text().splitlines()
    assert "subfunctionalization" in md
    assert tsv[0].startswith("subfamily\tverdict")
    assert tsv[1].split("\t")[2] == "2"          # n_evidence_for
