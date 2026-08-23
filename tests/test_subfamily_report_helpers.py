"""subfamily_report.py evidence-collection helpers.

These were extracted from a 143-line main() so they could be tested at all
(the refactor was verified byte-identical on the real PEPC report first).
"""
import json
import sys
import types
from pathlib import Path

import pytest

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
                sites_alignment=None, region_alignment=None, coord_bridge=None,
                alignment=None, family_name="fam", focal_subfamily="SF")
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
    header = tsv[0].split("\t")
    row = tsv[1].split("\t")
    assert header[:2] == ["subfamily", "verdict"]
    assert row[header.index("n_evidence_for")] == "2"
    # a verdict that never showed its coordinate system says so in the table
    assert row[header.index("coordinates_verified")] == "False"


# --- coordinate verification (issue #42) ---------------------------------

_MEME = {"MLE": {"headers": [["p-value", ""]],
                 "content": {"0": [[0.01], [1.0], [0.02], [0.01], [0.01]]}}}


def _meme_json(tmp_path):
    j = tmp_path / "meme.json"
    j.write_text(json.dumps(_MEME))
    return str(j)


def _fasta(tmp_path, name, seqs):
    p = tmp_path / name
    p.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    return str(p)


def test_region_counts_are_unverified_unless_both_alignments_are_given(tmp_path):
    ev = collect_selection_evidence(
        _args(meme_json=_meme_json(tmp_path), meme_region="3-5"))

    assert ev["coordinates_verified"] is False


def test_the_same_alignment_on_both_sides_verifies_the_count(tmp_path):
    aln = _fasta(tmp_path, "a.fa", {"g1": "MK-WQ", "g2": "MKYWQ"})

    ev = collect_selection_evidence(
        _args(meme_json=_meme_json(tmp_path), meme_region="3-5",
              sites_alignment=aln, region_alignment=aln))

    assert ev["coordinates_verified"] is True
    assert ev["meme_sites_in_region"] == 3      # sites 3, 4 and 5


def test_different_alignments_are_translated_through_the_bridge(tmp_path):
    sites_aln = _fasta(tmp_path, "s.fa", {"g1": "MK-WQ", "g2": "MKYWQ"})
    region_aln = _fasta(tmp_path, "r.fa", {"g1": "--MKWQ", "g2": "XXMKWQ"})

    ev = collect_selection_evidence(
        _args(meme_json=_meme_json(tmp_path), meme_region="5-6",
              sites_alignment=sites_aln, region_alignment=region_aln,
              coord_bridge="g1"))

    # sites 4 and 5 (residues 3 and 4 of g1) become columns 5 and 6; site 3 is
    # a gap in the bridge and site 1 lands on column 3, outside the region
    assert ev["coordinates_verified"] is True
    assert ev["meme_sites_in_region"] == 2
    assert ev["meme_sites_untranslatable"] == 1


def test_the_region_alignment_defaults_to_the_report_alignment(tmp_path):
    aln = _fasta(tmp_path, "a.fa", {"g1": "MK-WQ", "g2": "MKYWQ"})

    ev = collect_selection_evidence(
        _args(meme_json=_meme_json(tmp_path), meme_region="3-5",
              sites_alignment=aln, alignment=aln))

    assert ev["coordinates_verified"] is True


# --- structural coherence, sequence-identity control (issue #39 item 3) ------

def _pairs_args(**kw):
    base = dict(pairs=None, pair_columns=None, pair_metric="bits",
                min_group=5, delimiter="_")
    base.update(kw)
    return types.SimpleNamespace(**base)


_WIDE_COLUMNS = "query,target,fident,alnlen,evalue,bits,qtmscore,alntmscore"

_GROUPS = {"SF_A": ["g1", "g2", "g3"], "SF_B": ["g4", "g5"]}


def _wide_pairs(tmp_path):
    """Six pairs whose score is exactly 1000 * fident — no signal beyond identity."""
    p = tmp_path / "pairs.tsv"
    rows = [("g1", "g2", 0.90), ("g1", "g3", 0.88), ("g2", "g3", 0.91),
            ("g1", "g4", 0.40), ("g2", "g5", 0.42), ("g4", "g5", 0.95)]
    p.write_text("".join(
        f"{a}\t{b}\t{f:.3f}\t900\t0.0\t{1000 * f:.3f}\t0.9\t0.9\n"
        for a, b, f in rows))
    return p


def test_the_confounding_is_measured_when_fident_is_present(tmp_path):
    from subfamily_report import build_structure_coherence

    rows, confounding = build_structure_coherence(
        _pairs_args(pairs=str(_wide_pairs(tmp_path)),
                    pair_columns=_WIDE_COLUMNS), _GROUPS)

    assert confounding["r"] == pytest.approx(1.0)
    assert rows


def test_fident_being_present_is_not_the_same_as_the_control_holding(tmp_path):
    # these six pairs put every within-pair at high identity and every
    # between-pair low, so there is no shared identity range to compare over —
    # the real PEPC layout in miniature (#39)
    from subfamily_report import build_structure_coherence

    rows, _ = build_structure_coherence(
        _pairs_args(pairs=str(_wide_pairs(tmp_path)),
                    pair_columns=_WIDE_COLUMNS), _GROUPS)

    for row in rows:
        assert row["verdict"] == "no_interpretation_available"
        assert row["sequence_controlled"] is False
        assert row["coherent_controlled"] is None


def test_coherence_rows_carry_a_warning_when_fident_is_absent(tmp_path):
    from subfamily_report import build_structure_coherence

    p = tmp_path / "pairs.tsv"
    p.write_text("g1\tg2\t0.0\t900\t0.9\n"
                 "g1\tg4\t0.0\t400\t0.4\n")

    rows, confounding = build_structure_coherence(
        _pairs_args(pairs=str(p)), _GROUPS)

    assert confounding is None
    assert all(r["sequence_controlled"] is False for r in rows)
    assert all("UNCONTROLLED" in r["warning"] for r in rows)


def test_the_written_table_keeps_the_warning_column(tmp_path):
    from subfamily_report import build_structure_coherence, write_tsv

    p = tmp_path / "pairs.tsv"
    p.write_text("g1\tg2\t0.0\t900\t0.9\ng1\tg4\t0.0\t400\t0.4\n")
    rows, _ = build_structure_coherence(_pairs_args(pairs=str(p)), _GROUPS)
    out = tmp_path / "structure_coherence.tsv"

    write_tsv(rows, out)

    header = out.read_text().splitlines()[0].split("\t")
    for column in ("metric", "sequence_controlled", "coherent_controlled",
                   "warning"):
        assert column in header


# --- anchor transferability: the grade has to reach the TSV ----------------

_TREE = "((((P1,q1)100/100,q2)100/76,((P2,r1)100/100,r2)100/99)90/90,out);"


def _anchor_args(tmp_path, tree=_TREE, **kw):
    (tmp_path / "tree.nwk").write_text(tree)
    (tmp_path / "anchors.tsv").write_text("P1\tppc-1E1\nP2\tppc-1E2\n")
    base = dict(family_tree=str(tmp_path / "tree.nwk"),
                anchors=str(tmp_path / "anchors.tsv"),
                anchor_min_support=None, query_prefix=None,
                cross_tree=None, tree_name="this tree")
    base.update(kw)
    return types.SimpleNamespace(**base)


def test_no_tree_means_no_transferability_table(tmp_path):
    from subfamily_report import build_transferability

    args = _anchor_args(tmp_path, family_tree=None)
    assert build_transferability(args) == []


def test_transferability_rows_carry_the_grade(tmp_path):
    from subfamily_report import build_transferability

    by = {r["label"]: r for r in build_transferability(_anchor_args(tmp_path))}

    assert by["ppc-1E1"]["grade"] == "PROVISIONAL"
    assert by["ppc-1E1"]["support"] == "100/76"
    assert "76" in by["ppc-1E1"]["grade_reason"]
    assert by["ppc-1E1"]["cross_tree_evaluated"] == "False"


def test_cross_tree_datasets_are_read_from_their_own_trees(tmp_path):
    """--cross-tree NAME=tree runs the same clade rule on an independent
    dataset; agreement is the membership coming back identical."""
    from subfamily_report import build_transferability

    other = tmp_path / "aa.nwk"
    other.write_text("((((P1,q1)100/100,q2)100/100,"
                     "((P2,r1)100/100,r2)100/100)90/90,out);")
    args = _anchor_args(tmp_path, tree_name="codon123",
                        cross_tree=[f"aa102={other}"])

    by = {r["label"]: r for r in build_transferability(args)}

    assert by["ppc-1E2"]["grade"] == "HIGH"
    assert by["ppc-1E2"]["consistent_datasets"] == "codon123;aa102"
    assert by["ppc-1E1"]["grade"] == "PROVISIONAL"    # UFboot 76 in this tree


def test_a_disagreeing_cross_tree_is_named_not_hidden(tmp_path):
    from subfamily_report import build_transferability

    other = tmp_path / "aa.nwk"
    other.write_text("(((P1,q1)100/100,(q2,(P2,r1)100/100)100/100)90/90,out);")
    args = _anchor_args(tmp_path, cross_tree=[f"aa102={other}"])

    by = {r["label"]: r for r in build_transferability(args)}

    assert by["ppc-1E2"]["grade"] == "UNRESOLVED"
    assert by["ppc-1E2"]["inconsistent_datasets"] == "aa102"


def test_the_written_table_keeps_the_grade_columns(tmp_path):
    from subfamily_report import build_transferability, write_tsv

    out = tmp_path / "anchor_transferability.tsv"
    write_tsv(build_transferability(_anchor_args(tmp_path)), out)

    header = out.read_text().splitlines()[0].split("\t")
    for column in ("grade", "grade_reason", "sh_alrt", "ufboot",
                   "n_consistent_datasets", "consistent_datasets",
                   "inconsistent_datasets", "transferable"):
        assert column in header
