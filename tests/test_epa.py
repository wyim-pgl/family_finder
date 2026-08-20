"""Tests for steps/epa.py (issue #23): jplace parsing (fields-order
independence), LWR acceptance rules, edge->family mapping on a synthetic
merged clan, and the thin adjudicate orchestration.

No epa-ng/raxml-ng/hmmer required: subprocess wrappers are monkeypatched
and the parser is fed synthetic jplace JSON fixtures (house convention).
"""

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from config import Config  # noqa: E402
from steps import epa  # noqa: E402
from steps.epa import (  # noqa: E402
    JointReference,
    Placement,
    accept_placement,
    adjudicate,
    align_query,
    edge_family_map,
    evaluate_model,
    parse_jplace,
    run_epa,
    split_alignment,
    write_adjudications_tsv,
)


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

DEFAULT_FIELDS = [
    "edge_num", "likelihood", "like_weight_ratio",
    "distal_length", "pendant_length",
]

# Two-family clan: famA = {Sp1_a1, Sp2_a2}, famB = {Sp1_b1, Sp2_b2}.
CLAN_TREE = (
    "((Sp1_a1:0.1{0},Sp2_a2:0.1{1}):0.05{2},"
    "(Sp1_b1:0.1{3},Sp2_b2:0.1{4}):0.05{5}){6};"
)
CLAN_FAMILIES = {
    "Sp1_a1": "famA", "Sp2_a2": "famA",
    "Sp1_b1": "famB", "Sp2_b2": "famB",
}


def _row(fields, *, edge, lwr, lh=-100.0, distal=0.01, pendant=0.02):
    values = {
        "edge_num": edge, "likelihood": lh, "like_weight_ratio": lwr,
        "distal_length": distal, "pendant_length": pendant,
    }
    return [values[name] for name in fields]


def _write_jplace(path, tree, placements, fields=None):
    fields = list(DEFAULT_FIELDS if fields is None else fields)
    data = {
        "version": 3,
        "tree": tree,
        "placements": placements,
        "fields": fields,
        "metadata": {"invocation": "test"},
    }
    path = Path(path)
    path.write_text(json.dumps(data))
    return path


def _placement(edge, lwr, pendant=0.02):
    return Placement(
        edge_num=edge, likelihood=-100.0, like_weight_ratio=lwr,
        distal_length=0.01, pendant_length=pendant,
    )


# ---------------------------------------------------------------------------
# parse_jplace
# ---------------------------------------------------------------------------

def test_parse_jplace_reads_fields_ordering(tmp_path):
    entry = {"p": [_row(DEFAULT_FIELDS, edge=3, lwr=0.9)], "n": ["Sp1_q1"]}
    path = _write_jplace(tmp_path / "a.jplace", CLAN_TREE, [entry])
    parsed = parse_jplace(path)
    assert parsed["Sp1_q1"] == [_placement(3, 0.9)]


def test_parse_jplace_fields_order_independence(tmp_path):
    # Same placement serialized under a SHUFFLED fields array must parse
    # to the identical Placement — column positions are never hardcoded.
    shuffled = [
        "like_weight_ratio", "pendant_length", "edge_num",
        "distal_length", "likelihood",
    ]
    entry_default = {"p": [_row(DEFAULT_FIELDS, edge=3, lwr=0.9)], "n": ["q"]}
    entry_shuffled = {"p": [_row(shuffled, edge=3, lwr=0.9)], "n": ["q"]}
    p1 = _write_jplace(tmp_path / "a.jplace", CLAN_TREE, [entry_default])
    p2 = _write_jplace(tmp_path / "b.jplace", CLAN_TREE, [entry_shuffled],
                       fields=shuffled)
    assert parse_jplace(p1) == parse_jplace(p2)


def test_parse_jplace_multi_placement_lwr_ordering(tmp_path):
    rows = [
        _row(DEFAULT_FIELDS, edge=0, lwr=0.2),
        _row(DEFAULT_FIELDS, edge=5, lwr=0.7),
        _row(DEFAULT_FIELDS, edge=2, lwr=0.1),
    ]
    path = _write_jplace(tmp_path / "a.jplace", CLAN_TREE,
                         [{"p": rows, "n": ["q"]}])
    lwrs = [p.like_weight_ratio for p in parse_jplace(path)["q"]]
    assert lwrs == [0.7, 0.2, 0.1]


def test_parse_jplace_nm_name_convention(tmp_path):
    entry = {"p": [_row(DEFAULT_FIELDS, edge=1, lwr=0.8)],
             "nm": [["q1", 1], ["q2", 2]]}
    path = _write_jplace(tmp_path / "a.jplace", CLAN_TREE, [entry])
    parsed = parse_jplace(path)
    assert set(parsed) == {"q1", "q2"}
    assert parsed["q1"] == parsed["q2"]


def test_parse_jplace_missing_required_field_raises(tmp_path):
    fields = ["edge_num", "likelihood", "like_weight_ratio"]  # no lengths
    entry = {"p": [[3, -100.0, 0.9]], "n": ["q"]}
    path = _write_jplace(tmp_path / "a.jplace", CLAN_TREE, [entry],
                         fields=fields)
    with pytest.raises(ValueError, match="missing required fields"):
        parse_jplace(path)


# ---------------------------------------------------------------------------
# accept_placement
# ---------------------------------------------------------------------------

def test_accept_single_high_lwr_placement():
    config = Config()
    assert accept_placement([_placement(0, 0.95)], config) == _placement(0, 0.95)


def test_accept_rejects_low_lwr():
    config = Config()  # epa_min_lwr = 0.8
    assert accept_placement([_placement(0, 0.79)], config) is None


def test_accept_rejects_failing_margin():
    config = Config()  # epa_lwr_margin = 0.3
    placements = [_placement(0, 0.85), _placement(5, 0.7)]  # margin 0.15
    assert accept_placement(placements, config) is None


def test_accept_passes_with_margin():
    config = Config()
    placements = [_placement(5, 0.9), _placement(0, 0.05)]
    accepted = accept_placement(placements, config)
    assert accepted is not None and accepted.edge_num == 5


def test_accept_exact_tie_is_never_broken():
    # Even with a ZERO margin threshold an exact LWR tie stays ambiguous.
    config = Config(epa_lwr_margin=0.0)
    placements = [_placement(0, 0.9), _placement(5, 0.9)]
    assert accept_placement(placements, config) is None


def test_accept_empty_placements():
    assert accept_placement([], Config()) is None


def test_accept_pendant_gate_rejects_long_pendant_branch():
    # High LWR only says "least-bad edge" — a huge pendant branch means
    # the query belongs to no reference family (the abstention rule).
    config = Config(epa_max_pendant=1.0)
    assert accept_placement([_placement(0, 0.95, pendant=2.5)], config) is None


def test_accept_pendant_gate_disabled_by_default():
    config = Config()  # epa_max_pendant = 0.0 -> disabled
    accepted = accept_placement([_placement(0, 0.95, pendant=2.5)], config)
    assert accepted is not None and accepted.edge_num == 0


# ---------------------------------------------------------------------------
# edge_family_map
# ---------------------------------------------------------------------------

def test_edge_family_map_two_family_clan():
    mapping = edge_family_map(CLAN_TREE, CLAN_FAMILIES)
    assert mapping[0] == "famA"
    assert mapping[1] == "famA"
    assert mapping[2] == "famA"   # famA-side internal edge
    assert mapping[3] == "famB"
    assert mapping[4] == "famB"
    assert mapping[5] == "famB"   # famB-side internal edge
    # Root edge sees 2 famA vs 2 famB leaves: tied -> omitted (ambiguous).
    assert 6 not in mapping


def test_edge_family_map_majority_wins_in_mixed_clade():
    tree = "((Sp1_a1:0.1{0},Sp1_a2:0.1{1},Sp1_b1:0.1{2}):0.05{3},Sp1_c1:0.2{4});"
    families = {"Sp1_a1": "famA", "Sp1_a2": "famA",
                "Sp1_b1": "famB", "Sp1_c1": "famC"}
    mapping = edge_family_map(tree, families)
    assert mapping[3] == "famA"   # 2 famA vs 1 famB below edge 3
    assert mapping[2] == "famB"
    assert mapping[4] == "famC"


def test_edge_family_map_uninformative_leaves_ignored():
    mapping = edge_family_map(CLAN_TREE, {"Sp1_a1": "famA"})
    assert mapping[0] == "famA"
    assert mapping[2] == "famA"
    assert 3 not in mapping       # no informative leaf below edge 3
    assert mapping[6] == "famA"   # sole informative family -> no tie


# ---------------------------------------------------------------------------
# Subprocess wrappers (mocked _run_cmd; commands captured)
# ---------------------------------------------------------------------------

def _capture_run_cmd(calls):
    def fake(cmd, step, stdout_path=None):
        calls.append((list(cmd), step, stdout_path))
    return fake


def test_evaluate_model_returns_bestmodel_and_besttree(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(epa, "_run_cmd", _capture_run_cmd(calls))
    model_file, eval_tree = evaluate_model(tmp_path / "ref.fa",
                                           tmp_path / "ref.nwk",
                                           tmp_path / "out", Config())
    cmd = calls[0][0]
    assert cmd[0] == "raxml-ng" and "--evaluate" in cmd
    assert str(model_file).endswith(".raxml.bestModel")
    # Placement must run on the re-optimized tree, so it is returned too.
    assert str(eval_tree).endswith(".raxml.bestTree")


def test_align_query_hmmalign_backend(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(epa, "_run_cmd", _capture_run_cmd(calls))
    combined = align_query(tmp_path / "q.fa", tmp_path / "ref.fa",
                           tmp_path / "out", Config())
    assert [step for _, step, _ in calls] == ["hmmbuild (EPA reference)", "hmmalign"]
    hmmalign_cmd = calls[1][0]
    assert hmmalign_cmd[0] == "hmmalign" and "--mapali" in hmmalign_cmd
    assert combined.name == "combined.fa"


def test_align_query_mafft_add_backend(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(epa, "_run_cmd", _capture_run_cmd(calls))
    config = Config(epa_query_align="mafft-add")
    combined = align_query(tmp_path / "q.fa", tmp_path / "ref.fa",
                           tmp_path / "out", config)
    cmd, _, stdout_path = calls[0]
    assert cmd[0] == "mafft" and "--add" in cmd and "--keeplength" in cmd
    assert stdout_path == combined  # mafft writes the alignment to stdout


def test_align_query_unknown_backend_raises(tmp_path):
    config = Config(epa_query_align="muscle")
    with pytest.raises(ValueError, match="epa_query_align"):
        align_query(tmp_path / "q.fa", tmp_path / "ref.fa",
                    tmp_path / "out", config)


def test_split_and_run_epa_paths(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(epa, "_run_cmd", _capture_run_cmd(calls))
    config = Config()
    ref_msa, query_msa = split_alignment(tmp_path / "combined.fa",
                                         tmp_path / "ref.fa",
                                         tmp_path / "out", config)
    assert ref_msa == tmp_path / "ref.fa"       # epa-ng writes only query.fasta
    assert query_msa == tmp_path / "out" / "query.fasta"
    jplace = run_epa(tmp_path / "ref.nwk", ref_msa, query_msa,
                     tmp_path / "eval.raxml.bestModel", tmp_path / "out", config)
    assert jplace == tmp_path / "out" / "epa_result.jplace"
    epa_cmd = calls[1][0]
    assert "--model" in epa_cmd and "--redo" in epa_cmd


# ---------------------------------------------------------------------------
# adjudicate (thin orchestration on synthetic jplace + small newick)
# ---------------------------------------------------------------------------

def _patch_adjudicate_stack(monkeypatch, tmp_path, placement_entries):
    """Mock the four wrappers; run_epa writes a synthetic clan jplace."""
    monkeypatch.setattr(epa, "evaluate_model",
                        lambda *a, **k: (tmp_path / "eval.raxml.bestModel",
                                         tmp_path / "eval.raxml.bestTree"))
    monkeypatch.setattr(epa, "align_query",
                        lambda *a, **k: tmp_path / "combined.fa")
    monkeypatch.setattr(epa, "split_alignment",
                        lambda *a, **k: (tmp_path / "ref.fa",
                                         tmp_path / "query.fasta"))

    def fake_run_epa(ref_tree, ref_msa, query_msa, model_file, outdir, config):
        return _write_jplace(Path(outdir) / "epa_result.jplace",
                             CLAN_TREE, placement_entries)

    monkeypatch.setattr(epa, "run_epa", fake_run_epa)


def _clan_ref(tmp_path):
    return JointReference(
        alignment=tmp_path / "joint.fa",
        tree=tmp_path / "joint.nwk",
        family_of_leaf=dict(CLAN_FAMILIES),
    )


def test_adjudicate_places_query_on_famA_side(tmp_path, monkeypatch):
    entry = {"p": [_row(DEFAULT_FIELDS, edge=2, lwr=0.95)], "n": ["Sp1_q"]}
    _patch_adjudicate_stack(monkeypatch, tmp_path, [entry])
    family = adjudicate("Sp1_q", ["famA", "famB"], _clan_ref(tmp_path),
                        tmp_path / "q.fa", tmp_path / "adj", Config())
    assert family == "famA"


def test_adjudicate_places_query_on_famB_side(tmp_path, monkeypatch):
    entry = {"p": [_row(DEFAULT_FIELDS, edge=4, lwr=0.95)], "n": ["Sp1_q"]}
    _patch_adjudicate_stack(monkeypatch, tmp_path, [entry])
    family = adjudicate("Sp1_q", ["famA", "famB"], _clan_ref(tmp_path),
                        tmp_path / "q.fa", tmp_path / "adj", Config())
    assert family == "famB"


def test_adjudicate_ambiguous_on_low_lwr_or_margin(tmp_path, monkeypatch):
    entry = {"p": [_row(DEFAULT_FIELDS, edge=2, lwr=0.55),
                   _row(DEFAULT_FIELDS, edge=4, lwr=0.45)], "n": ["Sp1_q"]}
    _patch_adjudicate_stack(monkeypatch, tmp_path, [entry])
    assert adjudicate("Sp1_q", ["famA", "famB"], _clan_ref(tmp_path),
                      tmp_path / "q.fa", tmp_path / "adj", Config()) is None


def test_adjudicate_tied_root_edge_is_ambiguous(tmp_path, monkeypatch):
    # High-confidence placement on the ROOT edge of the balanced clan:
    # the edge has no family majority, so the gene stays ambiguous.
    entry = {"p": [_row(DEFAULT_FIELDS, edge=6, lwr=0.95)], "n": ["Sp1_q"]}
    _patch_adjudicate_stack(monkeypatch, tmp_path, [entry])
    assert adjudicate("Sp1_q", ["famA", "famB"], _clan_ref(tmp_path),
                      tmp_path / "q.fa", tmp_path / "adj", Config()) is None


def test_adjudicate_outside_candidates_is_rejected(tmp_path, monkeypatch):
    entry = {"p": [_row(DEFAULT_FIELDS, edge=2, lwr=0.95)], "n": ["Sp1_q"]}
    _patch_adjudicate_stack(monkeypatch, tmp_path, [entry])
    # famA wins the edge but was not among the nominated candidates.
    assert adjudicate("Sp1_q", ["famB", "famC"], _clan_ref(tmp_path),
                      tmp_path / "q.fa", tmp_path / "adj", Config()) is None


def test_adjudicate_places_on_the_evaluated_tree(tmp_path, monkeypatch):
    # jplace edge numbers refer to the tree given to epa-ng: adjudicate
    # must pass the raxml-ng bestTree, not the raw input reference tree.
    entry = {"p": [_row(DEFAULT_FIELDS, edge=2, lwr=0.95)], "n": ["Sp1_q"]}
    _patch_adjudicate_stack(monkeypatch, tmp_path, [entry])
    seen = {}
    orig_run_epa = epa.run_epa

    def spying_run_epa(ref_tree, *args, **kwargs):
        seen["tree"] = ref_tree
        return orig_run_epa(ref_tree, *args, **kwargs)

    monkeypatch.setattr(epa, "run_epa", spying_run_epa)
    adjudicate("Sp1_q", ["famA", "famB"], _clan_ref(tmp_path),
               tmp_path / "q.fa", tmp_path / "adj", Config())
    assert str(seen["tree"]).endswith(".raxml.bestTree")


def test_adjudicate_query_absent_from_jplace(tmp_path, monkeypatch):
    entry = {"p": [_row(DEFAULT_FIELDS, edge=2, lwr=0.95)], "n": ["Sp1_other"]}
    _patch_adjudicate_stack(monkeypatch, tmp_path, [entry])
    assert adjudicate("Sp1_q", ["famA", "famB"], _clan_ref(tmp_path),
                      tmp_path / "q.fa", tmp_path / "adj", Config()) is None


# ---------------------------------------------------------------------------
# TSV output + config defaults
# ---------------------------------------------------------------------------

def test_write_adjudications_tsv(tmp_path):
    path = tmp_path / "adjudications.tsv"
    write_adjudications_tsv({"Sp1_q2": None, "Sp1_q1": "famA"}, path)
    lines = path.read_text().splitlines()
    assert lines[0] == "gene_id\tfamily_id\tstatus"
    assert lines[1] == "Sp1_q1\tfamA\tplaced"
    assert lines[2] == "Sp1_q2\t\tambiguous"


def test_epa_config_defaults():
    config = Config()
    assert config.epa_ng_bin == "epa-ng"
    assert config.raxml_ng_bin == "raxml-ng"
    assert config.hmmalign_bin == "hmmalign"
    assert config.epa_query_align == "hmmalign"
    assert config.epa_min_lwr == 0.8
    assert config.epa_lwr_margin == 0.3
    assert config.epa_max_pendant == 0.0  # abstention gate off by default
