"""Tests for steps/plant_glm.py (issue #18) and steps/ecforest.py (issue #20),
plus the annotate_families.py helpers.

No gLM/ECForest binaries required: subprocess wrappers are monkeypatched,
parsers eat synthetic text fixtures, and tree mapping reuses the NWK fixture
style of tests/test_retargeting.py.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from config import Config  # noqa: E402
from steps import ecforest, plant_glm  # noqa: E402
from steps.ecforest import (  # noqa: E402
    consensus_disagreements,
    ec_switch_events,
    family_ec_consensus,
    load_cache,
    parse_ecforest,
    run_ecforest,
    save_cache,
)
from steps.plant_glm import (  # noqa: E402
    parse_glm_scores,
    run_glm_scores,
    splice_plausibility_report,
    write_report_tsv,
)

# FastTree-style family tree (same fixture style as test_retargeting.py)
NWK = "((Mcry_g1:0.5,Cgig_g1:0.1)0.99:0.2,(Ococ_g1:0.1,Obas_g1:0.1)0.95:0.2);"


class _FakeCompleted:
    def __init__(self, returncode=0, stderr=""):
        self.returncode = returncode
        self.stderr = stderr
        self.stdout = ""


def _pred(ec, is_enzyme=True, conf=0.9):
    return {"is_enzyme": is_enzyme, "ec": ec, "confidence": conf}


# ---------------------------------------------------------------------------
# plant_glm: parser
# ---------------------------------------------------------------------------

def test_parse_glm_scores_with_header(tmp_path):
    tsv = tmp_path / "glm_scores.tsv"
    tsv.write_text(
        "seq_id\tscore\n"
        "Cgig_g1.t1\t-1.25\n"
        "Cgig_g1.t2\t-3.75\n"
        "broken row\n"
    )
    scores = parse_glm_scores(tsv)
    assert scores == {"Cgig_g1.t1": -1.25, "Cgig_g1.t2": -3.75}


def test_parse_glm_scores_tolerates_missing_header(tmp_path):
    tsv = tmp_path / "scores.tsv"
    tsv.write_text("iso1\t0.5\niso2\t0.25\n")
    assert parse_glm_scores(tsv) == {"iso1": 0.5, "iso2": 0.25}


# ---------------------------------------------------------------------------
# plant_glm: splice plausibility (pure function)
# ---------------------------------------------------------------------------

def test_splice_plausibility_ranks_and_flags():
    scores = {"g1.t1": -1.0, "g1.t2": -1.4, "g1.t3": -9.0, "g2.t1": -2.0}
    gene_map = {"g1": ["g1.t1", "g1.t2", "g1.t3"], "g2": ["g2.t1"]}
    rows = splice_plausibility_report(scores, gene_map, min_gap=2.0)

    g1 = [r for r in rows if r["gene"] == "g1"]
    assert [r["isoform"] for r in g1] == ["g1.t1", "g1.t2", "g1.t3"]  # score desc
    assert [r["rank"] for r in g1] == [1, 2, 3]
    # t2 is only 0.4 below best -> plausible; t3 is 8.0 below -> artifact candidate
    assert [r["low_plausibility"] for r in g1] == [False, False, True]
    assert g1[2]["gap"] == pytest.approx(8.0)

    # single-isoform gene: one unflagged row (no within-gene comparison)
    g2 = [r for r in rows if r["gene"] == "g2"]
    assert len(g2) == 1 and not g2[0]["low_plausibility"]


def test_splice_plausibility_skips_unscored_isoforms_and_genes():
    rows = splice_plausibility_report(
        {"g1.t1": 1.0}, {"g1": ["g1.t1", "g1.t2"], "g3": ["g3.t1"]}, min_gap=1.0
    )
    assert [(r["gene"], r["isoform"]) for r in rows] == [("g1", "g1.t1")]


def test_write_report_tsv(tmp_path):
    rows = splice_plausibility_report(
        {"g1.t1": 0.0, "g1.t2": -5.0}, {"g1": ["g1.t1", "g1.t2"]}, min_gap=2.0
    )
    path = tmp_path / "report.tsv"
    write_report_tsv(rows, path)
    lines = path.read_text().splitlines()
    assert lines[0].startswith("gene\tisoform\trank")
    assert lines[2].endswith("\t1")  # flagged isoform


# ---------------------------------------------------------------------------
# plant_glm: wrapper
# ---------------------------------------------------------------------------

def test_run_glm_scores_requires_configured_template(tmp_path):
    with pytest.raises(ValueError, match="not configured"):
        run_glm_scores(tmp_path / "iso.fa", tmp_path / "out", Config())


def test_run_glm_scores_returns_interchange_tsv(tmp_path, monkeypatch):
    outdir = tmp_path / "out"

    def fake_run(cmd, **kwargs):
        # the adapter's contract: write {outdir}/glm_scores.tsv
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / plant_glm.GLM_SCORES_TSV).write_text("seq_id\tscore\niso1\t0.1\n")
        return _FakeCompleted(0)

    monkeypatch.setattr(plant_glm.subprocess, "run", fake_run)
    config = Config(glm_score_cmd="python adapter.py --fasta {fasta} --outdir {outdir}")
    tsv = run_glm_scores(tmp_path / "iso.fa", outdir, config)
    assert tsv.name == "glm_scores.tsv"
    assert parse_glm_scores(tsv) == {"iso1": 0.1}


def test_run_glm_scores_raises_when_adapter_breaks_contract(tmp_path, monkeypatch):
    monkeypatch.setattr(
        plant_glm.subprocess, "run", lambda cmd, **kw: _FakeCompleted(0)
    )
    config = Config(glm_score_cmd="adapter {fasta} {outdir}")
    with pytest.raises(FileNotFoundError, match="interchange"):
        run_glm_scores(tmp_path / "iso.fa", tmp_path / "out", config)


# ---------------------------------------------------------------------------
# ecforest: parser + cache
# ---------------------------------------------------------------------------

def test_parse_ecforest_synthetic_csv(tmp_path):
    csv_path = tmp_path / "predictions.csv"
    csv_path.write_text(
        "gene_id,is_enzyme,EC,confidence\n"
        "Cgig_g1,True,4.1.1.31,0.97\n"
        "Mcry_g1,1,4.1.1.31,0.88\n"
        "Ococ_g1,False,,0.91\n"
        "Obas_g1,no,,\n"
    )
    preds = parse_ecforest(csv_path)
    assert preds["Cgig_g1"] == {"is_enzyme": True, "ec": "4.1.1.31", "confidence": 0.97}
    assert preds["Mcry_g1"]["is_enzyme"] is True
    assert preds["Ococ_g1"]["is_enzyme"] is False and preds["Ococ_g1"]["ec"] == ""
    assert preds["Obas_g1"]["confidence"] == 0.0  # unparseable -> 0.0


def test_parse_ecforest_without_stage1_column(tmp_path):
    # No is_enzyme column: a non-empty EC implies enzyme
    csv_path = tmp_path / "p.csv"
    csv_path.write_text("id,EC_number,conf\ng1,1.1.1.40,0.8\ng2,,0.9\n")
    preds = parse_ecforest(csv_path)
    assert preds["g1"]["is_enzyme"] is True and preds["g1"]["ec"] == "1.1.1.40"
    assert preds["g2"]["is_enzyme"] is False


def test_ecforest_cache_roundtrip(tmp_path):
    preds = {
        "g1": {"is_enzyme": True, "ec": "2.7.9.1", "confidence": 0.75},
        "g2": {"is_enzyme": False, "ec": "", "confidence": 0.6},
    }
    path = tmp_path / "ecforest_cache.tsv"
    save_cache(preds, path)
    assert path.read_text().startswith("gene_id\tis_enzyme\tec\tconfidence\n")
    assert load_cache(path) == preds


def test_run_ecforest_requires_configured_template(tmp_path):
    with pytest.raises(ValueError, match="not configured"):
        run_ecforest(tmp_path / "embeddings.tsv", tmp_path / "out", Config())


def test_run_ecforest_returns_newest_csv(tmp_path, monkeypatch):
    outdir = tmp_path / "out"

    def fake_run(cmd, **kwargs):
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "predictions.csv").write_text("gene_id,EC\ng1,4.1.1.31\n")
        return _FakeCompleted(0)

    monkeypatch.setattr(ecforest.subprocess, "run", fake_run)
    config = Config(
        ecforest_cmd="conda run -n enzyllm_env predict --input {input} --outdir {outdir}"
    )
    csv_path = run_ecforest(tmp_path / "embeddings.tsv", outdir, config)
    assert csv_path.name == "predictions.csv"


# ---------------------------------------------------------------------------
# ecforest: family consensus + advisory disagreement flags
# ---------------------------------------------------------------------------

def test_family_ec_consensus_majority():
    preds = {
        "g1": _pred("4.1.1.31"), "g2": _pred("4.1.1.31"),
        "g3": _pred("1.1.1.40"),
        "g4": _pred("", is_enzyme=False),  # abstains, does not dilute
    }
    assert family_ec_consensus({"g1", "g2", "g3", "g4"}, preds) == "4.1.1.31"


def test_family_ec_consensus_below_min_agree_and_tie():
    preds = {"g1": _pred("4.1.1.31"), "g2": _pred("1.1.1.40")}
    # 50/50 top tie -> never picked by dict order
    assert family_ec_consensus({"g1", "g2"}, preds) is None
    # majority present but below a stricter min_agree
    preds3 = {"g1": _pred("4.1.1.31"), "g2": _pred("4.1.1.31"), "g3": _pred("1.1.1.40")}
    assert family_ec_consensus({"g1", "g2", "g3"}, preds3, min_agree=0.8) is None
    assert family_ec_consensus({"g1", "g2", "g3"}, preds3, min_agree=0.5) == "4.1.1.31"


def test_family_ec_consensus_non_enzyme_family_is_none():
    preds = {"g1": _pred("", is_enzyme=False), "g2": _pred("", is_enzyme=False)}
    assert family_ec_consensus({"g1", "g2"}, preds) is None


def test_consensus_disagreements_are_advisory_flags():
    preds = {
        "g1": _pred("4.1.1.31"), "g2": _pred("4.1.1.31"),
        "g3": _pred("1.1.1.40"),           # disagrees -> flagged
        "g4": _pred("", is_enzyme=False),  # no usable EC -> never flagged
    }
    members = {"g1", "g2", "g3", "g4"}
    consensus = family_ec_consensus(members, preds)
    assert consensus_disagreements(members, preds, consensus) == ["g3"]
    assert consensus_disagreements(members, preds, None) == []


# ---------------------------------------------------------------------------
# ecforest: EC-switch events on the family tree (Fitch reuse)
# ---------------------------------------------------------------------------

EC_PREDS = {
    "Mcry_g1": _pred("1.1.1.40"),
    "Cgig_g1": _pred("4.1.1.31"),
    "Ococ_g1": _pred("4.1.1.31"),
    "Obas_g1": _pred("4.1.1.31"),
}


def test_ec_switch_events_detects_mcry_switch():
    events = ec_switch_events("R1_OG1", NWK, EC_PREDS)
    assert len(events) == 1
    ev = events[0]
    assert ev.child_ec == "1.1.1.40" and ev.parent_ec == "4.1.1.31"
    assert ev.clade_leaves == frozenset({"Mcry_g1"})
    assert ev.clade_id == "Mcry_g1"


def test_ec_switch_exclude_suppresses_fake_event():
    # Truncation safeguard, exactly as retargeting: excluding the switching
    # gene removes the event instead of reporting an artifact
    assert ec_switch_events("R1_OG1", NWK, EC_PREDS, exclude={"Mcry_g1"}) == []


def test_ec_switch_uniform_or_uninformative_family_has_no_events():
    uniform = {g: _pred("4.1.1.31") for g in
               ("Mcry_g1", "Cgig_g1", "Ococ_g1", "Obas_g1")}
    assert ec_switch_events("R1_OG1", NWK, uniform) == []
    non_enzyme = {g: _pred("", is_enzyme=False) for g in uniform}
    assert ec_switch_events("R1_OG1", NWK, non_enzyme) == []


# ---------------------------------------------------------------------------
# annotate_families.py helpers
# ---------------------------------------------------------------------------

def _write_run_dir(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "summary.tsv").write_text(
        "family_id\tround\tn_genes\tn_species\tgene_list\n"
        "R1_OG0000001\t1\t4\t4\tCgig_g1,Mcry_g1,Obas_g1,Ococ_g1\n"
        "R2_OG0000002\t2\t2\t2\tObas_g2,Ococ_g2\n"
    )
    (run_dir / "unplaced_proteins.fa").write_text(">Sp1_u1\nMKV\n>Sp2_u2\nMLA\n")
    return run_dir


def test_read_families_and_prep_fold_list(tmp_path):
    from annotate_families import prep_fold_list, read_families

    run_dir = _write_run_dir(tmp_path)
    families = read_families(run_dir)
    assert families["R1_OG0000001"] == {"Cgig_g1", "Mcry_g1", "Obas_g1", "Ococ_g1"}

    outdir = tmp_path / "annot"
    outdir.mkdir()
    prep_fold_list(run_dir, families, outdir)
    lines = (outdir / "fold_list.tsv").read_text().splitlines()
    assert lines[0] == "gene_id\trole\tfamily_id"
    unplaced = [l for l in lines[1:] if "\tunplaced\t" in l]
    reps = [l for l in lines[1:] if "\trepresentative\t" in l]
    assert len(unplaced) == 2 and len(reps) == 2
    # deterministic representative: first sorted member
    assert "Cgig_g1\trepresentative\tR1_OG0000001" in lines


def test_ec_layer_writes_consensus_and_events(tmp_path):
    from annotate_families import ec_layer, read_families

    run_dir = _write_run_dir(tmp_path)
    # family tree for R1_OG0000001 at the pipeline's on-disk location
    base = run_dir / "round_01" / "orthogroups" / "OG0000001"
    base.mkdir(parents=True)
    # balanced tree so Fitch resolves the switch direction unambiguously
    (base / "confirmed_tree.nwk").write_text(
        "((Cgig_g1:0.1,Ococ_g1:0.1):0.1,(Mcry_g1:0.4,Obas_g1:0.1):0.1);"
    )

    preds = dict(EC_PREDS)
    preds.update({"Obas_g2": _pred("2.7.9.1"), "Ococ_g2": _pred("2.7.9.1")})
    outdir = tmp_path / "annot"
    outdir.mkdir()
    families = read_families(run_dir)
    ec_layer(run_dir, families, preds, exclude=set(), min_agree=0.5, outdir=outdir)

    consensus_lines = (outdir / "ec_consensus.tsv").read_text().splitlines()
    row1 = dict(zip(consensus_lines[0].split("\t"), consensus_lines[1].split("\t")))
    assert row1["family_id"] == "R1_OG0000001"
    assert row1["consensus_ec"] == "4.1.1.31"
    assert row1["disagreeing_genes"] == "Mcry_g1"  # advisory flag, not a filter

    event_lines = (outdir / "ec_switch_events.tsv").read_text().splitlines()
    assert len(event_lines) == 2  # header + the Mcry switch
    assert event_lines[1].startswith("R1_OG0000001\t4.1.1.31\t1.1.1.40\t1\tMcry_g1")
