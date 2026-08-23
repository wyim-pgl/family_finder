"""Run manifest + config-hash resume guard.

The pipeline resumes from per-round confirmed_families.tsv plus the newest
completed checkpoint, and nothing used to check that the configuration was
the one that produced them: change a threshold, resume, and the output
directory silently holds two configurations' results. These tests pin the
manifest that records what produced a run and the guard that refuses to
resume across a configuration change.

No conda env and no external tools: ete4 is stubbed, and the only binary any
test resolves is sys.executable.
"""

import json
import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# steps.prune imports ete4 at module level; stub it before importing anything
# that pulls in pipeline.
if "ete4" not in sys.modules:
    ete4_stub = types.ModuleType("ete4")
    ete4_stub.Tree = object
    sys.modules["ete4"] = ete4_stub

from config import Config  # noqa: E402
from utils import manifest as mf  # noqa: E402


def _tree(tmp_path, text="((A:0.1,B:0.1):0.1,C:0.2);", name="species_tree.nwk"):
    p = tmp_path / name
    p.write_text(text)
    return p


def _inputs(tmp_path):
    """Minimal protein/cds dirs with one FASTA each."""
    pep = tmp_path / "pep"
    cds = tmp_path / "cds"
    pep.mkdir()
    cds.mkdir()
    (pep / "Aaa.fa").write_text(">Aaa_g1\nMKV\n")
    (cds / "Aaa.fa").write_text(">Aaa_g1\nATGAAAGTT\n")
    return str(pep), str(cds)


# --------------------------------------------------------------------------
# 1. the hash split
# --------------------------------------------------------------------------

def test_hash_ignores_resource_and_location_knobs(tmp_path):
    """Threads, workers, sbatch extras, tool paths, device: same hash."""
    tree = _tree(tmp_path)
    a = Config()
    b = Config(
        n_workers=64,
        orthofinder_threads=64,
        hmmer_chunk_size=500,
        hmmer_chunk_concurrent=40,
        hmmer_chunk_sbatch_extra="--partition=cpu-s2-core-0 --account=cpu-s2-pgl-0",
        mafft_bin="/some/other/env/bin/mafft",
        fasttree_bin="/some/other/env/bin/FastTree",
        codeml_bin="/opt/paml4.10.10/bin/codeml",
        deeploc_device="cuda",
        orthofinder_reuse_blast_dir="/scratch/cache/WorkingDirectory",
    )
    assert mf.config_hash(a, str(tree)) == mf.config_hash(b, str(tree))


@pytest.mark.parametrize("field,value", [
    ("prune_criterion", "absolute"),
    ("prune_relative_threshold", 2.5),
    ("prune_score_floor", 1.0),
    ("distance_ratio_threshold", 4.0),
    ("min_species_for_pruning", 4),
    ("orthofinder_inflation", 1.5),
    ("orthofinder_search_program", "diamond_ultra_sens"),
    ("clustering_method", "sonicparanoid"),
    ("min_orthogroup_size", 5),
    ("min_family_size", 3),
    ("max_rounds", 12),
    ("convergence_threshold", 10),
    ("profile_assign_per_round", True),
    ("profile_margin_bits", 5.0),
    ("profile_min_coverage", 0.7),
    ("hmmer_rescue", True),
    ("hmmer_evalue", 1e-10),
    ("tree_builder", "iqtree"),
    ("mafft_strategy", "linsi"),
    ("species_delimiter", "|"),
])
def test_hash_changes_with_result_affecting_setting(tmp_path, field, value):
    tree = _tree(tmp_path)
    a = Config()
    b = Config(**{field: value})
    assert getattr(b, field) != getattr(a, field), "test would be vacuous"
    assert mf.config_hash(a, str(tree)) != mf.config_hash(b, str(tree)), field


def test_hash_covers_species_tree_contents_not_its_path(tmp_path):
    same_a = _tree(tmp_path, name="a.nwk")
    same_b = _tree(tmp_path, name="b.nwk")
    other = _tree(tmp_path, text="((A:0.9,B:0.1):0.1,C:0.2);", name="c.nwk")
    cfg = Config()
    assert mf.config_hash(cfg, str(same_a)) == mf.config_hash(cfg, str(same_b))
    assert mf.config_hash(cfg, str(same_a)) != mf.config_hash(cfg, str(other))


def test_list_field_order_does_not_change_the_hash(tmp_path):
    tree = _tree(tmp_path)
    a = Config(clustering_species_exclude=["CgigH", "Mcry"])
    b = Config(clustering_species_exclude=["Mcry", "CgigH"])
    assert mf.config_hash(a, str(tree)) == mf.config_hash(b, str(tree))
    c = Config(clustering_species_exclude=["CgigH"])
    assert mf.config_hash(a, str(tree)) != mf.config_hash(c, str(tree))


def test_every_config_field_is_hashed_unless_explicitly_excluded():
    """Fail-safe policy: the hash is a denylist, not an allowlist.

    A parameter added to Config later must be covered by the guard by
    default; forgetting to add it to an allowlist would make the guard
    quietly decorative for exactly the newest knob.
    """
    all_fields = set(Config.__dataclass_fields__)
    hashed = set(mf.hashed_config_fields(Config()))
    assert hashed | mf.EXCLUDED_FROM_HASH == all_fields
    assert not (hashed & mf.EXCLUDED_FROM_HASH)
    # the excluded set must not have drifted into covering real parameters
    assert mf.EXCLUDED_FROM_HASH <= all_fields


# --------------------------------------------------------------------------
# 2. the manifest
# --------------------------------------------------------------------------

def test_start_manifest_records_provenance(tmp_path):
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    cfg = Config(mafft_bin=sys.executable)

    data = mf.start_manifest(out, cfg, str(tree), pep, cds, resume=False)

    stored = json.loads((out / mf.MANIFEST_FILE).read_text())
    assert stored == data
    assert stored["status"] == "running"
    assert stored["resume_count"] == 0
    assert stored["run_id"]
    assert stored["config"]["prune_relative_threshold"] == 3.0
    assert stored["config_hash"] == mf.config_hash(cfg, str(tree))
    assert stored["species_tree"]["sha256"] == mf.sha256_file(tree)
    assert set(stored["git"]) >= {"commit", "dirty"}
    # input checksums, keyed by role then file name
    assert stored["inputs"]["protein_dir"]["files"]["Aaa.fa"] == mf.sha256_file(
        Path(pep) / "Aaa.fa"
    )
    assert stored["inputs"]["cds_dir"]["files"]["Aaa.fa"]
    # tool records go through config.resolve_tool
    assert stored["tools"]["mafft_bin"]["resolved"] == sys.executable
    assert stored["history"][0]["event"] == "start"


def test_unresolvable_tool_is_recorded_not_fatal(tmp_path):
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    cfg = Config(mafft_bin="/definitely/not/here/mafft")

    data = mf.start_manifest(out, cfg, str(tree), pep, cds, resume=False)
    rec = data["tools"]["mafft_bin"]
    assert rec["resolved"] is None
    assert "mafft_bin" in rec["error"]


def test_finish_manifest_sets_terminal_status(tmp_path):
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    mf.start_manifest(out, Config(), str(tree), pep, cds, resume=False)

    mf.finish_manifest(out, "completed")
    stored = json.loads((out / mf.MANIFEST_FILE).read_text())
    assert stored["status"] == "completed"
    assert stored["finished_at"]

    mf.finish_manifest(out, "failed", error="boom")
    stored = json.loads((out / mf.MANIFEST_FILE).read_text())
    assert stored["status"] == "failed"
    assert stored["error"] == "boom"


# --------------------------------------------------------------------------
# 3. the guard
# --------------------------------------------------------------------------

def test_resume_with_same_config_increments_resume_count(tmp_path):
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    cfg = Config()
    first = mf.start_manifest(out, cfg, str(tree), pep, cds, resume=False)

    second = mf.start_manifest(out, cfg, str(tree), pep, cds, resume=True)
    assert second["resume_count"] == 1
    assert second["run_id"] == first["run_id"]      # a resume continues the run
    assert second["status"] == "running"
    assert second["history"][-1]["event"] == "resume"

    third = mf.start_manifest(out, cfg, str(tree), pep, cds, resume=True)
    assert third["resume_count"] == 2


def test_resume_with_changed_config_names_the_settings(tmp_path):
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    mf.start_manifest(out, Config(), str(tree), pep, cds, resume=False)

    changed = Config(prune_relative_threshold=2.5, orthofinder_inflation=1.5)
    with pytest.raises(mf.ConfigMismatchError) as exc:
        mf.start_manifest(out, changed, str(tree), pep, cds, resume=True)

    msg = str(exc.value)
    assert "prune_relative_threshold" in msg
    assert "3.0" in msg and "2.5" in msg          # old value -> new value
    assert "orthofinder_inflation" in msg
    assert "1.2" in msg and "1.5" in msg
    assert "--allow-config-change" in msg
    # untouched settings are not listed as differences
    assert "min_family_size" not in msg
    # the refused resume must not have overwritten the stored manifest
    stored = json.loads((out / mf.MANIFEST_FILE).read_text())
    assert stored["config"]["prune_relative_threshold"] == 3.0
    assert stored["resume_count"] == 0


def test_resume_after_species_tree_edit_is_refused(tmp_path):
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    mf.start_manifest(out, Config(), str(tree), pep, cds, resume=False)

    tree.write_text("((A:0.5,B:0.1):0.1,C:0.2);")
    with pytest.raises(mf.ConfigMismatchError) as exc:
        mf.start_manifest(out, Config(), str(tree), pep, cds, resume=True)
    assert "species_tree" in str(exc.value)


def test_resume_after_input_edit_is_refused(tmp_path):
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    mf.start_manifest(out, Config(), str(tree), pep, cds, resume=False)

    (Path(pep) / "Bbb.fa").write_text(">Bbb_g1\nMKV\n")
    with pytest.raises(mf.ConfigMismatchError) as exc:
        mf.start_manifest(out, Config(), str(tree), pep, cds, resume=True)
    assert "protein_dir" in str(exc.value)
    assert "Bbb.fa" in str(exc.value)


def test_override_flag_allows_the_change_and_records_it(tmp_path):
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    mf.start_manifest(out, Config(), str(tree), pep, cds, resume=False)

    changed = Config(prune_relative_threshold=2.5)
    data = mf.start_manifest(
        out, changed, str(tree), pep, cds, resume=True, allow_config_change=True
    )
    assert data["config_override_used"] is True
    assert data["resume_count"] == 1
    assert data["config"]["prune_relative_threshold"] == 2.5
    assert data["config_hash"] == mf.config_hash(changed, str(tree))
    last = data["history"][-1]
    assert last["event"] == "resume"
    diffs = {d["setting"]: d for d in last["config_changes"]}
    assert diffs["prune_relative_threshold"]["old"] == 3.0
    assert diffs["prune_relative_threshold"]["new"] == 2.5
    # the flag stays recorded on later clean resumes: this output dir is mixed
    later = mf.start_manifest(out, changed, str(tree), pep, cds, resume=True)
    assert later["config_override_used"] is True


def test_resume_without_a_manifest_warns_but_proceeds(tmp_path, caplog):
    """Output dirs from before manifests existed must stay resumable."""
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    (out / "round_01").mkdir()

    with caplog.at_level("WARNING", logger="family_finder"):
        data = mf.start_manifest(out, Config(), str(tree), pep, cds, resume=True)
    assert data["unverified_resume"] is True
    assert any("no run_manifest.json" in r.message.lower()
               or "no run manifest" in r.message.lower()
               for r in caplog.records)


def test_changed_git_commit_is_recorded_but_does_not_refuse(tmp_path, monkeypatch):
    pep, cds = _inputs(tmp_path)
    tree = _tree(tmp_path)
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr(mf, "git_state", lambda: {"commit": "aaa", "dirty": False})
    mf.start_manifest(out, Config(), str(tree), pep, cds, resume=False)

    monkeypatch.setattr(mf, "git_state", lambda: {"commit": "bbb", "dirty": True})
    data = mf.start_manifest(out, Config(), str(tree), pep, cds, resume=True)
    assert data["git"]["commit"] == "bbb"
    assert data["history"][-1]["git"]["commit"] == "bbb"


# --------------------------------------------------------------------------
# 4. wiring
# --------------------------------------------------------------------------

def test_cli_exposes_the_opt_out_flag(monkeypatch):
    import family_finder

    base = ["family_finder.py", "--protein-dir", "p", "--cds-dir", "c",
            "--species-tree", "t", "--outdir", "o"]
    monkeypatch.setattr(sys, "argv", base)
    assert family_finder.parse_args().allow_config_change is False
    monkeypatch.setattr(sys, "argv", base + ["--allow-config-change"])
    assert family_finder.parse_args().allow_config_change is True


def test_pipeline_run_accepts_and_forwards_the_flag(monkeypatch, tmp_path):
    """pipeline.run must start the manifest before doing any round work."""
    import inspect
    import pipeline

    sig = inspect.signature(pipeline.run)
    assert "allow_config_change" in sig.parameters
    assert sig.parameters["allow_config_change"].default is False

    calls = {}

    def fake_start(outdir, config, species_tree_path, protein_dir, cds_dir,
                   resume=False, allow_config_change=False):
        calls["allow"] = allow_config_change
        calls["resume"] = resume
        raise RuntimeError("stop here")

    monkeypatch.setattr(pipeline, "start_manifest", fake_start)
    with pytest.raises(RuntimeError, match="stop here"):
        pipeline.run(
            protein_dir=str(tmp_path / "pep"), cds_dir=str(tmp_path / "cds"),
            species_tree_path=str(tmp_path / "t.nwk"), outdir=str(tmp_path / "o"),
            config=Config(), resume=True, allow_config_change=True,
        )
    assert calls == {"allow": True, "resume": True}
