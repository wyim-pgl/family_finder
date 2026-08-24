"""Integrity gates around the final fragmentation rescan chain."""

import sys
import types
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

if "ete4" not in sys.modules:
    ete4_stub = types.ModuleType("ete4")
    ete4_stub.Tree = object
    sys.modules["ete4"] = ete4_stub

import pipeline  # noqa: E402
from config import Config  # noqa: E402


def test_final_output_refuses_a_gene_missing_from_both_sides(tmp_path):
    # Arrange: Sp2_b was in the original input, but is neither placed nor in
    # the remaining pool. Without the final equation this writes a plausible
    # one-gene summary and silently loses Sp2_b.
    families = {"R1_OG0000001": {"Sp1_a"}}

    # Act / Assert
    with pytest.raises(RuntimeError, match="Sp2_b"):
        pipeline._write_final_output(
            families, {}, {}, tmp_path, Config(),
            expected_gene_ids={"Sp1_a", "Sp2_b"},
        )
    assert not (tmp_path / "summary.tsv").exists()


def test_final_output_refuses_a_gene_in_two_families(tmp_path):
    families = {
        "R1_OG0000001": {"Sp1_a"},
        "R1_OG0000002": {"Sp1_a"},
    }

    with pytest.raises(RuntimeError, match="multiple families"):
        pipeline._write_final_output(
            families, {}, {}, tmp_path, Config(),
            expected_gene_ids={"Sp1_a"},
        )
    assert not (tmp_path / "summary.tsv").exists()


def _run_zero_work_resume(tmp_path, monkeypatch, scan_impl):
    """Resume after round 1 so only final-output/merge-scan wiring executes."""
    proteins = {"Sp1_a": "MPEP", "Sp2_b": "MPEP"}
    cds = {"Sp1_a": "ATG", "Sp2_b": "ATG"}
    pep_dir = tmp_path / "pep"
    cds_dir = tmp_path / "cds"
    outdir = tmp_path / "out"
    round_dir = outdir / "round_01"
    round_dir.mkdir(parents=True)
    (round_dir / "outlier_pool.fa").write_text(">Sp2_b\nMPEP\n")

    monkeypatch.setattr(pipeline, "start_manifest", lambda *a, **k: None)
    monkeypatch.setattr(
        pipeline, "load_species_tree",
        lambda path: SimpleNamespace(leaves=lambda: [
            SimpleNamespace(name="Sp1"), SimpleNamespace(name="Sp2")
        ]),
    )
    monkeypatch.setattr(pipeline, "compute_pairwise_distances", lambda tree: {})
    monkeypatch.setattr(pipeline, "validate_species_tree", lambda *a: [])
    monkeypatch.setattr(pipeline, "_check_id_agreement", lambda *a: {})
    monkeypatch.setattr(
        pipeline, "find_latest_checkpoint",
        lambda out: {"status": "completed", "round_number": 1},
    )
    monkeypatch.setattr(
        pipeline, "_load_round_families",
        lambda out, maximum: {"R1_OG0000001": {"Sp1_a"}},
    )

    from utils import seqio

    def fake_seq_map(path):
        return dict(cds if Path(path) == cds_dir else proteins)

    monkeypatch.setattr(seqio, "build_seq_map", fake_seq_map)
    written = []
    monkeypatch.setattr(
        pipeline, "_write_final_output",
        lambda *a, **k: written.append(k.get("expected_gene_ids")),
    )
    monkeypatch.setattr(pipeline, "scan_for_fragmented_families", scan_impl)
    finished = []
    monkeypatch.setattr(
        pipeline, "finish_manifest", lambda *a, **k: finished.append((a, k)),
    )

    config = Config(max_rounds=1, merge_scan=True, hmmer_rescue=False)
    pipeline.run(
        protein_dir=str(pep_dir), cds_dir=str(cds_dir),
        species_tree_path=str(tmp_path / "tree.nwk"), outdir=str(outdir),
        config=config, resume=True,
    )
    return proteins, written, finished


def test_merge_scan_without_rescue_receives_the_full_protein_pool(
        tmp_path, monkeypatch):
    # Before the fix, full_protein_pool was assigned only inside the rescue
    # branch; the resulting UnboundLocalError was swallowed and the manifest
    # still ended completed.
    seen = []

    def fake_scan(families, protein_pool, outdir, config):
        seen.append(dict(protein_pool))

    proteins, written, finished = _run_zero_work_resume(
        tmp_path, monkeypatch, fake_scan,
    )

    assert seen == [proteins]
    assert written == [set(proteins)]
    assert finished


def test_configured_merge_scan_failure_propagates_after_base_output(
        tmp_path, monkeypatch):
    # A failed optional analysis may leave the base summary recoverable, but it
    # must not be swallowed and followed by finish_manifest(..., completed).
    def failed_scan(*args, **kwargs):
        raise RuntimeError("rescan failed")

    with pytest.raises(RuntimeError, match="rescan failed"):
        _run_zero_work_resume(tmp_path, monkeypatch, failed_scan)
