"""Regression tests for the family-emission size gate (issue #11).

process_single_orthogroup must emit a family when >= min_family_size genes
survive pruning, instead of dissolving the whole orthogroup (old behavior:
compared against min_orthogroup_size=4 and recycled pruning survivors too).

External tools (MAFFT/FastTree/ete4) are stubbed so the gate logic runs
in isolation.
"""

import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# steps.prune imports ete4 at module level; stub it before importing pipeline
# so these tests run without the conda environment.
if "ete4" not in sys.modules:
    ete4_stub = types.ModuleType("ete4")
    ete4_stub.Tree = object
    sys.modules["ete4"] = ete4_stub

import pipeline  # noqa: E402
from config import Config  # noqa: E402


GENES = [f"Sp{i}_g{i}" for i in range(1, 6)]  # 5 genes, 5 "species"


def _run_og(tmp_path, monkeypatch, confirmed, outliers):
    """Drive process_single_orthogroup with stubbed align/tree/prune."""
    monkeypatch.setattr(
        pipeline, "align_protein",
        lambda seqs, out, cfg: Path(out).write_text("stub") or Path(out),
    )
    monkeypatch.setattr(pipeline, "codon_align", lambda *a, **k: (None, set()))
    monkeypatch.setattr(
        pipeline, "build_tree",
        lambda aln, out, cfg: Path(out).write_text("(stub);") or Path(out),
    )
    monkeypatch.setattr(
        pipeline, "prune_orthogroup", lambda *a, **k: (confirmed, outliers)
    )

    protein_pool = {g: "MSEQ" for g in GENES}
    cds_pool = {g: "ATGTCC" for g in GENES}
    args = ("OG0000001", list(GENES), protein_pool, cds_pool, {}, Config(), tmp_path)
    return pipeline.process_single_orthogroup(args)


def test_config_default_min_family_size():
    # Arrange / Act
    cfg = Config()
    # Assert
    assert cfg.min_family_size == 2
    assert cfg.min_orthogroup_size == 4  # processing floor unchanged


def test_emits_small_family_when_two_survive_pruning(tmp_path, monkeypatch):
    # Arrange: pruning keeps 2 of 5 (below old min_orthogroup_size=4)
    confirmed = {GENES[0], GENES[1]}
    outliers = set(GENES[2:])

    # Act
    og_id, got_confirmed, got_outliers = _run_og(
        tmp_path, monkeypatch, confirmed, outliers
    )

    # Assert: family emitted; only true outliers recycled (old code returned
    # None and recycled all 5, including the 2 that passed pruning)
    assert got_confirmed == confirmed
    assert got_outliers == outliers
    assert (tmp_path / "orthogroups" / og_id / "confirmed_proteins.fa").exists()


def test_dissolves_when_below_min_family_size(tmp_path, monkeypatch):
    # Arrange: only 1 survivor — cannot form a family
    confirmed = {GENES[0]}
    outliers = set(GENES[1:])

    # Act
    og_id, got_confirmed, got_outliers = _run_og(
        tmp_path, monkeypatch, confirmed, outliers
    )

    # Assert: whole OG recycled
    assert got_confirmed is None
    assert got_outliers == set(GENES)


def test_processing_floor_still_uses_min_orthogroup_size(tmp_path, monkeypatch):
    # Arrange: 3-gene OG is below the processing floor (4) — skipped before
    # align/tree, regardless of min_family_size
    small = GENES[:3]
    protein_pool = {g: "MSEQ" for g in small}
    cds_pool = {g: "ATGTCC" for g in small}
    args = ("OG0000002", small, protein_pool, cds_pool, {}, Config(), tmp_path)

    # Act
    og_id, got_confirmed, got_outliers = pipeline.process_single_orthogroup(args)

    # Assert
    assert got_confirmed is None
    assert got_outliers == set(small)
