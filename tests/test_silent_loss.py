"""Regression tests for silent gene-loss paths (issues #15, #2, #4).

Covers:
  - stop-codon filter: terminal stop followed by alignment gaps is NOT internal
  - codon_align returns (path_or_None, removed_ids) so callers can recycle
  - process_single_orthogroup recycles stop-removed genes into outliers
  - round-level gene-conservation audit detects leaked genes
  - find_latest_checkpoint returns the newest COMPLETED checkpoint
  - per-round confirmed_families.tsv write/reload round-trip (resume)
  - pep/CDS ID-agreement report (issue #2)

External tools (MAFFT/pal2nal/FastTree/ete4) are stubbed following the
conventions in test_size_gate.py so tests run without the conda environment.
"""

import logging
import sys
import types
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# steps.prune imports ete4 at module level; stub it before importing pipeline
if "ete4" not in sys.modules:
    ete4_stub = types.ModuleType("ete4")
    ete4_stub.Tree = object
    sys.modules["ete4"] = ete4_stub

import pipeline  # noqa: E402
from config import Config  # noqa: E402
from steps import align  # noqa: E402
from utils.checkpoint import find_latest_checkpoint, save_checkpoint  # noqa: E402
from utils.seqio import write_fasta  # noqa: E402


# ---------------------------------------------------------------------------
# Stop-codon filter (steps/align.py)
# ---------------------------------------------------------------------------

def test_terminal_stop_before_gaps_is_not_internal(tmp_path):
    # Arrange: aligned protein where the terminal stop is followed by gaps
    aln = tmp_path / "proteins.afa"
    write_fasta({"SpA_g1": "MPROT*----", "SpB_g2": "MPROTMAAAA"}, str(aln))
    cds = {"SpA_g1": "ATG" * 6, "SpB_g2": "ATG" * 10}

    # Act
    kept_path, kept_cds, removed = align._filter_internal_stops(aln, cds)

    # Assert: nothing removed — the stop is terminal, not internal
    assert removed == set()
    assert kept_path == aln
    assert kept_cds == cds


def test_genuine_internal_stop_is_removed(tmp_path):
    # Arrange
    aln = tmp_path / "proteins.afa"
    write_fasta(
        {"SpA_g1": "MP*ROT----", "SpB_g2": "MPROTMAAAA", "SpC_g3": "MPROTMAAAA"},
        str(aln),
    )
    cds = {g: "ATG" * 10 for g in ("SpA_g1", "SpB_g2", "SpC_g3")}

    # Act
    kept_path, kept_cds, removed = align._filter_internal_stops(aln, cds)

    # Assert
    assert removed == {"SpA_g1"}
    assert "SpA_g1" not in kept_cds
    assert kept_path is not None and kept_path != aln


# ---------------------------------------------------------------------------
# codon_align return contract: (path_or_None, removed_ids)
# ---------------------------------------------------------------------------

def test_codon_align_returns_removed_when_too_few_remain(tmp_path):
    # Arrange: 2 sequences, one with an internal stop → 1 remains (< 2)
    aln = tmp_path / "proteins.afa"
    write_fasta({"SpA_g1": "MP*ROT", "SpB_g2": "MPROTM"}, str(aln))
    cds = {"SpA_g1": "ATG" * 6, "SpB_g2": "ATG" * 6}

    # Act
    result = align.codon_align(aln, cds, tmp_path / "codon.afa", Config())

    # Assert: tuple contract, failure path still reports removed ids
    assert isinstance(result, tuple) and len(result) == 2
    codon_path, removed = result
    assert codon_path is None
    assert removed == {"SpA_g1"}


def test_codon_align_success_returns_path_and_removed(tmp_path, monkeypatch):
    # Arrange: 3 sequences, one with an internal stop; pal2nal stubbed
    aln = tmp_path / "proteins.afa"
    write_fasta(
        {"SpA_g1": "MP*ROT", "SpB_g2": "MPROTM", "SpC_g3": "MPROTM"}, str(aln)
    )
    cds = {g: "ATG" * 6 for g in ("SpA_g1", "SpB_g2", "SpC_g3")}
    monkeypatch.setattr(
        align.subprocess, "run",
        lambda *a, **k: SimpleNamespace(
            returncode=0, stdout=">SpB_g2\nATGATG\n>SpC_g3\nATGATG\n", stderr=""
        ),
    )

    # Act
    codon_path, removed = align.codon_align(
        aln, cds, tmp_path / "codon.afa", Config()
    )

    # Assert
    assert codon_path == tmp_path / "codon.afa"
    assert codon_path.exists()
    assert removed == {"SpA_g1"}


# ---------------------------------------------------------------------------
# process_single_orthogroup recycles stop-removed genes (pipeline.py)
# ---------------------------------------------------------------------------

GENES = [f"Sp{i}_g{i}" for i in range(1, 6)]  # 5 genes, 5 "species"


def test_stop_removed_genes_are_recycled_as_outliers(tmp_path, monkeypatch):
    # Arrange: codon_align drops GENES[4] (internal stop); pruning never sees
    # it (not a tree leaf) and flags nothing
    monkeypatch.setattr(
        pipeline, "align_protein",
        lambda seqs, out, cfg: Path(out).write_text("stub") or Path(out),
    )
    monkeypatch.setattr(
        pipeline, "codon_align", lambda *a, **k: (None, {GENES[4]})
    )
    monkeypatch.setattr(
        pipeline, "build_tree",
        lambda aln, out, cfg: Path(out).write_text("(stub);") or Path(out),
    )
    confirmed = set(GENES[:4])
    monkeypatch.setattr(
        pipeline, "prune_orthogroup", lambda *a, **k: (confirmed, set())
    )
    protein_pool = {g: "MSEQ" for g in GENES}
    cds_pool = {g: "ATGTCC" for g in GENES}
    args = ("OG0000001", list(GENES), protein_pool, cds_pool, {}, Config(), tmp_path)

    # Act
    og_id, got_confirmed, got_outliers = pipeline.process_single_orthogroup(args)

    # Assert: the stop-removed gene is recycled instead of vanishing
    assert got_confirmed == confirmed
    assert GENES[4] in got_outliers


# ---------------------------------------------------------------------------
# Round-level gene-conservation audit (pipeline.py)
# ---------------------------------------------------------------------------

def test_conservation_audit_detects_leak(caplog):
    # Arrange
    pool_ids = {"SpA_g1", "SpA_g2", "SpB_g3"}
    new_families = {"R1_OG0000000": {"SpA_g1"}}
    outliers = {"SpB_g3"}

    # Act
    with caplog.at_level(logging.ERROR, logger="family_finder"):
        leaked = pipeline._audit_gene_conservation(
            pool_ids, new_families, outliers, "_"
        )

    # Assert: leak detected, logged at ERROR with per-species breakdown
    assert leaked == {"SpA_g2"}
    assert any(r.levelno == logging.ERROR for r in caplog.records)
    assert "SpA=1" in caplog.text
    assert "SpA_g2" in caplog.text


def test_conservation_audit_clean_round(caplog):
    # Arrange
    pool_ids = {"SpA_g1", "SpB_g3"}
    new_families = {"R1_OG0000000": {"SpA_g1"}}
    outliers = {"SpB_g3"}

    # Act
    with caplog.at_level(logging.ERROR, logger="family_finder"):
        leaked = pipeline._audit_gene_conservation(
            pool_ids, new_families, outliers, "_"
        )

    # Assert
    assert leaked == set()
    assert not any(r.levelno == logging.ERROR for r in caplog.records)


# ---------------------------------------------------------------------------
# Checkpoint: newest COMPLETED checkpoint wins (utils/checkpoint.py)
# ---------------------------------------------------------------------------

def test_find_latest_checkpoint_skips_incomplete(tmp_path):
    # Arrange: round 1 completed, round 2 died mid-processing
    save_checkpoint(tmp_path / "round_01", 1, 100, "completed")
    save_checkpoint(tmp_path / "round_02", 2, 50, "processing")

    # Act
    cp = find_latest_checkpoint(tmp_path)

    # Assert: resume must restart from round 1, not see round 2
    assert cp is not None
    assert cp["round_number"] == 1
    assert cp["status"] == "completed"


def test_find_latest_checkpoint_returns_newest_completed(tmp_path):
    # Arrange
    save_checkpoint(tmp_path / "round_01", 1, 100, "completed")
    save_checkpoint(tmp_path / "round_02", 2, 50, "completed")

    # Act
    cp = find_latest_checkpoint(tmp_path)

    # Assert
    assert cp["round_number"] == 2


def test_find_latest_checkpoint_none_when_no_completed(tmp_path):
    # Arrange
    save_checkpoint(tmp_path / "round_01", 1, 100, "started")

    # Act / Assert
    assert find_latest_checkpoint(tmp_path) is None


# ---------------------------------------------------------------------------
# Per-round confirmed_families.tsv write / reload (resume path)
# ---------------------------------------------------------------------------

def test_round_families_roundtrip(tmp_path):
    # Arrange
    r1 = tmp_path / "round_01"
    r2 = tmp_path / "round_02"
    fams_r1 = {"R1_OG0000000": {"SpA_g1", "SpB_g2"}}
    fams_r2 = {"R2_OG0000000": {"SpC_g3", "SpD_g4"}}
    pipeline._write_round_families(r1, fams_r1)
    pipeline._write_round_families(r2, fams_r2)

    # Act
    loaded_all = pipeline._load_round_families(tmp_path, max_round=2)
    loaded_r1_only = pipeline._load_round_families(tmp_path, max_round=1)

    # Assert
    assert loaded_all == {**fams_r1, **fams_r2}
    assert loaded_r1_only == fams_r1


def test_resume_replays_profile_assignment_into_an_earlier_round_family(tmp_path):
    # Arrange: round 2 assigns SpC_u into an R1 family and therefore removes it
    # from round_02/outlier_pool.fa. confirmed_families.tsv for round 2 cannot
    # carry that R1 family without duplicating the family across round files.
    r1 = tmp_path / "round_01"
    r2 = tmp_path / "round_02"
    pipeline._write_round_families(r1, {"R1_OG0000000": {"SpA_g1"}})
    pipeline._write_round_families(r2, {"R2_OG0000000": {"SpB_g2"}})
    pa = r2 / "profile_assign"
    pa.mkdir()
    pa.joinpath("profile_assignments.tsv").write_text(
        "gene_id\tfamily_id\tfull_evalue\tfull_bits\tprofile_cov\tquery_cov\n"
        "SpC_u\tR1_OG0000000\t1e-50\t200.0\t0.8\t0.9\n"
    )

    # Act
    completed_r2 = pipeline._load_round_families(tmp_path, max_round=2)
    completed_r1 = pipeline._load_round_families(tmp_path, max_round=1)

    # Assert: the assignment appears only once its own round is completed.
    assert completed_r2["R1_OG0000000"] == {"SpA_g1", "SpC_u"}
    assert completed_r1["R1_OG0000000"] == {"SpA_g1"}


def test_resume_refuses_profile_assignment_to_an_unknown_or_future_family(tmp_path):
    # Arrange: replaying this row by silently skipping it would lose SpC_u,
    # because the completed round's outlier pool no longer contains the gene.
    r1 = tmp_path / "round_01"
    pipeline._write_round_families(r1, {"R1_OG0000000": {"SpA_g1"}})
    pa = r1 / "profile_assign"
    pa.mkdir()
    pa.joinpath("profile_assignments.tsv").write_text(
        "gene_id\tfamily_id\nSpC_u\tR2_OG0000000\n"
    )

    # Act / Assert
    with pytest.raises(ValueError, match="unknown/future family"):
        pipeline._load_round_families(tmp_path, max_round=1)


def test_resume_refuses_duplicate_family_rows_instead_of_overwriting(tmp_path):
    # Arrange: the old dict assignment retained only SpB_g2 from the second
    # row, so SpA_g1 disappeared before any conservation check saw it.
    r1 = tmp_path / "round_01"
    r1.mkdir()
    r1.joinpath("confirmed_families.tsv").write_text(
        "R1_OG0000000\tSpA_g1\nR1_OG0000000\tSpB_g2\n"
    )

    # Act / Assert
    with pytest.raises(ValueError, match="duplicate family id"):
        pipeline._load_round_families(tmp_path, max_round=1)


def test_resume_refuses_a_completed_round_without_its_outlier_pool(tmp_path):
    # Treating this as an empty/no checkpoint artifact used to advance to R2
    # with the original full input, reprocessing genes already placed in R1.
    (tmp_path / "round_01").mkdir()

    with pytest.raises(ValueError, match="cannot reconstruct"):
        pipeline._load_completed_round_pool(tmp_path, round_num=1)


def test_resume_accepts_an_empty_completed_outlier_pool(tmp_path):
    # All genes placed is a legitimate completed state and is distinct from a
    # missing artifact.
    round_dir = tmp_path / "round_01"
    round_dir.mkdir()
    (round_dir / "outlier_pool.fa").write_text("")

    assert pipeline._load_completed_round_pool(tmp_path, round_num=1) == {}


def test_a_duplicate_family_id_aborts_before_round_outputs_are_written(
        tmp_path, monkeypatch):
    # Arrange: a broken worker/result path returns the same OG id twice with
    # different members. Last-write-wins would manufacture a plausible round
    # while silently moving the first occurrence's genes back to the pool.
    proteins = {
        "Sp1_a": "MPEP", "Sp2_b": "MPEP",
        "Sp3_c": "MPEP", "Sp4_d": "MPEP",
        "Sp5_e": "MPEP", "Sp6_f": "MPEP",
    }
    cds = {gene: "ATGTAA" for gene in proteins}
    pep_dir = tmp_path / "pep"
    cds_dir = tmp_path / "cds"
    outdir = tmp_path / "out"

    monkeypatch.setattr(pipeline, "start_manifest", lambda *a, **k: None)
    monkeypatch.setattr(
        pipeline, "load_species_tree",
        lambda path: SimpleNamespace(
            leaves=lambda: [SimpleNamespace(name=f"Sp{i}")
                            for i in range(1, 7)]),
    )
    monkeypatch.setattr(pipeline, "compute_pairwise_distances", lambda tree: {})
    monkeypatch.setattr(pipeline, "validate_species_tree", lambda *a: [])
    monkeypatch.setattr(pipeline, "_check_id_agreement", lambda *a: {})
    monkeypatch.setattr(pipeline, "split_by_species", lambda *a, **k: None)
    monkeypatch.setattr(pipeline, "run_orthofinder", lambda *a, **k: tmp_path)
    monkeypatch.setattr(
        pipeline, "parse_orthogroups",
        lambda path: {"OG0000001": list(proteins)},
    )
    monkeypatch.setattr(
        pipeline, "parallel_map",
        lambda *a, **k: [
            ("OG0000001", {"Sp1_a", "Sp2_b"}, set()),
            ("OG0000001", {"Sp3_c", "Sp4_d"}, set()),
            ("OG0000001", {"Sp5_e", "Sp6_f"}, set()),
        ],
    )
    monkeypatch.setattr(pipeline, "finish_manifest", lambda *a, **k: None)

    from utils import seqio

    def fake_seq_map(path):
        return dict(cds if Path(path) == cds_dir else proteins)

    monkeypatch.setattr(seqio, "build_seq_map", fake_seq_map)
    config = Config(
        max_rounds=1, n_workers=1, hmmer_rescue=False,
        pseudogene_detection=False, merge_scan=False,
    )

    # Act / Assert: identify the exact family and collision count, and stop
    # before any durable round/final result can make the run look completed.
    with pytest.raises(
            RuntimeError,
            match=r"R1_OG0000001.*collision_count=2.*occurrence_count=3"):
        pipeline.run(
            protein_dir=str(pep_dir), cds_dir=str(cds_dir),
            species_tree_path=str(tmp_path / "tree.nwk"),
            outdir=str(outdir), config=config,
        )

    round_dir = outdir / "round_01"
    assert not (round_dir / "confirmed_families.tsv").exists()
    assert not (round_dir / "outlier_pool.fa").exists()
    assert not (round_dir / "round_stats.json").exists()
    assert not (outdir / "summary.tsv").exists()
    # The diagnostic checkpoint remains "processing", so resume ignores this
    # partial round and returns to the last genuinely completed checkpoint.
    assert find_latest_checkpoint(outdir) is None


# ---------------------------------------------------------------------------
# pep/CDS ID-agreement report (issue #2)
# ---------------------------------------------------------------------------

def test_id_agreement_report(tmp_path, caplog):
    # Arrange: SpA has 4 proteins but only 2 with matching CDS (50% < 95%);
    # SpB agrees fully
    protein_pool = {
        "SpA_g1": "M", "SpA_g2": "M", "SpA_g3": "M", "SpA_g4": "M",
        "SpB_g1": "M",
    }
    cds_pool = {"SpA_g1": "ATG", "SpA_g2": "ATG", "SpB_g1": "ATG"}

    # Act
    with caplog.at_level(logging.WARNING, logger="family_finder"):
        stats = pipeline._check_id_agreement(protein_pool, cds_pool, tmp_path, "_")

    # Assert: stats correct
    assert stats["SpA"][1:4] == (4, 2, 2)   # n_pep, n_cds, n_shared
    assert stats["SpA"][4] == pytest.approx(50.0)
    assert stats["SpB"][1:4] == (1, 1, 1)
    assert stats["SpB"][4] == pytest.approx(100.0)

    # Warning emitted only for the disagreeing species
    assert "SpA" in caplog.text
    assert "SpB" not in caplog.text

    # TSV written
    tsv = tmp_path / "id_agreement.tsv"
    assert tsv.exists()
    lines = tsv.read_text().strip().split("\n")
    assert lines[0] == "species\tn_pep\tn_cds\tn_shared\tpct_shared"
    assert len(lines) == 3
