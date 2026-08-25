"""Tests for steps/merge_candidates.py.

No foldseek installation is required. The subprocess wrappers are monkeypatched
and the mapping logic is exercised with synthetic cluster TSV fixtures.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from config import Config  # noqa: E402
from steps import merge_candidates  # noqa: E402
from steps.merge_candidates import (  # noqa: E402
    cluster_pairs_to_clusters,
    gene_to_family_index,
    load_all_round_families,
    load_confirmed_families,
    parse_foldseek_cluster_tsv,
    run_merge_candidates,
    summarize_clusters,
)


def test_load_confirmed_families_parses_round_artifact(tmp_path):
    path = tmp_path / "confirmed_families.tsv"
    path.write_text(
        "R1_OG0000001\tSpA_g1,SpB_g2\n"
        "R1_OG0000002\tSpC_g3\n"
    )

    got = load_confirmed_families(path)

    assert got == {
        "R1_OG0000001": {"SpA_g1", "SpB_g2"},
        "R1_OG0000002": {"SpC_g3"},
    }


def test_parse_foldseek_cluster_tsv_uses_first_two_columns_only(tmp_path):
    path = tmp_path / "clusters.tsv"
    path.write_text(
        "SpA_g1\tSpA_g1\t100\n"
        "SpA_g1\tSpB_g2\t200\n"
    )

    assert parse_foldseek_cluster_tsv(path) == [
        ("SpA_g1", "SpA_g1"),
        ("SpA_g1", "SpB_g2"),
    ]


def test_parse_foldseek_cluster_tsv_raises_on_truncated_row(tmp_path):
    # A partial final line (disk full, preemption mid-write) must be loud:
    # silently dropping it computes statuses from incomplete membership.
    path = tmp_path / "clusters.tsv"
    path.write_text("SpA_g1\tSpA_g1\nSpA_g1\n")

    with pytest.raises(ValueError, match="malformed foldseek cluster row"):
        parse_foldseek_cluster_tsv(path)


def test_load_all_round_families_unions_across_rounds(tmp_path):
    (tmp_path / "round_01").mkdir()
    (tmp_path / "round_02").mkdir()
    (tmp_path / "round_01" / "confirmed_families.tsv").write_text(
        "R1_OG0000001\tSpA_g1,SpB_g2\n"
    )
    (tmp_path / "round_02" / "confirmed_families.tsv").write_text(
        "R2_OG0000001\tSpC_g3\n"
    )

    got = load_all_round_families(tmp_path)

    assert got == {
        "R1_OG0000001": {"SpA_g1", "SpB_g2"},
        "R2_OG0000001": {"SpC_g3"},
    }


def test_load_all_round_families_rejects_cross_round_gene_reuse(tmp_path):
    (tmp_path / "round_01").mkdir()
    (tmp_path / "round_02").mkdir()
    (tmp_path / "round_01" / "confirmed_families.tsv").write_text(
        "R1_OG0000001\tSpA_g1\n"
    )
    (tmp_path / "round_02" / "confirmed_families.tsv").write_text(
        "R2_OG0000001\tSpA_g1\n"
    )

    with pytest.raises(ValueError, match="occurs in both"):
        load_all_round_families(tmp_path)


def test_load_all_round_families_requires_round_dirs(tmp_path):
    with pytest.raises(FileNotFoundError, match="round_"):
        load_all_round_families(tmp_path)


def test_cluster_pairs_to_clusters_rejects_member_in_two_clusters():
    with pytest.raises(ValueError, match="appears in both"):
        cluster_pairs_to_clusters([
            ("rep1", "gene1"),
            ("rep2", "gene1"),
        ])


def test_summarize_clusters_emits_candidate_and_tracks_unplaced_and_species():
    gene_to_family = gene_to_family_index({
        "R1_OG0000001": {"SpA_g1", "SpB_g2"},
        "R1_OG0000002": {"SpC_g3"},
    })
    clusters = [
        {"representative": "SpA_g1", "members": ["SpA_g1", "SpB_g2", "SpC_g3", "SpA_g9"]},
    ]

    records, summary = summarize_clusters(
        clusters, gene_to_family, max_cluster_size=10, species_delimiter="_"
    )

    assert records[0]["status"] == "candidate"
    assert records[0]["family_gene_counts"] == {
        "R1_OG0000001": 2,
        "R1_OG0000002": 1,
    }
    assert records[0]["unplaced_gene_ids"] == ["SpA_g9"]
    assert records[0]["species_counts"] == {"SpA": 2, "SpB": 1, "SpC": 1}
    assert summary["candidate_clusters"] == 1
    assert summary["candidate_unplaced_genes"] == 1


def test_summarize_clusters_overflows_runaway_cluster():
    gene_to_family = gene_to_family_index({"R1_OG0000001": {"SpA_g1"}})
    clusters = [{
        "representative": "SpA_g1",
        "members": ["SpA_g1", "SpB_g2", "SpC_g3"],
    }]

    records, summary = summarize_clusters(
        clusters, gene_to_family, max_cluster_size=2, species_delimiter="_"
    )

    assert records[0]["status"] == "overflow"
    assert summary["overflow_clusters"] == 1
    assert summary["candidate_clusters"] == 0


def test_run_merge_candidates_writes_candidate_and_overflow_outputs(tmp_path, monkeypatch):
    # Families deliberately split across two rounds: the cluster spanning
    # R1_OG0000001 (round 1) and R2_OG0000001 (round 2) is a candidate only
    # if the CUMULATIVE table is loaded — a single round's file would count
    # the other round's genes as unplaced and drop it.
    run_dir = tmp_path / "run"
    (run_dir / "round_01").mkdir(parents=True)
    (run_dir / "round_02").mkdir(parents=True)
    (run_dir / "round_01" / "confirmed_families.tsv").write_text(
        "R1_OG0000001\tSpA_g1,SpB_g2\n"
        "R1_OG0000003\tSpD_g4\n"
    )
    (run_dir / "round_02" / "confirmed_families.tsv").write_text(
        "R2_OG0000001\tSpC_g3\n"
    )

    out_dir = tmp_path / "merge_candidates"
    db_path = tmp_path / "all_proteins_db"

    monkeypatch.setattr(merge_candidates, "verify_3di_db", lambda path: 6)
    monkeypatch.setattr(
        merge_candidates, "run_foldseek_cluster",
        lambda db, cluster_db, tmp_dir, cfg: cluster_db,
    )

    def fake_createtsv(db, cluster_db, out_tsv, cfg):
        out_tsv.write_text(
            "SpA_g1\tSpA_g1\n"
            "SpA_g1\tSpB_g2\n"
            "SpA_g1\tSpC_g3\n"
            "SpD_g4\tSpD_g4\n"
            "SpD_g4\tSpE_g5\n"
            "SpD_g4\tSpF_g6\n"
            "SpD_g4\tSpG_g7\n"
        )
        return out_tsv

    monkeypatch.setattr(merge_candidates, "run_foldseek_createtsv", fake_createtsv)

    outputs = run_merge_candidates(
        run_dir,
        db_path,
        out_dir,
        Config(merge_candidate_max_cluster_size=3),
    )

    assert outputs["candidates"].name == "structural_merge_candidates.tsv"
    assert outputs["overflow"].name == "structural_merge_overflow.tsv"
    assert outputs["summary"].name == "structural_merge_summary.tsv"

    candidate_lines = outputs["candidates"].read_text().splitlines()
    overflow_lines = outputs["overflow"].read_text().splitlines()
    summary_lines = outputs["summary"].read_text().splitlines()

    assert candidate_lines[0].startswith("cluster_id\trepresentative_gene_id")
    assert candidate_lines[1].split("\t")[0] == "FSCLUST_000001"
    assert "R1_OG0000001,R2_OG0000001" in candidate_lines[1]
    assert "\t0\t-\t" in candidate_lines[1]

    assert overflow_lines[0].startswith("cluster_id\trepresentative_gene_id")
    assert overflow_lines[1].split("\t")[0] == "FSCLUST_000002"
    assert summary_lines[0] == "metric\tvalue"
    assert "verified_3di_entries\t6" in summary_lines
    assert "candidate_clusters\t1" in summary_lines
    assert "overflow_clusters\t1" in summary_lines
