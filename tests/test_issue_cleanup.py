"""Four residual defects from #49, #50 and #51, each with a failing scenario.

They are grouped because they share one shape: a wrong thing happens and the
run still looks fine. AGENTS.md calls that this repository's only failure mode.

The clustering one is the dangerous member of the set. #48 has just been fixed
so resume no longer loses genes - but if a run whose clustering DIED still
finishes "completed", that partial output becomes the baseline the repaired
resume reads from, and the loss returns through the front door. #49 says so in
its own body.
"""
import sys
import types
from pathlib import Path

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


# ---------------------------------------------------------------------------
# #49-1: clustering failure must not finish "completed"
# ---------------------------------------------------------------------------

def test_clustering_failure_is_a_terminal_failure():
    # Arrange: round 3 clustering dies. The command must not exit 0 with a
    # two-round summary.tsv and a "completed" manifest - that output is what
    # the next resume trusts. Asserting the class exists is not enough: the
    # defect was a handler that BROKE out of the loop, and the class can sit
    # unused beside a break.
    import inspect
    import re

    import pipeline
    source = inspect.getsource(pipeline.run)
    handlers = re.findall(r"except Exception[^\n]*:\n((?:[ \t]+[^\n]*\n)+)",
                          source)
    clustering = [h for h in handlers
                  if "ClusteringFailed" in h or "clustering" in h.lower()]

    # Act / Assert
    assert clustering, "no clustering exception handler found in pipeline.run"
    for handler in clustering:
        assert "raise ClusteringFailed" in handler
        assert not any(line.strip() == "break"
                       for line in handler.splitlines())


def test_clustering_failed_names_the_round_and_cause():
    # Arrange
    from pipeline import ClusteringFailed

    # Act
    exc = ClusteringFailed("orthofinder", 3, "disk full")

    # Assert: a batch log must say which round and why, not just "failed"
    text = str(exc)
    assert "orthofinder" in text and "3" in text and "disk full" in text


# ---------------------------------------------------------------------------
# #51-2: the low-yield convergence branch must hand rescue the CURRENT pool
# ---------------------------------------------------------------------------

def test_every_convergence_branch_updates_the_pool():
    # Arrange: three branches leave the loop; two assigned current_pool =
    # new_outlier_pool before breaking and one did not, so a low-yield
    # convergence handed rescue the PREVIOUS round's pool.
    import inspect

    import pipeline
    source = inspect.getsource(pipeline.run)
    start = source.index("if should_stop_iterating(")
    branch = source[start:source.index("break", start)]

    # Act / Assert
    assert "current_pool = new_outlier_pool" in branch


# ---------------------------------------------------------------------------
# #51-3: duplicate input ids are silently overwritten before any audit
# ---------------------------------------------------------------------------

def test_read_fasta_refuses_a_duplicate_id(tmp_path):
    # Arrange: two different sequences under one id. dict last-write-wins
    # loses one BEFORE clustering, so no downstream audit can see it - the
    # parse_orthogroups guard only watches OrthoFinder's output.
    from utils.seqio import read_fasta
    p = tmp_path / "dup.fa"
    p.write_text(">Sp1_a\nMPEP\n>Sp1_a\nMKKK\n")

    # Act / Assert
    try:
        read_fasta(str(p))
    except ValueError as e:
        assert "Sp1_a" in str(e)
    else:
        raise AssertionError("expected a duplicate id to be refused")


def test_read_fasta_dir_catches_ids_duplicated_across_files(tmp_path):
    # Arrange: the same id in two proteomes - the collision this repository
    # actually risks, since Cgig and CgigH annotate one genome twice
    from utils.seqio import read_fasta_dir
    (tmp_path / "a.fa").write_text(">Sp1_a\nMPEP\n")
    (tmp_path / "b.fa").write_text(">Sp1_a\nMKKK\n")

    # Act / Assert
    try:
        read_fasta_dir(str(tmp_path))
    except ValueError as e:
        assert "Sp1_a" in str(e)
    else:
        raise AssertionError("expected a cross-file duplicate to be refused")


def test_read_fasta_allows_an_identical_repeated_record(tmp_path):
    # Arrange: the same id carrying the SAME sequence loses nothing, so
    # refusing it would only break legitimate concatenated inputs
    from utils.seqio import read_fasta
    p = tmp_path / "same.fa"
    p.write_text(">Sp1_a\nMPEP\n>Sp1_a\nMPEP\n")

    # Act / Assert
    assert read_fasta(str(p)) == {"Sp1_a": "MPEP"}


# ---------------------------------------------------------------------------
# #50-3: summary_v3 has seven columns and build_supermatrix unpacks five
# ---------------------------------------------------------------------------

def test_supermatrix_reads_a_merged_summary(tmp_path):
    # Arrange: summary_v3.tsv adds merged_from and cluster_id, so an exact
    # five-way unpack raises ValueError on the table that now ships
    from build_supermatrix import single_copy_families
    p = tmp_path / "summary_v3.tsv"
    p.write_text(
        "family_id\tround\tn_genes\tn_species\tgene_list\tmerged_from\tcluster_id\n"
        "R1_OG0000001\t1\t2\t2\tSp1_a,Sp2_b\tR1_OG0000001+R1_OG0000009\tC0001\n"
    )

    # Act
    markers = single_copy_families(str(p), ["Sp1", "Sp2"])

    # Assert
    assert [m.family for m in markers] == ["R1_OG0000001"]


def test_supermatrix_still_reads_the_five_column_summary(tmp_path):
    # Arrange
    from build_supermatrix import single_copy_families
    p = tmp_path / "summary.tsv"
    p.write_text("family_id\tround\tn_genes\tn_species\tgene_list\n"
                 "R1_OG0000001\t1\t2\t2\tSp1_a,Sp2_b\n")

    # Act / Assert
    assert [m.family for m in single_copy_families(str(p), ["Sp1", "Sp2"])] \
        == ["R1_OG0000001"]
