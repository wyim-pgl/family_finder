"""One-way vote edges and their clusters, in the repo instead of a lost script.

The v3 family table was built from two files - `vote_edges.tsv` (one-way
from/to/votes/from_size/frac) and `fragmentation_clusters.tsv` (connected
components of those edges) - and **nothing in this repository wrote either
one**. `detect_merge_candidates` returns reciprocal PAIRS; the one-way edge
dump and the clustering came from a script that is gone, so the campaign that
merged 3,611 families could not be reproduced.

The numbers were already inside the detector. `detect_merge_candidates` counts
votes per ordered pair and then throws the counts away, keeping only pairs that
clear `merge_min_reciprocal`. Exporting the counts is an extra exit, not a new
algorithm.

Why the un-cut edges matter: the 15-species run's `vote_edges.tsv` has a
minimum frac of exactly 0.600, and 8,902 of 23,744 families have no outgoing
edge at all. The Mcry PEPC flagship's family votes for `R1_OG0019668` - a
member of the very cluster the PEPC pieces merged into - at frac 2/9 = 0.22,
so the threshold, not the evidence, is what kept it out. The tree is the
arbiter (see steps/cluster_validate.py); a nomination cut this high decides
before the arbiter ever sees the case.
"""
import sys
import types
from pathlib import Path

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import Config
from steps.profile_assign import (ProfileHit, detect_merge_candidates,
                                  fragmentation_clusters, vote_edges)


def _hit(gene, fam, *, full_e=1e-40, full_bits=200.0, hmm_len=100,
         seq_len=100, domains=None):
    if domains is None:
        domains = [(1, hmm_len, 1, seq_len, full_bits)]
    return ProfileHit(
        gene_id=gene, family_id=fam, full_evalue=full_e, full_bits=full_bits,
        domains=domains, hmm_len=hmm_len, seq_len=seq_len,
    )


# ---------------------------------------------------------------------------
# (A2) the profile-coverage floor is a config field, not a module constant
# ---------------------------------------------------------------------------

def test_merge_coverage_floor_is_read_from_config():
    # Arrange: cross-hits cover 0.15 of the profile - below the 0.2 default
    families = {"FamA": {"Sp1_a1"}, "FamB": {"Sp1_b1"}}
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB", domains=[(1, 15, 1, 100, 200.0)])],
        "Sp1_b1": [_hit("Sp1_b1", "FamA", domains=[(1, 15, 1, 100, 200.0)])],
    }

    # Act
    at_default = detect_merge_candidates(hits, families, Config())
    at_relaxed = detect_merge_candidates(
        hits, families, Config(merge_min_profile_cov=0.1))

    # Assert
    assert at_default == []
    assert at_relaxed == [("FamA", "FamB", 1.0)]


def test_merge_coverage_floor_default_is_the_measured_one():
    # Arrange / Act / Assert: 0.2, not the 0.5 the module constant held.
    # Reproduction sweep against the shipped 15sp vote_edges.tsv (200 sampled
    # families, 4,182 genes, all 23,744 profiles): exact agreement peaks at
    # 0.2 (191/200, 0 missing) and collapses at 0.5 (84/200, 19 missing).
    # Independently, rescue measured healthy full-length members covering a
    # median 0.22 of their profile. See config.py for the full rationale.
    assert Config().merge_min_profile_cov == 0.2


def test_merge_coverage_floor_gate_redirects_votes_not_just_drops_them():
    # Arrange: the failure mode the sweep exposed, with the numbers of the
    # gene that showed it - Sole_Spov3_chr6.01383, whose strongest hit was
    # R1_OG0001440 at 44 bits and cov 0.43. Too high a floor does not merely
    # lose that vote, it hands it to a weaker but better-covered family,
    # which is worse than silence.
    families = {"FamA": {"Sp1_a1"}, "TrueNeighbour": {"Sp1_t1"},
                "Weak": {"Sp1_w1"}}
    hits = {
        "Sp1_a1": [
            _hit("Sp1_a1", "TrueNeighbour", full_bits=44.0,
                 domains=[(1, 43, 1, 100, 44.0)]),      # strong, cov 0.43
            _hit("Sp1_a1", "Weak", full_bits=28.0,
                 domains=[(1, 90, 1, 100, 28.0)]),      # weak, cov 0.90
        ],
    }

    # Act
    at_high = vote_edges(hits, families, Config(merge_min_profile_cov=0.5))
    at_default = vote_edges(hits, families, Config())

    # Assert
    assert at_high == [("FamA", "Weak", 1, 1, 1.0)]
    assert at_default == [("FamA", "TrueNeighbour", 1, 1, 1.0)]


def test_merge_coverage_floor_enters_the_config_hash():
    # Arrange: a resume must refuse when this knob changed, the way it does
    # for every other scientific parameter (utils/manifest.py). The 15-species
    # run shipped merge_candidates at cov 0.2 AND 0.5; with the floor as a
    # module constant, neither file could be tied to a config hash.
    from utils.manifest import hashed_config_fields

    # Act
    base = hashed_config_fields(Config())
    changed = hashed_config_fields(Config(merge_min_profile_cov=0.5))

    # Assert
    assert "merge_min_profile_cov" in base
    assert base != changed


# ---------------------------------------------------------------------------
# (A3) one-way vote edges
# ---------------------------------------------------------------------------

def test_vote_edges_report_votes_size_and_fraction():
    # Arrange: both FamA members' best non-self hit is FamB
    families = {"FamA": {"Sp1_a1", "Sp2_a2"}, "FamB": {"Sp1_b1"}}
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=150.0)],
        "Sp2_a2": [_hit("Sp2_a2", "FamB", full_bits=140.0)],
    }

    # Act
    edges = vote_edges(hits, families, Config())

    # Assert
    assert edges == [("FamA", "FamB", 2, 2, 1.0)]


def test_vote_edges_keep_the_fraction_the_reciprocity_cut_would_discard():
    # Arrange: the shape that lost the Mcry PEPC flagship - 2 of 9 members
    # vote, frac 0.22, far below merge_min_reciprocal 0.6.
    members = {f"Sp{i}_a{i}" for i in range(9)}
    families = {"FamA": members, "FamB": {"Sp1_b1"}}
    hits = {
        "Sp0_a0": [_hit("Sp0_a0", "FamB", full_bits=150.0)],
        "Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=140.0)],
    }

    # Act
    edges = vote_edges(hits, families, Config())
    candidates = detect_merge_candidates(hits, families, Config())

    # Assert: the nomination cut drops it, the edge dump keeps it
    assert candidates == []
    assert edges == [("FamA", "FamB", 2, 9, round(2 / 9, 3))]


def test_vote_edges_can_be_filtered_by_fraction():
    # Arrange
    members = {f"Sp{i}_a{i}" for i in range(9)}
    families = {"FamA": members, "FamB": {"Sp1_b1"}}
    hits = {"Sp0_a0": [_hit("Sp0_a0", "FamB")],
            "Sp1_a1": [_hit("Sp1_a1", "FamB")]}

    # Act
    kept = vote_edges(hits, families, Config(), min_frac=0.2)
    dropped = vote_edges(hits, families, Config(), min_frac=0.6)

    # Assert
    assert len(kept) == 1
    assert dropped == []


def test_vote_edges_count_one_vote_per_gene():
    # Arrange: a gene hitting two other families votes only for the best
    families = {"FamA": {"Sp1_a1"}, "FamB": {"Sp1_b1"}, "FamC": {"Sp1_c1"}}
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=300.0),
                   _hit("Sp1_a1", "FamC", full_bits=200.0)],
    }

    # Act
    edges = vote_edges(hits, families, Config())

    # Assert
    assert edges == [("FamA", "FamB", 1, 1, 1.0)]


def test_vote_edges_stay_within_one_vote_per_gene():
    # Arrange: the invariant the vote-inflation trap violated. domtblout is
    # ordered by PROFILE, not by gene, so a gene-contiguous assumption
    # double-counts. The bound is votes <= members, NOT edges <= families -
    # see test_a_family_may_emit_several_edges for why the latter is false.
    families = {"FamA": {"Sp1_a1"}, "FamB": {"Sp1_b1"}, "FamC": {"Sp1_c1"}}
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB"), _hit("Sp1_a1", "FamC")],
        "Sp1_b1": [_hit("Sp1_b1", "FamA"), _hit("Sp1_b1", "FamC")],
        "Sp1_c1": [_hit("Sp1_c1", "FamA")],
    }

    # Act
    edges = vote_edges(hits, families, Config())

    # Assert
    assert sum(v for _, _, v, _, _ in edges) <= sum(len(m) for m in families.values())


def test_vote_edges_ignore_self_hits_and_unknown_genes():
    # Arrange
    families = {"FamA": {"Sp1_a1"}}
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamA", full_bits=400.0)],
        "Sp9_zz": [_hit("Sp9_zz", "FamA", full_bits=400.0)],  # not a member
    }

    # Act
    edges = vote_edges(hits, families, Config())

    # Assert
    assert edges == []


# ---------------------------------------------------------------------------
# (A3) connected components over the edges
# ---------------------------------------------------------------------------

def test_fragmentation_clusters_are_connected_components():
    # Arrange: A->B, B->C is one component; D->E is another
    edges = [("A", "B", 5, 5, 1.0), ("B", "C", 4, 4, 1.0),
             ("D", "E", 3, 3, 1.0)]

    # Act
    clusters = fragmentation_clusters(edges)

    # Assert
    assert clusters == [["A", "B", "C"], ["D", "E"]]


def test_fragmentation_clusters_drop_singletons():
    # Arrange: an edge onto itself leaves nothing to merge
    edges = [("A", "A", 5, 5, 1.0)]

    # Act
    clusters = fragmentation_clusters(edges)

    # Assert
    assert clusters == []


def test_fragmentation_clusters_are_ordered_largest_first():
    # Arrange
    edges = [("D", "E", 1, 1, 1.0), ("A", "B", 1, 1, 1.0),
             ("B", "C", 1, 1, 1.0)]

    # Act
    clusters = fragmentation_clusters(edges)

    # Assert: biggest component first, members sorted - stable cluster ids
    assert clusters == [["A", "B", "C"], ["D", "E"]]


def test_fragmentation_clusters_merge_through_direction():
    # Arrange: the relation is undirected - B->A and C->A tie A, B, C together
    edges = [("B", "A", 1, 1, 1.0), ("C", "A", 1, 1, 1.0)]

    # Act
    clusters = fragmentation_clusters(edges)

    # Assert
    assert clusters == [["A", "B", "C"]]


# ---------------------------------------------------------------------------
# (A3) the two files the v3 campaign consumed, now written by this repository
# ---------------------------------------------------------------------------

def test_write_vote_edges_schema(tmp_path):
    # Arrange
    from pipeline import write_vote_edges

    # Act
    write_vote_edges([("R2_OG0000359", "R1_OG0019668", 2, 9, 0.222)], tmp_path)
    rows = (tmp_path / "vote_edges.tsv").read_text().splitlines()

    # Assert
    assert rows[0].split("\t") == ["from", "to", "votes", "from_size", "frac"]
    assert rows[1] == "R2_OG0000359\tR1_OG0019668\t2\t9\t0.222"


def test_write_vote_edges_writes_the_file_when_empty(tmp_path):
    # Arrange: a missing file cannot be told apart from a scan that never ran
    from pipeline import write_vote_edges

    # Act
    path = write_vote_edges([], tmp_path)

    # Assert
    assert path.exists()
    assert path.read_text().splitlines() == ["from\tto\tvotes\tfrom_size\tfrac"]


def test_write_fragmentation_clusters_numbers_and_sizes(tmp_path):
    # Arrange
    from pipeline import write_fragmentation_clusters
    families = {"A": {"g1", "g2"}, "B": {"g3"}, "D": {"g4"}, "E": {"g5"}}

    # Act
    write_fragmentation_clusters([["A", "B"], ["D", "E"]], families, tmp_path)
    rows = (tmp_path / "fragmentation_clusters.tsv").read_text().splitlines()

    # Assert
    assert rows[0].split("\t") == ["cluster_id", "n_families", "n_genes",
                                   "families"]
    assert rows[1] == "C0001\t2\t3\tA,B"
    assert rows[2] == "C0002\t2\t2\tD,E"


def test_write_fragmentation_clusters_writes_the_file_when_empty(tmp_path):
    # Arrange
    from pipeline import write_fragmentation_clusters

    # Act
    path = write_fragmentation_clusters([], {}, tmp_path)

    # Assert
    assert path.exists()


# ---------------------------------------------------------------------------
# (A3) the scan actually emits both files - the wiring, not just the writers
# ---------------------------------------------------------------------------

def _stub_scan(monkeypatch, tmp_path, hits):
    """Point scan_for_fragmented_families at synthetic hits, no hmmsearch."""
    import steps.profile_assign as pa

    (tmp_path / "merge_scan").mkdir(parents=True, exist_ok=True)
    (tmp_path / "merge_scan" / "all_families.hmm").write_text("stub")
    monkeypatch.setattr(pa, "_hmmsearch_domtblout",
                        lambda db, q, out, cfg: Path(out))
    monkeypatch.setattr(pa, "parse_domtblout", lambda path: hits)


def test_scan_emits_edges_and_clusters(tmp_path, monkeypatch):
    # Arrange: FamA's two members both vote FamB - one edge, one cluster
    from pipeline import scan_for_fragmented_families
    families = {"FamA": {"Sp1_a1", "Sp2_a2"}, "FamB": {"Sp1_b1"}}
    pool = {g: "MPEP" for fam in families.values() for g in fam}
    hits = {"Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=150.0)],
            "Sp2_a2": [_hit("Sp2_a2", "FamB", full_bits=140.0)]}
    _stub_scan(monkeypatch, tmp_path, hits)

    # Act
    scan_for_fragmented_families(families, pool, tmp_path, Config())

    # Assert
    edges = (tmp_path / "vote_edges.tsv").read_text().splitlines()
    clusters = (tmp_path / "fragmentation_clusters.tsv").read_text().splitlines()
    assert edges[1] == "FamA\tFamB\t2\t2\t1.000"
    assert clusters[1] == "C0001\t2\t3\tFamA,FamB"


def test_scan_writes_both_files_when_nothing_votes(tmp_path, monkeypatch):
    # Arrange: no cross-family hits at all
    from pipeline import scan_for_fragmented_families
    families = {"FamA": {"Sp1_a1"}}
    _stub_scan(monkeypatch, tmp_path, {})

    # Act
    scan_for_fragmented_families(families, {"Sp1_a1": "MPEP"}, tmp_path,
                                 Config())

    # Assert: "no fragmentation" is a result, not a missing file
    assert (tmp_path / "vote_edges.tsv").exists()
    assert (tmp_path / "fragmentation_clusters.tsv").exists()


def test_scan_refuses_to_write_inflated_edges(tmp_path, monkeypatch):
    # Arrange: inflation looks like a family casting more votes than it has
    # members - frac above 1.0. That, not the edge count, is the symptom of
    # reading profile-ordered domtblout as gene-contiguous.
    import steps.profile_assign as pa
    from pipeline import scan_for_fragmented_families
    families = {"FamA": {"Sp1_a1"}, "FamB": {"Sp1_b1"}}
    _stub_scan(monkeypatch, tmp_path, {})
    # scan_for_fragmented_families imports vote_edges from the module at call
    # time, so the module attribute is the binding that matters.
    monkeypatch.setattr(
        pa, "vote_edges",
        lambda h, f, c, min_frac=0.0: [("FamA", "FamB", 4, 1, 4.0)])

    # Act
    scan_for_fragmented_families(families, {"Sp1_a1": "MPEP"}, tmp_path,
                                 Config())

    # Assert: the file exists and is empty rather than plausibly wrong
    rows = (tmp_path / "vote_edges.tsv").read_text().splitlines()
    assert rows == ["from\tto\tvotes\tfrom_size\tfrac"]


def test_scan_writes_more_edges_than_families(tmp_path, monkeypatch):
    # Arrange: the real-data ratio. 200 sampled families produced 353 edges,
    # so a guard keyed on "edges > families" would have emptied the output of
    # every production run. Two families, three edges, all legitimate.
    import steps.profile_assign as pa
    from pipeline import scan_for_fragmented_families
    families = {"FamA": {"Sp1_a1", "Sp2_a2", "Sp3_a3"}, "FamB": {"Sp1_b1"}}
    pool = {g: "MPEP" for fam in families.values() for g in fam}
    _stub_scan(monkeypatch, tmp_path, {})
    monkeypatch.setattr(
        pa, "vote_edges",
        lambda h, f, c, min_frac=0.0: [("FamA", "FamB", 2, 3, 0.667),
                                       ("FamA", "FamC", 1, 3, 0.333),
                                       ("FamB", "FamA", 1, 1, 1.0)])

    # Act
    scan_for_fragmented_families(families, pool, tmp_path, Config())

    # Assert
    rows = (tmp_path / "vote_edges.tsv").read_text().splitlines()
    assert len(rows) == 4          # header + 3 edges, none refused
    assert len(rows) - 1 > len(families)


# ---------------------------------------------------------------------------
# the inflation invariant, corrected
# ---------------------------------------------------------------------------

def test_a_family_may_emit_several_edges():
    # Arrange: the shape that broke the first guard. A family's members do
    # not have to agree on a neighbour - measured on the 15sp sample, 200
    # families produced 353 edges (1.77 each). "more edges than families" is
    # therefore NORMAL, and a guard using it empties the output on real data.
    families = {"FamA": {"Sp1_a1", "Sp2_a2", "Sp3_a3"},
                "FamB": {"Sp1_b1"}, "FamC": {"Sp1_c1"}}
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=150.0)],
        "Sp2_a2": [_hit("Sp2_a2", "FamB", full_bits=140.0)],
        "Sp3_a3": [_hit("Sp3_a3", "FamC", full_bits=130.0)],
    }

    # Act
    edges = vote_edges(hits, families, Config())

    # Assert
    assert len(edges) == 2 > 0
    assert {(a, b) for a, b, *_ in edges} == {("FamA", "FamB"),
                                              ("FamA", "FamC")}


def test_total_votes_never_exceed_the_voting_population():
    # Arrange: the invariant that actually holds - one gene, one vote. This
    # is what catches treating profile-ordered domtblout as gene-contiguous.
    families = {"FamA": {"Sp1_a1", "Sp2_a2"}, "FamB": {"Sp1_b1"}}
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=300.0),
                   _hit("Sp1_a1", "FamB", full_bits=200.0)],   # 2 rows, 1 gene
        "Sp2_a2": [_hit("Sp2_a2", "FamB", full_bits=140.0)],
    }

    # Act
    edges = vote_edges(hits, families, Config())

    # Assert
    n_members = sum(len(m) for m in families.values())
    assert sum(v for _, _, v, _, _ in edges) <= n_members


def test_no_source_family_casts_more_votes_than_it_has_members():
    # Arrange: frac > 1.0 is the visible symptom of inflation
    families = {"FamA": {f"Sp{i}_a{i}" for i in range(5)}, "FamB": {"Sp1_b1"}}
    hits = {f"Sp{i}_a{i}": [_hit(f"Sp{i}_a{i}", "FamB")] for i in range(5)}

    # Act
    edges = vote_edges(hits, families, Config())

    # Assert
    for _, _, votes, from_size, frac in edges:
        assert votes <= from_size
        assert frac <= 1.0
