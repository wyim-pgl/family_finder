"""Tests for steps/profile_assign.py (issue #13).

No HMMER binaries required: hmmsearch output is synthesized as domtblout
text, and the hmmbuild/hmmpress/hmmsearch wrappers are monkeypatched.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from config import Config  # noqa: E402
from steps import profile_assign as pa  # noqa: E402
from steps.profile_assign import (  # noqa: E402
    ProfileHit,
    assign,
    build_profiles,
    detect_merge_candidates,
    parse_domtblout,
    run_profile_assignment,
)


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

def _dom_row(gene, fam, *, tlen, qlen, full_e, full_bits,
             dom_bits, hf, ht, af, at, i_evalue="1e-20"):
    """One synthetic domtblout data row (23 whitespace-delimited columns)."""
    return (
        f"{gene} - {tlen} {fam} - {qlen} {full_e} {full_bits} 0.1 "
        f"1 2 {i_evalue} {i_evalue} {dom_bits} 0.1 "
        f"{hf} {ht} {af} {at} {af} {at} 0.95 -\n"
    )


def _write_domtbl(tmp_path, rows, name="hits.domtblout"):
    path = tmp_path / name
    header = "# target name  ...\n#---\n"
    path.write_text(header + "".join(rows))
    return path


def _hit(gene, fam, *, full_e=1e-40, full_bits=200.0, hmm_len=100,
         seq_len=100, domains=None):
    if domains is None:
        # Single domain covering the full profile and query
        domains = [(1, hmm_len, 1, seq_len, full_bits)]
    return ProfileHit(
        gene_id=gene, family_id=fam, full_evalue=full_e, full_bits=full_bits,
        domains=domains, hmm_len=hmm_len, seq_len=seq_len,
    )


# ---------------------------------------------------------------------------
# parse_domtblout
# ---------------------------------------------------------------------------

def test_parse_domtblout_merges_overlapping_domains(tmp_path):
    # Arrange: two overlapping domains — hmm 1-60 and 41-100 (union = 100),
    # ali 1-50 and 45-80 (union = 80 of tlen 200). Naive summing would give
    # profile_cov 1.21 and query_cov 0.43.
    rows = [
        _dom_row("Sp1_g1", "R1_OG0000001", tlen=200, qlen=100,
                 full_e="1e-50", full_bits=180.5, dom_bits=100.0,
                 hf=1, ht=60, af=1, at=50),
        _dom_row("Sp1_g1", "R1_OG0000001", tlen=200, qlen=100,
                 full_e="1e-50", full_bits=180.5, dom_bits=80.0,
                 hf=41, ht=100, af=45, at=80),
    ]
    path = _write_domtbl(tmp_path, rows)

    # Act
    hits_by_gene = parse_domtblout(path)

    # Assert: one aggregated hit, coverages from merged intervals
    assert set(hits_by_gene) == {"Sp1_g1"}
    (hit,) = hits_by_gene["Sp1_g1"]
    assert hit.family_id == "R1_OG0000001"
    assert hit.full_bits == pytest.approx(180.5)
    assert hit.full_evalue == pytest.approx(1e-50)
    assert len(hit.domains) == 2
    assert hit.profile_cov == pytest.approx(1.0)
    assert hit.query_cov == pytest.approx(80 / 200)
    assert hit.nbits == pytest.approx(180.5 / 100)


def test_parse_domtblout_sorts_hits_by_full_bits_desc(tmp_path):
    # Arrange: weaker family appears FIRST in the file
    rows = [
        _dom_row("Sp1_g1", "R2_OG0000009", tlen=100, qlen=100,
                 full_e="1e-10", full_bits=80.0, dom_bits=80.0,
                 hf=1, ht=100, af=1, at=100),
        _dom_row("Sp1_g1", "R1_OG0000001", tlen=100, qlen=100,
                 full_e="1e-40", full_bits=250.0, dom_bits=250.0,
                 hf=1, ht=100, af=1, at=100),
    ]
    path = _write_domtbl(tmp_path, rows)

    # Act
    hits = parse_domtblout(path)["Sp1_g1"]

    # Assert
    assert [h.family_id for h in hits] == ["R1_OG0000001", "R2_OG0000009"]


# ---------------------------------------------------------------------------
# assign: new-gene acceptance criteria
# ---------------------------------------------------------------------------

def test_assign_accepts_best_hit_meeting_all_criteria():
    # Arrange: clear best (margin 100 bits), full coverage
    hits = {"Sp1_g1": [
        _hit("Sp1_g1", "FamA", full_bits=200.0),
        _hit("Sp1_g1", "FamB", full_bits=100.0),
    ]}

    # Act
    new_assignments, ambiguous, moves = assign(hits, {"Sp1_g1": None}, Config())

    # Assert
    assert new_assignments == {"Sp1_g1": "FamA"}
    assert ambiguous == []
    assert moves == {}


def test_assign_rejects_on_profile_coverage():
    # Arrange: only 30 of 100 HMM positions covered (< 0.5)
    hits = {"Sp1_g1": [
        _hit("Sp1_g1", "FamA", full_bits=200.0,
             domains=[(1, 30, 1, 90, 200.0)]),
    ]}

    # Act
    new_assignments, ambiguous, _ = assign(hits, {"Sp1_g1": None}, Config())

    # Assert: rejected outright — not assigned, not ambiguous
    assert new_assignments == {}
    assert ambiguous == []


def test_assign_rejects_on_query_coverage():
    # Arrange: 30 of 100 query residues covered (< 0.4)
    hits = {"Sp1_g1": [
        _hit("Sp1_g1", "FamA", full_bits=200.0,
             domains=[(1, 100, 1, 30, 200.0)]),
    ]}

    # Act
    new_assignments, ambiguous, _ = assign(hits, {"Sp1_g1": None}, Config())

    # Assert
    assert new_assignments == {}
    assert ambiguous == []


def test_assign_rejects_on_evalue_gate():
    # Arrange: weak significance despite coverage
    hits = {"Sp1_g1": [_hit("Sp1_g1", "FamA", full_e=1e-3, full_bits=20.0)]}

    # Act
    new_assignments, ambiguous, _ = assign(hits, {"Sp1_g1": None}, Config())

    # Assert
    assert new_assignments == {}
    assert ambiguous == []


def test_assign_margin_only_failure_goes_ambiguous():
    # Arrange: both candidates pass gates; margin 5 < max(10, 0.05*200=10)
    hits = {"Sp1_g1": [
        _hit("Sp1_g1", "FamA", full_bits=200.0),
        _hit("Sp1_g1", "FamB", full_bits=195.0),
    ]}

    # Act
    new_assignments, ambiguous, _ = assign(hits, {"Sp1_g1": None}, Config())

    # Assert
    assert new_assignments == {}
    assert len(ambiguous) == 1
    entry = ambiguous[0]
    assert entry["gene"] == "Sp1_g1"
    assert len(entry["candidates"]) == 2
    fams = {c["family_id"] for c in entry["candidates"]}
    assert fams == {"FamA", "FamB"}
    for cand in entry["candidates"]:
        assert "full_bits" in cand
        assert "profile_cov" in cand
        assert "query_cov" in cand


def test_assign_margin_frac_dominates_for_large_bits():
    # Arrange: margin 30 bits passes the absolute floor (10) but fails the
    # fractional one (0.05 * 1000 = 50)
    hits = {"Sp1_g1": [
        _hit("Sp1_g1", "FamA", full_bits=1000.0),
        _hit("Sp1_g1", "FamB", full_bits=970.0),
    ]}

    # Act
    new_assignments, ambiguous, _ = assign(hits, {"Sp1_g1": None}, Config())

    # Assert
    assert new_assignments == {}
    assert len(ambiguous) == 1


def test_assign_ambiguous_candidates_capped_at_five():
    # Arrange: 7 near-tied candidate families
    hits = {"Sp1_g1": [
        _hit("Sp1_g1", f"Fam{i}", full_bits=200.0 - i) for i in range(7)
    ]}

    # Act
    _, ambiguous, _ = assign(hits, {"Sp1_g1": None}, Config())

    # Assert
    assert len(ambiguous) == 1
    assert len(ambiguous[0]["candidates"]) == 5


@pytest.mark.parametrize("first_fam,second_fam", [
    ("R10_OG0000001", "R1_OG0000002"),   # lexicographically first wins in legacy code
    ("R1_OG0000002", "R10_OG0000001"),   # ... and reversed order
])
def test_equal_full_bits_tie_never_assigned(tmp_path, first_fam, second_fam):
    """Regression for the R10/R1 lexicographic tie-break bug.

    hmmer_rescue broke exact E-value ties (underflow to 0 or 2-sig-fig
    rounding) by HMM-db line order, where sorted() puts R10_* before
    R1_*. Equal full_bits must go to ambiguous regardless of family id
    or file order — never to whichever came first.
    """
    # Arrange: identical bit scores, both passing every other gate
    rows = [
        _dom_row("Ococ_g1", first_fam, tlen=100, qlen=100,
                 full_e="0", full_bits=3135.2, dom_bits=3135.2,
                 hf=1, ht=100, af=1, at=100),
        _dom_row("Ococ_g1", second_fam, tlen=100, qlen=100,
                 full_e="0", full_bits=3135.2, dom_bits=3135.2,
                 hf=1, ht=100, af=1, at=100),
    ]
    hits_by_gene = parse_domtblout(_write_domtbl(tmp_path, rows))

    # Act
    new_assignments, ambiguous, _ = assign(
        hits_by_gene, {"Ococ_g1": None}, Config()
    )

    # Assert: NOT assigned; reported ambiguous with both families
    assert new_assignments == {}
    assert len(ambiguous) == 1
    fams = {c["family_id"] for c in ambiguous[0]["candidates"]}
    assert fams == {first_fam, second_fam}


# ---------------------------------------------------------------------------
# assign: reassignment screen
# ---------------------------------------------------------------------------

def test_reassignment_candidate_detected_via_nbits_margin():
    # Arrange: current FamA nbits 0.5; FamB nbits 0.7 (margin 0.2 >= 0.15)
    # with profile_cov 1.0 (>= 0.6)
    hits = {"Sp1_g1": [
        _hit("Sp1_g1", "FamB", full_bits=70.0, hmm_len=100),
        _hit("Sp1_g1", "FamA", full_bits=50.0, hmm_len=100),
    ]}

    # Act
    new_assignments, ambiguous, moves = assign(
        hits, {"Sp1_g1": "FamA"}, Config()
    )

    # Assert: returned as a move candidate, never auto-applied as assignment
    assert moves == {"Sp1_g1": ("FamA", "FamB")}
    assert new_assignments == {}
    assert ambiguous == []


def test_reassignment_not_flagged_below_nbits_margin():
    # Arrange: FamB nbits 0.6 vs current 0.5 — margin 0.1 < 0.15
    hits = {"Sp1_g1": [
        _hit("Sp1_g1", "FamB", full_bits=60.0, hmm_len=100),
        _hit("Sp1_g1", "FamA", full_bits=50.0, hmm_len=100),
    ]}

    # Act
    _, _, moves = assign(hits, {"Sp1_g1": "FamA"}, Config())

    # Assert
    assert moves == {}


def test_reassignment_requires_candidate_profile_coverage():
    # Arrange: FamB has a big nbits margin but covers only half the profile
    hits = {"Sp1_g1": [
        _hit("Sp1_g1", "FamB", full_bits=90.0, hmm_len=100,
             domains=[(1, 50, 1, 100, 90.0)]),
        _hit("Sp1_g1", "FamA", full_bits=50.0, hmm_len=100),
    ]}

    # Act
    _, _, moves = assign(hits, {"Sp1_g1": "FamA"}, Config())

    # Assert
    assert moves == {}


# ---------------------------------------------------------------------------
# detect_merge_candidates
# ---------------------------------------------------------------------------

def test_detect_merge_candidates_reciprocal_pair():
    # Arrange: every FamA member's best non-self hit is FamB and vice versa
    families = {
        "FamA": {"Sp1_a1", "Sp2_a2"},
        "FamB": {"Sp1_b1", "Sp2_b2"},
    }
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=150.0)],
        "Sp2_a2": [_hit("Sp2_a2", "FamB", full_bits=140.0)],
        "Sp1_b1": [_hit("Sp1_b1", "FamA", full_bits=130.0)],
        "Sp2_b2": [_hit("Sp2_b2", "FamA", full_bits=120.0)],
    }

    # Act
    candidates = detect_merge_candidates(hits, families, Config())

    # Assert
    assert candidates == [("FamA", "FamB", 1.0)]


def test_detect_merge_candidates_rejects_one_way_attraction():
    # Arrange: FamA members all hit FamB, but no FamB member hits FamA
    families = {
        "FamA": {"Sp1_a1", "Sp2_a2"},
        "FamB": {"Sp1_b1", "Sp2_b2"},
    }
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=150.0)],
        "Sp2_a2": [_hit("Sp2_a2", "FamB", full_bits=140.0)],
    }

    # Act
    candidates = detect_merge_candidates(hits, families, Config())

    # Assert
    assert candidates == []


def test_detect_merge_candidates_below_reciprocal_fraction():
    # Arrange: only 1 of 3 FamA members votes for FamB (frac 0.33 < 0.6)
    families = {
        "FamA": {"Sp1_a1", "Sp2_a2", "Sp3_a3"},
        "FamB": {"Sp1_b1"},
    }
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=150.0)],
        "Sp1_b1": [_hit("Sp1_b1", "FamA", full_bits=130.0)],
    }

    # Act
    candidates = detect_merge_candidates(hits, families, Config())

    # Assert
    assert candidates == []


def test_detect_merge_candidates_ignores_low_coverage_and_weak_hits():
    # Arrange: cross-hits fail profile_cov>=0.5 or the E-value cutoff
    families = {"FamA": {"Sp1_a1"}, "FamB": {"Sp1_b1"}}
    hits = {
        "Sp1_a1": [_hit("Sp1_a1", "FamB", full_bits=150.0,
                        domains=[(1, 40, 1, 100, 150.0)])],  # cov 0.4
        "Sp1_b1": [_hit("Sp1_b1", "FamA", full_bits=130.0, full_e=1e-2)],
    }

    # Act
    candidates = detect_merge_candidates(hits, families, Config())

    # Assert
    assert candidates == []


# ---------------------------------------------------------------------------
# build_profiles: caching semantics
# ---------------------------------------------------------------------------

def _fake_hmm_tools(monkeypatch, calls):
    def fake_hmmbuild(family_id, alignment_path, hmm_path, config):
        calls.append(family_id)
        Path(hmm_path).write_text(f"HMM {family_id}\n")
        return True

    monkeypatch.setattr(pa, "_hmmbuild", fake_hmmbuild)
    monkeypatch.setattr(pa, "_hmmpress", lambda db, cfg: True)


def test_build_profiles_reuses_valid_cache_and_rebuilds_on_member_change(
    tmp_path, monkeypatch
):
    # Arrange
    calls = []
    _fake_hmm_tools(monkeypatch, calls)
    aln = tmp_path / "fam.afa"
    aln.write_text(">Sp1_g1\nMSEQ\n")
    hmm_dir = tmp_path / "hmm_profiles"
    families = {"FamA": {"Sp1_g1", "Sp2_g2"}}
    cfg = Config(n_workers=1)

    # Act 1: cold build
    db1 = build_profiles(families, lambda fid: aln, hmm_dir, cfg)
    # Act 2: identical membership — cached
    db2 = build_profiles(families, lambda fid: aln, hmm_dir, cfg)
    # Act 3: membership changed — sha1 sidecar mismatch forces rebuild
    families_changed = {"FamA": {"Sp1_g1", "Sp2_g2", "Sp3_g3"}}
    db3 = build_profiles(families_changed, lambda fid: aln, hmm_dir, cfg)

    # Assert
    assert calls == ["FamA", "FamA"]  # built cold, skipped, rebuilt
    assert db1 == db2 == db3
    assert db1.exists()
    assert (hmm_dir / "FamA.members").exists()


def test_build_profiles_rebuilds_zero_byte_cached_hmm(tmp_path, monkeypatch):
    # Arrange: zero-byte .hmm with a matching sidecar (the hmmer_rescue
    # cache bug treated any existing file as built)
    calls = []
    _fake_hmm_tools(monkeypatch, calls)
    aln = tmp_path / "fam.afa"
    aln.write_text(">Sp1_g1\nMSEQ\n")
    hmm_dir = tmp_path / "hmm_profiles"
    hmm_dir.mkdir()
    families = {"FamA": {"Sp1_g1"}}
    (hmm_dir / "FamA.hmm").write_text("")  # zero bytes
    (hmm_dir / "FamA.members").write_text(pa._member_digest({"Sp1_g1"}) + "\n")

    # Act
    db = build_profiles(families, lambda fid: aln, hmm_dir, Config(n_workers=1))

    # Assert
    assert calls == ["FamA"]
    assert db is not None


def test_build_profiles_returns_none_without_alignments(tmp_path, monkeypatch):
    # Arrange
    calls = []
    _fake_hmm_tools(monkeypatch, calls)

    # Act
    db = build_profiles(
        {"FamA": {"Sp1_g1"}}, lambda fid: None,
        tmp_path / "hmm_profiles", Config(n_workers=1),
    )

    # Assert
    assert db is None
    assert calls == []


# ---------------------------------------------------------------------------
# run_profile_assignment: orchestration with mocked subprocess wrappers
# ---------------------------------------------------------------------------

def test_run_profile_assignment_end_to_end(tmp_path, monkeypatch):
    # Arrange: two families; one unplaced gene assigns cleanly, one is a
    # bit-score tie (ambiguous), one has no hits
    calls = []
    _fake_hmm_tools(monkeypatch, calls)

    aln = tmp_path / "fam.afa"
    aln.write_text(">Sp1_g1\nMSEQ\n")
    families = {
        "R1_OG0000001": {"Sp1_g1", "Sp2_g2"},
        "R10_OG0000005": {"Sp1_g5", "Sp2_g6"},
    }
    unplaced = {"Ococ_u1": "MSEQA", "Ococ_u2": "MSEQB", "Ococ_u3": "MSEQC"}

    rows = [
        # clean winner for u1
        _dom_row("Ococ_u1", "R1_OG0000001", tlen=100, qlen=100,
                 full_e="1e-60", full_bits=300.0, dom_bits=300.0,
                 hf=1, ht=100, af=1, at=100),
        _dom_row("Ococ_u1", "R10_OG0000005", tlen=100, qlen=100,
                 full_e="1e-20", full_bits=90.0, dom_bits=90.0,
                 hf=1, ht=100, af=1, at=100),
        # exact tie for u2 -> ambiguous
        _dom_row("Ococ_u2", "R10_OG0000005", tlen=100, qlen=100,
                 full_e="0", full_bits=250.0, dom_bits=250.0,
                 hf=1, ht=100, af=1, at=100),
        _dom_row("Ococ_u2", "R1_OG0000001", tlen=100, qlen=100,
                 full_e="0", full_bits=250.0, dom_bits=250.0,
                 hf=1, ht=100, af=1, at=100),
    ]

    def fake_hmmsearch(hmm_db, query_fasta, outpath, config):
        Path(outpath).write_text("#header\n" + "".join(rows))
        return Path(outpath)

    monkeypatch.setattr(pa, "_hmmsearch_domtblout", fake_hmmsearch)

    outdir = tmp_path / "profile_assign"

    # Act
    assigned, still_unplaced = run_profile_assignment(
        families, unplaced, lambda fid: aln, outdir, Config(n_workers=1)
    )

    # Assert
    assert assigned == {"R1_OG0000001": {"Ococ_u1"}}
    assert still_unplaced == {"Ococ_u2", "Ococ_u3"}
    assert (outdir / "profile_assignments.tsv").exists()
    assert (outdir / "ambiguous_assignments.tsv").exists()

    tsv = (outdir / "profile_assignments.tsv").read_text().splitlines()
    assert tsv[0].startswith("gene_id\tfamily_id")
    assert any(line.startswith("Ococ_u1\tR1_OG0000001") for line in tsv[1:])

    amb = (outdir / "ambiguous_assignments.tsv").read_text()
    assert "Ococ_u2" in amb


def test_run_profile_assignment_empty_pool_is_noop(tmp_path):
    # Act
    assigned, still_unplaced = run_profile_assignment(
        {"FamA": {"Sp1_g1"}}, {}, lambda fid: None,
        tmp_path / "pa", Config(),
    )

    # Assert
    assert assigned == {}
    assert still_unplaced == set()


# ---------------------------------------------------------------------------
# Config defaults
# ---------------------------------------------------------------------------

def test_config_profile_assignment_defaults():
    cfg = Config()
    assert cfg.profile_assign_per_round is False  # opt-in
    assert cfg.profile_min_coverage == 0.5
    assert cfg.profile_min_query_coverage == 0.4
    assert cfg.profile_margin_bits == 10.0
    assert cfg.profile_margin_frac == 0.05
    assert cfg.profile_reassign_margin_nbits == 0.15
    assert cfg.merge_min_reciprocal == 0.6
