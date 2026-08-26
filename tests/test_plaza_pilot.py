"""plaza_pilot — cohort determinism, the SAFE_ORTHOLOGY gate contract,
UniProt one-to-many policy, and the pre-committed decision rules.

Pure dict fixtures; no DIAMOND/foldseek run (house convention).
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from plaza_pilot import (  # noqa: E402
    Gates,
    afdb_accession,
    afdb_named_calls,
    conflict_rate,
    decision_rules,
    load_idmapping,
    parse_diamond_tsv,
    pick_representatives,
    resolve_uniprot,
    safe_orthology,
    sample_comparable_unplaced,
    wilson_upper,
)


# ---------------------------------------------------------------- cohort --

def test_representative_is_longest_member_with_lexicographic_tiebreak():
    families = {"F1": ["Mcry_b", "Mcry_a", "Mcry_c"]}
    lengths = {"Mcry_a": 300, "Mcry_b": 300, "Mcry_c": 100}

    assert pick_representatives(families, lengths) == {"F1": "Mcry_a"}


def test_representative_missing_from_panel_refuses():
    with pytest.raises(ValueError, match="missing from the panel"):
        pick_representatives({"F1": ["Mcry_a", "Mcry_gone"]}, {"Mcry_a": 10})


def test_unplaced_sample_is_deterministic_and_capped():
    eligible = {"Mcry": [f"Mcry_g{i}" for i in range(250)],
                "Cgig": ["Cgig_g1", "Cgig_g2"]}

    one = sample_comparable_unplaced(eligible, per_species=100)
    two = sample_comparable_unplaced(eligible, per_species=100)

    assert one == two                          # salted hash, not RNG
    assert sum(1 for g in one if g.startswith("Mcry_")) == 100
    assert sum(1 for g in one if g.startswith("Cgig_")) == 2


# ------------------------------------------------------------- the gate --

def _hit(sseqid, bits, pident=80.0, qcov=90.0, scov=90.0, evalue=1e-50):
    return {"qseqid": "Mcry_g1", "sseqid": sseqid, "pident": pident,
            "length": "100", "qlen": "100", "slen": "100", "qstart": "1",
            "qend": "100", "sstart": "1", "send": "100", "evalue": evalue,
            "bitscore": bits, "qcovhsp": qcov, "scovhsp": scov}


REVERSE = {("Mcry", "AT1G53310"): "Mcry_g1"}


def test_gate_accepts_a_clean_reciprocal_hit():
    v = safe_orthology("Mcry_g1", [_hit("AT1G53310.1", 500.0)],
                       REVERSE, Gates())
    assert v["status"] == "SAFE" and v["agi"] == "AT1G53310"


def test_isoforms_of_one_locus_are_not_margin_competition():
    hits = [_hit("AT1G53310.1", 500.0), _hit("AT1G53310.2", 499.0)]
    v = safe_orthology("Mcry_g1", hits, REVERSE, Gates())
    assert v["status"] == "SAFE"


def test_close_second_locus_fails_the_margin_gate():
    hits = [_hit("AT1G53310.1", 500.0), _hit("AT3G14940.1", 480.0)]
    v = safe_orthology("Mcry_g1", hits, REVERSE, Gates())
    assert v["status"] == "MARGIN"


def test_low_coverage_abstains_with_its_reason():
    v = safe_orthology("Mcry_g1", [_hit("AT1G53310.1", 500.0, qcov=40.0)],
                       REVERSE, Gates())
    assert v["status"] == "QCOV"


def test_non_reciprocal_best_hit_abstains():
    reverse = {("Mcry", "AT1G53310"): "Mcry_OTHER"}
    v = safe_orthology("Mcry_g1", [_hit("AT1G53310.1", 500.0)],
                       reverse, Gates())
    assert v["status"] == "NOT_RBH"
    assert v["reverse_partner"] == "Mcry_OTHER"


def test_no_reverse_hit_is_its_own_reason_not_a_pass():
    v = safe_orthology("Mcry_g1", [_hit("AT1G53310.1", 500.0)], {}, Gates())
    assert v["status"] == "NO_REVERSE_HIT"


def test_parse_diamond_refuses_truncated_rows():
    with pytest.raises(ValueError, match="columns"):
        parse_diamond_tsv(["a\tb\t50.0\n"])


# ------------------------------------------------- conflicts and wilson --

def test_conflict_rate_counts_exact_description_disagreement():
    stats = conflict_rate({"g1": "PEPC", "g2": "MYB", "g3": "PEPC"},
                          {"g1": "PEPC", "g2": "bHLH", "g4": "x"})
    assert (stats["n_overlap"], stats["n_conflict"]) == (2, 1)
    assert stats["conflicting_genes"] == ["g2"]


def test_wilson_upper_shrinks_with_evidence():
    assert wilson_upper(0, 10) > wilson_upper(0, 1000)
    assert wilson_upper(0, 0) == 1.0


# --------------------------------------------------------- uniprot side --

IDMAP = ("P10490\tTAIR\tAT1G53310\n"
         "P10490\tGene_Name\tPPC1\n"
         "Q99999\tTAIR\tAT1G00001\n"
         "Q99999\tTAIR\tAT2G00002\n"
         "A0A000\tGene_Name\tX\n")


def test_uniprot_one_to_many_policy():
    mapping = load_idmapping(IDMAP.splitlines())
    assert resolve_uniprot("P10490", mapping) == ("AT1G53310", "OK")
    assert resolve_uniprot("Q99999", mapping)[1] == "AMBIGUOUS_AGI"
    assert resolve_uniprot("A0A000", mapping)[1] == "NO_AGI_MAP"
    assert resolve_uniprot("MISSING", mapping)[1] == "NO_AGI_MAP"


def test_afdb_accession_extraction():
    assert afdb_accession("AF-P10490-F1-model_v6") == "P10490"
    assert afdb_accession("pdb_1abc") is None


def test_afdb_calls_gate_coverage_and_name():
    mapping = load_idmapping(IDMAP.splitlines())
    m8 = ("Mcry_g1\tAF-P10490-F1-model_v6\t0.9\t500\t0.95\t0.95\t1e-80\t900\n"
          "Mcry_g2\tAF-P10490-F1-model_v6\t0.9\t100\t0.30\t0.95\t1e-20\t200\n")
    out = afdb_named_calls(m8.splitlines(), mapping)
    assert out["calls"]["Mcry_g1"]["name"] == "PPC1"
    assert "Mcry_g2" not in out["calls"]
    assert out["abstained"] == {"LOW_COVERAGE": 1}


# --------------------------------------------------------- decision rules --

GOOD = {"safe_pct": 62.0, "named_rep_pct": 41.0, "conflict_rate": 0.01,
        "conflict_wilson": 0.03, "n_overlap": 220,
        "worst_species_conflict": 0.05, "swissprot_named_pct": 22.0,
        "arath_addressable_pct": 18.0}


def test_decision_rules_full_pass_configuration():
    d = decision_rules(GOOD)
    assert d["a_full_panel_diamond"]["go"] is True
    assert d["b_arath_download"]["go"] is True
    assert d["d_abandon_structural"]["go"] is False


def test_conflict_over_ten_percent_abandons_structural():
    d = decision_rules({**GOOD, "conflict_rate": 0.11})
    assert d["a_full_panel_diamond"]["go"] is False
    assert d["d_abandon_structural"]["go"] is True


def test_missing_metrics_are_undecided_not_passing():
    d = decision_rules({"safe_pct": 90.0})
    assert d["a_full_panel_diamond"]["go"] is None
    assert d["b_arath_download"]["go"] is None


def test_reverse_best_stores_the_panel_gene_not_the_ath_query():
    # Storing qseqid (the ATH id itself) made every RBH fail against itself
    # and produced SAFE=0 over 25,131 queries on the first real run.
    from plaza_pilot import reverse_best_from_hits

    hits = {"AT1G53310.1": [
        {"qseqid": "AT1G53310.1", "sseqid": "Mcry_g1",
         "bitscore": 500.0, "evalue": 1e-80},
        {"qseqid": "AT1G53310.1", "sseqid": "Mcry_g9",
         "bitscore": 100.0, "evalue": 1e-10},
    ], "AT1G53310.2": [
        {"qseqid": "AT1G53310.2", "sseqid": "Mcry_g1",
         "bitscore": 510.0, "evalue": 1e-82},
    ]}

    out = reverse_best_from_hits("Mcry", hits)

    assert out == {("Mcry", "AT1G53310"): "Mcry_g1"}


# ------------------------------------------------------- CLI file paths --

def test_direct_table_newlines_do_not_manufacture_conflicts(tmp_path):
    """qc-orthology read --direct values without rstrip: every curated
    description carried its newline, so conflict_rate reported 100%
    disagreement against the rstripped ATH descriptions (Mcry curated
    calibration, 2026-08-26: rate 1.0 on all 36 grid cells while the
    AGI-level agreement was 175/188)."""
    from plaza_pilot import main

    row = ("Mcry_g1\tAT1G01000.1\t80.0\t400\t400\t400\t1\t400\t1\t400"
           "\t1e-100\t500\t90.0\t90.0\n")
    rev = ("AT1G01000.1\tMcry_g1\t80.0\t400\t400\t400\t1\t400\t1\t400"
           "\t1e-100\t500\t90.0\t90.0\n")
    (tmp_path / "forward.tsv").write_text(row)
    (tmp_path / "Mcry.tsv").write_text(rev)
    (tmp_path / "cohort.tsv").write_text(
        "family_id\tgene_id\tspecies\n" "F1\tMcry_g1\tMcry\n")
    (tmp_path / "direct.tsv").write_text("Mcry_g1\tPPC1\n")
    (tmp_path / "desc.tsv").write_text("AT1G01000\tPPC1\n")

    out = tmp_path / "qc"
    main(["qc-orthology",
          "--cohort", str(tmp_path / "cohort.tsv"),
          "--forward", str(tmp_path / "forward.tsv"),
          "--reverse", str(tmp_path / "Mcry.tsv"),
          "--direct", str(tmp_path / "direct.tsv"),
          "--descriptions", str(tmp_path / "desc.tsv"),
          "--outdir", str(out)])

    lines = (out / "gate_calibration.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    i_over = header.index("n_overlap")
    i_conf = header.index("n_conflict")
    for line in lines[1:]:
        parts = line.split("\t")
        assert parts[i_over] == "1"
        assert parts[i_conf] == "0"
