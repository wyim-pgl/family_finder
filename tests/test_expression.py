"""steps/expression.py — expression share as an axis that declares its coverage.

Expression matrices exist for two of the fifteen species in the PEPC clan
(Mcry and Ococ). An axis that quietly reports only those two reads as though
the other thirteen had been measured and found silent, so every species
without a matrix must come back as `expression unavailable` and the coverage
fraction must be on the summary.

Shares are always WITHIN one species: summing TPM across species compares
numbers normalised against different transcriptomes.

Pure text fixtures, no external tools (repo convention).
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from steps.expression import (
    UNAVAILABLE,
    family_expression,
    gene_means,
    read_matrix,
    species_coverage,
    subfamily_shares,
)

# One high-expression copy and three quiet ones, the SF3 shape in miniature.
OCOC = """\
gene_id\tT00\tT04\tT08
OcoChr03G21370\t1000\t2000\t3000
OcoChr03G21380\t10\t20\t30
OcoChr03G21390\t5\t5\t5
OcoChr03G13300\t100\t100\t100
"""

MCRY = """\
Geneid\trep1\trep2
Mcr8G11630\t40\t60
Mcr7G08600\t10\t10
"""


def _matrix(tmp_path, name, text, species):
    p = tmp_path / name
    p.write_text(text)
    return read_matrix(p, species)


def test_read_matrix_keeps_sample_names_and_values(tmp_path):
    m = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    assert m.samples == ("T00", "T04", "T08")
    assert m.values["OcoChr03G21370"] == pytest.approx((1000.0, 2000.0, 3000.0))
    assert m.species == "Ococ"


def test_gene_means_join_through_the_transcript_suffix(tmp_path):
    m = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    means, report = gene_means(m, ["Ococ_OcoChr03G21370.t1"])
    assert means["Ococ_OcoChr03G21370.t1"] == pytest.approx(2000.0)
    assert report["n_unmatched"] == 0


def test_gene_means_refuse_when_too_many_ids_stay_unmatched(tmp_path):
    m = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    with pytest.raises(ValueError, match="unmatched"):
        gene_means(m, ["Ococ_OcoChr03G21370.t1", "Ococ_nope1", "Ococ_nope2"],
                   max_unmatched=0.1)


def test_a_gene_absent_from_the_matrix_is_reported_not_zeroed(tmp_path):
    m = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    means, report = gene_means(m, ["Ococ_OcoChr03G21370.t1", "Ococ_missing"],
                               max_unmatched=0.9)
    assert "Ococ_missing" not in means
    assert report["unmatched"] == ["Ococ_missing"]


def test_shares_are_computed_within_one_species(tmp_path):
    m = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    members = ["Ococ_OcoChr03G21370.t1", "Ococ_OcoChr03G21380.t1",
               "Ococ_OcoChr03G21390.t1", "Ococ_OcoChr03G13300.t1"]
    report = family_expression(members, {"Ococ": m})
    rows = {r["gene"]: r for r in report["rows"]}
    total = 2000.0 + 20.0 + 5.0 + 100.0
    assert rows["Ococ_OcoChr03G21370.t1"]["share"] == pytest.approx(2000.0 / total)
    assert rows["Ococ_OcoChr03G21370.t1"]["status"] == "measured"
    assert report["by_species"]["Ococ"]["family_total_tpm"] == pytest.approx(total)


def test_subfamily_shares_sum_the_member_means(tmp_path):
    means = {"a": 4776.8, "b": 100.0, "c": 77.7}
    groups = {"SF3": ["a"], "SF1": ["b", "c"]}
    out = {r["group"]: r for r in subfamily_shares(means, groups)}
    assert out["SF3"]["total_tpm"] == pytest.approx(4776.8)
    assert out["SF3"]["share"] == pytest.approx(4776.8 / 4954.5, rel=1e-6)
    assert out["SF1"]["n_members"] == 2


def test_a_species_with_no_matrix_is_unavailable_never_dropped(tmp_path):
    m = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    members = ["Ococ_OcoChr03G21370.t1", "Cgig_something", "Sund_other"]
    report = family_expression(members, {"Ococ": m})
    rows = {r["gene"]: r for r in report["rows"]}
    assert len(rows) == 3
    assert rows["Cgig_something"]["status"] == UNAVAILABLE
    assert rows["Cgig_something"]["mean_tpm"] is None
    assert rows["Sund_other"]["status"] == UNAVAILABLE


def test_coverage_states_how_many_species_were_measured(tmp_path):
    m = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    members = ["Ococ_OcoChr03G21370.t1", "Cgig_a", "Sund_b"]
    cov = species_coverage(members, {"Ococ": m})
    assert cov["n_species"] == 3
    assert cov["n_species_with_matrix"] == 1
    assert cov["species_without_matrix"] == ["Cgig", "Sund"]
    assert cov["coverage"] == "1 of 3"


def test_two_matrices_are_kept_apart(tmp_path):
    ococ = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    mcry = _matrix(tmp_path, "mcry.tsv", MCRY, "Mcry")
    members = ["Ococ_OcoChr03G21370.t1", "Mcry_Mcr8G11630", "Mcry_Mcr7G08600"]
    report = family_expression(members, {"Ococ": ococ, "Mcry": mcry})
    rows = {r["gene"]: r for r in report["rows"]}
    # Mcry shares are against the Mcry total (50+10), not a pooled one.
    assert rows["Mcry_Mcr8G11630"]["share"] == pytest.approx(50.0 / 60.0)
    assert rows["Ococ_OcoChr03G21370.t1"]["share"] == pytest.approx(1.0)
    assert report["coverage"]["n_species_with_matrix"] == 2


def test_the_report_never_offers_a_cross_species_total(tmp_path):
    ococ = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    mcry = _matrix(tmp_path, "mcry.tsv", MCRY, "Mcry")
    report = family_expression(
        ["Ococ_OcoChr03G21370.t1", "Mcry_Mcr8G11630"],
        {"Ococ": ococ, "Mcry": mcry})
    assert "family_total_tpm" not in report
    assert set(report["by_species"]) == {"Ococ", "Mcry"}


def test_subfamily_shares_are_reported_per_species(tmp_path):
    ococ = _matrix(tmp_path, "ococ.tsv", OCOC, "Ococ")
    members = ["Ococ_OcoChr03G21370.t1", "Ococ_OcoChr03G13300.t1",
               "Cgig_absent"]
    groups = {"SF3": ["Ococ_OcoChr03G21370.t1"],
              "SF1": ["Ococ_OcoChr03G13300.t1", "Cgig_absent"]}
    report = family_expression(members, {"Ococ": ococ}, groups=groups)
    sf = {r["group"]: r for r in report["by_species"]["Ococ"]["groups"]}
    assert sf["SF3"]["share"] == pytest.approx(2000.0 / 2100.0)
    # the Cgig member is in the group but has no matrix -- say so
    assert sf["SF1"]["n_unmeasured"] == 1
