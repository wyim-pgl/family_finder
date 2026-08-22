"""measure_v2.py — v2 acceptance metric (per-species unplaced rate).

Pure fixtures; no pipeline run needed.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _make_run(tmp_path, families, unplaced):
    (tmp_path / "summary.tsv").write_text(
        "family_id\ta\tb\tc\tgene_list\n"
        + "".join(f"F{i}\t.\t.\t.\t{','.join(m)}\n" for i, m in enumerate(families))
    )
    (tmp_path / "unplaced_proteins.fa").write_text(
        "".join(f">{g}\nMA\n" for g in unplaced)
    )
    return tmp_path


def _run(run_dir):
    r = subprocess.run([sys.executable, str(ROOT / "measure_v2.py"), str(run_dir)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


def test_reports_pass_when_mcry_beats_baseline(tmp_path):
    run = _make_run(tmp_path, [["Mcry_a", "Cgig_b"]], ["Mcry_x"])
    out = _run(run)
    assert "-> PASS" in out            # 1/25226 is far below 14.9%


def test_reports_fail_when_mcry_worse(tmp_path):
    unplaced = [f"Mcry_g{i}" for i in range(4000)]   # 15.9% of 25226
    run = _make_run(tmp_path, [["Mcry_a"]], unplaced)
    out = _run(run)
    assert "-> FAIL" in out


def test_counts_species_by_prefix(tmp_path):
    run = _make_run(tmp_path, [["Ococ_a", "Ococ_b", "Obas_c"]], ["Obas_z"])
    out = _run(run)
    ococ = [l for l in out.splitlines() if l.startswith("Ococ")][0]
    obas = [l for l in out.splitlines() if l.startswith("Obas")][0]
    assert ococ.split()[2] == "2"      # placed
    assert obas.split()[3] == "1"      # unplaced


def test_flags_gene_count_mismatch(tmp_path):
    """Accounted genes far below the known total means the run is partial —
    the metric must say so rather than report a flattering rate."""
    run = _make_run(tmp_path, [["Mcry_a"]], ["Mcry_x"])
    out = _run(run)
    assert "[!]" in out


def test_missing_unplaced_file_is_tolerated(tmp_path):
    (tmp_path / "summary.tsv").write_text(
        "family_id\ta\tb\tc\tgene_list\nF0\t.\t.\t.\tMcry_a\n")
    out = _run(tmp_path)
    assert "HEADLINE" in out


# ---------------------------------------------------------------------------
# Length stratification (issue #36)
#
# The raw per-species unplaced rate compares annotation policy as much as
# pipeline behaviour: Cgig.pep.fa has a hard floor at 151 aa (zero proteins
# under 100), Mcry's shortest is 3 aa with 3,374 under 100. The metric has to
# make that visible instead of burying it in one number.
# ---------------------------------------------------------------------------

def _make_pep_dir(tmp_path, lengths_by_species):
    """lengths_by_species: {species: {gene: protein_length}}"""
    pep = tmp_path / "pep"
    pep.mkdir(exist_ok=True)
    for sp, genes in lengths_by_species.items():
        (pep / f"{sp}.pep.fa").write_text(
            "".join(f">{g}\n{'M' * n}\n" for g, n in genes.items())
        )
    return pep


def _run_stratified(run_dir, pep_dir, *extra):
    r = subprocess.run(
        [sys.executable, str(ROOT / "measure_v2.py"), str(run_dir),
         "--pep-dir", str(pep_dir), *extra],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


def _row(out, prefix):
    for line in out.splitlines():
        if line.split()[:len(prefix.split())] == prefix.split():
            return line.split()
    raise AssertionError(f"no row starting {prefix!r} in:\n{out}")


def test_reports_minimum_protein_length_per_species(tmp_path):
    # Arrange: Cgig-like hard floor vs Mcry-like no floor
    pep = _make_pep_dir(tmp_path, {
        "Mcry": {"Mcry_a": 3, "Mcry_b": 300},
        "Cgig": {"Cgig_a": 151, "Cgig_b": 400},
    })
    run = _make_run(tmp_path, [["Mcry_b", "Cgig_b"]], ["Mcry_a", "Cgig_a"])

    # Act
    out = _run_stratified(run, pep)

    # Assert: the floor that causes the confound is on screen
    assert _row(out, "Mcry")[1] == "3" or "3" in out
    assert "151" in out


def test_warns_when_species_minimum_lengths_differ_materially(tmp_path):
    pep = _make_pep_dir(tmp_path, {
        "Mcry": {"Mcry_a": 3, "Mcry_b": 300},
        "Cgig": {"Cgig_a": 151, "Cgig_b": 400},
    })
    run = _make_run(tmp_path, [["Mcry_b", "Cgig_b"]], ["Mcry_a", "Cgig_a"])

    out = _run_stratified(run, pep)

    assert "CONFOUND" in out


def test_no_confound_warning_when_annotations_share_a_floor(tmp_path):
    pep = _make_pep_dir(tmp_path, {
        "Mcry": {"Mcry_a": 120, "Mcry_b": 300},
        "Cgig": {"Cgig_a": 130, "Cgig_b": 400},
    })
    run = _make_run(tmp_path, [["Mcry_b", "Cgig_b"]], ["Mcry_a", "Cgig_a"])

    out = _run_stratified(run, pep)

    assert "CONFOUND" not in out


def test_stratified_rates_use_per_stratum_denominators(tmp_path):
    # Arrange: Mcry has 4 short (2 unplaced) and 2 long (1 unplaced)
    pep = _make_pep_dir(tmp_path, {"Mcry": {
        "Mcry_s1": 10, "Mcry_s2": 20, "Mcry_s3": 30, "Mcry_s4": 40,
        "Mcry_l1": 500, "Mcry_l2": 600,
    }})
    run = _make_run(tmp_path,
                    [["Mcry_s3", "Mcry_s4", "Mcry_l2"]],
                    ["Mcry_s1", "Mcry_s2", "Mcry_l1"])

    # Act
    out = _run_stratified(run, pep)

    # Assert: 2/4 = 50.0% short, 1/2 = 50.0% long — denominators are per stratum
    short = _row(out, "Mcry short")
    long_ = _row(out, "Mcry long")
    assert short[2] == "4" and short[3] == "2" and "50.0%" in short[4]
    assert long_[2] == "2" and long_[3] == "1" and "50.0%" in long_[4]


def test_stratification_floor_is_configurable(tmp_path):
    pep = _make_pep_dir(tmp_path, {"Mcry": {
        "Mcry_a": 150, "Mcry_b": 250, "Mcry_c": 350,
    }})
    run = _make_run(tmp_path, [["Mcry_b", "Mcry_c"]], ["Mcry_a"])

    # With a 200 aa floor, Mcry_a (150) falls in the short stratum
    out = _run_stratified(run, pep, "--floor", "200")

    assert "200" in out
    short = _row(out, "Mcry short")
    assert short[2] == "1" and short[3] == "1"


def test_species_with_no_short_genes_reports_the_empty_stratum(tmp_path):
    # Cgig's real state: zero proteins under the floor. The row must exist and
    # say so, not vanish — its absence is the whole finding.
    pep = _make_pep_dir(tmp_path, {"Cgig": {"Cgig_a": 151, "Cgig_b": 400}})
    run = _make_run(tmp_path, [["Cgig_b"]], ["Cgig_a"])

    out = _run_stratified(run, pep)

    short = _row(out, "Cgig short")
    assert short[2] == "0"


def test_raw_headline_is_preserved_for_continuity_with_the_v1_baseline(tmp_path):
    pep = _make_pep_dir(tmp_path, {"Mcry": {"Mcry_a": 3, "Mcry_b": 300}})
    run = _make_run(tmp_path, [["Mcry_b"]], ["Mcry_a"])

    out = _run_stratified(run, pep)

    assert "HEADLINE" in out and "14.9%" in out


def test_stratified_headline_is_reported_alongside_the_raw_one(tmp_path):
    pep = _make_pep_dir(tmp_path, {"Mcry": {"Mcry_a": 3, "Mcry_b": 300}})
    run = _make_run(tmp_path, [["Mcry_b"]], ["Mcry_a"])

    out = _run_stratified(run, pep)

    assert "COMPARABLE" in out


def test_without_a_pep_dir_the_tool_still_runs_and_says_why_it_cannot_stratify(tmp_path):
    run = _make_run(tmp_path, [["Mcry_a"]], ["Mcry_x"])

    out = _run(run)

    assert "HEADLINE" in out
    assert "--pep-dir" in out
