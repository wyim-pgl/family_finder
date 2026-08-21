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
