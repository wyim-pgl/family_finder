"""region_disorder.py — is a subfamily's region disordered, and ONLY there?

Synthetic PDB/alignment fixtures; no structure prediction involved.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from region_disorder import (
    columns_to_residues,
    mann_whitney_p,
    measure,
    residue_plddt,
)


def _pdb(path, plddts):
    """One CA atom per residue, pLDDT in the B-factor column (cols 61-66)."""
    lines = []
    for i, b in enumerate(plddts, start=1):
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
            f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{b:6.2f}\n"
        )
    path.write_text("".join(lines))


def test_columns_to_residues_skips_gaps():
    # cols:      1234567890
    assert columns_to_residues("AB--CD", 1, 6) == [1, 2, 3, 4]
    assert columns_to_residues("AB--CD", 5, 6) == [3, 4]


def test_residue_plddt_reads_bfactor(tmp_path):
    p = tmp_path / "x.pdb"
    _pdb(p, [0.9, 0.5, 0.8])
    assert residue_plddt(p) == {1: 0.9, 2: 0.5, 3: 0.8}


def test_mann_whitney_separates_clearly():
    p = mann_whitney_p([1, 1, 1, 1, 1], [9, 9, 9, 9, 9])
    assert p < 0.05


def test_mann_whitney_identical_groups_not_significant():
    assert mann_whitney_p([1, 2, 3], [1, 2, 3]) > 0.05


def _setup(tmp_path, focal_low=True):
    """3 focal + 4 other; region = columns 3-7."""
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    seqs = {}
    for i in range(3):
        g = f"Focal_g{i}"
        seqs[g] = "A" * 10
        # residues 3-7 low, rest high
        _pdb(pdb_dir / f"{g}.pdb",
             [0.9, 0.9, 0.4, 0.4, 0.4, 0.4, 0.4, 0.9, 0.9, 0.9]
             if focal_low else [0.9] * 10)
    for i in range(4):
        g = f"Other_g{i}"
        seqs[g] = "A" * 10
        _pdb(pdb_dir / f"{g}.pdb", [0.9] * 10)
    aln = tmp_path / "aln.fa"
    aln.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    return pdb_dir, aln


def test_measure_detects_focal_specific_disorder(tmp_path):
    pdb_dir, aln = _setup(tmp_path)
    r = measure(pdb_dir, aln, 3, 7, {"Focal_g0", "Focal_g1", "Focal_g2"})
    assert r["n_focal"] == 3 and r["n_other"] == 4
    assert r["delta"] < -0.1
    assert abs(r["delta_other"]) < 1e-9
    assert r["all_below"] is True
    assert r["p"] < 0.05


def test_measure_reports_no_effect_when_region_is_ordinary(tmp_path):
    pdb_dir, aln = _setup(tmp_path, focal_low=False)
    r = measure(pdb_dir, aln, 3, 7, {"Focal_g0", "Focal_g1", "Focal_g2"})
    assert abs(r["delta"]) < 1e-9
    assert r["all_below"] is False
    assert r["p"] > 0.05


def test_measure_matches_dotted_ids_to_underscored_structures(tmp_path):
    """Structure files drop dots from gene ids — the member list uses the
    real ids, so a mismatch here would silently produce an empty focal set."""
    pdb_dir, aln = _setup(tmp_path)
    (pdb_dir / "Focal_g0.pdb").rename(pdb_dir / "Focal_g0.pdb")
    r = measure(pdb_dir, aln, 3, 7, {"Focal.g0", "Focal.g1", "Focal.g2"})
    assert r["n_focal"] == 3, "dotted member ids must match underscored files"


def test_measure_skips_sequences_with_too_little_region(tmp_path):
    pdb_dir, aln = _setup(tmp_path)
    # a gappy member contributes < min_residues in the window
    extra = "AA--------"
    txt = aln.read_text() + f">Gappy_g\n{extra}\n"
    aln.write_text(txt)
    _pdb(pdb_dir / "Gappy_g.pdb", [0.9, 0.9])
    r = measure(pdb_dir, aln, 3, 7, {"Focal_g0", "Focal_g1", "Focal_g2"})
    assert r["n_focal"] + r["n_other"] == 7   # Gappy_g excluded
