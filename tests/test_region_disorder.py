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


def test_measure_reports_genes_it_skipped_rather_than_dropping_them(tmp_path):
    # A structure whose name is absent from the alignment used to vanish
    # silently; shrinking coverage then looks like an absent signal (#42).
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    for gene in ("g1", "ghost"):
        (pdb_dir / f"{gene}.pdb").write_text(
            "\n".join(
                f"ATOM  {i:>5}  CA  ALA A{i:>4}       1.000   2.000   3.000  1.00 80.00           C"
                for i in range(1, 11)
            ) + "\nEND\n"
        )
    aln = tmp_path / "aln.fa"
    aln.write_text(">g1\n" + "M" * 10 + "\n>g2\n" + "M" * 10 + "\n")

    import region_disorder as rd
    result = rd.measure(pdb_dir, aln, 1, 10, {"g1"})

    assert "ghost" in result["skipped"]["not_in_alignment"]
    assert result["n_skipped"] == 1


def test_a_gene_outside_the_region_is_distinguishable_from_a_merely_short_one(tmp_path):
    # 'region_too_short' covers two different facts: a gene that reaches the
    # region and stops inside it, and a gene with no residue there at all.
    # On the PEPC clan all ten skips were the second kind — the block is
    # absent, not short — and the bare label cannot say so.
    pdb_dir, aln = _setup(tmp_path)
    aln.write_text(aln.read_text()
                   + ">Partial_g\nAAAA------\n"   # residues 3-4 inside 3-7
                   + ">Absent_g\nAA--------\n")   # nothing inside 3-7
    _pdb(pdb_dir / "Partial_g.pdb", [0.9] * 4)
    _pdb(pdb_dir / "Absent_g.pdb", [0.9] * 2)

    r = measure(pdb_dir, aln, 3, 7, {"Focal_g0", "Focal_g1", "Focal_g2"})
    detail = {d["gene"]: d for d in r["skipped_detail"]}

    assert detail["Partial_g"]["n_region_residues"] == 2
    assert detail["Absent_g"]["n_region_residues"] == 0
    assert detail["Absent_g"]["reason"] == "region_too_short"
    assert detail["Absent_g"]["n_residues"] == 2


def test_an_unmatched_structure_is_in_the_detail_with_nothing_invented(tmp_path):
    pdb_dir, aln = _setup(tmp_path)
    _pdb(pdb_dir / "ghost_g.pdb", [0.9] * 10)

    r = measure(pdb_dir, aln, 3, 7, {"Focal_g0", "Focal_g1", "Focal_g2"})
    detail = {d["gene"]: d for d in r["skipped_detail"]}

    assert detail["ghost_g"]["reason"] == "not_in_alignment"
    # It has no alignment row, so it has no residue counts. None, not zero:
    # zero would read as a gene that covers nothing.
    assert detail["ghost_g"]["n_region_residues"] is None
    assert detail["ghost_g"]["n_residues"] is None


def test_the_detail_accounts_for_every_skip(tmp_path):
    pdb_dir, aln = _setup(tmp_path)
    aln.write_text(aln.read_text() + ">Absent_g\nAA--------\n")
    _pdb(pdb_dir / "Absent_g.pdb", [0.9] * 2)
    _pdb(pdb_dir / "ghost_g.pdb", [0.9] * 10)

    r = measure(pdb_dir, aln, 3, 7, {"Focal_g0", "Focal_g1", "Focal_g2"})

    assert len(r["skipped_detail"]) == r["n_skipped"] == 2
    assert r["n_focal"] + r["n_other"] + r["n_skipped"] == r["n_structures"]


def test_alignment_members_with_no_structure_are_reported_too(tmp_path):
    # The other direction of the same mismatch. On the PEPC clan seven of the
    # 102 aligned sequences -- the ATH and Aco anchors -- were never folded,
    # and nothing in the output said so: 85 measured out of 102 looked like
    # loss, when 95 was the largest number the two inputs could ever share.
    pdb_dir, aln = _setup(tmp_path)
    aln.write_text(aln.read_text() + ">Unfolded_g\n" + "A" * 10 + "\n")

    r = measure(pdb_dir, aln, 3, 7, {"Focal_g0", "Focal_g1", "Focal_g2"})

    assert r["n_alignment_without_structure"] == 1
    assert r["alignment_without_structure"] == ["Unfolded_g"]


def test_measure_matches_transcript_suffixed_structures_through_the_matcher(tmp_path):
    # The PEPC pilot FASTA and #40's groups.json disagreed by a `.t1` suffix
    # (#34); region_disorder's own '.'->'_' rule cannot recover that, the
    # shared matcher can.
    pdb_dir, aln = _setup(tmp_path)
    (pdb_dir / "Focal_g0.pdb").rename(pdb_dir / "Focal_g0.t1.pdb")

    r = measure(pdb_dir, aln, 3, 7, {"Focal_g0", "Focal_g1", "Focal_g2"})

    assert r["n_focal"] == 3
    assert r["n_skipped"] == 0


def test_measure_refuses_when_told_a_ceiling_and_too_much_is_unmatched(tmp_path):
    import pytest as _pytest

    pdb_dir, aln = _setup(tmp_path)
    for name in ("ghost1", "ghost2", "ghost3", "ghost4"):
        _pdb(pdb_dir / f"{name}.pdb", [0.9] * 10)

    with _pytest.raises(ValueError, match="unmatched"):
        measure(pdb_dir, aln, 3, 7, {"Focal_g0"}, max_unmatched=0.1)


def test_measure_records_which_id_level_the_match_needed(tmp_path):
    pdb_dir, aln = _setup(tmp_path)

    r = measure(pdb_dir, aln, 3, 7, {"Focal_g0", "Focal_g1", "Focal_g2"})

    assert r["id_match_level"] == "exact"
