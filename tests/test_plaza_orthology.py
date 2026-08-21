"""Tests for PLAZA-orthology naming (issue #30).

Two layers under test:
  steps/family_naming.py — map_orthology / combine_sources / weighted
    name_groups (provenance direct vs orthology, direct outweighs).
  extract_plaza_orthology.py — pure parsers over DIAMOND outfmt-6 text:
    best hits, transcript->gene collapse, RBH, annotation join.

House convention: plain dicts in/out, no external tools, no network.
"""

import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from steps.family_naming import (  # noqa: E402
    combine_sources,
    map_annotations,
    map_orthology,
    name_groups,
)
from extract_plaza_orthology import (  # noqa: E402
    best_hits,
    best_hits_per_species,
    build_orthology_rows,
    reciprocal_best,
    strip_transcript,
)

FAMILIES = {
    "F1": {"Mcry_g1", "Mcry_g2", "Ococ_x1"},
    "F2": {"Mcry_g3.t1", "Ococ_x2"},
}
GROUPS = {"Mcry_g1": "SF_A", "Mcry_g2": "SF_A", "Ococ_x2": "SF_B"}


# ---------------------------------------------------------------- naming

def test_map_orthology_rows_carry_provenance_and_ortholog():
    orthology = {"Mcry_g1": ("AT1G53310", "PPC1"), "zzz": ("AT0", "X")}
    rows, unmatched = map_orthology(orthology, FAMILIES, groups=GROUPS)
    assert unmatched == ["zzz"]
    (row,) = rows
    assert row["family_id"] == "F1"
    assert row["subfamily"] == "SF_A"
    assert row["source"] == "orthology"
    assert row["ortholog"] == "AT1G53310"
    assert row["description"] == "PPC1"


def test_map_orthology_matches_id_variants():
    # pipeline stored the isoform, orthology table has the locus
    orthology = {"g3": ("AT2G42600", "PPC2")}
    rows, unmatched = map_orthology(orthology, FAMILIES, species="Mcry")
    assert unmatched == []
    assert rows[0]["pipeline_gene"] == "Mcry_g3.t1"


def test_direct_rows_get_direct_source():
    rows, _ = map_annotations({"Mcry_g1": "PEPC"}, FAMILIES)
    assert rows[0]["source"] == "direct"
    assert rows[0]["ortholog"] == ""


def test_combine_sources_direct_wins_per_gene():
    direct, _ = map_annotations({"Mcry_g1": "PEPC kinase"}, FAMILIES)
    ortho, _ = map_orthology(
        {"Mcry_g1": ("AT1G53310", "PPC1"), "Ococ_x1": ("AT1G53310", "PPC1")},
        FAMILIES,
    )
    combined = combine_sources(direct, ortho)
    by_gene = {r["pipeline_gene"]: r for r in combined}
    # Mcry_g1 keeps only its direct annotation; Ococ_x1 orthology survives
    assert by_gene["Mcry_g1"]["source"] == "direct"
    assert by_gene["Ococ_x1"]["source"] == "orthology"
    assert len(combined) == 2


def test_name_groups_weighted_direct_outweighs_orthology():
    # 1 direct "kinase" (w=1.0) vs 3 orthology "PPC1" (w=0.5 each = 1.5):
    # orthology wins on weight; flip the weight to 0.25 and direct wins.
    rows = [
        {"family_id": "F1", "description": "kinase", "source": "direct"},
        {"family_id": "F1", "description": "PPC1", "source": "orthology"},
        {"family_id": "F1", "description": "PPC1", "source": "orthology"},
        {"family_id": "F1", "description": "PPC1", "source": "orthology"},
    ]
    named = name_groups(rows, key="family_id")[0]
    assert named["name"] == "PPC1"
    assert named["provenance"] == "orthology"
    assert named["support"] == pytest.approx(1.5 / 2.5)
    assert named["n_direct"] == 1 and named["n_orthology"] == 3

    named = name_groups(
        rows, key="family_id", weights={"direct": 1.0, "orthology": 0.25}
    )[0]
    assert named["name"] == "kinase"
    assert named["provenance"] == "direct"
    assert named["support"] == pytest.approx(1.0 / 1.75)


def test_name_groups_tie_prefers_direct_then_lexicographic():
    rows = [
        {"family_id": "F1", "description": "beta", "source": "orthology"},
        {"family_id": "F1", "description": "beta", "source": "orthology"},
        {"family_id": "F1", "description": "alpha", "source": "direct"},
    ]
    # beta 2*0.5 = 1.0 ties alpha 1.0 -> direct wins the tie
    assert name_groups(rows)[0]["name"] == "alpha"
    rows_all_ortho = [
        {"family_id": "F1", "description": "beta", "source": "orthology"},
        {"family_id": "F1", "description": "alpha", "source": "orthology"},
    ]
    assert name_groups(rows_all_ortho)[0]["name"] == "alpha"


def test_name_groups_legacy_rows_without_source_are_direct():
    rows = [
        {"family_id": "F1", "description": "PEPC 1"},
        {"family_id": "F1", "description": "PEPC 1"},
        {"family_id": "F1", "description": "kinase"},
    ]
    named = name_groups(rows)[0]
    assert named["name"] == "PEPC 1"
    assert named["support"] == pytest.approx(2 / 3)
    assert named["provenance"] == "direct"


# ------------------------------------------------------- extraction script

DMND = """\
Mcry_g1\tAT1G53310.2\t85.0\t900\t100\t2\t1\t900\t1\t900\t0.0\t1500
Mcry_g1\tAT2G42600.1\t80.0\t900\t120\t2\t1\t900\t1\t900\t0.0\t1400
Ococ_x1\tAT2G42600.1\t82.0\t900\t100\t2\t1\t900\t1\t900\t1e-200\t1300
"""

# Reverse search runs against a MULTI-species DB: RBH must be judged
# within the query gene's own species (one ATH gene can be the ortholog
# of one gene PER species, not one gene overall).
REV = """\
AT1G53310.2\tMcry_g1\t85.0\t900\t100\t2\t1\t900\t1\t900\t0.0\t1500
AT1G53310.2\tOcoc_x7\t86.0\t900\t100\t2\t1\t900\t1\t900\t0.0\t1550
AT2G42600.1\tMcry_g99\t83.0\t900\t100\t2\t1\t900\t1\t900\t0.0\t1350
AT2G42600.1\tOcoc_x9\t83.0\t900\t100\t2\t1\t900\t1\t900\t0.0\t1340
"""


def test_best_hits_keeps_top_bitscore_per_query():
    hits = best_hits(DMND.splitlines())
    assert hits["Mcry_g1"]["subject"] == "AT1G53310.2"
    assert hits["Mcry_g1"]["bitscore"] == 1500.0
    assert hits["Ococ_x1"]["subject"] == "AT2G42600.1"


def test_strip_transcript():
    assert strip_transcript("AT1G53310.2") == "AT1G53310"
    assert strip_transcript("AT1G53310") == "AT1G53310"


def test_best_hits_per_species_splits_on_gene_prefix():
    rev = best_hits_per_species(REV.splitlines())
    assert rev["AT1G53310.2"]["Mcry"]["subject"] == "Mcry_g1"
    assert rev["AT1G53310.2"]["Ococ"]["subject"] == "Ococ_x7"


def test_reciprocal_best_is_species_aware():
    rev = best_hits_per_species(REV.splitlines())
    # Mcry_g1 is AT1G53310's best hit IN Mcry, even though the Ococ copy
    # scores higher overall.
    assert reciprocal_best("Mcry_g1", "AT1G53310.2", rev) is True
    assert reciprocal_best("Ococ_x1", "AT2G42600.1", rev) is False


def test_build_orthology_rows_joins_annotations_and_flags_rbh():
    fwd = best_hits(DMND.splitlines())
    rev = best_hits_per_species(REV.splitlines())
    ann = {"AT1G53310": "PPC1"}
    rows = build_orthology_rows(fwd, rev, ann)
    by_gene = {r["gene_id"]: r for r in rows}
    assert by_gene["Mcry_g1"]["ath_gene"] == "AT1G53310"
    assert by_gene["Mcry_g1"]["description"] == "PPC1"
    assert by_gene["Mcry_g1"]["rbh"] == "yes"
    # no annotation -> falls back to the ATH gene id, never blank
    assert by_gene["Ococ_x1"]["description"] == "AT2G42600"
    assert by_gene["Ococ_x1"]["rbh"] == "no"


def test_build_orthology_rows_rbh_only_filter():
    fwd = best_hits(DMND.splitlines())
    rev = best_hits_per_species(REV.splitlines())
    rows = build_orthology_rows(fwd, rev, {}, rbh_only=True)
    assert [r["gene_id"] for r in rows] == ["Mcry_g1"]


def test_build_orthology_rows_without_reverse_marks_unknown():
    fwd = best_hits(DMND.splitlines())
    rows = build_orthology_rows(fwd, None, {})
    assert all(r["rbh"] == "" for r in rows)


# ------------------------------------------------------------ CLI reader

def test_read_orthology_tsv(tmp_path):
    from name_families import read_orthology_tsv
    p = tmp_path / "ortho.tsv"
    p.write_text(
        "gene_id\tath_gene\tdescription\tpident\trbh\n"
        "Mcry_g1\tAT1G53310\tPPC1\t85.0\tyes\n"
        "# comment\n"
        "Ococ_x1\tAT2G42600\tPPC2\t82.0\tno\n"
    )
    ortho = read_orthology_tsv(p)
    assert ortho == {
        "Mcry_g1": ("AT1G53310", "PPC1"),
        "Ococ_x1": ("AT2G42600", "PPC2"),
    }
