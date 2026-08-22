"""Tests for build_supermatrix.py (issue #41).

Pure helpers only — no MAFFT/pal2nal/IQ-TREE is ever invoked. Marker
selection and codon-column filtering take plain text/dicts precisely so
they test without the cluster or any external tool.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import build_supermatrix as bs


SPECIES5 = ["Mcry", "Cgig", "CgigH", "Obas", "Ococ"]


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


# --------------------------------------------------------------------------
# strict single-copy selection from OrthoFinder GeneCount
# --------------------------------------------------------------------------

GENECOUNT = "\n".join([
    "Orthogroup\tAhyp\tMcry\tObas\tTotal",
    "OG0000001\t1\t1\t1\t3",     # strict single copy
    "OG0000002\t1\t2\t1\t4",     # Mcry duplicated
    "OG0000003\t1\t0\t1\t2",     # Mcry missing
    "OG0000004\t1\t1\t1\t3",     # strict single copy
    "OG0000005\t3\t3\t3\t9",     # multi-copy
]) + "\n"


def test_selects_only_orthogroups_with_exactly_one_gene_per_species(tmp_path):
    # Arrange
    gc = write(tmp_path, "gc.tsv", GENECOUNT)

    # Act
    ogs = bs.strict_single_copy_orthogroups(gc, ["Ahyp", "Mcry", "Obas"])

    # Assert
    assert ogs == ["OG0000001", "OG0000004"]


def test_ignores_the_total_column_when_matching_species(tmp_path):
    # "Total" must never be treated as a species, or nothing is ever single-copy
    gc = write(tmp_path, "gc.tsv", GENECOUNT)

    ogs = bs.strict_single_copy_orthogroups(gc, ["Ahyp", "Mcry", "Obas"])

    assert "OG0000001" in ogs  # Total=3, would fail a naive ==1 check


def test_raises_when_a_requested_species_is_absent_from_the_header(tmp_path):
    gc = write(tmp_path, "gc.tsv", GENECOUNT)

    with pytest.raises(ValueError, match="Nope"):
        bs.strict_single_copy_orthogroups(gc, ["Ahyp", "Nope"])


# --------------------------------------------------------------------------
# single-copy family selection from summary.tsv (legacy 5sp path)
# --------------------------------------------------------------------------

SUMMARY = "\n".join([
    "family\tround\tn_genes\tn_species\tgenes",
    "R1_OG0000001\t1\t5\t5\tMcry_a,Cgig_b,CgigH_c,Obas_d,Ococ_e",
    "R1_OG0000002\t1\t6\t5\tMcry_a,Mcry_a2,Cgig_b,CgigH_c,Obas_d,Ococ_e",
    "R2_OG0000003\t2\t5\t4\tMcry_a,Cgig_b,Cgig_b2,Obas_d,Ococ_e",
    "R2_OG0000004\t2\t5\t5\tMcry_f,Cgig_g,CgigH_h,Obas_i,Ococ_j",
]) + "\n"


def test_summary_markers_require_one_gene_per_species(tmp_path):
    s = write(tmp_path, "summary.tsv", SUMMARY)

    markers = bs.single_copy_families(s, SPECIES5)

    assert [m.family for m in markers] == ["R1_OG0000001", "R2_OG0000004"]
    assert markers[0].round == 1
    assert markers[1].round == 2
    assert markers[0].members == ["Mcry_a", "Cgig_b", "CgigH_c", "Obas_d", "Ococ_e"]


def test_family_id_maps_to_its_round_orthogroup_directory():
    assert bs.orthogroup_dir("out", "R2_OG0000004") == "out/round_02/orthogroups/OG0000004"
    assert bs.orthogroup_dir("out", "R10_OG0000004") == "out/round_10/orthogroups/OG0000004"


def test_orthogroup_dir_rejects_an_unparseable_family_id():
    with pytest.raises(ValueError):
        bs.orthogroup_dir("out", "OG0000004")


# --------------------------------------------------------------------------
# codon-column filtering
# --------------------------------------------------------------------------

def test_drops_codon_columns_that_are_gaps_in_every_species():
    # codon 2 is all gaps in both species and must disappear
    by_sp = {"Mcry": "AAA---CCC", "Obas": "AAT---CCG"}

    kept, n = bs.drop_all_gap_codons(by_sp)

    assert n == 2
    assert kept == {"Mcry": "AAACCC", "Obas": "AATCCG"}


def test_keeps_a_codon_column_present_in_any_single_species():
    by_sp = {"Mcry": "AAAGGG", "Obas": "AAT---"}

    kept, n = bs.drop_all_gap_codons(by_sp)

    assert n == 2
    assert kept == {"Mcry": "AAAGGG", "Obas": "AAT---"}


def test_trailing_partial_codon_is_discarded():
    by_sp = {"Mcry": "AAACC", "Obas": "AATCC"}  # 5 nt: one whole codon + 2

    kept, n = bs.drop_all_gap_codons(by_sp)

    assert n == 1
    assert kept == {"Mcry": "AAA", "Obas": "AAT"}


def test_rejects_ragged_alignments_rather_than_silently_truncating():
    by_sp = {"Mcry": "AAACCC", "Obas": "AAT"}

    with pytest.raises(ValueError, match="unequal"):
        bs.drop_all_gap_codons(by_sp)


def test_returns_nothing_when_every_column_is_gapped():
    by_sp = {"Mcry": "------", "Obas": "------"}

    kept, n = bs.drop_all_gap_codons(by_sp)

    assert n == 0
    assert kept == {"Mcry": "", "Obas": ""}


def test_sequences_are_upcased_so_pal2nal_soft_masking_does_not_split_states():
    by_sp = {"Mcry": "aaaCCC", "Obas": "AATccg"}

    kept, _ = bs.drop_all_gap_codons(by_sp)

    assert kept == {"Mcry": "AAACCC", "Obas": "AATCCG"}


# --------------------------------------------------------------------------
# concatenation + partitions
# --------------------------------------------------------------------------

def test_partition_offsets_are_contiguous_and_1_based_inclusive():
    concat, parts = bs.concatenate(
        [("locus_a", {"Mcry": "AAACCC", "Obas": "AATCCG"}),
         ("locus_b", {"Mcry": "GGG", "Obas": "GGT"})],
        ["Mcry", "Obas"],
    )

    assert concat == {"Mcry": "AAACCCGGG", "Obas": "AATCCGGGT"}
    assert parts == [("locus_a", 1, 6), ("locus_b", 7, 9)]


def test_concatenation_refuses_a_locus_missing_a_species():
    with pytest.raises(ValueError, match="Obas"):
        bs.concatenate([("locus_a", {"Mcry": "AAA"})], ["Mcry", "Obas"])


def test_partition_file_uses_iqtree_DNA_syntax():
    text = bs.format_partitions([("locus_a", 1, 6), ("locus_b", 7, 9)])

    assert text == "DNA, locus_a = 1-6\nDNA, locus_b = 7-9\n"


# --------------------------------------------------------------------------
# fasta round-trip
# --------------------------------------------------------------------------

def test_read_fasta_keeps_only_the_id_before_the_first_space(tmp_path):
    p = write(tmp_path, "x.fa", ">Mcry_a some description\nAAA\nCCC\n>Obas_b\nGGG\n")

    seqs = bs.read_fasta(p)

    assert seqs == {"Mcry_a": "AAACCC", "Obas_b": "GGG"}
