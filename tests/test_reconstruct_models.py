"""Reconstructing full-length gene models from anchors (issue #40).

Six sequences entered the PEPC clan set with an `_M` suffix, replacing
truncated released annotations. No script produced them and none of the
project's documentation mentions the convention, so the models behind a
published claim could not be regenerated. Four of the six turned out to be
genuine — encoded in their own genome at 98-99.7% — and two were not:
Ococ at 85.5% and Pami at 79.0%.

The defect was never the naming. It was that nothing checked the output
against the genome it claimed to come from. These tests pin that gate.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from reconstruct_models import (
    Prediction,
    accept,
    best_per_locus,
    parse_miniprot,
)


GFF = """##gff-version 3
##PAF\tMcry_Mcr8G11630\t966\t0\t966\t+\tChr3\t91527810\t19252376\t19259269\t2514
##STA\tMATAKLEKLASIDAHLR
Chr3\tminiprot\tmRNA\t19252377\t19259272\t4458\t+\t.\tID=MP000001;Rank=1;Identity=0.8675;Target=Mcry_Mcr8G11630 1 966
Chr3\tminiprot\tCDS\t19252377\t19252544\t238\t+\t0\tParent=MP000001;Rank=1
##STA\tMARNLEKMASIDAQLRLL
Chr10\tminiprot\tmRNA\t36209236\t36227802\t5000\t+\t.\tID=MP000002;Rank=1;Identity=0.9950;Target=Mcry_Mcr7G08600 1 968
"""


def test_translated_protein_is_taken_from_the_sta_line_above_its_mrna():
    preds = parse_miniprot(GFF)

    assert [p.protein for p in preds] == ["MATAKLEKLASIDAHLR", "MARNLEKMASIDAQLRLL"]


def test_locus_identity_and_anchor_are_carried_with_the_prediction():
    first = parse_miniprot(GFF)[0]

    assert first.contig == "Chr3"
    assert (first.start, first.end) == (19252377, 19259272)
    assert first.identity == pytest.approx(0.8675)
    assert first.anchor == "Mcry_Mcr8G11630"
    assert first.rank == 1


def test_cds_rows_do_not_become_predictions():
    assert len(parse_miniprot(GFF)) == 2


def test_an_mrna_without_a_preceding_sta_line_is_skipped_not_guessed():
    text = GFF.replace("##STA\tMATAKLEKLASIDAHLR\n", "")

    preds = parse_miniprot(text)

    assert [p.anchor for p in preds] == ["Mcry_Mcr7G08600"]


# ------------------------------------------------------------- the gate ----

def make(identity=0.99, protein="M" * 900, contig="Chr1"):
    return Prediction(contig=contig, start=1, end=9000, strand="+",
                      identity=identity, rank=1, anchor="a", protein=protein)


def test_a_model_well_supported_by_its_own_genome_is_accepted():
    ok, reason = accept(make(identity=0.99))

    assert ok and reason == "accepted"


def test_a_model_the_genome_does_not_support_is_rejected():
    # Ococ sat at 0.8550 and Pami at 0.7898 against their own assemblies
    ok, reason = accept(make(identity=0.855))

    assert not ok
    assert "identity" in reason


def test_the_gate_boundary_is_stated_not_implied():
    assert accept(make(identity=0.95))[0] is True
    assert accept(make(identity=0.9499))[0] is False


def test_a_truncated_prediction_is_rejected_even_at_high_identity():
    # the whole point is replacing truncated models, so emitting one is a bug
    ok, reason = accept(make(protein="M" * 100), min_length=300)

    assert not ok
    assert "length" in reason


def test_a_prediction_carrying_a_stop_codon_is_rejected():
    ok, reason = accept(make(protein="M" * 400 + "*" + "M" * 400))

    assert not ok
    assert "stop" in reason


# --------------------------------------------------------- best per locus --

def test_overlapping_predictions_collapse_to_the_highest_identity():
    a = Prediction("Chr1", 100, 200, "+", 0.90, 1, "x", "MMM")
    b = Prediction("Chr1", 150, 250, "+", 0.99, 2, "y", "MMMM")
    c = Prediction("Chr2", 100, 200, "+", 0.80, 1, "z", "MM")

    best = best_per_locus([a, b, c])

    assert {p.identity for p in best} == {0.99, 0.80}


def test_predictions_on_the_same_contig_that_do_not_overlap_are_both_kept():
    a = Prediction("Chr1", 100, 200, "+", 0.90, 1, "x", "MMM")
    b = Prediction("Chr1", 5000, 6000, "+", 0.95, 1, "y", "MMM")

    assert len(best_per_locus([a, b])) == 2
