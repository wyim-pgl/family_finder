"""steps/gene_structure.py — intron structure as a GENE-MODEL QUALITY axis.

The measured PEPC result is why this axis is not a subfamily one: all nine
conserved intron positions are shared by ppc-1E1 and ppc-1E2 at 92-100% and
none is diagnostic, while the deviation rate tracks the annotation programme
(Helixer 0.59 off-set introns per gene, AUGUSTUS/UNR 0.00). These tests pin
both halves: the geometry has to be right, and no intron number may leave the
module without the programme covariate beside it.

Pure text fixtures, no GFF library and no external tools (repo convention).
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from steps.gene_structure import (
    GeneModel,
    UNKNOWN_PROGRAMME,
    build_rows,
    cds_structure,
    codon_to_column,
    conserved_columns,
    deviation_by_programme,
    deviations,
    diagnostic_columns,
    intron_columns,
    join_by_species,
    join_gff_models,
    parse_gff_cds,
    structure_report,
)

# Two genes on opposite strands plus a gene carrying a short second isoform.
GFF = """\
##gff-version 3
chr1\tHelixer\tgene\t1\t1000\t.\t+\t.\tID=g1
chr1\tHelixer\tmRNA\t1\t1000\t.\t+\t.\tID=g1.t1;Parent=g1
chr1\tHelixer\tCDS\t1\t168\t.\t+\t0\tParent=g1.t1
chr1\tHelixer\tCDS\t300\t500\t.\t+\t0\tParent=g1.t1
chr1\tHelixer\texon\t1\t168\t.\t+\t.\tParent=g1.t1
chr2\tAUGUSTUS\tgene\t1\t1000\t.\t-\t.\tID=g2
chr2\tAUGUSTUS\tmRNA\t1\t1000\t.\t-\t.\tID=g2.t1;Parent=g2
chr2\tAUGUSTUS\tCDS\t1\t168\t.\t-\t0\tParent=g2.t1
chr2\tAUGUSTUS\tCDS\t300\t500\t.\t-\t0\tParent=g2.t1
chr2\tAUGUSTUS\tmRNA\t1\t400\t.\t-\t.\tID=g2.t2;Parent=g2
chr2\tAUGUSTUS\tCDS\t300\t359\t.\t-\t0\tParent=g2.t2
"""


def _write(tmp_path, text, name="test.gff3"):
    p = tmp_path / name
    p.write_text(text)
    return p


# ---------------------------------------------------------------------------
# geometry: genomic blocks -> intron position inside the coding sequence
# ---------------------------------------------------------------------------

def test_intron_position_is_the_cumulative_cds_length_in_codons():
    # Arrange: 168 nt then 201 nt, one intron between them.
    blocks = [(1, 168), (300, 500)]

    # Act
    n_cds, cds_len, positions = cds_structure(blocks, "+")

    # Assert: 168 nt = 56 codons exactly, so the intron sits at codon 56.0.
    assert (n_cds, cds_len) == (2, 369)
    assert positions == pytest.approx([56.0])


def test_a_phase_shifted_intron_keeps_its_fraction():
    # 170 nt is 56 codons and 2 bases -- the intron interrupts codon 57.
    _n, _len, positions = cds_structure([(1, 170), (300, 500)], "+")
    assert positions == pytest.approx([56.6667], abs=1e-3)


def test_minus_strand_blocks_are_ordered_by_transcription_not_by_coordinate():
    # On the minus strand the FIRST coding block is the highest-coordinate one,
    # so ordering by start would put the intron at the wrong codon.
    _n, _len, positions = cds_structure([(1, 168), (300, 500)], "-")
    assert positions == pytest.approx([67.0])


def test_codon_position_maps_through_alignment_gaps():
    # residue 3 of MACDEF is C, which sits in column 5 of this alignment
    assert codon_to_column("MA--CDEF", 3.0) == 5
    # a phase-shifted position lands on the codon it interrupts, residue 4
    assert codon_to_column("MA--CDEF", 3.33) == 6


def test_a_codon_position_past_the_end_of_the_sequence_is_not_invented():
    assert codon_to_column("MACD", 99.0) is None


def test_intron_columns_drops_nothing_silently():
    cols = intron_columns("MA--CDEF", [3.0, 99.0])
    assert cols == [5]


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def test_parse_gff_cds_keeps_the_annotation_programme_from_column_two(tmp_path):
    models = {m.gene_id: m for m in parse_gff_cds(_write(tmp_path, GFF))}
    assert models["g1"].source == "Helixer"
    assert models["g2"].source == "AUGUSTUS"


def test_parse_gff_cds_keeps_the_longest_coding_transcript_per_gene(tmp_path):
    models = {m.gene_id: m for m in parse_gff_cds(_write(tmp_path, GFF))}
    assert models["g2"].transcript_id == "g2.t1"
    assert models["g2"].blocks == ((1, 168), (300, 500))


def test_parse_gff_cds_ignores_non_cds_features(tmp_path):
    models = {m.gene_id: m for m in parse_gff_cds(_write(tmp_path, GFF))}
    assert len(models["g1"].blocks) == 2   # the exon line must not be counted


# ---------------------------------------------------------------------------
# conserved set and per-gene deviation
# ---------------------------------------------------------------------------

def test_conserved_columns_tolerate_a_two_column_window():
    per_gene = {
        "a": [100, 500], "b": [101, 500], "c": [102, 500],
        "d": [100], "e": [100],
    }
    conserved = conserved_columns(per_gene, min_fraction=0.6, window=2)
    columns = [c["column"] for c in conserved]
    # 100/101/102 are one position seen in 5/5 genes; 500 is in 3/5 = 0.6
    assert columns == [100, 500]
    assert conserved[0]["fraction"] == pytest.approx(1.0)


def test_a_position_below_the_fraction_floor_is_not_conserved():
    per_gene = {"a": [100], "b": [100], "c": [700], "d": [], "e": []}
    columns = [c["column"] for c in
               conserved_columns(per_gene, min_fraction=0.6, window=2)]
    assert columns == []


def test_one_intron_split_by_an_alignment_indel_is_one_position():
    # The measured PEPC case: the first intron is in column 118 for some
    # species and 145 for others because an N-terminal indel is not aligned.
    # Neither half clears 60% alone; together they are 100% of the family.
    per_gene = {"a": [118, 400], "b": [118, 400], "c": [145, 400],
                "d": [145, 400], "e": [145, 400]}
    conserved = conserved_columns(per_gene, min_fraction=0.6, window=2)
    first = [c for c in conserved if c["column"] < 200]
    assert len(first) == 1
    assert first[0]["columns"] == [118, 145]
    assert first[0]["split_by_alignment"] is True
    # and no gene is then charged with an extra intron it does not have
    assert deviations(per_gene["a"], conserved, window=2)["extra"] == []


def test_one_aberrant_model_does_not_veto_a_split_position():
    # 1 gene of 21 carries both columns -- a genuine extra intron in that
    # model. The other twenty still share one position and must keep it.
    per_gene = {f"a{i}": [118] for i in range(10)}
    per_gene.update({f"b{i}": [145] for i in range(10)})
    per_gene["c"] = [118, 145]
    conserved = conserved_columns(per_gene, min_fraction=0.6, window=2)
    assert [c["columns"] for c in conserved] == [[118, 145]]


def test_a_column_only_one_gene_has_is_not_absorbed_into_a_conserved_one():
    # Nineteen genes share column 300; one has 356 instead. That is a model
    # defect and has to stay in the extra count, not join the conserved set.
    per_gene = {f"a{i}": [300] for i in range(19)}
    per_gene["b"] = [356]
    conserved = conserved_columns(per_gene, min_fraction=0.6, window=2)
    assert [c["columns"] for c in conserved] == [[300]]
    assert deviations([356], conserved, window=2)["extra"] == [356]


def test_two_positions_that_co_occur_in_a_gene_are_never_merged():
    # Gene "a" has both, so they are two introns, not one split one.
    per_gene = {"a": [118, 145], "b": [118], "c": [145], "d": [], "e": []}
    conserved = conserved_columns(per_gene, min_fraction=0.6, window=2)
    assert conserved == []


def test_a_third_intron_between_the_halves_blocks_the_merge():
    # "a" and "b" carry an intron at 130, between the two candidate halves --
    # so 118 and 145 are two positions with something in the middle, not one
    # position the alignment split.
    per_gene = {"a": [118, 130], "b": [118, 130],
                "c": [145], "d": [145], "e": [145]}
    conserved = conserved_columns(per_gene, min_fraction=0.6, window=2)
    assert not any(118 in entry["columns"] and 145 in entry["columns"]
                   for entry in conserved)


def test_two_distant_rare_introns_are_not_merged_into_a_conserved_one():
    per_gene = {"a": [100], "b": [100], "c": [700], "d": [], "e": []}
    assert conserved_columns(per_gene, min_fraction=0.6, window=2) == []


def test_deviations_separate_extra_introns_from_missing_ones():
    conserved = [{"column": 100}, {"column": 500}]
    dev = deviations([101, 900], conserved, window=2)
    assert dev["extra"] == [900]
    assert dev["missing"] == [500]


# ---------------------------------------------------------------------------
# the redefinition: diagnostic power, measured rather than assumed
# ---------------------------------------------------------------------------

def test_a_position_present_in_only_one_group_is_diagnostic():
    per_gene = {"a1": [200], "a2": [200], "a3": [201],
                "b1": [], "b2": [], "b3": []}
    groups = {"A": ["a1", "a2", "a3"], "B": ["b1", "b2", "b3"]}
    conserved = [{"column": 200}]
    rows = diagnostic_columns(per_gene, groups, conserved, window=2)
    assert rows[0]["diagnostic"] is True
    assert rows[0]["by_group"]["A"] == pytest.approx(1.0)
    assert rows[0]["by_group"]["B"] == pytest.approx(0.0)


def test_a_position_shared_by_both_groups_is_not_diagnostic():
    # This is the PEPC case in miniature: conserved everywhere, therefore
    # useless for telling the subfamilies apart.
    per_gene = {"a1": [200], "a2": [201], "b1": [200], "b2": [202]}
    groups = {"A": ["a1", "a2"], "B": ["b1", "b2"]}
    rows = diagnostic_columns(per_gene, groups, [{"column": 200}], window=2)
    assert rows[0]["diagnostic"] is False
    assert rows[0]["by_group"] == {"A": pytest.approx(1.0), "B": pytest.approx(1.0)}


# ---------------------------------------------------------------------------
# the covariate: no intron number without its annotation programme
# ---------------------------------------------------------------------------

def _model(gene, source, blocks=((1, 168), (300, 500))):
    return GeneModel(gene_id=gene, transcript_id=gene + ".t1", source=source,
                     strand="+", blocks=tuple(blocks))


def test_every_emitted_row_carries_the_annotation_programme():
    rows = build_rows(
        {"Xx_g1": "MACDEFGHIK"},
        {"Xx_g1": _model("g1", "Helixer")},
        conserved=[],
    )
    assert rows[0]["annotation_source"] == "Helixer"
    assert set(rows[0]) >= {"gene", "species", "annotation_source",
                            "n_introns", "n_extra", "n_missing", "status"}


def test_a_gene_with_no_annotation_programme_gets_no_intron_count():
    # An uncontrolled count is exactly what must not be emitted: without the
    # programme there is no way to tell a real intron gain from an annotator's.
    rows = build_rows(
        {"Xx_g1": "MACDEFGHIK"},
        {"Xx_g1": _model("g1", "")},
        conserved=[],
    )
    row = rows[0]
    assert row["annotation_source"] == UNKNOWN_PROGRAMME
    assert row["n_introns"] is None
    assert row["n_extra"] is None and row["n_missing"] is None
    assert "uncontrolled" in row["status"]


def test_a_dot_in_gff_column_two_is_not_an_annotation_programme():
    # GFF3 writes "." for an undefined source. Two of the PEPC clan's genomes
    # do, and taking it as a programme name would invent a control group.
    rows = build_rows({"Xx_g1": "MACDEFGHIK"},
                      {"Xx_g1": _model("g1", ".")}, conserved=[])
    assert rows[0]["annotation_source"] == UNKNOWN_PROGRAMME
    assert rows[0]["n_introns"] is None


def test_deviation_rate_is_reported_per_annotation_programme():
    rows = [
        {"annotation_source": "Helixer", "n_extra": 1, "n_missing": 0},
        {"annotation_source": "Helixer", "n_extra": 0, "n_missing": 1},
        {"annotation_source": "AUGUSTUS", "n_extra": 0, "n_missing": 0},
        {"annotation_source": UNKNOWN_PROGRAMME, "n_extra": None,
         "n_missing": None},
    ]
    by = {r["annotation_source"]: r for r in deviation_by_programme(rows)}
    assert by["Helixer"]["n_genes"] == 2
    assert by["Helixer"]["extra_per_gene"] == pytest.approx(0.5)
    assert by["AUGUSTUS"]["extra_per_gene"] == pytest.approx(0.0)
    # the uncontrolled gene is counted but contributes no rate
    assert by[UNKNOWN_PROGRAMME]["n_genes"] == 1
    assert by[UNKNOWN_PROGRAMME]["extra_per_gene"] is None


def test_the_report_says_the_deviation_is_confounded_when_programmes_differ():
    alignment = {f"Xx_g{i}": "M" * 40 for i in range(1, 5)}
    models = {"Xx_g1": _model("g1", "Helixer", [(1, 60), (100, 160)]),
              "Xx_g2": _model("g2", "Helixer", [(1, 60), (100, 160)]),
              "Xx_g3": _model("g3", "AUGUSTUS", [(1, 60), (100, 160)]),
              "Xx_g4": _model("g4", "AUGUSTUS", [(1, 60), (100, 160)])}
    report = structure_report(alignment, models)
    assert "by_programme" in report["summary"]
    assert "confounded" in report["summary"]


# ---------------------------------------------------------------------------
# joining GFF ids onto pep ids -- through match_ids, with a cap
# ---------------------------------------------------------------------------

def test_join_matches_a_transcript_suffix_the_gff_does_not_carry():
    models = [_model("OcoChr03G21370", "EVM")]
    mapping, report = join_gff_models(["Ococ_OcoChr03G21370.t1"], models)
    assert mapping["Ococ_OcoChr03G21370.t1"].gene_id == "OcoChr03G21370"
    assert report["n_unmatched"] == 0


def test_join_recovers_an_assembly_version_suffix():
    # `Bvul_EL10Ac4g08333.EL10_1.0` against a GFF that knows only the locus.
    models = [_model("EL10Ac4g08333", "phytozomev13")]
    mapping, report = join_gff_models(["Bvul_EL10Ac4g08333.EL10_1.0"], models)
    assert mapping["Bvul_EL10Ac4g08333.EL10_1.0"].gene_id == "EL10Ac4g08333"
    assert report["level"] == "dot-trimmed"


def test_join_refuses_when_too_many_ids_stay_unmatched():
    models = [_model("known", "EVM")]
    with pytest.raises(ValueError, match="unmatched"):
        join_gff_models(["Xx_known", "Xx_absent1", "Xx_absent2"],
                        models, max_unmatched=0.15)


def test_join_reports_the_unmatched_ids_when_under_the_cap():
    models = [_model("known", "EVM")]
    mapping, report = join_gff_models(["Xx_known", "Xx_absent"], models,
                                      max_unmatched=0.9)
    assert set(mapping) == {"Xx_known"}
    assert report["unmatched"] == ["Xx_absent"]
    assert report["unmatched_fraction"] == pytest.approx(0.5)


def test_species_are_joined_separately_so_a_shared_locus_name_survives():
    # `g1` is an AUGUSTUS-style name two genomes can both use. Pooling the
    # GFFs would make it ambiguous and drop BOTH genes.
    by_species = {"Aa": [_model("g1", "AUGUSTUS")],
                  "Bb": [_model("g1", "Helixer")]}
    mapping, report = join_by_species(["Aa_g1", "Bb_g1"], by_species)
    assert mapping["Aa_g1"].source == "AUGUSTUS"
    assert mapping["Bb_g1"].source == "Helixer"
    assert report["n_unmatched"] == 0


def test_a_species_with_no_gff_is_named_not_merely_missing():
    mapping, report = join_by_species(["Aa_g1", "Zz_g9"],
                                      {"Aa": [_model("g1", "EVM")]})
    assert set(mapping) == {"Aa_g1"}
    assert report["species_without_gff"] == ["Zz"]
    assert report["unmatched"] == ["Zz_g9"]


def test_join_by_species_refuses_over_the_cap():
    with pytest.raises(ValueError, match="unmatched"):
        join_by_species(["Aa_g1", "Zz_g9"], {"Aa": [_model("g1", "EVM")]},
                        max_unmatched=0.1)


def test_an_unrelated_gene_that_squashes_alike_does_not_cost_the_match():
    # `evm.model.Hap1Chr1.244` and `evm.model.Hap1Chr12.44` squash to the same
    # string. Treating the whole genome as the reference set made that pair
    # ambiguous and lost the Dcar gene the family actually asked for.
    models = [_model("evm.model.Hap1Chr1.244", "EVM"),
              _model("evm.model.Hap1Chr12.44", "EVM")]
    mapping, report = join_gff_models(["Dcar_evm.model.Hap1Chr1.244"], models)
    assert mapping["Dcar_evm.model.Hap1Chr1.244"].gene_id == \
        "evm.model.Hap1Chr1.244"
    assert report["n_unmatched"] == 0
    assert report["n_ambiguous_aliases"] == 0


def test_an_alias_two_entries_claim_is_given_up_not_guessed():
    # `shared` and `shared.t1` cannot be separated once the transcript suffix
    # is stripped. One spelling is given up and counted, so the join reports a
    # reduced reference rather than handing the family an arbitrary model.
    models = [_model("shared", "EVM"), _model("shared.t1", "Helixer")]
    mapping, report = join_gff_models(["Xx_shared.t2"], models,
                                      max_unmatched=1.0)
    assert report["n_ambiguous_aliases"] >= 1
    assert mapping["Xx_shared.t2"].gene_id == "shared"
