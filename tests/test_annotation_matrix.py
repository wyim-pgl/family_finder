"""annotation_matrix.py — merge every annotation axis into one per-gene TSV
(project goal: reusable multi-tool stack; distant members judged by axis
agreement, never one tool).

Pure text fixtures, no external tools (repo convention).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from annotation_matrix import build_matrix, load_axes, membership_verdict

EMAPPER = """\
#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
CORE_gene\tseed1\t0.0\t1500.0\tOG@1\tlvl\tE\tphosphoenolpyruvate carboxylase\tPPC1\tGO:1\t4.1.1.31\t-\t-\t-\t-\t-\t-\t-\t-\t-\tPEPcase
INTRUDER_gene\tseed2\t1e-50\t200.0\tOG@2\tlvl\tI\tSDR family member\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tadh_short
"""

CLEAN = """\
CORE_gene,EC:4.1.1.31/0.9915
INTRUDER_gene,EC:1.1.1.100/0.0021
REMOTE_gene,EC:3.5.99.5/0.0011
"""

FOLDSEEK = """\
gene_id\tbest_hit\tuniprot\tprotein_name\tgene_symbol\tec\tgo_ids\tqtmscore\tfident\tevalue\tbits\ttarget_taxon
CORE_gene\tAF-P10490-F1-model_v6\tP10490\tPhosphoenolpyruvate carboxylase 1\tPPC1\t4.1.1.31\tGO:1\t0.95\t1.000\t1e-100\t2000\tMesembryanthemum crystallinum
INTRUDER_gene\tAF-A0PJE2-2-F1-model_v6\tA0PJE2\tDehydrogenase/reductase SDR\tDHRS12\t1.1.-.-\tGO:3\t0.80\t0.380\t1e-20\t300\tHomo sapiens
"""

DEEPLOC = """\
Protein_ID,Localizations,Signals,Membrane types,Cytoplasm,Nucleus
CORE_gene,Cytoplasm,Nuclear export signal,Peripheral|Soluble,0.72,0.24
INTRUDER_gene,Cytoplasm,,Soluble,0.9,0.1
REMOTE_gene,Cytoplasm|Mitochondrion,Mitochondrial transit peptide,Soluble,0.5,0.1
"""

SIGNALP = """\
# SignalP-6.0\tOrganism: Eukarya
# ID\tPrediction\tOTHER\tSP(Sec/SPI)\tCS Position
CORE_gene\tOTHER\t1.000000\t0.000002
INTRUDER_gene\tOTHER\t1.000000\t0.000001
REMOTE_gene\tSP\t0.000189\t0.999789\tCS pos: 18-19. Pr: 0.9806
"""


def _axes(tmp_path):
    files = {}
    for name, text in [
        ("a.emapper.annotations", EMAPPER),
        ("clean.csv", CLEAN),
        ("foldseek.tsv", FOLDSEEK),
        ("deeploc.csv", DEEPLOC),
        ("signalp.txt", SIGNALP),
    ]:
        p = tmp_path / name
        p.write_text(text)
        files[name.split(".")[0]] = p
    return load_axes(
        emapper=files["a"],
        clean_csv=files["clean"],
        foldseek_tsv=files["foldseek"],
        deeploc_csv=files["deeploc"],
        signalp_txt=files["signalp"],
    )


# --- load_axes -----------------------------------------------------------

def test_load_axes_collects_all_genes(tmp_path):
    axes = _axes(tmp_path)
    matrix = build_matrix(axes)
    assert set(matrix) == {"CORE_gene", "INTRUDER_gene", "REMOTE_gene"}


def test_axes_are_optional(tmp_path):
    p = tmp_path / "clean.csv"
    p.write_text(CLEAN)
    axes = load_axes(clean_csv=p)
    matrix = build_matrix(axes)
    assert matrix["CORE_gene"]["clean_ec"] == "4.1.1.31"
    assert matrix["CORE_gene"]["emapper_ec"] == ""


# --- per-axis columns ----------------------------------------------------

def test_matrix_row_core_gene(tmp_path):
    row = build_matrix(_axes(tmp_path))["CORE_gene"]
    assert row["emapper_ec"] == "4.1.1.31"
    assert row["emapper_symbol"] == "PPC1"
    assert row["clean_ec"] == "4.1.1.31"
    assert abs(row["clean_conf"] - 0.9915) < 1e-9
    assert row["fs_uniprot"] == "P10490"
    assert abs(row["fs_qtm"] - 0.95) < 1e-9
    assert row["deeploc_signals"] == "Nuclear export signal"
    assert row["signalp"] == "OTHER"


def test_matrix_missing_axis_values_empty(tmp_path):
    row = build_matrix(_axes(tmp_path))["REMOTE_gene"]
    assert row["emapper_ec"] == ""       # not in emapper file
    assert row["fs_uniprot"] == ""       # foldseek no-hit
    assert row["signalp"] == "SP"


# --- membership verdict --------------------------------------------------

def test_verdict_core_member(tmp_path):
    matrix = build_matrix(_axes(tmp_path))
    v = membership_verdict(matrix["CORE_gene"], expected_ec="4.1.1.31")
    assert v["verdict"] == "member"
    assert v["n_support"] >= 3  # emapper + clean + foldseek all agree


def test_verdict_intruder(tmp_path):
    matrix = build_matrix(_axes(tmp_path))
    v = membership_verdict(matrix["INTRUDER_gene"], expected_ec="4.1.1.31")
    # emapper: no EC/SDR desc; clean: wrong EC at ~0 conf; foldseek: SDR hit
    assert v["verdict"] == "intruder"
    assert v["n_support"] == 0


def test_verdict_remote_is_review_not_intruder(tmp_path):
    matrix = build_matrix(_axes(tmp_path))
    v = membership_verdict(matrix["REMOTE_gene"], expected_ec="4.1.1.31")
    # no axis supports, but no axis CONTRADICTS with confidence either
    # (clean conf ~0 is abstention, foldseek absent) -> review, not intruder
    assert v["verdict"] == "review"


def test_dot_underscore_ids_merge_to_one_row(tmp_path):
    """foldseek PDB filenames replace '.' with '_' — Gene.1 (emapper) and
    Gene_1 (foldseek) are the same gene and must land in ONE row."""
    em = tmp_path / "b.emapper.annotations"
    em.write_text(EMAPPER.replace("CORE_gene", "Bvul_EL10.1"))
    fs = tmp_path / "fs.tsv"
    fs.write_text(FOLDSEEK.replace("CORE_gene", "Bvul_EL10_1"))
    matrix = build_matrix(load_axes(emapper=em, foldseek_tsv=fs))
    merged = [g for g in matrix if g.startswith("Bvul_EL10")]
    assert len(merged) == 1
    row = matrix[merged[0]]
    assert row["emapper_ec"] == "4.1.1.31" and row["fs_uniprot"] == "P10490"
    # display id keeps the dotted (original annotation) spelling
    assert merged[0] == "Bvul_EL10.1"


def test_verdict_evidence_strings_reusable(tmp_path):
    matrix = build_matrix(_axes(tmp_path))
    v = membership_verdict(matrix["CORE_gene"], expected_ec="4.1.1.31")
    assert "emapper" in v["support_axes"]
    assert "clean" in v["support_axes"]
    assert "foldseek" in v["support_axes"]


# --- gene-structure and expression axes (issue #38) ----------------------

GENE_STRUCTURE = """\
gene\tspecies\tsubfamily\tannotation_source\tgff_gene_id\tn_cds\tcds_len_nt\tn_introns\tintron_columns\tn_extra\tn_missing\textra_columns\tmissing_columns\tstatus
CORE_gene\tOcoc\tSF1\tHelixer\tg1\t10\t2900\t9\t118,298\t1\t0\t900\t\tdeviates from the family's conserved intron set
INTRUDER_gene\tOcoc\t\tunknown\tg2\t3\t900\tNA\t\tNA\tNA\t\t\tuncontrolled: the annotation programme is unknown
"""

EXPRESSION = """\
gene\tspecies\tsubfamily\tmean_tpm\tshare\tn_samples\tmatrix\tstatus
CORE_gene\tOcoc\tSF1\t4776.8\t0.741\t7\tRSEM_TPM.average.tsv\tmeasured
INTRUDER_gene\tCgig\t\tNA\tNA\tNA\t\texpression unavailable
"""


def _extra_axes(tmp_path):
    gs = tmp_path / "gene_structure.tsv"
    gs.write_text(GENE_STRUCTURE)
    ex = tmp_path / "expression.tsv"
    ex.write_text(EXPRESSION)
    return gs, ex


def test_gene_structure_columns_arrive_with_their_programme(tmp_path):
    gs, _ex = _extra_axes(tmp_path)
    row = build_matrix(load_axes(gene_structure_tsv=gs))["CORE_gene"]
    assert row["gs_source"] == "Helixer"
    assert row["gs_n_introns"] == "9"
    assert row["gs_extra"] == "1"


def test_an_uncontrolled_gene_carries_no_intron_count_into_the_matrix(tmp_path):
    gs, _ex = _extra_axes(tmp_path)
    row = build_matrix(load_axes(gene_structure_tsv=gs))["INTRUDER_gene"]
    assert row["gs_source"] == "unknown"
    assert row["gs_n_introns"] == "NA"
    assert "uncontrolled" in row["gs_status"]


def test_expression_columns_join(tmp_path):
    _gs, ex = _extra_axes(tmp_path)
    row = build_matrix(load_axes(expression_tsv=ex))["CORE_gene"]
    assert row["expr_mean_tpm"] == "4776.8"
    assert row["expr_share"] == "0.741"
    assert row["expr_status"] == "measured"


def test_a_species_with_no_matrix_says_so_in_the_matrix(tmp_path):
    _gs, ex = _extra_axes(tmp_path)
    row = build_matrix(load_axes(expression_tsv=ex))["INTRUDER_gene"]
    assert row["expr_status"] == "expression unavailable"
    assert row["expr_mean_tpm"] == "NA"


def test_a_gene_absent_from_the_expression_axis_is_not_read_as_silent(tmp_path):
    """The axis ran but never saw this gene. Blank would read as zero TPM."""
    _gs, ex = _extra_axes(tmp_path)
    matrix = build_matrix(load_axes(expression_tsv=ex, clean_csv=_clean(tmp_path)))
    assert matrix["REMOTE_gene"]["expr_status"] == "expression unavailable"
    assert matrix["REMOTE_gene"]["expr_mean_tpm"] == ""


def _clean(tmp_path):
    p = tmp_path / "clean2.csv"
    p.write_text(CLEAN)
    return p


def test_both_new_axes_merge_alongside_the_existing_ones(tmp_path):
    gs, ex = _extra_axes(tmp_path)
    axes = _axes(tmp_path)
    axes.update(load_axes(gene_structure_tsv=gs, expression_tsv=ex))
    row = build_matrix(axes)["CORE_gene"]
    assert row["emapper_ec"] == "4.1.1.31"
    assert row["gs_source"] == "Helixer"
    assert row["expr_mean_tpm"] == "4776.8"
