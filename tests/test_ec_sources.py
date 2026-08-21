"""steps/ec_sources.py — eggNOG-mapper / CLEAN parsers into the EC-layer
prediction dict (issue #28: ECForest replacement after its known-answer
rejection).

Pure text-fixture tests, no external tools (repo convention).
"""

import sys
import types
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from steps.ec_sources import merge_ec_predictions, parse_clean, parse_emapper

EMAPPER_TEXT = """\
## some emapper preamble
#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs
ATH_AT1G53310.2\t3702.AT1G53310.2\t0.0\t1500.0\tOG@1\tlvl\tE\tphosphoenolpyruvate carboxylase\tPPC1\tGO:1,GO:2\t4.1.1.31\tko:K01595\t-\t-\t-\t-\t-\t-\t-\t-\tPEPcase
Ccac_g10054\t123.XP_1\t1e-50\t200.0\tOG@2\tlvl\tI\tSDR family member\t-\tGO:3\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tadh_short
Mcry_MULTI\t3702.AT2\t0.0\t900.0\tOG@3\tlvl\tE\tbifunctional thing\tBIF1\t-\t1.1.1.1,2.2.2.2\t-\t-\t-\t-\t-\t-\t-\t-\t-\tDom
## 3 queries scanned
"""

CLEAN_TEXT = """\
ATH_AT1G53310.2,EC:4.1.1.31/0.9915
Bvul_remote,EC:3.5.99.5/0.0011
Mcry_MULTI,EC:1.1.1.1/0.4;EC:2.2.2.2/0.9
"""


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


# --- emapper -------------------------------------------------------------

def test_parse_emapper_basic(tmp_path):
    preds = parse_emapper(_write(tmp_path, "a.emapper.annotations", EMAPPER_TEXT))
    assert preds["ATH_AT1G53310.2"] == {
        "is_enzyme": True,
        "ec": "4.1.1.31",
        "confidence": 0.0,
    }


def test_parse_emapper_dash_ec_is_non_enzyme(tmp_path):
    preds = parse_emapper(_write(tmp_path, "a.annotations", EMAPPER_TEXT))
    assert preds["Ccac_g10054"] == {
        "is_enzyme": False,
        "ec": "",
        "confidence": 0.0,
    }


def test_parse_emapper_multi_ec_keeps_first(tmp_path):
    preds = parse_emapper(_write(tmp_path, "a.annotations", EMAPPER_TEXT))
    assert preds["Mcry_MULTI"]["ec"] == "1.1.1.1"
    assert preds["Mcry_MULTI"]["is_enzyme"] is True


def test_parse_emapper_skips_comment_lines(tmp_path):
    preds = parse_emapper(_write(tmp_path, "a.annotations", EMAPPER_TEXT))
    assert len(preds) == 3


# --- CLEAN ---------------------------------------------------------------

def test_parse_clean_basic(tmp_path):
    preds = parse_clean(_write(tmp_path, "maxsep.csv", CLEAN_TEXT))
    p = preds["ATH_AT1G53310.2"]
    assert p["ec"] == "4.1.1.31"
    assert p["is_enzyme"] is True
    assert abs(p["confidence"] - 0.9915) < 1e-9


def test_parse_clean_low_confidence_still_reported(tmp_path):
    # CLEAN always emits an EC; the confidence carries the signal (ECForest
    # rejection showed uniform noise — here near-zero conf is the flag).
    preds = parse_clean(_write(tmp_path, "maxsep.csv", CLEAN_TEXT))
    assert preds["Bvul_remote"]["ec"] == "3.5.99.5"
    assert preds["Bvul_remote"]["confidence"] < 0.01


def test_parse_clean_multi_ec_takes_highest_confidence(tmp_path):
    preds = parse_clean(_write(tmp_path, "maxsep.csv", CLEAN_TEXT))
    assert preds["Mcry_MULTI"]["ec"] == "2.2.2.2"
    assert abs(preds["Mcry_MULTI"]["confidence"] - 0.9) < 1e-9


# --- merge ---------------------------------------------------------------

def test_merge_emapper_wins_on_conflict(tmp_path):
    emapper = parse_emapper(_write(tmp_path, "a.annotations", EMAPPER_TEXT))
    clean = parse_clean(_write(tmp_path, "maxsep.csv", CLEAN_TEXT))
    merged = merge_ec_predictions(emapper, clean)
    # both sources: emapper's EC is authoritative (orthology > embedding)
    assert merged["Mcry_MULTI"]["ec"] == "1.1.1.1"
    assert merged["Mcry_MULTI"]["source"] == "emapper+clean"
    assert merged["Mcry_MULTI"]["agree"] is False


def test_merge_agreement_flag_and_confidence_from_clean(tmp_path):
    emapper = parse_emapper(_write(tmp_path, "a.annotations", EMAPPER_TEXT))
    clean = parse_clean(_write(tmp_path, "maxsep.csv", CLEAN_TEXT))
    merged = merge_ec_predictions(emapper, clean)
    p = merged["ATH_AT1G53310.2"]
    assert p["agree"] is True
    assert abs(p["confidence"] - 0.9915) < 1e-9  # CLEAN conf attached


def test_merge_single_source_genes_pass_through(tmp_path):
    emapper = parse_emapper(_write(tmp_path, "a.annotations", EMAPPER_TEXT))
    clean = parse_clean(_write(tmp_path, "maxsep.csv", CLEAN_TEXT))
    merged = merge_ec_predictions(emapper, clean)
    assert merged["Ccac_g10054"]["source"] == "emapper"
    assert merged["Bvul_remote"]["source"] == "clean"
    assert merged["Bvul_remote"]["ec"] == "3.5.99.5"


def test_merge_output_compatible_with_ecforest_cache(tmp_path):
    """merged dicts must round-trip through steps.ecforest.save_cache/load_cache
    (the ec_layer downstream contract)."""
    from steps.ecforest import load_cache, save_cache

    emapper = parse_emapper(_write(tmp_path, "a.annotations", EMAPPER_TEXT))
    clean = parse_clean(_write(tmp_path, "maxsep.csv", CLEAN_TEXT))
    merged = merge_ec_predictions(emapper, clean)
    cache = tmp_path / "cache.tsv"
    save_cache(merged, cache)
    loaded = load_cache(cache)
    assert loaded["ATH_AT1G53310.2"]["ec"] == "4.1.1.31"
    assert loaded["Ccac_g10054"]["is_enzyme"] is False
