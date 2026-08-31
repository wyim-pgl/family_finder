"""Tests for steps/label_symbols.py — multi-symbol family adjudication.

Written against MULTI_SYMBOL_LABEL_DESIGN.md plus the review findings of
2026-08-31 (fail-open booleans, UNPARSED-as-abstention, calibrated-only
coverage, catalog member leaks, roman space boundary, collision guard).
The module is pure (plain dicts/dataclasses, no ete4), so no stubs.

Deliberate deviations from the design doc are documented in the module
docstring; the invariant tests here enforce the resolved versions.
"""

import json
import random
from pathlib import Path

import pytest

from steps.label_symbols import (
    AUTO_VERDICTS,
    NormalizationCollision,
    OverrideConflict,
    Policy,
    adjudicate_family,
    different_stem_margin,
    load_overrides,
    normalize_symbol,
    parse_symbol,
    stale_overrides,
    symbol_key,
)

OVERRIDE_PATH = Path(__file__).resolve().parents[1] / "data" / "symbol_stem_overrides_v1.tsv"


@pytest.fixture(scope="module")
def overrides():
    return load_overrides(OVERRIDE_PATH.read_text().splitlines())


def policy(**kw):
    kw.setdefault("min_label_coverage", 0.2)
    return Policy(**kw)


def call(source, gene, symbol, gate_passed=True, gate_reason="", calibrated=None):
    if calibrated is None:
        calibrated = source != "afdb_swissprot"
    return {
        "source": source,
        "gene": gene,
        "symbol": symbol,
        "gate_passed": gate_passed,
        "gate_reason": gate_reason,
        "calibrated": calibrated,
    }


# ---------------------------------------------------------------- §12.1 normalization


def test_normalization_is_idempotent():
    for raw in ["  CRY  2 ", "Cipk10 (Sip1)", "Fba2;1", "PHY‐C", "Snrk2.3"]:
        once = normalize_symbol(raw)
        twice = normalize_symbol(once.display)
        assert once.key == twice.key
        assert once.display == twice.display


def test_cry_space_and_case_variants_share_a_key():
    keys = {normalize_symbol(r).key for r in ["CRY2", "Cry 2", "Cry2", "cry2"]}
    assert len(keys) == 1


def test_unicode_dash_variants_unify():
    assert normalize_symbol("Vha‑a3").key == normalize_symbol("Vha-a3").key


def test_two_independent_word_tokens_are_not_glued():
    norm = normalize_symbol("NAC domain")
    assert norm.status == "AMBIGUOUS"


def test_parenthetical_alias_splits_into_primary_and_alias():
    norm = normalize_symbol("Cipk10 (Sip1)")
    assert norm.status == "OK"
    assert norm.primary == "Cipk10"
    assert norm.aliases == ("Sip1",)


def test_unbalanced_parentheses_are_ambiguous():
    assert normalize_symbol("Cipk10 (Sip1").status == "AMBIGUOUS"
    assert normalize_symbol("Cipk10 (Sip1))").status == "AMBIGUOUS"


def test_semicolon_allele_is_preserved_in_key():
    assert ";1" in normalize_symbol("Fba2;1").key


def test_greek_letters_named_in_key_only():
    norm = normalize_symbol("PPCα")
    assert "alpha" in norm.key
    assert "α" in norm.display


# ------------------------------------------------------------------- §12.1 overrides


def test_override_file_loads(overrides):
    assert overrides[(symbol_key("PHYC"), "*")].stem == "PHY"
    assert overrides[(symbol_key("Vha-a3"), "*")].stem == "Vha-a"


def test_conflicting_duplicate_override_raises():
    lines = [
        "symbol_key\tstem\tsuffix\tscope_source\treason\tadded_in",
        "phyc\tPHY\tC\t*\tr\tv1",
        "phyc\tPHYC\t\t*\tr\tv1",
    ]
    with pytest.raises(OverrideConflict):
        load_overrides(lines)


def test_identical_duplicate_override_is_tolerated():
    lines = [
        "symbol_key\tstem\tsuffix\tscope_source\treason\tadded_in",
        "phyc\tPHY\tC\t*\tr\tv1",
        "phyc\tPHY\tC\t*\tr\tv1",
    ]
    assert load_overrides(lines)[(symbol_key("PHYC"), "*")].stem == "PHY"


def test_unnormalized_override_key_raises_instead_of_never_matching():
    lines = [
        "symbol_key\tstem\tsuffix\tscope_source\treason\tadded_in",
        "PHYC\tPHY\tC\t*\tr\tv1",
    ]
    with pytest.raises(ValueError, match="normalized"):
        load_overrides(lines)


def test_stale_overrides_are_reported(overrides):
    used = {(symbol_key("PHYC"), "*"), (symbol_key("PHYD"), "*")}
    stale = stale_overrides(overrides, used)
    assert (symbol_key("Vha-a3"), "*") in stale
    assert (symbol_key("PHYC"), "*") not in stale


def test_source_scoped_override_applies_only_to_its_source():
    lines = [
        "symbol_key\tstem\tsuffix\tscope_source\treason\tadded_in",
        "atppa1\tPPa\t1\tmcry_curated\tAt prefix alias\tv1",
    ]
    ov = load_overrides(lines)
    scoped = parse_symbol(normalize_symbol("AtPPa1"), ov, source="mcry_curated")
    assert (scoped.stem, scoped.status) == ("PPa", "OVERRIDE")
    other = parse_symbol(normalize_symbol("AtPPa1"), ov, source="ath_diamond")
    assert (other.stem, other.status) == ("AtPPa", "PARSED")


# -------------------------------------------------------------------------- aliases


ALIAS_PATH = Path(__file__).resolve().parents[1] / "data" / "symbol_alias_v1.tsv"


def alias_lines(pairs):
    lines = ["stem_key_a\tstem_key_b\tn_evidence\texample_family\tsource\tadded_in"]
    lines += ["%s\t%s\t1\tF\ttest\tv1" % (a, b) for a, b in pairs]
    return lines


def test_alias_table_builds_transitive_classes():
    from steps.label_symbols import load_aliases
    aliases = load_aliases(alias_lines([("dnm", "drp"), ("drp", "adl")]))
    assert aliases["dnm"] == aliases["drp"] == aliases["adl"] == "adl"


def test_alias_keys_must_be_normalized():
    from steps.label_symbols import load_aliases
    with pytest.raises(ValueError, match="normalized"):
        load_aliases(alias_lines([("DNM", "drp")]))


def test_repo_alias_table_loads():
    from steps.label_symbols import load_aliases
    aliases = load_aliases(ALIAS_PATH.read_text().splitlines())
    assert aliases["dnm"] == aliases["drp"]  # Homo sapiens dynamin evidence
    assert aliases["apg"] == aliases["atg"]


def test_alias_equivalence_collapses_a_cross_axis_conflict(overrides):
    """TT/CHS-type case: without the alias it is a stem conflict, with the
    alias it is one class and the curated spelling wins the display."""
    from steps.label_symbols import load_aliases
    members = {"g1", "g2"}
    calls = [call("mcry_curated", "g1", "Drp3"), call("ath_diamond", "g2", "DNM1")]
    without = adjudicate_family("F1", members, calls, overrides=overrides,
                                policy=policy())
    assert without["verdict"] == "NEEDS_REVIEW_MULTI_STEM"
    aliases = load_aliases(alias_lines([("dnm", "drp")]))
    with_alias = adjudicate_family("F1", members, calls, overrides=overrides,
                                   policy=policy(), aliases=aliases)
    assert with_alias["verdict"] == "SINGLE_STEM_MULTI_SUFFIX"
    fam = [r for r in with_alias["rows"] if r["transfer_tier"] == "stem_auto"][0]
    assert fam["stem"] == "Drp"  # curated spelling of the class


# ---------------------------------------------------------------------- §4.2 parser


KNOWN_ANSWERS = [
    ("Tpt1", "Tpt", "1"),
    ("Tpt2", "Tpt", "2"),
    ("Bam5", "Bam", "5"),
    ("Pth1;4", "Pth", "1;4"),
    ("CKB1;2", "CKB", "1;2"),
    ("Fba2;1", "Fba", "2;1"),
    ("Ppck1a", "Ppck", "1a"),
    ("Ppck1b", "Ppck", "1b"),
    ("Snrk2.3", "Snrk", "2.3"),
    ("Cry 2", "Cry", "2"),
    ("CRY2", "CRY", "2"),
    ("Almt2;2", "Almt", "2;2"),
    ("ZFN1", "ZFN", "1"),
    ("JMJ14", "JMJ", "14"),
    ("AtPPa1", "AtPPa", "1"),
]


@pytest.mark.parametrize("raw,stem,suffix", KNOWN_ANSWERS)
def test_parser_known_answers(overrides, raw, stem, suffix):
    norm = normalize_symbol(raw)
    parsed = parse_symbol(norm, overrides)
    assert parsed.status in ("PARSED", "OVERRIDE")
    assert (parsed.stem, parsed.suffix) == (stem, suffix)


def test_cipk_alias_is_evidence_not_a_vote(overrides):
    norm = normalize_symbol("Cipk10 (Sip1)")
    parsed = parse_symbol(norm, overrides)
    assert (parsed.stem, parsed.suffix) == ("Cipk", "10")
    assert norm.aliases == ("Sip1",)


def test_vha_a_is_ambiguous_without_override():
    parsed = parse_symbol(normalize_symbol("Vha-a3"), {})
    assert parsed.status == "AMBIGUOUS"


def test_vha_a_override_wins(overrides):
    parsed = parse_symbol(normalize_symbol("Vha-a3"), overrides)
    assert parsed.status == "OVERRIDE"
    assert (parsed.stem, parsed.suffix) == ("Vha-a", "3")


def test_phyc_is_bare_stem_without_override_and_split_with(overrides):
    bare = parse_symbol(normalize_symbol("PHYC"), {})
    assert (bare.stem, bare.suffix, bare.status) == ("PHYC", "", "PARSED")
    over = parse_symbol(normalize_symbol("PHYC"), overrides)
    assert (over.stem, over.suffix, over.status) == ("PHY", "C", "OVERRIDE")


def test_terminal_letters_are_not_stripped(overrides):
    parsed = parse_symbol(normalize_symbol("Acl"), overrides)
    assert (parsed.stem, parsed.suffix, parsed.status) == ("Acl", "", "PARSED")


def test_attached_roman_terminal_is_not_parsed(overrides):
    parsed = parse_symbol(normalize_symbol("SODI"), overrides)
    assert (parsed.stem, parsed.suffix) == ("SODI", "")


@pytest.mark.parametrize("raw", ["Fer-II", "Fer_II", "Fer II"])
def test_bounded_roman_terminal_is_a_suffix(overrides, raw):
    parsed = parse_symbol(normalize_symbol(raw), overrides)
    assert (parsed.stem, parsed.suffix) == ("Fer", "II")


def test_digits_only_or_empty_is_unparsed(overrides):
    assert parse_symbol(normalize_symbol("123"), overrides).status == "UNPARSED"
    assert parse_symbol(normalize_symbol(""), overrides).status == "UNPARSED"


@pytest.mark.parametrize("raw", ["Bam1-3", "Cry-2", "Hsp70-1", "Nhx1.2.3",
                                 "Ppck1/2", "Bam1;1;2"])
def test_no_grammar_match_with_separators_is_unparsed_not_a_stem(overrides, raw):
    """Review C2: zero grammar candidates must abstain, not assert a stem."""
    parsed = parse_symbol(normalize_symbol(raw), overrides)
    assert parsed.status == "UNPARSED"
    assert parsed.stem == ""


# ------------------------------------------------------------------ §12.2 axis gates


def test_margin_is_a_ratio_against_a_different_stem_only():
    hits = [(500.0, "bam"), (480.0, "bam"), (300.0, "fba")]
    margin, status = different_stem_margin(hits)
    assert status == "OK"
    assert margin == pytest.approx((500.0 - 300.0) / 500.0)


def test_no_rival_stem_is_na_not_infinite():
    margin, status = different_stem_margin([(500.0, "bam"), (480.0, "bam")])
    assert margin is None
    assert status == "NA_single_stem_candidate"


def test_equal_best_bits_abstains_regardless_of_order():
    hits = [(500.0, "bam"), (500.0, "fba")]
    for ordering in (hits, hits[::-1]):
        margin, status = different_stem_margin(ordering)
        assert margin is None
        assert status == "ABSTAIN_tied_best"


def test_nonpositive_best_bits_abstains():
    margin, status = different_stem_margin([(0.0, "bam"), (-5.0, "fba")])
    assert margin is None
    assert status == "ABSTAIN_nonpositive_best"


def test_unknown_call_source_raises(overrides):
    with pytest.raises(ValueError, match="unknown call source"):
        adjudicate_family("F1", {"g1"}, [call("afdb_swissport", "g1", "Bam1")],
                          overrides=overrides, policy=policy())


def test_gate_failed_call_is_abstention_not_a_vote(overrides):
    members = {"g1", "g2"}
    calls = [
        call("mcry_curated", "g1", "Bam1"),
        call("ath_diamond", "g2", "Fba2", gate_passed=False, gate_reason="low margin"),
    ]
    res = adjudicate_family("F1", members, calls, overrides=overrides, policy=policy())
    assert res["verdict"] in ("SINGLE_STEM", "SINGLE_STEM_MULTI_SUFFIX")
    fam_rows = [r for r in res["rows"] if r["transfer_tier"] == "stem_auto"]
    assert fam_rows[0]["stem"] == "Bam"
    assert res["coverage"] == 0.5  # only g1 is an accepted evidence gene


def test_duplicate_hits_collapse_to_one_evidence_unit(overrides):
    members = {"g1", "g2"}
    calls = [
        call("ath_diamond", "g1", "Bam1"),
        call("ath_diamond", "g1", "Bam1"),
    ]
    res = adjudicate_family("F1", members, calls, overrides=overrides, policy=policy())
    assert res["coverage"] == 0.5
    assert res["support"] == 1.0


def test_second_distinct_symbol_keeps_no_vote_but_stays_visible(overrides):
    """Review M7: the dropped symbol's suffix still reaches the tree gate."""
    members = {"g1", "g2"}
    calls = [call("ath_diamond", "g1", "Bam1"), call("ath_diamond", "g1", "Bam7")]
    res = adjudicate_family("F1", members, calls, overrides=overrides, policy=policy())
    assert res["coverage"] == 0.5  # still one voting gene
    suffixes = sorted(r["suffix"] for r in res["rows"]
                      if r["transfer_tier"] == "suffix_needs_tree_gate")
    assert suffixes == ["1", "7"]
    assert any(a.startswith("VARIANT:") for a in res["audit"])


def test_high_support_low_coverage_does_not_pass_the_gate(overrides):
    members = {f"g{i}" for i in range(20)}
    calls = [call("ath_diamond", "g1", "Bam1"), call("ath_diamond", "g2", "Bam2")]
    res = adjudicate_family("F1", members, calls, overrides=overrides, policy=policy())
    assert res["support"] == 1.0
    assert res["coverage"] == pytest.approx(0.1)
    assert res["verdict"] == "NEEDS_REVIEW_LOW_COVERAGE"
    assert all(r["transfer_tier"] == "needs_review" for r in res["rows"])


def test_uncalibrated_genes_do_not_unlock_the_coverage_gate(overrides):
    """Review H3: AFDB pre-calibration may corroborate, never unlock."""
    members = {f"g{i}" for i in range(10)}
    calls = [call("ath_diamond", "g1", "Bam1")]
    calls += [call("afdb_swissprot", f"g{i}", "Bam1", calibrated=False)
              for i in range(2, 10)]
    res = adjudicate_family("F1", members, calls, overrides=overrides,
                            policy=policy(min_label_coverage=0.5))
    assert res["verdict"] == "NEEDS_REVIEW_LOW_COVERAGE"
    assert res["coverage"] == pytest.approx(0.1)
    assert res["coverage_all_sources"] == pytest.approx(0.9)


def test_normalization_collision_raises_on_divergent_parse(overrides):
    """Review M1/§3.3: same key, different stems, one side not overridden."""
    lines = [
        "symbol_key\tstem\tsuffix\tscope_source\treason\tadded_in",
        "atppa1\tPPa\t1\tath_diamond\tAt prefix alias\tv1",
    ]
    ov = load_overrides(lines)
    calls = [
        call("ath_diamond", "g1", "AtPPa1"),     # override -> PPa
        call("mcry_curated", "g2", "AtPPa1"),    # grammar -> AtPPa
    ]
    with pytest.raises(NormalizationCollision):
        adjudicate_family("F1", {"g1", "g2"}, calls, overrides=ov, policy=policy())


# --------------------------------------------------------------- §12.3 state machine


def test_no_calls_is_no_evidence(overrides):
    res = adjudicate_family("F1", {"g1"}, [], overrides=overrides, policy=policy())
    assert res["verdict"] == "NO_EVIDENCE"
    assert res["rows"] == []


def test_evidence_gene_outside_family_raises(overrides):
    with pytest.raises(ValueError):
        adjudicate_family("F1", {"g1"}, [call("ath_diamond", "gX", "Bam1")],
                          overrides=overrides, policy=policy())


def test_single_stem_multi_suffix_emits_stem_auto_plus_gated_suffixes(overrides):
    members = {"g1", "g2", "g3", "g4"}
    calls = [call("mcry_curated", f"g{i}", f"Bam{n}") for i, n in enumerate((1, 3, 5, 7), 1)]
    res = adjudicate_family("F1", members, calls, overrides=overrides, policy=policy())
    assert res["verdict"] == "SINGLE_STEM_MULTI_SUFFIX"
    fam = [r for r in res["rows"] if r["transfer_tier"] == "stem_auto"]
    assert len(fam) == 1 and fam[0]["stem"] == "Bam"
    assert fam[0]["symbol"] == "Bam"  # a stem row never displays a suffixed symbol
    suffixes = sorted(r["suffix"] for r in res["rows"] if r["transfer_tier"] == "suffix_needs_tree_gate")
    assert suffixes == ["1", "3", "5", "7"]


def test_parse_failure_in_an_accepted_call_forces_review(overrides):
    members = {"g1", "g2"}
    calls = [call("mcry_curated", "g1", "Bam1"), call("mcry_curated", "g2", "Vha-a9")]
    res = adjudicate_family("F1", members, calls, overrides={}, policy=policy())
    assert res["verdict"] == "NEEDS_REVIEW_PARSE"
    assert any("PARSE_AMBIGUOUS" in r["reason_code"] for r in res["rows"])


def test_curation_vs_ath_stem_conflict_is_never_weight_resolved(overrides):
    members = {"g1", "g2", "g3"}
    calls = [
        call("mcry_curated", "g1", "Nhx2"),
        call("mcry_curated", "g2", "Nhx1"),
        call("ath_diamond", "g3", "Acl"),
    ]
    res = adjudicate_family("F1", members, calls, overrides=overrides, policy=policy())
    assert res["verdict"] == "NEEDS_REVIEW_MULTI_STEM"
    assert any("STEM_CONFLICT" in r["reason_code"] for r in res["rows"])
    assert all(r["transfer_tier"] == "needs_review" for r in res["rows"])


def test_uncalibrated_afdb_conflict_is_audit_not_veto(overrides):
    members = {"g1", "g2", "g3"}
    calls = [
        call("mcry_curated", "g1", "Bam1"),
        call("ath_diamond", "g2", "Bam3"),
        call("afdb_swissprot", "g3", "Fba1", calibrated=False),
    ]
    res = adjudicate_family("F1", members, calls, overrides=overrides, policy=policy())
    assert res["verdict"] == "SINGLE_STEM_MULTI_SUFFIX"
    assert any("AFDB_CONFLICT" in a for a in res["audit"])
    fam = [r for r in res["rows"] if r["transfer_tier"] == "stem_auto"][0]
    assert fam["stem"] == "Bam" and fam["reason_code"] == ""


def test_afdb_only_uncalibrated_is_needs_review(overrides):
    members = {"g1"}
    calls = [call("afdb_swissprot", "g1", "Bam1", calibrated=False)]
    res = adjudicate_family("F1", members, calls, overrides=overrides, policy=policy())
    assert res["verdict"] == "NEEDS_REVIEW_UNCALIBRATED_SOURCE"


def test_orthographic_variants_pick_display_by_source_priority(overrides):
    members = {"g1", "g2"}
    calls = [call("ath_diamond", "g2", "CRY2"), call("mcry_curated", "g1", "Cry2")]
    res = adjudicate_family("F1", members, calls, overrides=overrides, policy=policy())
    assert res["verdict"] in ("SINGLE_STEM", "ORTHOGRAPHIC_VARIANT")
    suffix_row = [r for r in res["rows"] if r["suffix"] == "2"][0]
    assert suffix_row["symbol"] == "Cry2"  # curated spelling wins the display


HIGH_CATALOG = {
    "OG_A": {"grade": "HIGH", "members": {"g1", "g2"}},
    "OG_B": {"grade": "HIGH", "members": {"g3", "g4"}},
}

TREES_OK = {
    ("nhx", "OG_A"): {"transferable": True, "grade": "HIGH", "coverage": 0.9},
    ("acl", "OG_B"): {"transferable": True, "grade": "HIGH", "coverage": 0.9},
}

MULTI_CALLS = [
    call("mcry_curated", "g1", "Nhx2"),
    call("mcry_curated", "g2", "Nhx1"),
    call("ath_diamond", "g3", "Acl"),
    call("ath_diamond", "g4", "Acl"),
]


def test_disjoint_high_partition_auto_confirms_subfamily_rows_only(overrides):
    res = adjudicate_family(
        "F1", {"g1", "g2", "g3", "g4"}, MULTI_CALLS, overrides=overrides,
        catalog=HIGH_CATALOG, tree_results=TREES_OK, policy=policy())
    assert res["verdict"] == "MULTI_STEM_PARTITIONED"
    sub = [r for r in res["rows"] if r["level"] == "subfamily"]
    assert {r["target_id"] for r in sub} == {"OG_A", "OG_B"}
    assert {r["stem"] for r in sub} == {"Nhx", "Acl"}
    assert all(r["transfer_tier"] == "stem_auto" for r in sub)
    fam = [r for r in res["rows"] if r["level"] == "family"]
    assert len(fam) == 1
    assert fam[0]["transfer_tier"] == "needs_review"
    assert fam[0]["stem"] == ""  # no forced family-level single stem


def test_two_stems_in_one_high_og_is_review(overrides):
    catalog = {"OG_A": {"grade": "HIGH", "members": {"g1", "g2", "g3", "g4"}}}
    res = adjudicate_family(
        "F1", {"g1", "g2", "g3", "g4"}, MULTI_CALLS, overrides=overrides,
        catalog=catalog, tree_results=TREES_OK, policy=policy())
    assert res["verdict"] == "NEEDS_REVIEW_MULTI_STEM"


def test_provisional_catalog_blocks_auto_confirmation(overrides):
    catalog = {
        "OG_A": {"grade": "HIGH", "members": {"g1", "g2"}},
        "OG_B": {"grade": "PROVISIONAL", "members": {"g3", "g4"}},
    }
    res = adjudicate_family(
        "F1", {"g1", "g2", "g3", "g4"}, MULTI_CALLS, overrides=overrides,
        catalog=catalog, tree_results=TREES_OK, policy=policy())
    sub = [r for r in res["rows"] if r["level"] == "subfamily"]
    assert {r["target_id"] for r in sub} == {"OG_A"}
    fam = [r for r in res["rows"] if r["level"] == "family"][0]
    assert "CATALOG_NOT_HIGH" in fam["reason_code"]
    assert any(w["reason"] == "CATALOG_NOT_HIGH" and w["og"] == "OG_B"
               for w in res["withheld"])


def test_tree_not_transferable_blocks_the_subfamily_row(overrides):
    trees = dict(TREES_OK)
    trees[("acl", "OG_B")] = {"transferable": False, "grade": "HIGH", "coverage": 0.9}
    res = adjudicate_family(
        "F1", {"g1", "g2", "g3", "g4"}, MULTI_CALLS, overrides=overrides,
        catalog=HIGH_CATALOG, tree_results=trees, policy=policy())
    sub = [r for r in res["rows"] if r["level"] == "subfamily"]
    assert {r["target_id"] for r in sub} == {"OG_A"}
    fam = [r for r in res["rows"] if r["level"] == "family"][0]
    assert "TREE_NOT_TRANSFERABLE" in fam["reason_code"]


def test_missing_tree_is_distinguished_from_a_negative_tree(overrides):
    trees = {("nhx", "OG_A"): TREES_OK[("nhx", "OG_A")]}  # no entry for Acl/OG_B
    res = adjudicate_family(
        "F1", {"g1", "g2", "g3", "g4"}, MULTI_CALLS, overrides=overrides,
        catalog=HIGH_CATALOG, tree_results=trees, policy=policy())
    fam = [r for r in res["rows"] if r["level"] == "family"][0]
    assert "TREE_MISSING" in fam["reason_code"]
    assert "TREE_NOT_TRANSFERABLE" not in fam["reason_code"]


def test_catalog_member_leak_withholds_the_row(overrides):
    """Review H4: a label never covers a gene outside the family."""
    catalog = {
        "OG_A": {"grade": "HIGH", "members": {"g1", "g2", "OUTSIDER"}},
        "OG_B": {"grade": "HIGH", "members": {"g3", "g4"}},
    }
    res = adjudicate_family(
        "F1", {"g1", "g2", "g3", "g4"}, MULTI_CALLS, overrides=overrides,
        catalog=catalog, tree_results=TREES_OK, policy=policy())
    sub = [r for r in res["rows"] if r["level"] == "subfamily"]
    assert {r["target_id"] for r in sub} == {"OG_B"}
    assert any(w["reason"] == "CATALOG_MEMBER_MISMATCH" and w["og"] == "OG_A"
               for w in res["withheld"])


def test_partial_high_coverage_does_not_spread_the_label(overrides):
    members = {"g1", "g2", "g3", "g4", "g5", "g6"}  # g5,g6 uncovered by HIGH OGs
    res = adjudicate_family(
        "F1", members, MULTI_CALLS, overrides=overrides,
        catalog=HIGH_CATALOG, tree_results=TREES_OK, policy=policy())
    assert res["verdict"] == "MULTI_STEM_PARTITIONED"
    fam = [r for r in res["rows"] if r["level"] == "family"][0]
    assert "PARTIAL_HIGH_COVERAGE" in fam["reason_code"]
    labelled = set()
    for r in res["rows"]:
        if r["level"] == "subfamily":
            labelled |= HIGH_CATALOG[r["target_id"]]["members"]
    assert "g5" not in labelled and "g6" not in labelled


def test_shuffled_input_order_is_byte_identical(overrides):
    members = {"g1", "g2", "g3", "g4"}
    baseline = None
    for seed in range(5):
        calls = list(MULTI_CALLS)
        random.Random(seed).shuffle(calls)
        res = adjudicate_family(
            "F1", members, calls, overrides=overrides,
            catalog=HIGH_CATALOG, tree_results=TREES_OK, policy=policy())
        blob = json.dumps(res, sort_keys=True, default=str)
        if baseline is None:
            baseline = blob
        assert blob == baseline


# ------------------------------------------------------------- §12.4 row invariants


def all_fixture_results(overrides):
    cases = [
        ("F_bam", {"g1", "g2", "g3", "g4"},
         [call("mcry_curated", f"g{i}", f"Bam{n}") for i, n in enumerate((1, 3, 5, 7), 1)],
         None, None),
        ("F_conflict", {"g1", "g2"},
         [call("mcry_curated", "g1", "Fbp1"), call("ath_diamond", "g2", "Sbp1")],
         None, None),
        ("F_part", {"g1", "g2", "g3", "g4"}, MULTI_CALLS, HIGH_CATALOG, TREES_OK),
        ("F_none", {"g1"}, [], None, None),
        ("F_afdb", {"g1"}, [call("afdb_swissprot", "g1", "Bam1", calibrated=False)], None, None),
    ]
    for fam, members, calls, catalog, trees in cases:
        yield adjudicate_family(fam, members, calls, overrides=overrides,
                                catalog=catalog, tree_results=trees, policy=policy())


def test_the_flags_and_verdicts_never_disagree(overrides):
    for res in all_fixture_results(overrides):
        for row in res["rows"]:
            needs_review = row["transfer_tier"] == "needs_review"
            assert needs_review is row["verdict"].startswith("NEEDS_REVIEW")
            auto = row["transfer_tier"] in ("stem_auto", "suffix_needs_tree_gate")
            assert auto is (row["verdict"] in AUTO_VERDICTS)
            assert (row["reason_code"] != "") is needs_review


def test_needs_review_rows_carry_reason_and_evidence(overrides):
    for res in all_fixture_results(overrides):
        for row in res["rows"]:
            if row["transfer_tier"] == "needs_review":
                assert row["reason_code"] != ""
                assert row["evidence_id"] != ""


def test_suffix_gate_rows_have_a_nonempty_suffix(overrides):
    for res in all_fixture_results(overrides):
        for row in res["rows"]:
            if row["transfer_tier"] == "suffix_needs_tree_gate":
                assert row["suffix"] != ""


def test_stem_auto_rows_have_positive_support(overrides):
    for res in all_fixture_results(overrides):
        for row in res["rows"]:
            if row["transfer_tier"] == "stem_auto":
                assert row["support"] > 0.0


def test_support_and_coverage_are_bounded(overrides):
    for res in all_fixture_results(overrides):
        assert 0.0 <= res["support"] <= 1.0
        assert 0.0 <= res["coverage"] <= 1.0
        assert res["coverage"] <= res["coverage_all_sources"] <= 1.0


# --------------------------------------------------- §12.5 measured known-answer set


CONFLICT_SET = [
    (["Tpt1", "Tpt2"], "SINGLE_STEM_MULTI_SUFFIX", {"Tpt"}),
    (["Bam1", "Bam3", "Bam5", "Bam7"], "SINGLE_STEM_MULTI_SUFFIX", {"Bam"}),
    (["Vha-a1", "Vha-a2", "Vha-a3"], "SINGLE_STEM_MULTI_SUFFIX", {"Vha-a"}),
    (["Pth1;1", "Pth1;2", "Pth1;3", "Pth1;4"], "SINGLE_STEM_MULTI_SUFFIX", {"Pth"}),
    (["CKB1;1", "CKB1;2"], "SINGLE_STEM_MULTI_SUFFIX", {"CKB"}),
    (["DiT1", "DiT2"], "SINGLE_STEM_MULTI_SUFFIX", {"DiT"}),
    (["PHYC", "PHYD"], "SINGLE_STEM_MULTI_SUFFIX", {"PHY"}),
    (["CRY2", "Cry 2"], None, {"CRY", "Cry"}),  # single stem either way
    (["Cipk10 (Sip1)", "Cipk11 (Sip4)"], "SINGLE_STEM_MULTI_SUFFIX", {"Cipk"}),
    (["Fba2;1", "Fba3"], "SINGLE_STEM_MULTI_SUFFIX", {"Fba"}),
    (["Ppck1a", "Ppck1b"], "SINGLE_STEM_MULTI_SUFFIX", {"Ppck"}),
    (["Almt7", "Almt2;2"], "SINGLE_STEM_MULTI_SUFFIX", {"Almt"}),
    (["Nhx2", "Acl"], "NEEDS_REVIEW_MULTI_STEM", None),
    (["Cipk1", "Snrk2.3", "Cbl3;3"], "NEEDS_REVIEW_MULTI_STEM", None),
    (["ZFN1", "JMJ14"], "NEEDS_REVIEW_MULTI_STEM", None),
    (["Fbp1", "Sbp1"], "NEEDS_REVIEW_MULTI_STEM", None),
]


@pytest.mark.parametrize("symbols,expected_verdict,expected_stems", CONFLICT_SET)
def test_measured_conflict_set(overrides, symbols, expected_verdict, expected_stems):
    members = {f"g{i}" for i in range(1, len(symbols) + 1)}
    calls = [call("mcry_curated", f"g{i}", s) for i, s in enumerate(symbols, 1)]
    res = adjudicate_family("F1", members, calls, overrides=overrides,
                            policy=policy(min_label_coverage=0.0))
    if expected_verdict is not None:
        assert res["verdict"] == expected_verdict
    if expected_stems is not None:
        fam_stems = {r["stem"] for r in res["rows"]
                     if r["transfer_tier"] == "stem_auto" and r["stem"]}
        assert fam_stems <= expected_stems and fam_stems
