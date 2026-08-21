"""steps/subfunctionalization.py — evidence integration + verdict + narrative
for HOW a subfamily diverged (sub- vs neo-functionalization vs conserved).

Fixtures mirror the real HyPhy JSON shapes (relax_sf3.json etc. on gpu).
No external tools invoked (repo convention).
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from steps.subfunctionalization import (
    classify,
    narrative,
    parse_absrel,
    parse_meme,
    parse_relax,
)

RELAX_JSON = {
    "test results": {
        "LRT": 68.297,
        "p-value": 1.11e-16,
        "relaxation or intensification parameter": 0.6362,
    },
}

MEME_JSON = {
    "MLE": {
        "headers": [["&alpha;", ""], ["&beta;<sup>1</sup>", ""], ["p<sup>1</sup>", ""],
                    ["&beta;<sup>+</sup>", ""], ["p<sup>+</sup>", ""], ["LRT", ""],
                    ["p-value", ""], ["# branches under selection", ""]],
        "content": {"0": [
            [0.1, 0, 1, 0, 0, 0, 1.0, 0],      # site 1: p=1.0
            [0.2, 0, 1, 5, 0.1, 9.0, 0.03, 2], # site 2: p=0.03 (episodic)
            [0.3, 0, 1, 0, 0, 0, 0.5, 0],      # site 3
        ]},
    },
}

ABSREL_JSON = {
    "branch attributes": {"0": {
        "Tfru_tip": {"Corrected P-value": 6.0e-13},
        "Node12": {"Corrected P-value": 1.0},
        "Mcry_tip": {"Corrected P-value": 0.8},
    }},
}


def _dump(tmp_path, name, obj):
    p = tmp_path / name
    p.write_text(json.dumps(obj))
    return p


# --- parsers -------------------------------------------------------------

def test_parse_relax(tmp_path):
    r = parse_relax(_dump(tmp_path, "relax.json", RELAX_JSON))
    assert abs(r["k"] - 0.6362) < 1e-9
    assert r["p"] < 1e-15
    assert r["direction"] == "relaxed"        # K < 1


def test_parse_relax_intensified(tmp_path):
    j = {"test results": {"LRT": 10, "p-value": 0.001,
                          "relaxation or intensification parameter": 2.5}}
    r = parse_relax(_dump(tmp_path, "relax.json", j))
    assert r["direction"] == "intensified"    # K > 1


def test_parse_meme_sites(tmp_path):
    sites = parse_meme(_dump(tmp_path, "meme.json", MEME_JSON), max_p=0.05)
    assert sites == [(2, 0.03)]               # 1-based site, p


def test_parse_absrel(tmp_path):
    br = parse_absrel(_dump(tmp_path, "absrel.json", ABSREL_JSON), max_p=0.05)
    assert br == [("Tfru_tip", 6.0e-13)]


# --- classification ------------------------------------------------------

def _sf3_evidence():
    """The real PEPC SF3 evidence pattern."""
    return {
        "relax": {"k": 0.636, "p": 1.1e-16, "direction": "relaxed"},
        "meme_sites_in_region": 0,   # episodic sites inside the signal region
        "meme_sites_total": 2,
        "absrel_stem": [],           # stem/internal branches: none positive
        "absrel_terminal": [("Tfru_tip", 6e-13)],
        "expression_share": 0.74,    # focal subfamily's share of family expression
        "signal_partition": "SF3-specific N-terminal region (col 165-209), 7/10 NLS",
        "retargeting_events": 0,
    }


def test_classify_subfunctionalization():
    v = classify(_sf3_evidence())
    assert v["verdict"] == "subfunctionalization"
    # relaxation + no stem positive selection + partitioned expression/signal
    assert "relaxed" in ",".join(v["evidence_for"])


def test_classify_neofunctionalization_when_stem_positive():
    ev = _sf3_evidence()
    ev["absrel_stem"] = [("Node5", 1e-6)]
    ev["retargeting_events"] = 2
    v = classify(ev)
    assert v["verdict"] == "neofunctionalization"


def test_classify_conserved_when_no_signal():
    ev = {
        "relax": {"k": 1.02, "p": 0.4, "direction": "relaxed"},
        "meme_sites_in_region": 0, "meme_sites_total": 0,
        "absrel_stem": [], "absrel_terminal": [],
        "expression_share": None, "signal_partition": "",
        "retargeting_events": 0,
    }
    assert classify(ev)["verdict"] == "conserved/inconclusive"


def test_classify_insufficient_without_selection_tests():
    assert classify({})["verdict"] == "insufficient evidence"


# --- narrative -----------------------------------------------------------

def test_narrative_contains_all_evidence_axes():
    ev = _sf3_evidence()
    text = narrative("PEPC clan", "SF3", classify(ev), ev)
    for token in ("K = 0.636", "1.1e-16", "74%", "165-209", "Tfru_tip",
                  "subfunctionalization"):
        assert token in text, token


def test_narrative_explains_the_mechanism_not_just_the_verdict():
    ev = _sf3_evidence()
    text = narrative("PEPC clan", "SF3", classify(ev), ev)
    # the HOW: relaxation under partitioned function, not positive selection
    assert "relaxed" in text.lower()
    assert "positive selection" in text.lower()


def test_apply_branch_names():
    from steps.subfunctionalization import apply_branch_names
    renamed = apply_branch_names(
        [("g101", 6e-13), ("Node12", 0.01)],
        {"g101": "Tfru_contig_062_000131"},
    )
    assert renamed == [("Tfru_contig_062_000131", 6e-13), ("Node12", 0.01)]


# --- branch-site axis ----------------------------------------------------

def test_classify_branchsite_significant_but_region_clean_is_still_subfunc():
    """Clade-wide branch-site significance driven by a TERMINAL branch does
    not make the split neofunctionalization — the stem is what defines it."""
    ev = _sf3_evidence()
    ev["branchsite"] = {"lrt": 15.363, "p": 8.87e-05,
                        "beb_sites_total": 14, "beb_sites_in_region": 0}
    v = classify(ev)
    assert v["verdict"] == "subfunctionalization"
    assert any("branch-site" in e for e in v["evidence_for"] + v["evidence_against"])


def test_narrative_reports_branchsite_and_its_localization():
    ev = _sf3_evidence()
    ev["branchsite"] = {"lrt": 15.363, "p": 8.87e-05,
                        "beb_sites_total": 14, "beb_sites_in_region": 0}
    text = narrative("PEPC clan", "SF3", classify(ev), ev)
    assert "8.87e-05" in text
    assert "14" in text
    # the region stays clean on BOTH site-level methods
    assert "MEME" in text and "BEB" in text


# --- structural disorder axis (pLDDT in the signal region) ---------------

def test_classify_uses_disorder_as_support():
    """A subfamily-specific region that is structurally DISORDERED supports
    divergence under relaxed constraint rather than adaptive remodelling."""
    ev = _sf3_evidence()
    ev["disorder"] = {"delta": -0.304, "n_focal": 9, "n_other": 76,
                      "p": 1.95e-05, "all_below": True}
    v = classify(ev)
    assert v["verdict"] == "subfunctionalization"
    assert any("disorder" in e.lower() or "pLDDT" in e for e in v["evidence_for"])


def test_disorder_without_significance_is_not_counted():
    """An unremarkable region must not be dressed up as evidence."""
    ev = _sf3_evidence()
    ev["disorder"] = {"delta": -0.02, "n_focal": 9, "n_other": 76,
                      "p": 0.6, "all_below": False}
    v = classify(ev)
    assert not any("pLDDT" in e for e in v["evidence_for"])


def test_narrative_reports_the_disorder_control():
    """The control (non-focal members at the same columns) is what makes the
    number meaningful — the narrative must carry it, not just the delta."""
    ev = _sf3_evidence()
    ev["disorder"] = {"delta": -0.304, "n_focal": 9, "n_other": 76,
                      "p": 1.95e-05, "all_below": True}
    text = narrative("PEPC clan", "SF3", classify(ev), ev)
    assert "76" in text and "9" in text
    assert "1.95e-05" in text or "1.9e-05" in text
    assert "disorder" in text.lower()
