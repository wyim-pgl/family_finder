"""steps/subfunctionalization.narrative() — which paragraphs get written.

narrative() was one 93->131-line linear string assembly; it is now a
heading, one body per verdict, and a trailing ledger, with each numbered
point its own builder. These lock what that extraction had to preserve:
the exact set of paragraphs emitted for a given evidence dict, their
order, and the conditions under which each is withheld.

Withholding is the part that matters scientifically. A point that is
silently dropped reads to a human as an absent signal rather than as an
unmeasured one, so every omission below is a deliberate rule, not a
formatting detail.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from steps.subfunctionalization import narrative


RELAX = {"k": 0.312, "p": 0.0041, "lrt": 8.24, "direction": "relaxed"}
SIGNAL = "cols 812-847"
BRANCHSITE = {"lrt": 12.4, "p": 0.00043,
              "beb_sites_total": 7, "beb_sites_in_region": 0}
DISORDER = {"p": 0.003, "delta": -0.1421, "n_focal": 4, "n_other": 17}


def full_evidence(**over):
    ev = {"relax": RELAX, "signal_partition": SIGNAL,
          "meme_sites_in_region": 0, "meme_sites_total": 11,
          "expression_share": 0.6667, "branchsite": BRANCHSITE,
          "disorder": DISORDER,
          "absrel_stem": [], "absrel_terminal": [("Cgig_g1", 0.0421)]}
    ev.update(over)
    return ev


def sub_verdict(coords=True, **over):
    v = {"verdict": "subfunctionalization",
         "evidence_for": ["relaxed selection (K = 0.312, p = 0.0041)"],
         "evidence_against": [],
         "coordinates_verified": coords}
    v.update(over)
    return v


def points(text):
    """The numbered-point prefixes present, in the order they appear."""
    return [ln.split(".")[0] for ln in text.splitlines()
            if ln[:2] in ("1.", "2.", "3.", "4.", "5.")]


# --- heading and ledger --------------------------------------------------

def test_heading_names_the_family_the_subfamily_and_the_call():
    text = narrative("OG42", "SF3", sub_verdict(), full_evidence())

    assert text.splitlines()[0] == "## OG42 / SF3: subfunctionalization"
    assert text.splitlines()[1] == ""


def test_every_verdict_ends_with_the_evidence_ledger():
    for call in ("subfunctionalization", "neofunctionalization",
                 "conserved/inconclusive", "insufficient evidence"):
        verdict = sub_verdict(verdict=call, evidence_for=["A"],
                              evidence_against=["B"])

        text = narrative("OG42", "SF3", verdict, full_evidence())

        assert text.endswith("\nEvidence summary:\n"
                             "- [supports] A\n"
                             "- [context]  B")


def test_evidence_that_could_not_be_judged_is_listed_not_dropped():
    verdict = sub_verdict(evidence_for=[], evidence_against=[],
                          cannot_judge=["MEME sites inside the signal region"])

    text = narrative("OG42", "SF3", verdict, full_evidence())

    assert text.endswith(
        "- [cannot judge] MEME sites inside the signal region")


# --- the subfunctionalization body ---------------------------------------

def test_complete_evidence_emits_all_five_points_in_order():
    text = narrative("OG42", "SF3", sub_verdict(), full_evidence())

    assert points(text) == ["1", "2", "3", "4", "5"]


def test_a_missing_axis_is_skipped_and_the_others_keep_their_numbers():
    # Arrange: no RELAX run, no expression table
    ev = full_evidence(relax=None, expression_share=None)

    # Act
    text = narrative("OG42", "SF3", sub_verdict(), ev)

    # Assert: the numbers label the axis, not the position in the list
    assert points(text) == ["2", "4", "5"]


def test_the_closing_terminal_branch_paragraph_needs_a_clean_stem():
    with_stem = narrative("OG42", "SF3", sub_verdict(),
                          full_evidence(absrel_stem=[("Node7", 0.001)]))
    clean_stem = narrative("OG42", "SF3", sub_verdict(), full_evidence())

    assert "The only positive selection detected is on terminal" not in with_stem
    assert "The only positive selection detected is on terminal" in clean_stem
    assert "Cgig_g1, p = 0.042" in clean_stem


def test_terminal_paragraph_is_withheld_when_no_branch_is_positive():
    text = narrative("OG42", "SF3", sub_verdict(),
                     full_evidence(absrel_terminal=[]))

    assert "The only positive selection detected" not in text


# --- point 2: MEME, and the coordinate gate ------------------------------

def test_a_clean_region_is_only_reported_when_the_coordinates_were_verified():
    verified = narrative("OG42", "SF3", sub_verdict(coords=True),
                         full_evidence())
    unverified = narrative("OG42", "SF3", sub_verdict(coords=False),
                           full_evidence())

    assert "2. **No positive selection where it would matter.**" in verified
    assert "2. **No positive selection where it would matter.**" not in unverified


def test_an_absent_coordinate_flag_counts_as_unverified():
    verdict = sub_verdict()
    del verdict["coordinates_verified"]

    text = narrative("OG42", "SF3", verdict, full_evidence())

    assert "2. **No positive selection" not in text


def test_no_signal_region_means_nothing_to_say_about_sites_inside_it():
    text = narrative("OG42", "SF3", sub_verdict(),
                     full_evidence(signal_partition=""))

    assert "2. **No positive selection" not in text


def test_a_region_with_sites_in_it_is_not_reported_as_clean():
    text = narrative("OG42", "SF3", sub_verdict(),
                     full_evidence(meme_sites_in_region=3))

    assert "2. **No positive selection" not in text


def test_the_family_wide_count_is_mentioned_only_when_there_is_one():
    with_total = narrative("OG42", "SF3", sub_verdict(), full_evidence())
    no_total = narrative("OG42", "SF3", sub_verdict(),
                         full_evidence(meme_sites_total=0))

    assert "11 episodic sites exist family-wide" in with_total
    assert "episodic sites exist family-wide" not in no_total


# --- point 3: expression -------------------------------------------------

def test_a_zero_expression_share_is_reported_rather_than_treated_as_missing():
    # 0.0 is a measurement; None is the absence of one. They must not collapse.
    measured = narrative("OG42", "SF3", sub_verdict(),
                         full_evidence(expression_share=0.0))
    absent = narrative("OG42", "SF3", sub_verdict(),
                       full_evidence(expression_share=None))

    assert "SF3 carries 0% of the family's expression" in measured
    assert "3. **The partition itself.**" not in absent


# --- point 4: branch-site ------------------------------------------------

def test_branch_site_localization_clause_flips_on_coordinate_verification():
    verified = narrative("OG42", "SF3", sub_verdict(coords=True),
                         full_evidence())
    unverified = narrative("OG42", "SF3", sub_verdict(coords=False),
                           full_evidence())

    # the test result itself survives either way
    for text in (verified, unverified):
        assert "LRT = 12.4, p = 0.00043" in text
        assert "BEB identifies 7 sites" in text
    assert "0 of those BEB sites fall" in verified
    assert "CANNOT BE JUDGED" not in verified
    assert "CANNOT BE JUDGED" in unverified
    assert "0 of those BEB sites fall" not in unverified


def test_branch_site_paragraph_is_absent_without_a_branch_site_run():
    text = narrative("OG42", "SF3", sub_verdict(), full_evidence(branchsite=None))

    assert "4. **What the branch-site test adds" not in text


# --- point 5: disorder ---------------------------------------------------

def test_disorder_is_reported_only_when_significant_and_downward():
    cases = {
        "significant drop": (DISORDER, True),
        "drop but p above alpha": ({**DISORDER, "p": 0.42}, False),
        "significant but upward": ({**DISORDER, "delta": 0.09}, False),
        "no p value": ({**DISORDER, "p": None}, False),
    }
    for label, (dis, expected) in cases.items():
        text = narrative("OG42", "SF3", sub_verdict(), full_evidence(disorder=dis))

        assert ("5. **The region is disordered" in text) is expected, label


def test_the_disorder_paragraph_names_the_control_it_rests_on():
    text = narrative("OG42", "SF3", sub_verdict(), full_evidence())

    assert "the same alignment columns in 17 non-focal" in text
    assert "Mann-Whitney p = 0.003" in text


# --- the other two bodies ------------------------------------------------

def test_neofunctionalization_body_bullets_the_supporting_evidence():
    verdict = {"verdict": "neofunctionalization",
               "evidence_for": ["stem positive selection: Node7 (p = 0.0012)",
                                "retargeting events on the family tree: 2"],
               "evidence_against": []}

    text = narrative("OG42", "SF3", verdict, full_evidence())

    assert "SF3 shows a neofunctionalization signature: a NEW state" in text
    assert "- stem positive selection: Node7 (p = 0.0012)" in text
    assert "1. **Selection regime.**" not in text      # no numbered case here


def test_an_undecided_verdict_labels_both_sides_of_the_evidence():
    verdict = {"verdict": "conserved/inconclusive",
               "evidence_for": ["relaxed selection"],
               "evidence_against": ["terminal-branch positive selection only"]}

    text = narrative("OG42", "SF3", verdict, full_evidence())

    assert ("The available evidence does not distinguish a functional shift "
            "(conserved/inconclusive)." in text)
    assert "- for: relaxed selection" in text
    assert "- against: terminal-branch positive selection only" in text


def test_insufficient_evidence_produces_a_report_rather_than_an_error():
    verdict = {"verdict": "insufficient evidence",
               "evidence_for": [], "evidence_against": []}

    text = narrative("OG42", "SF3", verdict, {})

    assert text == ("## OG42 / SF3: insufficient evidence\n"
                    "\n"
                    "The available evidence does not distinguish a functional "
                    "shift (insufficient evidence).\n"
                    "\n"
                    "Evidence summary:")
