"""steps/subfamily.anchor_transferability — can a reference species' subfamily
names be transferred to the query species at all?

The Musa-MYB style method (PLOS ONE 10.1371/journal.pone.0239275) classifies
query genes by clade membership with Arabidopsis anchors. That is sound only
when the reference subfamilies are ANCIENT — older than the split with the
query lineage. When two reference labels are each other's nearest relatives
with no query gene between them, they are a reference-lineage-specific
duplication and their names name nothing in the query species.

Measured motivation: Arabidopsis PPC1 (AT1G53310) and PPC3 (AT3G14940) are
sisters at support 100 in the PEPC clan tree, so "PPC1"/"PPC3" cannot label
Caryophyllales subfamilies — and naming by symbol duly produced six different
subfamilies all called PPC1. MYB subfamilies do not have this problem.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from steps.subfamily import anchor_transferability


def test_reference_specific_duplication_is_flagged():
    """(A1,A3) sisters, no query leaf between them -> names not transferable."""
    tree = "(((A1,A3),(q1,q2)),out);"
    labels = {"A1": "PPC1", "A3": "PPC3"}
    res = anchor_transferability(tree, labels)
    by = {r["label"]: r for r in res}
    assert by["PPC1"]["transferable"] is False
    assert by["PPC3"]["transferable"] is False
    assert "PPC3" in by["PPC1"]["blocked_by"]
    assert "PPC1" in by["PPC3"]["blocked_by"]


def test_ancient_subfamilies_are_transferable():
    """Each anchor sits with query genes -> the name designates a real group."""
    tree = "(((A_S1,(q1,q2)),(A_S4,(q3,q4))),out);"
    labels = {"A_S1": "S1", "A_S4": "S4"}
    res = anchor_transferability(tree, labels)
    by = {r["label"]: r for r in res}
    assert by["S1"]["transferable"] is True
    assert by["S4"]["transferable"] is True
    assert by["S1"]["n_query_in_clade"] == 2


def test_query_only_clade_is_reported_as_lineage_specific():
    """Musa-only clades: the paper called these lineage-specific expansions.
    They carry no anchor, so no reference name can apply."""
    tree = "(((A_S1,q1),(q2,q3)),out);"
    labels = {"A_S1": "S1"}
    res = anchor_transferability(tree, labels)
    unlabelled = [r for r in res if r["label"] is None]
    assert unlabelled, "a query-only clade must be reported"
    assert unlabelled[0]["n_query_in_clade"] >= 2


def test_single_anchor_label_with_no_query_neighbour_is_not_transferable():
    tree = "((A1,A2),(q1,q2));"
    labels = {"A1": "X", "A2": "Y"}
    by = {r["label"]: r for r in anchor_transferability(tree, labels)}
    assert by["X"]["transferable"] is False


def test_multiple_anchors_sharing_one_label_are_one_group():
    """Two reference genes both called S1 are not a 'duplication' problem —
    they are one named group."""
    tree = "(((A1a,A1b),(q1,q2)),(A_S4,q3));"
    labels = {"A1a": "S1", "A1b": "S1", "A_S4": "S4"}
    by = {r["label"]: r for r in anchor_transferability(tree, labels)}
    assert by["S1"]["transferable"] is True
    assert by["S1"]["blocked_by"] == []


def test_query_prefixes_restrict_what_counts_as_query():
    """With explicit query prefixes, an unlabelled REFERENCE leaf must not be
    mistaken for a query gene (that would hide a duplication)."""
    tree = "(((A1,A3),(ATH_other,q1)),out);"
    labels = {"A1": "PPC1", "A3": "PPC3"}
    by = {r["label"]: r
          for r in anchor_transferability(tree, labels, query_prefixes=["q"])}
    assert by["PPC1"]["transferable"] is False


def test_verdict_strings_are_actionable():
    tree = "(((A1,A3),(q1,q2)),out);"
    by = {r["label"]: r
          for r in anchor_transferability(tree, {"A1": "PPC1", "A3": "PPC3"})}
    assert "duplication" in by["PPC1"]["verdict"].lower()


# --- the subfamily is the MAXIMAL exclusive clade, not the nearest neighbour --

def test_returns_the_subfamily_not_the_nearest_neighbour():
    """first_informative_clade() stops at whatever joins first, which is the
    anchor's nearest relatives — not its subfamily. Measured consequence on the
    rebuilt PEPC tree: it reported Ppc2 as a 10-leaf all-Amaranthaceae clade and
    therefore said 'no query gene' the moment queries were restricted to cacti,
    even though the real 17-leaf subfamily does contain cactus genes."""
    #        (((P1,q1),q2), ((P2,r1),r2))
    tree = "((((P1,q1),q2),((P2,r1),r2)),out);"
    labels = {"P1": "Ppc1", "P2": "Ppc2"}
    by = {r["label"]: r for r in anchor_transferability(tree, labels)}
    # Ppc1's subfamily is everything on its side that excludes Ppc2: P1,q1,q2
    assert by["Ppc1"]["clade_size"] == 3
    assert by["Ppc2"]["clade_size"] == 3
    assert by["Ppc1"]["n_query_in_clade"] == 2


def test_subfamilies_are_disjoint():
    tree = "((((P1,q1),q2),((P2,r1),r2)),out);"
    rows = anchor_transferability(tree, {"P1": "Ppc1", "P2": "Ppc2"})
    m = {r["label"]: set(r["members"]) for r in rows if r["label"]}
    assert not (m["Ppc1"] & m["Ppc2"])


def test_support_threshold_rejects_a_weak_clade():
    """A clade below the support bar must not be reported as transferable."""
    tree = "((((P1,q1)30,q2)25,((P2,r1)99,r2)98)90,out);"
    by = {r["label"]: r
          for r in anchor_transferability(tree, {"P1": "Ppc1", "P2": "Ppc2"},
                                          min_support=80)}
    assert by["Ppc1"]["transferable"] is False
    assert "support" in by["Ppc1"]["verdict"].lower()
    assert by["Ppc2"]["transferable"] is True


def test_support_is_reported_even_without_a_threshold():
    tree = "((((P1,q1)30,q2)25,((P2,r1)99,r2)98)90,out);"
    by = {r["label"]: r
          for r in anchor_transferability(tree, {"P1": "Ppc1", "P2": "Ppc2"})}
    assert by["Ppc1"]["support"] == "25"


def test_slash_separated_support_uses_the_weaker_half():
    """IQ-TREE writes 'SH-aLRT/UFBoot'. Both must clear the bar."""
    tree = "((((P1,q1)96.2/94,q2)86.2/73,((P2,r1)99/99,r2)98/98)90/90,out);"
    by = {r["label"]: r
          for r in anchor_transferability(tree, {"P1": "Ppc1", "P2": "Ppc2"},
                                          min_support=80)}
    assert by["Ppc1"]["transferable"] is False    # 73 < 80


# --- the grade: a boolean cannot say 100/76 apart from 100/98 ---------------

# The two-anchor shape used below: Ppc1's subfamily is (P1,q1),q2 and Ppc2's is
# (P2,r1),r2, so the support written on the two inner-most-but-one nodes is the
# support of each subfamily.
def _graded(sh_uf_1, sh_uf_2="100/99", **kw):
    tree = (f"((((P1,q1)100/100,q2){sh_uf_1},"
            f"((P2,r1)100/100,r2){sh_uf_2})90/90,out);")
    return {r["label"]: r
            for r in anchor_transferability(
                tree, {"P1": "Ppc1", "P2": "Ppc2"}, **kw)}


def _agreeing(*names):
    """Cross-tree evidence in which every dataset recovers the same members."""
    members = {"Ppc1": ["P1", "q1", "q2"], "Ppc2": ["P2", "r1", "r2"]}
    return {n: members for n in names}


def test_high_needs_both_supports_and_recovered_membership():
    by = _graded("100/98", cross_tree_members=_agreeing("codon12", "aa"))
    assert by["Ppc1"]["grade"] == "HIGH"
    assert by["Ppc1"]["n_consistent_datasets"] == 3


def test_the_weaker_support_gates_the_grade():
    """100/76 is not the same claim as 100/98 — UFboot decides, not SH-aLRT."""
    by = _graded("100/76", cross_tree_members=_agreeing("codon12", "aa"))
    assert by["Ppc1"]["grade"] == "PROVISIONAL"
    assert "76" in by["Ppc1"]["grade_reason"]


def test_the_average_of_the_two_supports_must_not_decide():
    """mean(100, 90) = 95 would clear the HIGH bar; the weaker one does not."""
    by = _graded("100/90", cross_tree_members=_agreeing("codon12", "aa"))
    assert by["Ppc1"]["grade"] == "PROVISIONAL"


def test_sh_alrt_below_its_own_bar_is_unresolved():
    by = _graded("70/99", cross_tree_members=_agreeing("codon12", "aa"))
    assert by["Ppc1"]["grade"] == "UNRESOLVED"
    assert "SH-aLRT" in by["Ppc1"]["grade_reason"]


def test_ufboot_below_the_provisional_bar_is_unresolved():
    by = _graded("100/60", cross_tree_members=_agreeing("codon12", "aa"))
    assert by["Ppc1"]["grade"] == "UNRESOLVED"


def test_unevaluated_cross_tree_consistency_cannot_reach_high():
    """No independent dataset was supplied. That is not agreement, and the
    reason has to say so — 'we did not check' is not 'it checked out'."""
    by = _graded("100/98")
    assert by["Ppc1"]["grade"] == "PROVISIONAL"
    assert by["Ppc1"]["cross_tree_evaluated"] is False
    assert by["Ppc1"]["n_consistent_datasets"] is None
    reason = by["Ppc1"]["grade_reason"].lower()
    assert "not evaluated" in reason
    assert "disagree" not in reason and "failed" not in reason


def test_cross_tree_disagreement_is_a_measured_failure():
    """A dataset that recovered different members is evidence AGAINST, and is
    reported differently from nobody having looked."""
    other = {"codon12": {"Ppc1": ["P1", "q1"], "Ppc2": ["P2", "r1", "r2"]}}
    by = _graded("100/98", cross_tree_members=other)
    assert by["Ppc1"]["grade"] == "UNRESOLVED"
    assert by["Ppc1"]["cross_tree_evaluated"] is True
    assert by["Ppc1"]["inconsistent_datasets"] == ["codon12"]
    assert "codon12" in by["Ppc1"]["grade_reason"]
    assert by["Ppc2"]["grade"] == "HIGH"          # this one did agree


def test_the_agreeing_datasets_are_named():
    by = _graded("100/98", tree_name="codon123",
                 cross_tree_members=_agreeing("codon12", "aa102"))
    assert by["Ppc1"]["consistent_datasets"] == ["codon123", "codon12", "aa102"]


def test_a_missing_label_in_another_dataset_is_not_agreement():
    by = _graded("100/98",
                 cross_tree_members={"aa": {"Ppc2": ["P2", "r1", "r2"]}})
    assert by["Ppc1"]["grade"] == "UNRESOLVED"
    assert by["Ppc1"]["inconsistent_datasets"] == ["aa"]


def test_an_unsupported_clade_cannot_be_graded():
    tree = "((((P1,q1),q2),((P2,r1)99/99,r2)98/98)90/90,out);"
    by = {r["label"]: r
          for r in anchor_transferability(tree, {"P1": "Ppc1", "P2": "Ppc2"})}
    assert by["Ppc1"]["grade"] == "UNRESOLVED"
    assert "no support" in by["Ppc1"]["grade_reason"].lower()


def test_a_blocked_label_is_never_graded_above_unresolved():
    """Reference-lineage-specific duplication: there is no assignment to grade,
    however strong the node holding the two anchors is."""
    tree = ("((((A1,x1)100/100,(A3,x2)100/100)100/100,"
            "(q1,q2)100/100)100/100,out);")
    by = {r["label"]: r for r in anchor_transferability(
        tree, {"A1": "PPC1", "A3": "PPC3"}, query_prefixes=["q"],
        cross_tree_members={"aa": {"PPC1": ["A1", "x1"],
                                   "PPC3": ["A3", "x2"]}})}
    assert by["PPC1"]["blocked_by"] == ["PPC3"]       # and its clade is 100/100
    assert by["PPC1"]["grade"] == "UNRESOLVED"


def test_query_only_clades_carry_a_grade_column_too():
    tree = "(((A_S1,q1),(q2,q3)),out);"
    rows = anchor_transferability(tree, {"A_S1": "S1"})
    assert all(r["grade"] == "UNRESOLVED" for r in rows if r["label"] is None)


def test_the_grade_does_not_change_the_existing_boolean():
    """Callers that only read `transferable` must see what they saw before."""
    by = _graded("100/76")
    assert by["Ppc1"]["transferable"] is True     # unchanged, no min_support
    assert by["Ppc1"]["grade"] == "PROVISIONAL"


def test_the_thresholds_are_named_constants():
    from steps.subfamily import (GRADE_HIGH_MIN_UFBOOT, GRADE_MIN_SH_ALRT,
                                 GRADE_MIN_CONSISTENT_DATASETS,
                                 GRADE_PROVISIONAL_MIN_UFBOOT)
    assert (GRADE_MIN_SH_ALRT, GRADE_HIGH_MIN_UFBOOT,
            GRADE_PROVISIONAL_MIN_UFBOOT, GRADE_MIN_CONSISTENT_DATASETS) \
        == (80, 95, 70, 2)


# --- known answer on the real corrected PEPC tree --------------------------

REAL_TREE = (Path(__file__).resolve().parent / "data"
             / "clan_corrected_tree.treefile")


def test_corrected_pepc_subfamilies_are_both_provisional():
    """`mM_repair/clan_corrected_tree.treefile`, the tree #40 was decided on.

    ppc-1E1 is 95/91 and ppc-1E2 is 100/76: one of them looks far stronger
    than the other and neither reaches UFboot 95, so both are PROVISIONAL.
    This is the outcome to preserve — the bar is not to be moved until
    ppc-1E1 clears it.
    """
    rows = anchor_transferability(
        REAL_TREE.read_text(),
        {"Mcry_Mcr8G11630": "ppc-1E1", "Mcry_Mcr7G08600": "ppc-1E2"},
    )
    by = {r["label"]: r for r in rows if r["label"]}
    assert by["ppc-1E1"]["support"] == "95/91"
    assert by["ppc-1E2"]["support"] == "100/76"
    assert by["ppc-1E1"]["grade"] == "PROVISIONAL"
    assert by["ppc-1E2"]["grade"] == "PROVISIONAL"
    # and the grade is not the boolean in disguise
    assert by["ppc-1E1"]["transferable"] is True
    assert by["ppc-1E2"]["transferable"] is True
    # nothing here was checked against a second dataset, and the reason says
    # so rather than counting the missing check as agreement
    assert by["ppc-1E1"]["cross_tree_evaluated"] is False
    assert "not evaluated" in by["ppc-1E1"]["grade_reason"].lower()
    # ...and had all three #40 datasets agreed, these supports would still not
    # reach HIGH. That is what the bar being 95 means.
    from steps.subfamily import grade_assignment
    for label in ("ppc-1E1", "ppc-1E2"):
        assert grade_assignment(by[label]["sh_alrt"], by[label]["ufboot"],
                                3)[0] == "PROVISIONAL"
