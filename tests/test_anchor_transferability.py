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
