"""Merge fragments only when the tree says they are one family (issue #47).

The one-way fragmentation clusters over-merge by construction (nearest-neighbor
vote concentration), so cluster membership alone must never merge anything.
The tree decides, and it can say three things about a fragment:

  INTERLEAVED   its members mix with another fragment's - lineage-axis
                fragmentation, the PEPC signature. Merge.
  MONOPHYLETIC  its members form their own clade. Topology alone CANNOT
                distinguish "subfamily of the same family" from "distinct
                neighbouring family" - PPC4 is monophyletic inside the PEPC
                cluster and was its own family in the 5-species run. Undecided,
                with the metrics reported.
  (missing)     members absent from the tree are counted, never assumed.
"""
import sys
import types
from pathlib import Path

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from steps.cluster_validate import fragment_verdict, write_verdicts


def test_interleaved_fragments_form_one_merge_group():
    """A(a1,a2) and B(b1,b2) alternate: neither is a clade -> merge."""
    v = fragment_verdict("((a1,b1),(a2,b2));",
                         {"A": ["a1", "a2"], "B": ["b1", "b2"]})
    assert v["merge_groups"] == [["A", "B"]]
    assert v["fragments"]["A"]["status"] == "INTERLEAVED"


def test_a_monophyletic_fragment_is_undecided_not_kept_and_not_merged():
    """C sits as its own clade beside the mixed pair: PPC4's shape."""
    v = fragment_verdict("(((a1,b1),(a2,b2)),(c1,c2));",
                         {"A": ["a1", "a2"], "B": ["b1", "b2"],
                          "C": ["c1", "c2"]})
    assert v["merge_groups"] == [["A", "B"]]
    assert v["fragments"]["C"]["status"] == "MONOPHYLETIC"
    assert "topology alone" in v["fragments"]["C"]["reason"]


def test_all_monophyletic_means_no_merges_at_all():
    v = fragment_verdict("((a1,a2),(b1,b2));",
                         {"A": ["a1", "a2"], "B": ["b1", "b2"]})
    assert v["merge_groups"] == []


def test_three_way_interleaving_merges_all_three():
    """The PEPC lineage chain: no pair is clean, all mix."""
    v = fragment_verdict("((a1,(b1,c1)),(b2,(a2,c2)));",
                         {"A": ["a1", "a2"], "B": ["b1", "b2"],
                          "C": ["c1", "c2"]})
    assert v["merge_groups"] == [["A", "B", "C"]]


def test_missing_members_are_counted_not_assumed():
    v = fragment_verdict("((a1,b1),(a2,b2));",
                         {"A": ["a1", "a2", "a3_gone"], "B": ["b1", "b2"]})
    assert v["fragments"]["A"]["n_missing"] == 1
    assert v["n_tree_leaves_unclaimed"] == 0


def test_leaves_no_fragment_claims_are_counted():
    v = fragment_verdict("((a1,x1),(a2,b1));",
                         {"A": ["a1", "a2"], "B": ["b1"]})
    assert v["n_tree_leaves_unclaimed"] == 1


def test_a_single_fragment_cluster_is_trivially_monophyletic():
    v = fragment_verdict("((a1,a2),a3);", {"A": ["a1", "a2", "a3"]})
    assert v["merge_groups"] == []
    assert v["fragments"]["A"]["status"] == "MONOPHYLETIC"


def test_verdict_report_carries_the_cluster_and_the_decision(tmp_path):
    rows = [{"cluster_id": "C0297", "verdict": fragment_verdict(
        "((a1,b1),(a2,b2));", {"A": ["a1", "a2"], "B": ["b1", "b2"]})}]
    path = write_verdicts(rows, tmp_path)
    lines = path.read_text().splitlines()
    assert lines[0] == ("cluster_id\tfragment\tstatus\tn_members\tn_missing\t"
                       "merge_group\treason")
    a = [l for l in lines if l.startswith("C0297\tA")][0].split("\t")
    assert a[2] == "INTERLEAVED" and a[5] == "A+B"
