"""Merge fragments only when the tree says they are one family (issue #47).

The one-way fragmentation clusters over-merge by construction (nearest-neighbor
vote concentration), so cluster membership alone must never merge anything.
The tree decides, and it can say three things about a fragment:

  INTERLEAVED   its members do not form an edge-defined clade. Pairwise
                inseparable fragments are lineage-axis fragmentation, the PEPC
                signature, and merge. An isolated internal strip does not.
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


# ---------------------------------------------------------------------------
# unrooted trees (issue #51) - the verdict must not depend on where the
# Newick happens to be rooted
# ---------------------------------------------------------------------------

# Same topology, two root placements. FastTree emits unrooted trees and puts
# the root wherever it likes, so any verdict that differs between these two is
# an artifact of notation, not of the data.
_SAME_TOPOLOGY = [
    "((a1,a2),((b1,c1),(b2,c2)));",
    "(a1,(a2,((b1,c1),(b2,c2))));",
]
_THREE_FRAGS = {"A": ["a1", "a2"], "B": ["b1", "b2"], "C": ["c1", "c2"]}


def _statuses(newick, frags):
    v = fragment_verdict(newick, frags)
    return ({f: v["fragments"][f]["status"] for f in v["fragments"]},
            [sorted(g) for g in v["merge_groups"]])


def test_verdict_is_invariant_to_rerooting():
    # Arrange / Act
    first = _statuses(_SAME_TOPOLOGY[0], _THREE_FRAGS)
    second = _statuses(_SAME_TOPOLOGY[1], _THREE_FRAGS)

    # Assert
    assert first == second


def test_a_clade_split_across_the_root_is_still_monophyletic():
    # Arrange: A's members sit on either side of the root, but one EDGE still
    # separates {a1,a2} from everything else, which is what monophyly means in
    # an unrooted tree.
    frags = {"A": ["a1", "a2"], "B": ["b1", "b2"]}

    # Act
    statuses, groups = _statuses("(a1,(a2,(b1,b2)));", frags)

    # Assert
    assert statuses["A"] == "MONOPHYLETIC"
    assert groups == []


def test_genuinely_interleaved_stays_interleaved_under_rerooting():
    # Arrange: no edge separates A from B in either notation
    frags = {"A": ["a1", "a2"], "B": ["b1", "b2"]}

    # Act
    first = _statuses("((a1,b1),(a2,b2));", frags)
    second = _statuses("(a1,(b1,(a2,b2)));", frags)

    # Assert
    assert first[0]["A"] == first[0]["B"] == "INTERLEAVED"
    assert first == second


def test_all_members_present_is_monophyletic_regardless_of_root():
    # Arrange: a fragment holding every leaf is trivially its own clade -
    # the complement is empty, which must not be read as "mixes with nothing
    # in particular" and flipped to INTERLEAVED.
    frags = {"A": ["a1", "a2", "b1"]}

    # Act
    statuses, groups = _statuses("((a1,a2),b1);", frags)

    # Assert
    assert statuses["A"] == "MONOPHYLETIC"
    assert groups == []


def test_tied_smallest_sides_do_not_make_merge_groups_root_dependent():
    # Arrange: one unrooted caterpillar, printed with the root at two internal
    # nodes. B={e,f} occupies an internal strip.  The equally small sides
    # containing B are A+B and B+C+D; picking whichever is visited first used
    # to call B mixed with A in one notation and C,D in the other.
    trees = [
        "(a,b,(c,(d,(e,(f,(g,(h,(i,j))))))));",
        "((((a,b),c),d),e,(f,(g,(h,(i,j)))));",
    ]
    frags = {
        "A": ["a", "b", "c", "d"],
        "B": ["e", "f"],
        "C": ["g", "i"],
        "D": ["h", "j"],
    }

    # Act
    first = _statuses(trees[0], frags)
    second = _statuses(trees[1], frags)
    first_reasons = {
        fam: row["reason"]
        for fam, row in fragment_verdict(trees[0], frags)["fragments"].items()
    }
    second_reasons = {
        fam: row["reason"]
        for fam, row in fragment_verdict(trees[1], frags)["fragments"].items()
    }

    # Assert: C and D really interleave; an edge separates B from every other
    # individual fragment, so B supplies no evidence for a merge.
    assert first == second
    assert first_reasons == second_reasons
    assert first[1] == [["C", "D"]]


def test_an_internal_strip_does_not_name_a_separable_fragment_as_a_partner():
    # Arrange: B is not itself one side of an edge, but the edge immediately
    # to either side separates it from A and from C respectively. Non-clade is
    # therefore not enough to claim that B mixes with either neighbour.
    tree = "((a1,a2),(b1,(b2,(c1,c2))));"
    frags = {
        "A": ["a1", "a2"],
        "B": ["b1", "b2"],
        "C": ["c1", "c2"],
    }

    # Act
    verdict = fragment_verdict(tree, frags)

    # Assert
    assert verdict["fragments"]["B"]["status"] == "INTERLEAVED"
    assert (
        "every other fragment is separable"
        in verdict["fragments"]["B"]["reason"]
    )
    assert verdict["merge_groups"] == []


def test_polytomy_interleaving_is_invariant_to_putting_root_on_an_edge():
    # Arrange: in a four-way star, both two-leaf convex hulls meet at the
    # central node, so no edge separates A from B. The second Newick inserts a
    # degree-two printed root on a pendant edge of the same unrooted star.
    frags = {"A": ["a1", "a2"], "B": ["b1", "b2"]}
    at_node = "(a1,a2,b1,b2);"
    on_edge = "(a1,(a2,b1,b2));"

    # Act / Assert
    assert _statuses(at_node, frags) == _statuses(on_edge, frags)
    assert _statuses(at_node, frags)[1] == [["A", "B"]]


def test_unclaimed_leaf_can_break_a_clade_but_cannot_invent_a_partner():
    # Arrange: x lies inside B's smallest rooted-looking span, but belongs to
    # no fragment. A and B are still separated from one another by an edge.
    tree = "((a1,a2),(b1,(x,b2)));"
    frags = {"A": ["a1", "a2"], "B": ["b1", "b2"]}

    # Act
    verdict = fragment_verdict(tree, frags)

    # Assert: reporting B as a non-clade is intentional; silently turning the
    # unclaimed leaf into evidence to merge A+B would not be.
    assert verdict["fragments"]["A"]["status"] == "MONOPHYLETIC"
    assert verdict["fragments"]["B"]["status"] == "INTERLEAVED"
    assert (
        "every other fragment is separable"
        in verdict["fragments"]["B"]["reason"]
    )
    assert verdict["merge_groups"] == []
    assert verdict["n_tree_leaves_unclaimed"] == 1


def test_complement_side_cannot_add_an_edge_separable_fragment_to_a_group():
    # Arrange: this is a merge-level counterexample to the claim that merely
    # adding complement sides made the old rooted grouping monotone. The old
    # code grouped A+D+E and left B alone; the first unrooted patch chose a
    # complement span for E and grouped A+B+D+E, even though the edge before i
    # separates every E member (g,h) from every B member (j,k).
    tree = "(((((((((a,b),c),d),e),f),g),h),i),(j,(k,l)));"
    frags = {
        "A": ["a", "b", "c", "i"],
        "B": ["j", "k"],
        "C": ["l"],
        "D": ["d", "e", "f"],
        "E": ["g", "h"],
    }

    # Act
    verdict = fragment_verdict(tree, frags)

    # Assert: the status monotonicity argument was valid, but the former
    # smallest-containing-side partner attribution was not.
    assert verdict["merge_groups"] == [["A", "D", "E"]]
    assert (
        "every other fragment is separable"
        in verdict["fragments"]["B"]["reason"]
    )
