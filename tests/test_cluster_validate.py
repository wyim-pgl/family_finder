"""Merge fragments only when the tree says they are one family (issue #47).

The one-way fragmentation clusters over-merge by construction (nearest-neighbor
vote concentration), so cluster membership alone must never merge anything.
The tree decides, and it can say three things about a fragment:

  INTERLEAVED   its members do not form an edge-defined clade. Pairwise
                inseparable fragments mix, and merge. An isolated internal
                strip does not.
  MONOPHYLETIC  its members form their own clade. Topology alone CANNOT
                distinguish "subfamily of the same family" from "distinct
                neighbouring family" - PPC4 is monophyletic inside the PEPC
                cluster and was its own family in the 5-species run. Undecided,
                with the metrics reported.
  (missing)     members absent from the tree are counted, never assumed.

⚠️ These fixtures are synthetic. PEPC is NOT the worked example for the
INTERLEAVED shape - this file used to call it "the PEPC signature". Measured
on the 15-species run the split runs largely along the SUBFAMILY axis (all 15
species appear in 3-5 of the six fragments), and a subfamily is a clade, so
the rule merges only two of the six. See steps/cluster_validate.py's module
docstring for the numbers.
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
    """Three fragments, no pair separable by an edge: all mix.

    Synthetic, not PEPC - see the module docstring.
    """
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


# ---------------------------------------------------------------------------
# outgroup-anchored verdicts
# ---------------------------------------------------------------------------
# Without an outgroup the tree cannot answer "are these one family". A
# subfamily is a clade by construction, so MONOPHYLETIC is permanent silence -
# which is why the PEPC clan sat undecided across every table version while
# non-topological evidence (structure, external anchors, the clan tree) said
# it was one family.
#
# An outgroup - a sequence known to lie OUTSIDE the family - makes the
# question decidable, and it is what clan_merge actually used: bacterial-type
# PEPC against plant-type, 0/95 intrusions. The test becomes "do these
# fragments group together to the exclusion of the outgroup", and it yields
# the verdict the rule could never reach before: SEPARATE.

def _og(newick, frags, outgroup):
    v = fragment_verdict(newick, frags, outgroup=outgroup)
    return ({f: v["fragments"][f]["status"] for f in v["fragments"]},
            [sorted(g) for g in v["merge_groups"]])


def test_outgroup_merges_fragments_that_exclude_it():
    # Arrange: A and B are both inside; the outgroup hangs off outside them.
    # This is the PEPC shape - two monophyletic subfamilies of one family.
    frags = {"A": ["a1", "a2"], "B": ["b1", "b2"]}

    # Act
    statuses, groups = _og("(((a1,a2),(b1,b2)),o1);", frags, ["o1"])

    # Assert: undecidable without the outgroup, decided with it
    assert statuses == {"A": "SAME_FAMILY", "B": "SAME_FAMILY"}
    assert groups == [["A", "B"]]


def test_outgroup_between_fragments_makes_them_separate():
    # Arrange: the outgroup nests on BOTH sides, so no edge puts A and B
    # together against it. Each fragment is closer to an outgroup member than
    # to the other fragment - by the outgroup's definition, distinct families.
    frags = {"A": ["a1", "a2"], "B": ["b1", "b2"]}

    # Act
    statuses, groups = _og("(((a1,a2),o1),((b1,b2),o2));", frags,
                           ["o1", "o2"])

    # Assert: the verdict the rule could not previously produce
    assert statuses == {"A": "SEPARATE", "B": "SEPARATE"}
    assert groups == []


def test_a_single_leaf_outgroup_cannot_separate_anything():
    # Arrange: one outgroup leaf hangs off its own pendant edge, so an edge
    # always separates it from everything else and every fragment comes back
    # SAME_FAMILY. That is honest - a lone outgroup says only "all of this is
    # inside relative to me" - but it means a single sequence carries no
    # discriminating power. Pin it so nobody reads such a merge as evidence.
    frags = {"A": ["a1", "a2"], "B": ["b1", "b2"]}

    # Act
    nested, _ = _og("((a1,a2),o1,(b1,b2));", frags, ["o1"])
    sister, _ = _og("(((a1,a2),(b1,b2)),o1);", frags, ["o1"])

    # Assert
    assert nested == sister == {"A": "SAME_FAMILY", "B": "SAME_FAMILY"}


def test_outgroup_splits_one_group_from_another():
    # Arrange: an outgroup member sits between C and A+B, so C is nearer to
    # something outside the family than A and B are to it. A+B merge; C does
    # not.
    frags = {"A": ["a1"], "B": ["b1"], "C": ["c1"]}

    # Act
    statuses, groups = _og("(((a1,b1),o1),(o2,c1));", frags, ["o1", "o2"])

    # Assert
    assert statuses["A"] == statuses["B"] == "SAME_FAMILY"
    assert statuses["C"] == "SEPARATE"
    assert groups == [["A", "B"]]


def test_outgroup_verdict_is_invariant_to_rerooting():
    # Arrange: same topology, two root placements
    frags = {"A": ["a1", "a2"], "B": ["b1", "b2"]}

    # Act
    first = _og("(((a1,a2),(b1,b2)),o1);", frags, ["o1"])
    second = _og("(a1,(a2,((b1,b2),o1)));", frags, ["o1"])

    # Assert
    assert first == second


def test_a_fragment_straddling_the_outgroup_is_reported_not_merged():
    # Arrange: A's members sit on both sides of the outgroup. The fragment
    # itself is not a coherent unit relative to this outgroup, and saying so
    # is more honest than merging or silently dropping it.
    frags = {"A": ["a1", "a2"], "B": ["b1"]}

    # Act
    statuses, groups = _og("((a1,o1),(o2,(a2,b1)));", frags, ["o1", "o2"])

    # Assert
    assert statuses["A"] == "OUTGROUP_SPLIT"
    assert "A" not in [f for g in groups for f in g]


def test_an_absent_outgroup_falls_back_to_the_topology_rule():
    # Arrange: the outgroup was requested but is not in the tree - the caller
    # must not be told SAME_FAMILY on the strength of a sequence that never
    # made it into the alignment.
    frags = {"A": ["a1", "a2"], "B": ["b1", "b2"]}

    # Act
    with_missing = _og("((a1,a2),(b1,b2));", frags, ["not_in_tree"])
    without = _og("((a1,a2),(b1,b2));", frags, [])

    # Assert
    assert with_missing == without
    assert with_missing[0]["A"] == "MONOPHYLETIC"


def test_outgroup_members_are_not_counted_as_unclaimed_leaves():
    # Arrange
    frags = {"A": ["a1"], "B": ["b1"]}

    # Act
    v = fragment_verdict("((a1,b1),o1);", frags, outgroup=["o1"])

    # Assert
    assert v["n_tree_leaves_unclaimed"] == 0


def test_a_fragment_claiming_an_outgroup_member_is_refused():
    # Arrange: the outgroup must lie outside the family by construction; a
    # fragment that owns it makes the whole test circular.
    frags = {"A": ["a1", "o1"]}

    # Act / Assert
    try:
        fragment_verdict("((a1,o1),b1);", frags, outgroup=["o1"])
    except ValueError as e:
        assert "o1" in str(e)
    else:
        raise AssertionError("expected a fragment claiming the outgroup to be refused")


def test_a_monophyletic_outgroup_always_yields_one_family():
    # Arrange: when the outgroup is itself a clade, its complement IS the
    # ingroup, so every fragment inside is one family by construction. That is
    # correct and it is also the normal case for a well-chosen outgroup - which
    # means SEPARATE only ever arises when the outgroup interleaves with the
    # fragments, i.e. something outside the family sits BETWEEN them. Pin the
    # property so nobody reads a clean merge as a strong result.
    frags = {"A": ["a1"], "B": ["b1"], "C": ["c1"]}

    # Act
    statuses, groups = _og("((a1,b1),((o1,o2),c1));", frags, ["o1", "o2"])

    # Assert
    assert set(statuses.values()) == {"SAME_FAMILY"}
    assert groups == [["A", "B", "C"]]
