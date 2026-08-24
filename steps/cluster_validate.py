"""Tree validation of fragmentation clusters (issue #47, family table v3).

The one-way clusters in fragmentation_clusters.tsv over-merge by construction:
a frac >= 0.6 vote edge marks the nearest-neighbour family, so cluster
membership NOMINATES; only the cluster's own tree merges.

They also UNDER-COVER, and that half was missing here until it was measured.
"every family has a nearest neighbour" is true of the data and false of the
file: the fifteen-species edge dump has a minimum frac of exactly 0.600 - the
cut, not the evidence, is what its floor records - and 8,902 of 23,744
families have no outgoing edge at all. 7,683 families therefore reached no
cluster and no tree, 32.4% of the table. The family holding the Mcry PEPC
flagship is one of them: it votes for R1_OG0019668, a member of C0297, the
very cluster the other PEPC pieces merged into, at frac 2/9 = 0.22.

So a verdict here is evidence about the fragments it was given, and silence
about a third of the table. Over-merging is corrected by the tree; this is
not, because the tree never sees the case. steps/profile_assign.vote_edges
exports the edges uncut for exactly this reason.

The tree can say three things about a fragment, and only one of them merges:

  INTERLEAVED   the fragment's members do not form their own edge-defined
                clade. Pairwise-inseparable fragments mix; that is lineage-axis
                fragmentation, the shape MCL fragmentation produces (the PEPC
                pieces split by species, not by subfamily). Connected mixing
                fragments merge. An isolated internal strip does not.
  MONOPHYLETIC  the fragment is its own clade. Topology alone cannot say
                whether that is a subfamily of the same family or a distinct
                neighbouring family - PPC4 is monophyletic inside the PEPC
                cluster and was its own family in the five-species run.
                Undecided: reported with its numbers, merged by nobody.
  missing       members absent from the tree are counted, never assumed in.

Merge groups are connected components of the "mixes with" relation.  That
relation is pairwise and unrooted: two fragments mix only when no tree edge
separates all members of one from all members of the other.  A fragment can
fail to be an edge-defined clade without mixing with a particular neighbour
(for example, when it occupies an internal strip of the tree); that shape is
reported as INTERLEAVED but produces no merge edge by itself.
"""
from pathlib import Path
from typing import Dict, List, Sequence

from utils.newick import parse_newick

_TOPOLOGY_CAVEAT = ("monophyletic — topology alone cannot separate a subfamily "
                    "of the same family from a distinct neighbouring family")
_NO_MIXING_PARTNER = ("non-clade, but every other fragment is separable by "
                      "an edge")


def _leaves(node) -> List[str]:
    if not node.children:
        return [node.name]
    return [x for c in node.children for x in _leaves(c)]


def _edge_splits(root) -> List[tuple]:
    """Return both leaf sides of every edge in the parsed tree.

    ``parse_newick`` represents the printed root as a node, but an unrooted
    tree has no edge above that node.  Every other node has exactly one parent
    edge.  A degree-two printed root is harmless: its two child edges describe
    the same unrooted split, and duplicate splits do not change either test.
    """
    descendants = {}
    for node in root.postorder():
        if not node.children:
            side = {node.name}
        else:
            side = set().union(*(descendants[id(c)] for c in node.children))
        descendants[id(node)] = side

    leaf_set = descendants[id(root)]
    return [
        (descendants[id(node)], leaf_set - descendants[id(node)])
        for node in root.postorder()
        if node is not root
    ]


def _separated_by_edge(a: set, b: set, splits: Sequence[tuple]) -> bool:
    """Whether some unrooted edge puts all of ``a`` opposite all of ``b``.

    This is equivalent to asking whether the two sets' convex hulls (their
    minimal connecting subtrees) are disjoint.  If the hulls are disjoint, an
    edge on the unique path between them separates the sets.  Conversely, a
    separating edge puts each entire hull in a different component, so the
    hulls cannot overlap.  Thus *no* separating edge is exactly pairwise
    interleaving, without choosing a root or an arbitrary containing side.
    """
    return any(
        (a <= left and b <= right) or (a <= right and b <= left)
        for left, right in splits
    )


def fragment_verdict(newick: str, fragments: Dict[str, Sequence[str]]) -> dict:
    """Judge one cluster: which fragments does its tree actually mix?

    Returns fragment statuses, merge groups (sorted lists of fragment ids),
    per-fragment missing-member counts, and the count of tree leaves no
    fragment claims - every mismatch is a number, not a silence.
    """
    root = parse_newick(newick.strip())
    tree_leaves = _leaves(root)
    leaf_set = set(tree_leaves)
    owner = {}
    for fam, members in fragments.items():
        for g in members:
            owner[g] = fam

    present: Dict[str, set] = {
        fam: {g for g in members if g in leaf_set}
        for fam, members in fragments.items()
    }

    # For each fragment present in the tree: is there an EDGE separating its
    # members from everything else?
    #
    # FastTree emits UNROOTED trees and puts the root wherever it likes, so
    # "clade" cannot mean "descendants of some node" (issue #51): that reads
    # two Newick strings of the same topology differently, and the campaign's
    # verdicts would then be an artifact of notation. In an unrooted tree
    # every non-root node defines a SPLIT - its descendants on one side, the
    # rest on the other - and both sides are candidate clades.
    splits = _edge_splits(root)
    mixes: Dict[str, set] = {fam: set() for fam in fragments}
    status: Dict[str, str] = {}
    for fam, mem in present.items():
        if not mem:
            status[fam] = "ABSENT"
            continue
        # The full leaf set is the one intentional exception to the literal
        # edge rule: a fragment holding the whole tree is trivially its own
        # clade even though no edge has an empty side.
        if mem == leaf_set or any(mem == side for split in splits
                                  for side in split):
            status[fam] = "MONOPHYLETIC"
            continue
        status[fam] = "INTERLEAVED"

    # "Smallest split side containing mem" is not well-defined in an
    # unrooted tree: equally small sides can name different neighbours, and
    # which one is visited first changes after rerooting.  Use the invariant
    # pairwise statement instead.  Equivalently, the two fragments' minimal
    # connecting subtrees overlap exactly when no edge separates them.
    interleaved = {f for f, s in status.items() if s == "INTERLEAVED"}
    ordered = sorted(interleaved)
    for i, fam in enumerate(ordered):
        for other in ordered[i + 1:]:
            if not _separated_by_edge(present[fam], present[other], splits):
                mixes[fam].add(other)
                mixes[other].add(fam)

    # Merge groups are connected components over the symmetric mixing
    # relation. Only fragments that failed monophyly participate. For disjoint
    # fragments that restriction follows from the definition anyway: a
    # monophyletic fragment's defining edge separates it from every other
    # fragment. Keeping the guard makes the merge policy explicit.
    par = {}
    def find(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for fam in interleaved:
        for other in mixes[fam]:
            if other in interleaved:
                a, b = find(fam), find(other)
                if a != b:
                    par[a] = b
    comps: Dict[str, List[str]] = {}
    for fam in interleaved:
        comps.setdefault(find(fam), []).append(fam)
    merge_groups = sorted(sorted(v) for v in comps.values() if len(v) > 1)

    out_frags = {}
    for fam, members in fragments.items():
        out_frags[fam] = {
            "status": status[fam],
            "n_members": len(members),
            "n_present": len(present[fam]),
            "n_missing": len(members) - len(present[fam]),
            "reason": (_TOPOLOGY_CAVEAT if status[fam] == "MONOPHYLETIC"
                       else "members mix with: " + ",".join(sorted(mixes[fam]))
                       if mixes[fam]
                       else _NO_MIXING_PARTNER
                       if status[fam] == "INTERLEAVED" else "no members in tree"),
        }
    return {
        "fragments": out_frags,
        "merge_groups": merge_groups,
        "n_tree_leaves_unclaimed": sum(1 for g in tree_leaves
                                       if g not in owner),
    }


def write_verdicts(rows: Sequence[dict], outdir: Path) -> Path:
    """One TSV row per (cluster, fragment), with the merge group named."""
    path = Path(outdir) / "cluster_verdicts.tsv"
    with open(path, "w") as f:
        f.write("cluster_id\tfragment\tstatus\tn_members\tn_missing\t"
                "merge_group\treason\n")
        for row in rows:
            v = row["verdict"]
            group_of = {}
            for grp in v["merge_groups"]:
                for fam in grp:
                    group_of[fam] = "+".join(grp)
            for fam in sorted(v["fragments"]):
                d = v["fragments"][fam]
                f.write(f"{row['cluster_id']}\t{fam}\t{d['status']}\t"
                        f"{d['n_members']}\t{d['n_missing']}\t"
                        f"{group_of.get(fam, '')}\t{d['reason']}\n")
    return path
