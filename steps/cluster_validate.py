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

  INTERLEAVED   the fragment's members do not form their own clade - they mix
                with another fragment's. That is lineage-axis fragmentation,
                the shape MCL fragmentation produces (the PEPC pieces split by
                species, not by subfamily), and the two fragments are one
                family. Merge.
  MONOPHYLETIC  the fragment is its own clade. Topology alone cannot say
                whether that is a subfamily of the same family or a distinct
                neighbouring family - PPC4 is monophyletic inside the PEPC
                cluster and was its own family in the five-species run.
                Undecided: reported with its numbers, merged by nobody.
  missing       members absent from the tree are counted, never assumed in.

Merge groups are connected components of the "mixes with" relation. A
monophyletic fragment joins no group even when an interleaved one mixes into
its clade span - mixing is measured on the fragment that fails monophyly.
"""
from pathlib import Path
from typing import Dict, List, Sequence

from utils.newick import parse_newick

_TOPOLOGY_CAVEAT = ("monophyletic — topology alone cannot separate a subfamily "
                    "of the same family from a distinct neighbouring family")


def _leaves(node) -> List[str]:
    if not node.children:
        return [node.name]
    return [x for c in node.children for x in _leaves(c)]


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

    # For each fragment present in the tree: is its MRCA's leaf set exactly
    # its own members? If not, the extra leaves say who it mixes with.
    mixes: Dict[str, set] = {fam: set() for fam in fragments}
    status: Dict[str, str] = {}
    for fam, mem in present.items():
        if not mem:
            status[fam] = "ABSENT"
            continue
        # smallest clade containing mem
        best = None
        stack = [root]
        while stack:
            n = stack.pop()
            L = set(_leaves(n))
            if mem <= L and (best is None or len(L) < len(best)):
                best = L
            stack.extend(n.children)
        extra = best - mem
        if not extra:
            status[fam] = "MONOPHYLETIC"
            continue
        status[fam] = "INTERLEAVED"
        for g in extra:
            other = owner.get(g)
            if other and other != fam:
                mixes[fam].add(other)

    # Merge groups: connected components over the mixing relation, but ONLY
    # fragments that themselves failed monophyly participate. A monophyletic
    # fragment nested inside a messy span stays undecided.
    par = {}
    def find(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    interleaved = {f for f, s in status.items() if s == "INTERLEAVED"}
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
