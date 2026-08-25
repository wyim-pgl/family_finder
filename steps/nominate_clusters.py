"""Budget-capped nomination clustering for fragmentation/merge candidates.

This module turns uncut one-way vote edges into deterministic nomination units
for the tree arbiter. It does not merge families and it does not choose a
global threshold: every input link is accounted for as either inside a
cluster, a pairwise nomination, or a deferred edge.
"""

from itertools import groupby
from typing import Dict, Iterable, List, Sequence, Tuple

from config import Config

VoteEdge = Tuple[str, str, int, int, float]
EdgeKey = Tuple[str, str]


def collapse_vote_edges(edges: Sequence[VoteEdge]) -> Dict[EdgeKey, dict]:
    """Collapse directed vote edges to undirected links keyed by sorted family ids.

    The link weight is the max directed frac, recomputed exactly as
    votes/from_size: vote_edges.tsv serializes frac at 3 decimals, so a
    disk round-trip would otherwise zero small weights (1/2500 -> 0.000)
    and fabricate ties (0.9374 and 0.9366 -> 0.937) that change the
    decomposition. Both directional vote/frac values are retained, and
    reciprocity is recorded (by votes, never by rounded frac) but never
    required. A directed edge appearing twice is an error: vote_edges
    emits each (from, to) once, so a duplicate means concatenated scans.
    """
    links: Dict[EdgeKey, dict] = {}
    seen_directed = set()
    for fam_from, fam_to, votes, from_size, _frac in edges:
        if fam_from == fam_to:
            raise ValueError(f"self-edge is invalid for nomination: {fam_from!r}")
        if from_size <= 0:
            raise ValueError(f"edge {fam_from!r} -> {fam_to!r} has from_size {from_size}")
        if votes <= 0:
            raise ValueError(f"edge {fam_from!r} -> {fam_to!r} has votes {votes}")
        if (fam_from, fam_to) in seen_directed:
            raise ValueError(f"duplicate directed edge {fam_from!r} -> {fam_to!r}")
        seen_directed.add((fam_from, fam_to))
        frac = votes / from_size
        key = tuple(sorted((fam_from, fam_to)))
        rec = links.setdefault(key, {
            "edge_key": key,
            "families": key,
            "fam_a": key[0],
            "fam_b": key[1],
            "votes_ab": 0,
            "votes_ba": 0,
            "frac_ab": 0.0,
            "frac_ba": 0.0,
            "weight": 0.0,
            "reciprocal": False,
        })
        if fam_from == rec["fam_a"] and fam_to == rec["fam_b"]:
            rec["votes_ab"] = votes
            rec["frac_ab"] = frac
        elif fam_from == rec["fam_b"] and fam_to == rec["fam_a"]:
            rec["votes_ba"] = votes
            rec["frac_ba"] = frac
        else:
            raise ValueError(f"edge key mismatch for {fam_from!r} -> {fam_to!r}")
        rec["weight"] = max(rec["frac_ab"], rec["frac_ba"])
        rec["reciprocal"] = rec["votes_ab"] > 0 and rec["votes_ba"] > 0
    return links


def _validate_family_sizes(links: Dict[EdgeKey, dict], family_sizes: Dict[str, int]) -> None:
    for link in links.values():
        for family in link["families"]:
            if family not in family_sizes:
                raise ValueError(f"family size missing for {family!r}")
            if family_sizes[family] <= 0:
                raise ValueError(f"family {family!r} has non-positive size {family_sizes[family]}")


def _validate_edge_sizes(edges: Sequence[VoteEdge], family_sizes: Dict[str, int]) -> None:
    """An edge's from_size must match the family table it is nominated against.

    A mismatch means the edge file and the family table come from different
    runs (or the table changed since the scan) — every frac in the file is
    then wrong for this table, so refuse loudly instead of nominating from
    stale fractions.
    """
    for fam_from, fam_to, _votes, from_size, _frac in edges:
        expected = family_sizes.get(fam_from)
        if expected is not None and expected != from_size:
            raise ValueError(
                f"edge {fam_from!r} -> {fam_to!r} carries from_size {from_size} "
                f"but the family table says {expected}; edge file and family "
                "table are from different runs"
            )


def _sum_genes(families: Iterable[str], family_sizes: Dict[str, int]) -> int:
    return sum(family_sizes[fam] for fam in families)


def _components_from_links(links: Dict[EdgeKey, dict]) -> List[Tuple[str, ...]]:
    graph: Dict[str, set] = {}
    for fam_a, fam_b in links:
        graph.setdefault(fam_a, set()).add(fam_b)
        graph.setdefault(fam_b, set()).add(fam_a)

    seen = set()
    components: List[Tuple[str, ...]] = []
    for start in sorted(graph):
        if start in seen:
            continue
        stack = [start]
        comp = []
        seen.add(start)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nxt in sorted(graph[cur]):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        comp_t = tuple(sorted(comp))
        if len(comp_t) > 1:
            components.append(comp_t)
    components.sort(key=lambda c: (-len(c), c[0]))
    return components


def _global_cluster(component: Tuple[str, ...], comp_links: List[dict],
                    family_sizes: Dict[str, int]) -> dict:
    edge_keys = tuple(sorted(link["edge_key"] for link in comp_links))
    return {
        "families": component,
        "total_genes": _sum_genes(component, family_sizes),
        "n_families": len(component),
        "min_admitted_weight": None,
        "decomposed": False,
        "edge_keys": edge_keys,
    }


def _pairwise_record(link: dict, family_sizes: Dict[str, int], reason: str) -> dict:
    fam_a, fam_b = link["families"]
    return {
        "families": link["families"],
        "edge_key": link["edge_key"],
        "weight": link["weight"],
        "votes_ab": link["votes_ab"],
        "votes_ba": link["votes_ba"],
        "frac_ab": link["frac_ab"],
        "frac_ba": link["frac_ba"],
        "reciprocal": link["reciprocal"],
        "total_genes": family_sizes[fam_a] + family_sizes[fam_b],
        "reason": reason,
    }


def _deferred_record(link: dict, family_sizes: Dict[str, int], reason: str) -> dict:
    rec = _pairwise_record(link, family_sizes, reason)
    rec["retry_count"] = 0
    rec["terminal"] = False
    return rec


def _decompose_component(component: Tuple[str, ...], comp_links: List[dict],
                         family_sizes: Dict[str, int], cap: int) -> Tuple[List[dict], List[dict], List[dict], List[dict]]:
    parent = {fam: fam for fam in component}
    comp_size = {fam: family_sizes[fam] for fam in component}
    admitted_weights: Dict[str, List[float]] = {fam: [] for fam in component}
    tie_groups: List[dict] = []
    pairwise: List[dict] = []
    deferred: List[dict] = []
    classified_edges: Dict[EdgeKey, str] = {}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str, weight: float) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if (comp_size[ra], ra) < (comp_size[rb], rb):
            ra, rb = rb, ra
        parent[rb] = ra
        comp_size[ra] += comp_size[rb]
        admitted_weights[ra].extend(admitted_weights[rb])
        admitted_weights[ra].append(weight)
        del comp_size[rb]
        del admitted_weights[rb]

    ordered = sorted(comp_links, key=lambda rec: (-rec["weight"], rec["edge_key"]))
    for weight, group_iter in groupby(ordered, key=lambda rec: rec["weight"]):
        group = list(group_iter)
        if len(group) > 1:
            tie_groups.append({
                "weight": weight,
                "edges": tuple(link["edge_key"] for link in group),
            })
        for link in group:
            fam_a, fam_b = link["families"]
            root_a, root_b = find(fam_a), find(fam_b)
            if root_a == root_b:
                continue
            if comp_size[root_a] + comp_size[root_b] <= cap:
                union(fam_a, fam_b, weight)
                continue
            pair_size = family_sizes[fam_a] + family_sizes[fam_b]
            if pair_size <= cap:
                pairwise.append(_pairwise_record(link, family_sizes, "refused_by_cap"))
                classified_edges[link["edge_key"]] = "pairwise"
            else:
                # The condition is on the pair SUM: neither endpoint need
                # exceed the cap on its own (two 300-gene families defer
                # under cap 500). The reason string must not claim more.
                deferred.append(_deferred_record(link, family_sizes, "pair_exceeds_cap"))
                classified_edges[link["edge_key"]] = "deferred"

    comp_members: Dict[str, List[str]] = {}
    for fam in component:
        comp_members.setdefault(find(fam), []).append(fam)

    clusters: List[dict] = []
    for root, families in comp_members.items():
        families_t = tuple(sorted(families))
        if len(families_t) < 2:
            continue
        edge_keys = tuple(sorted(
            link["edge_key"]
            for link in comp_links
            if link["edge_key"] not in classified_edges
            and find(link["fam_a"]) == root
            and find(link["fam_b"]) == root
        ))
        # min over ADMITTED spanning-edge weights. This is NOT a cut
        # threshold: budget refusal is not monotone in weight, so a refused
        # edge incident to this cluster may be heavier, and a cycle-closing
        # edge lighter than this value can still sit inside edge_keys.
        min_admitted = min(admitted_weights[root]) if admitted_weights[root] else None
        clusters.append({
            "families": families_t,
            "total_genes": _sum_genes(families_t, family_sizes),
            "n_families": len(families_t),
            "min_admitted_weight": min_admitted,
            "decomposed": True,
            "edge_keys": edge_keys,
        })

    # Edges not explicitly refused/deferred are inside one final cluster.
    final_root = {fam: find(fam) for fam in component}
    cluster_roots = {find(fams[0]) for fams in comp_members.values() if len(fams) > 1}
    for link in comp_links:
        key = link["edge_key"]
        if key in classified_edges:
            continue
        root_a = final_root[link["fam_a"]]
        root_b = final_root[link["fam_b"]]
        if root_a != root_b:
            raise AssertionError(f"edge {key} was neither classified nor kept inside a cluster")
        if root_a not in cluster_roots:
            raise AssertionError(f"edge {key} ended in a singleton component")

    clusters.sort(key=lambda c: (-c["total_genes"], -c["n_families"], c["families"]))
    pairwise.sort(key=lambda r: (-r["weight"], r["families"]))
    deferred.sort(key=lambda r: (-r["weight"], r["families"]))
    return clusters, pairwise, deferred, tie_groups


def nominate_clusters(edges: Sequence[VoteEdge], family_sizes: Dict[str, int],
                      config: Config) -> dict:
    """Nominate under-budget clusters, pairwise cut edges, and deferred edges."""
    links = collapse_vote_edges(edges)
    _validate_family_sizes(links, family_sizes)
    _validate_edge_sizes(edges, family_sizes)
    cap = config.max_cluster_genes
    if cap <= 0:
        raise ValueError(f"max_cluster_genes must be positive, got {cap}")

    clusters: List[dict] = []
    pairwise: List[dict] = []
    deferred: List[dict] = []
    tie_groups: List[dict] = []

    # Group links by component in one pass: both endpoints of a link are in
    # the same component by construction, so a family -> component index
    # replaces the per-component membership scan that was O(components x
    # links x component size) — 122 s on the 15-species hairball, seconds
    # this way.
    components = _components_from_links(links)
    comp_index: Dict[str, int] = {}
    for idx, component in enumerate(components):
        for fam in component:
            comp_index[fam] = idx
    links_by_comp: List[List[dict]] = [[] for _ in components]
    for key in sorted(links):
        links_by_comp[comp_index[key[0]]].append(links[key])

    for idx, component in enumerate(components):
        comp_links = links_by_comp[idx]
        total_genes = _sum_genes(component, family_sizes)
        if total_genes <= cap:
            clusters.append(_global_cluster(component, comp_links, family_sizes))
            continue
        dec_clusters, dec_pairwise, dec_deferred, dec_ties = _decompose_component(
            component, comp_links, family_sizes, cap
        )
        clusters.extend(dec_clusters)
        pairwise.extend(dec_pairwise)
        deferred.extend(dec_deferred)
        tie_groups.extend(dec_ties)

    clusters.sort(key=lambda c: (-c["total_genes"], -c["n_families"], c["families"]))
    pairwise.sort(key=lambda r: (-r["weight"], r["families"]))
    deferred.sort(key=lambda r: (-r["weight"], r["families"]))
    tie_groups.sort(key=lambda g: (-g["weight"], g["edges"]))

    accounted = set()
    for cluster in clusters:
        for edge_key in cluster["edge_keys"]:
            if edge_key in accounted:
                raise AssertionError(f"edge {edge_key} accounted twice")
            accounted.add(edge_key)
    for recs in (pairwise, deferred):
        for rec in recs:
            edge_key = tuple(rec["edge_key"])
            if edge_key in accounted:
                raise AssertionError(f"edge {edge_key} accounted twice")
            accounted.add(edge_key)
    missing = set(links) - accounted
    if missing:
        raise AssertionError(f"unaccounted nomination edges: {sorted(missing)}")

    return {
        "clusters": clusters,
        "pairwise_nominations": pairwise,
        "deferred_edges": deferred,
        "tie_groups": tie_groups,
    }


def advance_deferred_edges(deferred_edges: Sequence[dict],
                           changed_families: Iterable[str]) -> List[dict]:
    """Advance deferred-edge retry state after a caller attempted one retry.

    ``changed_families`` is the set of family ids whose definitions changed
    since the deferral. The check is per edge, not global: an edge is worth
    retrying only when one of ITS OWN endpoints changed — with a global flag,
    any unrelated family changing every cycle would let the same deferred
    edge be re-attempted forever and 'unchanged_after_retry' would be
    unreachable.
    """
    changed = set(changed_families)
    advanced: List[dict] = []
    for rec in deferred_edges:
        new = dict(rec)
        if new.get("terminal"):
            advanced.append(new)
            continue
        if changed.intersection(new["families"]):
            advanced.append(new)
            continue
        new["retry_count"] = int(new.get("retry_count", 0)) + 1
        new["terminal"] = True
        new["reason"] = "unchanged_after_retry"
        advanced.append(new)
    advanced.sort(key=lambda r: (-r["weight"], r["families"]))
    return advanced
