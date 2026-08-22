"""Subfamily diagnostics for ANY family (generalized from the PEPC pilot).

Three independent, pure layers — callers supply plain dicts so everything
tests without ete4 or external tools:

  * sdp_scan: specificity-determining positions — alignment columns where
    one subfamily is conserved on a residue the rest of the family rarely
    uses. This is the "WHY does it split" evidence at residue level.
  * taxonomic_composition: is each subfamily a lineage-specific block
    or an ancient paralog split that exists "by itself"? Default evidence
    (issue #27) is the species tree the pipeline already requires: a
    subfamily whose species set is a clade strictly inside the family's
    span is lineage-specific; a non-monophyletic set (interleaved with
    the other subfamilies' species) or one spanning the family's own
    MRCA is a paralog split. A taxonomy TSV is an optional LABEL layer
    on top (rank names like "Cactaceae (family)"); without a species
    tree the legacy rank-purity verdict is used.
  * structure_coherence: given foldseek all-vs-all scores, are subfamilies
    also structurally coherent (mean within-score > mean between-score)?

Foldseek's role: the pipeline already uses foldseek for tier-3 assignment
(gene vs family representative, steps/esm.py) and clan-merge evidence.
Here it is reused a third way — all-vs-all WITHIN a family to ask whether
the sequence-defined subfamilies are visible in structure space at all.
"""

import logging
import random
import statistics
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

from utils.newick import Node

logger = logging.getLogger("family_finder")

# Ranks checked lowest-to-highest for the lineage-specific verdict.
TAXONOMY_RANKS = ["species", "genus", "family", "order"]


# ---------------------------------------------------------------------------
# SDP scan
# ---------------------------------------------------------------------------

def sdp_scan(
    alignment: Dict[str, str],
    groups: Dict[str, List[str]],
    *,
    min_group: int = 5,
    in_cons: float = 0.8,
    out_max: float = 0.2,
    min_cover: float = 0.7,
    ref_seq_id: Optional[str] = None,
) -> List[dict]:
    """Subfamily-diagnostic alignment columns.

    A column is diagnostic for group G when >= in_cons of G's non-gap
    members share one residue, G's non-gap coverage is >= min_cover, and
    that residue's frequency among the REST of the alignment's non-gap
    residues is <= out_max. Groups smaller than min_group are skipped
    (their "conservation" is not meaningful). ref_seq_id adds a ref_pos
    column with that sequence's ungapped numbering (e.g. ATH_PPC1
    coordinates) so residues can be cited against published positions.
    """
    lengths = {len(s) for s in alignment.values()}
    if len(lengths) > 1:
        raise ValueError(f"ragged alignment: lengths {sorted(lengths)}")
    alen = lengths.pop() if lengths else 0

    ref_pos: Dict[int, int] = {}
    if ref_seq_id is not None:
        ref = alignment[ref_seq_id]
        n = 0
        for col, ch in enumerate(ref, start=1):
            if ch != "-":
                n += 1
                ref_pos[col] = n

    rows: List[dict] = []
    for group_id, members in sorted(groups.items()):
        members = [m for m in members if m in alignment]
        if len(members) < min_group:
            continue
        member_set = set(members)
        others = [g for g in alignment if g not in member_set]
        for col in range(alen):
            inside = [alignment[g][col] for g in members]
            nongap = [c for c in inside if c != "-"]
            if not nongap or len(nongap) / len(inside) < min_cover:
                continue
            res, cnt = Counter(nongap).most_common(1)[0]
            in_freq = cnt / len(nongap)
            if in_freq < in_cons:
                continue
            out_nongap = [alignment[g][col] for g in others
                          if alignment[g][col] != "-"]
            if not out_nongap:
                continue
            out_freq = sum(1 for c in out_nongap if c == res) / len(out_nongap)
            if out_freq > out_max:
                continue
            rest_res = Counter(out_nongap).most_common(1)[0][0]
            rows.append({
                "subfamily": group_id,
                "n_members": len(members),
                "aln_col": col + 1,
                "ref_pos": ref_pos.get(col + 1),
                "sf_residue": res,
                "sf_freq": round(in_freq, 3),
                "rest_residue": rest_res,
                "rest_freq_of_sf_residue": round(out_freq, 3),
            })
    return rows


def coverage_suppressed(alignment: Dict[str, str],
                        groups: Dict[str, List[str]],
                        min_cover: float = 0.7,
                        min_group: int = 5) -> Dict[str, dict]:
    """Which columns `sdp_scan` never judged for each group, and why.

    `sdp_scan` skips a column for a group whose non-gap coverage falls below
    `min_cover`, and skips a group smaller than `min_group` outright. Both are
    correct — conservation over three of seventy sequences is not conservation
    — but both are silent, so "no diagnostic columns here" is indistinguishable
    from "this group was never examined here". Issue #40 hit exactly that: the
    ppc-1E1 subfamily covers the N-terminal regulatory region in only 33 of 70
    members, so every N-terminal column was filtered out and the region read as
    signal-free.

    Call this next to `sdp_scan` and report it alongside the hits.
    """
    lengths = {len(s) for s in alignment.values()}
    if len(lengths) > 1:
        raise ValueError(f"ragged alignment: lengths {sorted(lengths)}")
    alen = lengths.pop() if lengths else 0

    report: Dict[str, dict] = {}
    for group_id, members in sorted(groups.items()):
        present = [m for m in members if m in alignment]
        coverage = []
        for col in range(alen):
            if not present:
                coverage.append(0.0)
                continue
            nongap = sum(1 for m in present if alignment[m][col] != "-")
            coverage.append(nongap / len(present))
        suppressed = [c + 1 for c, cov in enumerate(coverage) if cov < min_cover]
        report[group_id] = {
            "n_members": len(present),
            "skipped_too_small": len(present) < min_group,
            "n_suppressed": len(suppressed),
            "columns": suppressed,
            "coverage": coverage,
            "min_coverage": min(coverage) if coverage else None,
        }
    return report


def sdp_core_relationship(
    alignment: Dict[str, str],
    groups: Dict[str, List[str]],
    *,
    min_group: int = 5,
    in_cons: float = 0.8,
    out_max: float = 0.2,
    min_cover: float = 0.7,
    invariant_threshold: float = 0.95,
    n_null: int = 1000,
    seed: int = 0,
) -> dict:
    """Where diagnostic residues sit relative to the reference-free core.

    Issue #32 concluded that the PEPC diagnostic residues avoid the catalytic
    core. Stating that required knowing which motif is catalytic, which is
    available for a handful of families and for none of the rest. The
    reference-free stand-in is `utils.alignment.invariant_columns`: the columns
    every subfamily independently agrees on.

    **The obvious comparison is worthless.** A diagnostic column means one group
    differs; an invariant column means none do. The two sets are disjoint by
    definition, so their overlap is always zero and says nothing — measured on
    the PEPC clan it is zero at every threshold from 0.80 to 1.00. `overlap` is
    reported with `overlap_is_by_construction` set so nobody mistakes it for
    evidence.

    The claim that carries information is spatial: are the diagnostic columns
    *further from* the core than columns picked at random from the same
    scannable pool? That is testable against a null and needs no reference.

    Returns a verdict of `avoids_core`, `at_core`, `indistinguishable`, or
    `no_interpretation_available` with a `reason` — silence and "no signal"
    have to be distinguishable (issue #40).
    """
    from utils.alignment import group_occupancy, invariant_columns

    lengths = {len(s) for s in alignment.values()}
    if len(lengths) > 1:
        raise ValueError(f"ragged alignment: lengths {sorted(lengths)}")
    alen = lengths.pop() if lengths else 0

    skipped = sorted(name for name, members in groups.items()
                     if len([m for m in members if m in alignment]) < min_group)
    hits = sdp_scan(alignment, groups, min_group=min_group, in_cons=in_cons,
                    out_max=out_max, min_cover=min_cover)
    sdp_cols = sorted({h["aln_col"] for h in hits})
    invariant = invariant_columns(alignment, groups,
                                  threshold=invariant_threshold,
                                  min_cover=min_cover)
    inv_cols = sorted(invariant)

    # Columns sdp_scan was able to judge at all: every group covers them.
    occ = group_occupancy(alignment, groups)
    candidates = [c + 1 for c in range(alen)
                  if occ and all(v[c] >= min_cover for v in occ.values())]

    result = {
        "n_sdp_columns": len(sdp_cols),
        "n_invariant_columns": len(inv_cols),
        "n_candidate_columns": len(candidates),
        "skipped_groups": skipped,
        "overlap": len(set(sdp_cols) & set(inv_cols)),
        "overlap_is_by_construction": True,
        "observed_median_distance": None,
        "null_median_distance": None,
        "null_p05": None,
        "null_p95": None,
        "p_value": None,
        "verdict": "no_interpretation_available",
        "reason": "",
    }

    if skipped and len(skipped) == len(groups):
        result["reason"] = (
            f"every group is smaller than min_group={min_group}: {skipped}"
        )
        return result
    if not sdp_cols:
        result["reason"] = "no diagnostic columns were found for any group"
        return result
    if not inv_cols:
        result["reason"] = (
            f"no invariant columns at threshold {invariant_threshold} — the "
            "family has no core the groups agree on, so there is nothing to "
            "measure distance from"
        )
        return result
    if len(candidates) < len(sdp_cols):
        result["reason"] = (
            f"only {len(candidates)} column(s) were scannable, fewer than the "
            f"{len(sdp_cols)} diagnostic column(s); no null can be built"
        )
        return result

    def nearest(col: int) -> int:
        return min(abs(col - j) for j in inv_cols)

    observed = statistics.median(nearest(c) for c in sdp_cols)
    rng = random.Random(seed)
    null = sorted(statistics.median(nearest(c) for c in
                                    rng.sample(candidates, len(sdp_cols)))
                  for _ in range(n_null))
    p_far = sum(1 for v in null if v >= observed) / len(null)
    p_near = sum(1 for v in null if v <= observed) / len(null)

    result["observed_median_distance"] = observed
    result["null_median_distance"] = statistics.median(null)
    result["null_p05"] = null[int(0.05 * len(null))]
    result["null_p95"] = null[int(0.95 * len(null)) - 1]
    if p_far <= 0.05:
        result["verdict"] = "avoids_core"
        result["p_value"] = p_far
    elif p_near <= 0.05:
        result["verdict"] = "at_core"
        result["p_value"] = p_near
    else:
        result["verdict"] = "indistinguishable"
        result["p_value"] = min(p_far, p_near)
    result["reason"] = (
        f"{len(sdp_cols)} diagnostic vs {len(inv_cols)} invariant columns, "
        f"null drawn from {len(candidates)} scannable columns"
    )
    return result


# ---------------------------------------------------------------------------
# Taxonomy attribution
# ---------------------------------------------------------------------------

def load_taxonomy(path: Path) -> Dict[str, Dict[str, str]]:
    """TSV with a header (species, then any of genus/family/order/...).

    Returns species_prefix -> {rank: taxon}.
    """
    taxonomy: Dict[str, Dict[str, str]] = {}
    lines = Path(path).read_text().splitlines()
    if not lines:
        return taxonomy
    header = lines[0].rstrip("\n").split("\t")
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        row = dict(zip(header, fields))
        species = row.pop(header[0])
        taxonomy[species] = row
    return taxonomy


def _mrca_with_depth(node: Node, targets: set, depth: int = 0):
    """Deepest node whose subtree holds every target leaf, with its depth
    (edges from the root). None when targets is not fully under node."""
    for child in node.children:
        hit = _mrca_with_depth(child, targets, depth + 1)
        if hit is not None:
            return hit
    if targets <= node.leaf_names():
        return node, depth
    return None


def species_monophyly(species: set, tree: Node) -> dict:
    """Is a species set a clade on the species tree?

    Judged on the species actually present in the tree; the rest are
    reported in "missing". "intruders" are tree species inside the MRCA
    but outside the set — empty iff the set is monophyletic (a gene-loss
    candidate list when it is not). monophyletic is None when no species
    is in the tree at all.
    """
    tree_leaves = tree.leaf_names()
    in_tree = species & tree_leaves
    missing = sorted(species - tree_leaves)
    if not in_tree:
        return {"monophyletic": None, "mrca_name": "", "mrca_depth": None,
                "n_in_tree": 0, "intruders": [], "missing": missing}
    node, depth = _mrca_with_depth(tree, in_tree)
    intruders = sorted(node.leaf_names() - in_tree)
    return {"monophyletic": not intruders, "mrca_name": node.name,
            "mrca_depth": depth, "n_in_tree": len(in_tree),
            "intruders": intruders, "missing": missing}


def clade_rank_label(species: set, taxonomy: Dict[str, Dict[str, str]]) -> str:
    """Optional label layer: the lowest rank shared by ALL species, e.g.
    "Cactaceae (family)". Empty when taxonomy is absent/incomplete for the
    set (a partial label would be misleading) or no rank is shared."""
    species = set(species)
    if len(species) == 1:
        return f"{next(iter(species))} (species)"
    if not taxonomy or any(s not in taxonomy for s in species):
        return ""
    for rank in TAXONOMY_RANKS[1:]:
        taxa = {taxonomy[s].get(rank) for s in species}
        if len(taxa) == 1 and None not in taxa:
            return f"{taxa.pop()} ({rank})"
    return ""


def _rank_purities(species: List[str],
                   taxonomy: Dict[str, Dict[str, str]]) -> dict:
    """Per-rank dominant taxon + purity (largest single-taxon fraction)."""
    cols: dict = {}
    for rank in TAXONOMY_RANKS:
        if rank == "species":
            taxa = species
        else:
            taxa = [taxonomy.get(s, {}).get(rank, f"unknown:{s}")
                    for s in species]
        dominant, dom_n = Counter(taxa).most_common(1)[0]
        cols[f"{rank}_purity"] = round(dom_n / len(taxa), 3)
        cols[f"{rank}_dominant"] = dominant
    return cols


def taxonomic_composition(
    groups: Dict[str, List[str]],
    taxonomy: Optional[Dict[str, Dict[str, str]]] = None,
    delimiter: str = "_",
    species_tree: Optional[Node] = None,
) -> List[dict]:
    """Per subfamily: does the split follow the species phylogeny, or cross it?

    Species-tree path (default when species_tree is given, issue #27):
    the subfamily's species set is tested for monophyly on the species
    tree. A clade whose MRCA sits strictly below the family's own span
    (the MRCA of ALL subfamilies' species) is "lineage-specific (clade)"
    — the split maps onto one lineage. A non-monophyletic set, or one
    spanning the family's span itself (e.g. two ancient paralogs each
    present in every species), is a paralog split: the divergence
    predates the speciations that could explain it. taxonomy, when
    given, only adds rank labels (clade_label, per-rank purities).

    Legacy path (species_tree=None): for each rank (species < genus <
    family < order) the purity is the largest fraction of members in one
    taxon; the verdict is "lineage-specific (<rank>)" for the LOWEST
    pure rank and "paralog-split (crosses order)" when even order is
    mixed. Species absent from the taxonomy table still count (identity
    comes from the gene-id prefix); higher ranks treat them as singleton
    taxa and the row notes them.
    """
    taxonomy = taxonomy or {}
    per_group = {
        group_id: [m.split(delimiter, 1)[0] for m in members]
        for group_id, members in groups.items()
    }

    family_depth = None
    if species_tree is not None:
        union = set().union(*per_group.values()) if per_group else set()
        span = species_monophyly(union, species_tree)
        family_depth = span["mrca_depth"]

    rows: List[dict] = []
    for group_id, members in sorted(groups.items()):
        species = per_group[group_id]
        unknown = sorted({s for s in species if s not in taxonomy})
        row: dict = {
            "subfamily": group_id,
            "n_members": len(members),
            "n_species": len(set(species)),
        }
        notes: List[str] = []

        if species_tree is not None:
            mono = species_monophyly(set(species), species_tree)
            row["n_in_tree"] = mono["n_in_tree"]
            row["monophyletic"] = mono["monophyletic"]
            row["mrca_name"] = mono["mrca_name"]
            row["mrca_depth"] = mono["mrca_depth"]
            row["clade_label"] = clade_rank_label(set(species), taxonomy)
            if mono["monophyletic"] is None:
                verdict = "unknown (no species in species tree)"
            elif not mono["monophyletic"]:
                verdict = "paralog-split (non-monophyletic)"
            elif mono["mrca_depth"] == family_depth:
                verdict = "paralog-split (spans family root)"
            else:
                verdict = "lineage-specific (clade)"
            if mono["missing"]:
                notes.append(
                    "not in species tree: " + ",".join(mono["missing"]))
            if mono["intruders"]:
                notes.append(
                    "interleaved species: " + ",".join(mono["intruders"]))
            if taxonomy:
                row.update(_rank_purities(species, taxonomy))
        else:
            row.update(_rank_purities(species, taxonomy))
            verdict = "paralog-split (crosses order)"
            for rank in TAXONOMY_RANKS:  # lowest pure rank wins
                if row[f"{rank}_purity"] == 1.0:
                    verdict = f"lineage-specific ({rank})"
                    break

        if unknown and (species_tree is None or taxonomy):
            notes.append("unknown taxonomy for: " + ",".join(unknown))
        row["verdict"] = verdict
        row["notes"] = "; ".join(notes)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Structural coherence (foldseek all-vs-all)
# ---------------------------------------------------------------------------

def parse_pairwise_scores(tsv: Path) -> Dict[frozenset, float]:
    """Foldseek easy-search all-vs-all TSV -> symmetric best bits per pair.

    Columns as in steps.esm.FOLDSEEK_FORMAT_OUTPUT: query, target, evalue,
    bits, alntmscore (no header). Self-hits are dropped; the two search
    directions keep the better bits.
    """
    scores: Dict[frozenset, float] = {}
    with open(tsv) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                logger.debug(f"Skipping malformed pair line: {line[:80]!r}")
                continue
            query, target = fields[0], fields[1]
            if query == target:
                continue
            try:
                bits = float(fields[3])
            except ValueError:
                logger.debug(f"Skipping malformed pair line: {line[:80]!r}")
                continue
            key = frozenset((query, target))
            if bits > scores.get(key, float("-inf")):
                scores[key] = bits
    return scores


def structure_coherence(
    groups: Dict[str, List[str]],
    pair_scores: Dict[frozenset, float],
) -> List[dict]:
    """Mean within-subfamily vs between-subfamily structural score.

    Only OBSERVED pairs count — foldseek omits pairs below its reporting
    threshold, and inventing zeros for them would fake incoherence.
    mean_between is None (and coherent is None) when no cross-subfamily
    pair was observed for the group.
    """
    member_of = {}
    for group_id, members in groups.items():
        for m in members:
            member_of[m] = group_id

    rows: List[dict] = []
    for group_id, members in sorted(groups.items()):
        within, between = [], []
        for pair, score in pair_scores.items():
            a, b = tuple(pair)
            ga, gb = member_of.get(a), member_of.get(b)
            if group_id not in (ga, gb):
                continue
            if ga == gb == group_id:
                within.append(score)
            elif ga is not None and gb is not None:
                between.append(score)
        mean_w = sum(within) / len(within) if within else None
        mean_b = sum(between) / len(between) if between else None
        coherent = (mean_w > mean_b) if (mean_w is not None
                                         and mean_b is not None) else None
        rows.append({
            "subfamily": group_id,
            "n_within_pairs": len(within),
            "n_between_pairs": len(between),
            "mean_within": mean_w,
            "mean_between": mean_b,
            "coherent": coherent,
        })
    return rows


# ---------------------------------------------------------------------------
# Can a reference species' subfamily NAMES be transferred at all?
# ---------------------------------------------------------------------------

def _min_support(raw: Optional[str]) -> Optional[float]:
    """Node label -> the WEAKEST support value in it.

    IQ-TREE writes "SH-aLRT/UFBoot" (e.g. "96.2/94"); both halves have to
    clear the bar, so the weaker one decides.
    """
    if not raw:
        return None
    vals = []
    for part in str(raw).split("/"):
        try:
            vals.append(float(part))
        except ValueError:
            continue
    return min(vals) if vals else None


def anchor_transferability(
    newick: str,
    anchor_labels: Dict[str, str],
    query_prefixes: Optional[List[str]] = None,
    min_support: Optional[float] = None,
) -> List[dict]:
    """Test whether reference subfamily labels designate anything in the query.

    The standard comparative-genomics recipe (e.g. Musa MYB classification,
    PLOS ONE 10.1371/journal.pone.0239275) names query subfamilies by clade
    membership with reference-species anchors. That is only valid when the
    reference subfamilies predate the split with the query lineage.

    Two reference labels whose anchors are each other's nearest relatives —
    with no query gene between them — are a REFERENCE-LINEAGE-SPECIFIC
    duplication: the names cannot designate distinct query groups. Measured
    case: Arabidopsis PPC1 and PPC3 are sisters at support 100 in the PEPC
    clan, and naming by symbol duly labelled six different Caryophyllales
    subfamilies "PPC1".

    Returns one row per label (plus rows with ``label=None`` for anchor-free
    clades, the "lineage-specific expansion" case the paper reports by hand).
    """
    from utils.newick import parse_newick

    root = parse_newick(newick)

    def leaves_of(node):
        if not node.children:
            return [node.name]
        out: List[str] = []
        for c in node.children:
            out += leaves_of(c)
        return out

    all_leaves = leaves_of(root)
    anchors = set(anchor_labels)

    parent_of: Dict[int, object] = {}
    def _index_parents(node):
        for c in node.children:
            parent_of[id(c)] = node
            _index_parents(c)
    _index_parents(root)

    def is_query(leaf: str) -> bool:
        if leaf in anchors:
            return False
        if query_prefixes is None:
            return True
        return any(leaf.startswith(p) for p in query_prefixes)

    def smallest_clade_with(targets: set):
        best = None
        def rec(node):
            nonlocal best
            L = set(leaves_of(node))
            if targets <= L and (best is None or len(L) < len(best)):
                best = L
            for c in node.children:
                rec(c)
        rec(root)
        return best or set(all_leaves)

    def subfamily_clade(own: set, foreign: set):
        """LARGEST clade containing `own` and NO other label's anchors.

        The subfamily is the biggest group the label can claim without
        swallowing a rival label — not the anchor's nearest neighbours.
        Measured why this matters: on the rebuilt PEPC tree the
        nearest-neighbour rule returned a 10-leaf all-Amaranthaceae clade
        for Ppc2 and so reported "no query gene" once queries were
        restricted to cacti, while the real 17-leaf subfamily does contain
        cactus genes.

        Returns (leaf set, node) — the node carries the support label.
        """
        best = None
        def rec(node):
            nonlocal best
            L = set(leaves_of(node))
            if own <= L and not (L & foreign):
                if best is None or len(L) > len(best[0]):
                    best = (L, node)
            for c in node.children:
                rec(c)
        rec(root)
        return best or (set(own), None)

    rows: List[dict] = []
    labels = sorted(set(anchor_labels.values()))
    for label in labels:
        own = {a for a, l in anchor_labels.items() if l == label}
        foreign = anchors - own
        clade, node = subfamily_clade(own, foreign)
        others = sorted({anchor_labels[a] for a in clade & anchors} - {label})
        n_query = sum(1 for x in clade if is_query(x))

        # What sits immediately OUTSIDE the subfamily? If the very next node
        # up adds a rival label's anchors and no query gene, the two labels
        # are a reference-lineage-specific duplication: nothing of the query
        # separates them, so neither name designates a query group.
        parent = parent_of.get(id(node)) if node is not None else None
        added = (set(leaves_of(parent)) - clade) if parent is not None else set()
        added_labels = sorted({anchor_labels[a] for a in added & anchors})
        added_query = sum(1 for x in added if is_query(x))
        dup_with = added_labels if (added_labels and added_query == 0) else []
        raw_support = getattr(node, "name", None) if node is not None else None
        support_val = _min_support(raw_support)
        weak = (min_support is not None
                and (support_val is None or support_val < min_support))
        transferable = ((not others) and (not dup_with)
                        and n_query > 0 and not weak)
        if dup_with:
            verdict = ("reference-lineage-specific duplication with "
                       + ", ".join(dup_with)
                       + " — no query gene separates them, so these names do "
                         "not designate distinct query subfamilies")
        elif weak and not others and n_query > 0:
            verdict = (f"clade support {raw_support} is below the "
                       f"min_support={min_support:g} bar — the grouping is not "
                       "established, do not transfer the name")
        elif others and n_query == 0:
            verdict = ("reference-lineage-specific duplication with "
                       + ", ".join(others)
                       + " — no query gene separates them, so these names do "
                         "not designate distinct query subfamilies")
        elif others:
            verdict = ("shares its smallest clade with " + ", ".join(others)
                       + " — the split is not resolved, treat the names as one group")
        elif n_query == 0:
            verdict = "no query gene in the anchor's clade — nothing to name"
        else:
            verdict = f"transferable — {n_query} query genes in the labelled clade"
        rows.append({
            "label": label,
            "n_anchors": len(own),
            "clade_size": len(clade),
            "n_query_in_clade": n_query,
            "support": raw_support,
            "blocked_by": sorted(set(others) | set(dup_with)),
            "transferable": transferable,
            "verdict": verdict,
            "members": sorted(clade),
        })

    # anchor-free clades: the lineage-specific-expansion case
    def maximal_anchor_free(node) -> List[set]:
        L = set(leaves_of(node))
        if L & anchors:
            out: List[set] = []
            for c in node.children:
                out += maximal_anchor_free(c)
            return out
        return [L] if len(L) > 1 else []

    for clade in maximal_anchor_free(root):
        n_query = sum(1 for x in clade if is_query(x))
        if n_query < 2:
            continue
        rows.append({
            "label": None,
            "n_anchors": 0,
            "clade_size": len(clade),
            "n_query_in_clade": n_query,
            "support": None,
            "blocked_by": [],
            "transferable": False,
            "verdict": ("query-only clade — lineage-specific expansion, no "
                        "reference name applies"),
            "members": sorted(clade),
        })
    return rows
