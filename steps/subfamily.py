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
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from utils.gene_ids import match_ids
from utils.newick import Node

logger = logging.getLogger("family_finder")

# Ranks checked lowest-to-highest for the lineage-specific verdict.
TAXONOMY_RANKS = ["species", "genus", "family", "order"]


# ---------------------------------------------------------------------------
# Group id resolution
# ---------------------------------------------------------------------------

def resolve_groups(
    alignment: Dict[str, str],
    groups: Dict[str, List[str]],
    *,
    max_unmatched: Optional[float] = None,
) -> Tuple[Dict[str, List[str]], dict]:
    """Re-express each group's members in the alignment's own identifiers.

    Every layer here used to keep members with `m in alignment` — an exact
    string test that drops anything spelled differently and says nothing. That
    is how `Ococ_OcoChr10G09070.t1` failed to meet `groups.json` (#34): the
    subfamily was scanned with fewer members than it has, or with none, and the
    result read as "no diagnostic residues" rather than "never examined".

    Matching goes through `utils.gene_ids.match_ids`, the single normaliser, so
    a fix there fixes every consumer at once — including its collision guard,
    which refuses to loosen far enough to merge `_1338_000001` with
    `_1338_000002`. `max_unmatched` (a fraction of all group members) makes a
    coverage loss raise instead of shrinking the analysis quietly.

    Returns the resolved groups and a report: the level `match_ids` settled on,
    the overall unmatched count and fraction, and the unmatched ids per group.
    """
    return resolve_group_ids(list(alignment), groups,
                             max_unmatched=max_unmatched)


def resolve_group_ids(
    reference_ids: Sequence[str],
    groups: Dict[str, List[str]],
    *,
    max_unmatched: Optional[float] = None,
) -> Tuple[Dict[str, List[str]], dict]:
    """`resolve_groups` against any id universe, not only an alignment.

    The other universe that matters is the foldseek pair table, which is keyed
    on structure filenames rather than alignment names.
    """
    members = [m for group in sorted(groups) for m in groups[group]]
    matched = match_ids(members, list(reference_ids), max_unmatched=max_unmatched)

    resolved: Dict[str, List[str]] = {}
    per_group: Dict[str, dict] = {}
    for name, group_members in groups.items():
        hits = [matched.mapping[m] for m in group_members if m in matched.mapping]
        # A group may legitimately list the same gene twice; the alignment
        # cannot, so de-duplicate while keeping the caller's order.
        seen = set()
        resolved[name] = [g for g in hits if not (g in seen or seen.add(g))]
        per_group[name] = {
            "n_members": len(resolved[name]),
            "n_requested": len(group_members),
            "unmatched": [m for m in group_members if m not in matched.mapping],
        }

    report = {
        "level": matched.level,
        "n_requested": len(members),
        "n_unmatched": len(matched.unmatched),
        "unmatched_fraction": matched.unmatched_fraction,
        "groups": per_group,
    }
    return resolved, report


# ---------------------------------------------------------------------------
# Coordinate reference
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ReferenceChoice:
    """Which sequence a set of residue positions is numbered against."""
    seq_id: str
    source: str  # "explicit" | "characterised" | "automatic"
    n_characterised: int = 0
    unmatched_characterised: List[str] = field(default_factory=list)
    reason: str = ""


def resolve_reference(
    alignment: Dict[str, str],
    *,
    ref_seq_id: Optional[str] = None,
    characterised: Optional[Sequence[str]] = None,
) -> ReferenceChoice:
    """Pick the sequence diagnostic positions are cited in, and say why.

    A position is only citable against something. `utils.alignment.
    choose_reference` already answers the reference-free case — the most
    complete sequence in the family — but when the caller knows which members
    are characterised (a SwissProt-backed anchor, an enzyme with published
    residue numbers) those coordinates are the ones a reader can look up, so
    they win. The automatic representative is the fallback, never a silent
    one: `source` records which of the three routes was taken.

    Ids are matched through `utils.gene_ids.match_ids`, so a characterised
    anchor spelled with a transcript suffix the alignment does not carry is
    still found instead of dropping to the fallback unnoticed.
    """
    from utils.alignment import choose_reference

    if not alignment:
        raise ValueError("Cannot choose a coordinate reference from an "
                         "empty alignment")
    names = list(alignment)

    if ref_seq_id is not None:
        hit = match_ids([ref_seq_id], names).mapping.get(ref_seq_id)
        if hit is None:
            raise ValueError(
                f"Reference sequence {ref_seq_id!r} is absent from the "
                "alignment even after id normalisation — refusing to number "
                "positions against a sequence that is not there"
            )
        return ReferenceChoice(hit, "explicit",
                               reason=f"caller named {ref_seq_id!r}")

    wanted = list(characterised or [])
    if wanted:
        matched = match_ids(wanted, names)
        pool = sorted(set(matched.mapping.values()))
        if pool:
            return ReferenceChoice(
                choose_reference(alignment, candidates=pool),
                "characterised",
                n_characterised=len(pool),
                unmatched_characterised=list(matched.unmatched),
                reason=(f"most complete of {len(pool)} characterised member(s) "
                        f"matched at the {matched.level!r} level"),
            )
        logger.warning(
            "None of the %d characterised id(s) matched an alignment name "
            "even at the %r level — numbering falls back to the automatic "
            "representative: %s",
            len(wanted), matched.level, matched.unmatched[:5],
        )
        return ReferenceChoice(
            choose_reference(alignment), "automatic",
            unmatched_characterised=list(matched.unmatched),
            reason=(f"none of {len(wanted)} characterised id(s) is in the "
                    "alignment; most complete sequence used instead"),
        )

    return ReferenceChoice(
        choose_reference(alignment), "automatic",
        reason="no characterised member supplied; most complete sequence used",
    )


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
    characterised: Optional[Sequence[str]] = None,
    max_unmatched: Optional[float] = None,
) -> List[dict]:
    """Subfamily-diagnostic alignment columns.

    A column is diagnostic for group G when >= in_cons of G's non-gap
    members share one residue, G's non-gap coverage is >= min_cover, and
    that residue's frequency among the REST of the alignment's non-gap
    residues is <= out_max. Groups smaller than min_group are skipped
    (their "conservation" is not meaningful).

    Every row carries `ref_pos` — an ungapped residue number — and `ref_seq`,
    the sequence that number is in. `ref_pos` alone was never citable: the
    same integer means a different residue in every sequence of the family.
    The reference is `ref_seq_id` when given, otherwise the most complete of
    `characterised` (e.g. ATH_PPC1 and other SwissProt-backed anchors, so
    positions can be compared with published ones), otherwise the family's own
    representative — see `resolve_reference`.

    Member ids are resolved against the alignment through `resolve_groups`,
    not by exact string equality — pass `max_unmatched` to refuse rather than
    scan a group the id forms have quietly shrunk.
    """
    from utils.alignment import reference_positions

    lengths = {len(s) for s in alignment.values()}
    if len(lengths) > 1:
        raise ValueError(f"ragged alignment: lengths {sorted(lengths)}")
    alen = lengths.pop() if lengths else 0

    groups, id_report = resolve_groups(alignment, groups,
                                       max_unmatched=max_unmatched)
    if id_report["n_unmatched"]:
        logger.warning(
            "sdp_scan: %d of %d group members did not match an alignment id "
            "even at the %r level — those genes were not scanned",
            id_report["n_unmatched"], id_report["n_requested"],
            id_report["level"],
        )

    reference = resolve_reference(alignment, ref_seq_id=ref_seq_id,
                                  characterised=characterised)
    ref_pos = reference_positions(alignment[reference.seq_id])
    logger.info("sdp_scan: positions numbered against %s (%s — %s)",
                reference.seq_id, reference.source, reference.reason)

    rows: List[dict] = []
    for group_id, members in sorted(groups.items()):
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
                "ref_seq": reference.seq_id,
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
                        min_group: int = 5,
                        max_unmatched: Optional[float] = None) -> Dict[str, dict]:
    """Which columns `sdp_scan` never judged for each group, and why.

    `sdp_scan` skips a column for a group whose non-gap coverage falls below
    `min_cover`, and skips a group smaller than `min_group` outright. Both are
    correct — conservation over three of seventy sequences is not conservation
    — but both are silent, so "no diagnostic columns here" is indistinguishable
    from "this group was never examined here". Issue #40 hit exactly that: the
    ppc-1E1 subfamily covers the N-terminal regulatory region in only 33 of 70
    members, so every N-terminal column was filtered out and the region read as
    signal-free.

    Call this next to `sdp_scan` and report it alongside the hits. It resolves
    member ids the same way `sdp_scan` does, so the two agree on which genes a
    group actually has in the alignment.
    """
    lengths = {len(s) for s in alignment.values()}
    if len(lengths) > 1:
        raise ValueError(f"ragged alignment: lengths {sorted(lengths)}")
    alen = lengths.pop() if lengths else 0

    groups, _id_report = resolve_groups(alignment, groups,
                                        max_unmatched=max_unmatched)
    report: Dict[str, dict] = {}
    for group_id, present in sorted(groups.items()):
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
    ref_seq_id: Optional[str] = None,
    characterised: Optional[Sequence[str]] = None,
    distance_space: str = "residue",
    max_unmatched: Optional[float] = None,
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

    **Distance is in residues, not columns.** An alignment column is a
    bookkeeping slot: a 200-column insertion carried by one clade pushes the
    core and a diagnostic residue 200 columns apart in a protein where nothing
    moved. Both the observed distance and the null inherit that inflation, and
    they inherit it unevenly, because the null samples columns wherever the
    gaps happen to be. Distances are therefore measured in the ungapped
    numbering of one reference sequence (`resolve_reference`), and
    `distance_convention` records which sequence that is. `distance_space=
    "column"` restores the old convention, labelled as such. Columns the
    reference cannot place — it is gapped there — are counted in
    `n_*_untranslatable` rather than folded in at some arbitrary position.

    Returns a verdict of `avoids_core`, `at_core`, `indistinguishable`, or
    `no_interpretation_available` with a `reason` — silence and "no signal"
    have to be distinguishable (issue #40).
    """
    from utils.alignment import (group_occupancy, invariant_columns,
                                 invariant_suppressed, reference_positions)

    if distance_space not in ("residue", "column"):
        raise ValueError(
            f"unknown distance_space {distance_space!r}: use 'residue' (the "
            "reference's own numbering) or 'column' (alignment columns)"
        )

    lengths = {len(s) for s in alignment.values()}
    if len(lengths) > 1:
        raise ValueError(f"ragged alignment: lengths {sorted(lengths)}")
    alen = lengths.pop() if lengths else 0

    # Resolve ids once, so the scan, the invariant core and the candidate pool
    # all agree on which genes each group has. They used to disagree whenever
    # an id form differed (#34).
    groups, _id_report = resolve_groups(alignment, groups,
                                        max_unmatched=max_unmatched)
    skipped = sorted(name for name, members in groups.items()
                     if len(members) < min_group)

    # Which sequence the positions below are numbered in. Resolved before the
    # scan so the scan's ref_pos column and these distances agree.
    reference = (resolve_reference(alignment, ref_seq_id=ref_seq_id,
                                   characterised=characterised)
                 if distance_space == "residue" else None)

    hits = sdp_scan(alignment, groups, min_group=min_group, in_cons=in_cons,
                    out_max=out_max, min_cover=min_cover,
                    ref_seq_id=ref_seq_id, characterised=characterised)
    sdp_cols = sorted({h["aln_col"] for h in hits})
    invariant = invariant_columns(alignment, groups,
                                  threshold=invariant_threshold,
                                  min_cover=min_cover)
    inv_cols = sorted(invariant)

    # Columns min_cover kept the core from being judged on at all. Without
    # this a gappy region reads as "the groups do not agree here" when no
    # comparison was ever made (the same silence coverage_suppressed removed
    # from sdp_scan).
    core_suppressed = invariant_suppressed(alignment, groups,
                                           min_cover=min_cover)

    # Columns sdp_scan was able to judge at all: every group covers them.
    occ = group_occupancy(alignment, groups)
    candidates = [c + 1 for c in range(alen)
                  if occ and all(v[c] >= min_cover for v in occ.values())]

    result = {
        "n_sdp_columns": len(sdp_cols),
        "n_invariant_columns": len(inv_cols),
        "n_candidate_columns": len(candidates),
        "n_core_suppressed_columns": core_suppressed["n_suppressed"],
        "n_core_examined_columns": core_suppressed["n_examined"],
        "core_suppressed_columns": core_suppressed["columns"],
        "core_suppressed_by_group": {name: len(cols) for name, cols
                                     in core_suppressed["by_group"].items()},
        "skipped_groups": skipped,
        "overlap": len(set(sdp_cols) & set(inv_cols)),
        "overlap_is_by_construction": True,
        "distance_space": distance_space,
        "distance_convention": ("column" if reference is None
                                else f"residue:{reference.seq_id}"),
        "distance_reference": None if reference is None else reference.seq_id,
        "distance_reference_source": (None if reference is None
                                      else reference.source),
        "n_sdp_untranslatable": 0,
        "n_invariant_untranslatable": 0,
        "n_candidate_untranslatable": 0,
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
            f"measure distance from. {core_suppressed['n_suppressed']} of "
            f"{alen} column(s) were never examined for it at all: some group "
            f"covers them below min_cover={min_cover}"
        )
        return result

    ref_pos = ({} if reference is None
               else reference_positions(alignment[reference.seq_id]))

    def place(columns: Sequence[int]) -> Tuple[List[int], List[int]]:
        """Columns as reference residue numbers, plus the ones it cannot place."""
        if reference is None:
            return list(columns), []
        return ([ref_pos[c] for c in columns if c in ref_pos],
                [c for c in columns if c not in ref_pos])

    sdp_pos, sdp_lost = place(sdp_cols)
    inv_pos, inv_lost = place(inv_cols)
    cand_pos, cand_lost = place(candidates)
    result["n_sdp_untranslatable"] = len(sdp_lost)
    result["n_invariant_untranslatable"] = len(inv_lost)
    result["n_candidate_untranslatable"] = len(cand_lost)

    if not sdp_pos or not inv_pos:
        result["reason"] = (
            f"the reference {reference.seq_id!r} is gapped at "
            f"{len(sdp_lost)} of {len(sdp_cols)} diagnostic and "
            f"{len(inv_lost)} of {len(inv_cols)} invariant column(s), so they "
            "have no residue number in it; pick a reference that covers them "
            "or measure in alignment columns"
        )
        return result
    if len(cand_pos) < len(sdp_pos):
        result["reason"] = (
            f"only {len(cand_pos)} column(s) were scannable, fewer than the "
            f"{len(sdp_pos)} diagnostic column(s); no null can be built"
        )
        return result

    def nearest(pos: int) -> int:
        return min(abs(pos - j) for j in inv_pos)

    observed = statistics.median(nearest(c) for c in sdp_pos)
    rng = random.Random(seed)
    null = sorted(statistics.median(nearest(c) for c in
                                    rng.sample(cand_pos, len(sdp_pos)))
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
        f"null drawn from {len(candidates)} scannable columns, measured in "
        f"{result['distance_convention']}; "
        f"{core_suppressed['n_suppressed']} column(s) were never examined for "
        f"the core (min_cover={min_cover})"
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


def _family_span_depth(per_group: Dict[str, List[str]],
                       species_tree: Node) -> Optional[int]:
    """Depth of the MRCA of every subfamily's species, or None with no groups.

    This is the reference the root-span rule compares against: a subfamily
    whose own MRCA is that same node did not arise inside one lineage.
    """
    if not per_group:
        return None
    union = set().union(*per_group.values())
    return species_monophyly(union, species_tree)["mrca_depth"]


def _tree_attribution(species: List[str], species_tree: Node,
                      family_depth: Optional[int],
                      taxonomy: Dict[str, Dict[str, str]]) -> Tuple[dict, str, List[str]]:
    """Species-tree columns, verdict and notes for one subfamily."""
    mono = species_monophyly(set(species), species_tree)
    columns = {
        "n_in_tree": mono["n_in_tree"],
        "monophyletic": mono["monophyletic"],
        "mrca_name": mono["mrca_name"],
        "mrca_depth": mono["mrca_depth"],
        "clade_label": clade_rank_label(set(species), taxonomy),
    }
    if mono["monophyletic"] is None:
        verdict = "unknown (no species in species tree)"
    elif not mono["monophyletic"]:
        verdict = "paralog-split (non-monophyletic)"
    elif mono["mrca_depth"] == family_depth:
        # Ancient paralogs kept in every sampled species are trivially
        # monophyletic; the root-span rule is what stops them being read as
        # lineage-specific.
        verdict = "paralog-split (spans family root)"
    else:
        verdict = "lineage-specific (clade)"

    notes = []
    if mono["missing"]:
        notes.append("not in species tree: " + ",".join(mono["missing"]))
    if mono["intruders"]:
        notes.append("interleaved species: " + ",".join(mono["intruders"]))
    if taxonomy:
        columns.update(_rank_purities(species, taxonomy))
    return columns, verdict, notes


def _legacy_attribution(species: List[str],
                        taxonomy: Dict[str, Dict[str, str]]) -> Tuple[dict, str]:
    """Rank-purity columns and verdict when there is no species tree."""
    columns = _rank_purities(species, taxonomy)
    for rank in TAXONOMY_RANKS:  # lowest pure rank wins
        if columns[f"{rank}_purity"] == 1.0:
            return columns, f"lineage-specific ({rank})"
    return columns, "paralog-split (crosses order)"


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
    family_depth = (_family_span_depth(per_group, species_tree)
                    if species_tree is not None else None)

    rows: List[dict] = []
    for group_id, members in sorted(groups.items()):
        species = per_group[group_id]
        row = {"subfamily": group_id, "n_members": len(members),
               "n_species": len(set(species))}

        if species_tree is not None:
            columns, verdict, notes = _tree_attribution(
                species, species_tree, family_depth, taxonomy)
        else:
            columns, verdict = _legacy_attribution(species, taxonomy)
            notes = []
        row.update(columns)

        unknown = sorted({s for s in species if s not in taxonomy})
        if unknown and (species_tree is None or taxonomy):
            notes.append("unknown taxonomy for: " + ",".join(unknown))
        row["verdict"] = verdict
        row["notes"] = "; ".join(notes)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Structural coherence (foldseek all-vs-all)
# ---------------------------------------------------------------------------

# The layout steps.esm.run_foldseek asks for. A run that also wants the
# sequence-identity control must add `fident` (issue #39 item 3).
DEFAULT_PAIR_COLUMNS = "query,target,evalue,bits,alntmscore"


def parse_pairwise_table(tsv: Path,
                         columns: Optional[str] = None) -> Dict[frozenset, dict]:
    """Foldseek all-vs-all TSV -> one record per pair, keyed by the id pair.

    `columns` is the foldseek `--format-output` string the file was written
    with (default `DEFAULT_PAIR_COLUMNS`). It is checked against the actual
    field count rather than guessed at: a table read under the wrong layout
    silently returns numbers from the wrong column, which is the same class of
    failure as an unverified coordinate system.

    Self-hits are dropped. When both search directions report a pair, the
    higher-scoring alignment wins AND ITS OWN other fields are kept — mixing
    the bits of one alignment with the fident of another would break the very
    comparison the identity control needs.
    """
    names = [c.strip() for c in (columns or DEFAULT_PAIR_COLUMNS).split(",")]
    for required in ("query", "target"):
        if required not in names:
            raise ValueError(f"Pair-table columns must include {required!r}: {names}")
    rank_on = "bits" if "bits" in names else names[-1]

    table: Dict[frozenset, dict] = {}
    with open(tsv) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(names):
                raise ValueError(
                    f"{tsv}: expected {len(names)} columns ({','.join(names)}) "
                    f"but the table has {len(fields)} — pass the "
                    "--format-output string this file was written with rather "
                    "than reading numbers out of the wrong column"
                )
            record = dict(zip(names, fields))
            if record["query"] == record["target"]:
                continue
            try:
                for name in names[2:]:
                    record[name] = float(record[name])
            except ValueError:
                logger.debug(f"Skipping malformed pair line: {line[:80]!r}")
                continue
            key = frozenset((record["query"], record["target"]))
            prev = table.get(key)
            if prev is None or record[rank_on] > prev[rank_on]:
                table[key] = record
    return table


def _pair_column(tsv: Path, column: str,
                 columns: Optional[str] = None) -> Dict[frozenset, float]:
    names = [c.strip() for c in (columns or DEFAULT_PAIR_COLUMNS).split(",")]
    if column not in names:
        raise ValueError(
            f"{column!r} is not in the pair table's columns ({','.join(names)}). "
            f"Re-run foldseek with {column} in --format-output, or the number "
            "it would produce does not exist in this file."
        )
    return {key: rec[column]
            for key, rec in parse_pairwise_table(tsv, columns).items()}


def parse_pairwise_scores(tsv: Path, metric: str = "bits",
                          columns: Optional[str] = None) -> Dict[frozenset, float]:
    """Foldseek all-vs-all TSV -> symmetric best score per pair.

    `metric` names the column to score on. `bits` grows with both alignment
    length and identity; `alntmscore` is length-normalised and therefore the
    less confounded of the two, but neither is free of sequence similarity —
    see `structure_coherence`.
    """
    return _pair_column(tsv, metric, columns)


def parse_pair_identities(tsv: Path,
                          columns: Optional[str] = None) -> Dict[frozenset, float]:
    """Per-pair sequence identity (`fident`), for the identity control."""
    return _pair_column(tsv, "fident", columns)


def _fit_line(xs: List[float], ys: List[float]) -> Optional[Tuple[float, float]]:
    """Ordinary least squares `y = slope * x + intercept`, or None.

    None when there are too few points or x has no spread — with a constant
    x there is nothing to regress out, and pretending otherwise would report
    the raw difference as if it had been controlled.
    """
    n = len(xs)
    if n < 3:
        return None
    mean_x, mean_y = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    if sxx == 0.0:
        return None
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    slope = sxy / sxx
    return slope, mean_y - slope * mean_x


def sequence_confounding(pair_scores: Dict[frozenset, float],
                         identities: Dict[frozenset, float]) -> dict:
    """How tightly the structural score tracks sequence identity.

    This is the number that decides whether the coherence result says anything
    beyond "close relatives look alike". Report it next to the coherence table,
    not instead of it: a high `r` does not make the structural signal absent,
    it makes the UNCONTROLLED comparison uninformative.
    """
    shared = sorted(set(pair_scores) & set(identities), key=lambda k: sorted(k))
    xs = [identities[k] for k in shared]
    ys = [pair_scores[k] for k in shared]
    fit = _fit_line(xs, ys)
    out = {"n_pairs": len(shared), "r": None, "slope": None, "intercept": None}
    if fit is None:
        return out
    slope, intercept = fit
    n = len(xs)
    mean_x, mean_y = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    out.update(slope=slope, intercept=intercept,
               r=(sxy / (sxx * syy) ** 0.5) if sxx and syy else None)
    return out


IDENTITY_BIN_WIDTH = 0.05
MIN_SHARED_BINS = 3
MIN_WITHIN_IN_SHARED_BINS = 5
# A binned residual mean this many standard errors from zero is a systematic
# departure, not sampling noise. Three sigma is the conventional bar and was
# not chosen to make any particular dataset fail — on the real PEPC table the
# departure is an order of magnitude past it.
MAX_RESIDUAL_Z = 3.0


def _identity_bin(fident: float) -> int:
    return int(fident / IDENTITY_BIN_WIDTH)


def residual_diagnostics(pair_scores: Dict[frozenset, float],
                         identities: Dict[frozenset, float],
                         n_bins: int = 10) -> dict:
    """Is a straight line an adequate model of score-versus-identity?

    The identity control regresses the structural score on `fident` and
    compares residuals. That is only a control if the line actually fits: a
    misspecified model leaves residuals that depend on identity, and since
    within-subfamily pairs are concentrated at high identity, the
    within-versus-between residual difference then measures the
    misspecification rather than the structures.

    Measured on the real PEPC pair table (5,343 pairs) the line does NOT fit:
    r2 = 0.430 for bits and 0.497 for alntmscore, with binned residual means
    running +0.006 / -0.061 / -0.068 / +0.074 / +0.053 / +0.046 / +0.035 /
    +0.017 / -0.014 / -0.088 instead of hovering at zero. alntmscore saturates
    at 1.0 while fident keeps climbing, so every high-identity pair sits below
    the line whatever its structure.

    Residuals are binned on identity into `n_bins` equal-count bins; each
    bin's mean is scored against its own standard error. Any bin past
    `MAX_RESIDUAL_Z` marks the fit inadequate, and the caller must then refuse
    to draw a conclusion from the residuals rather than report the artifact.
    """
    shared = sorted(set(pair_scores) & set(identities), key=lambda k: sorted(k))
    xs = [identities[k] for k in shared]
    ys = [pair_scores[k] for k in shared]
    out = {"n_pairs": len(shared), "r2": None, "max_abs_z": None,
           "bin_residual_means": [], "bin_sizes": [],
           "linear_fit_adequate": False,
           "reason": "too few pairs, or no spread in identity, to fit a line"}
    fit = _fit_line(xs, ys)
    if fit is None:
        return out

    slope, intercept = fit
    residuals = [y - (slope * x + intercept) for x, y in zip(xs, ys)]
    mean_y = sum(ys) / len(ys)
    ss_tot = sum((y - mean_y) ** 2 for y in ys)
    ss_res = sum(r * r for r in residuals)
    out["r2"] = 1.0 - ss_res / ss_tot if ss_tot else None

    order = sorted(range(len(xs)), key=lambda i: xs[i])
    per_bin = max(1, len(order) // n_bins)
    bins = [order[i:i + per_bin] for i in range(0, len(order), per_bin)]
    # Chunking leaves a remainder; a three-point bin has a huge standard error
    # and would be judged on noise, so fold it into its neighbour.
    if len(bins) > n_bins:
        bins[-2].extend(bins.pop())
    sd = (ss_res / len(residuals)) ** 0.5
    # An exact fit leaves residuals at floating-point noise, and dividing that
    # noise by its own standard error produces meaningless z-scores. When the
    # residual spread is negligible against the spread of the response there
    # is no structure left to test for.
    y_sd = (ss_tot / len(ys)) ** 0.5 if ss_tot else 0.0
    if sd <= 1e-9 * y_sd:
        out.update(linear_fit_adequate=True, max_abs_z=0.0,
                   reason="the line fits exactly; no residual structure remains")
        return out
    means, sizes, zs = [], [], []
    for members in bins:
        if not members:
            continue
        m = sum(residuals[i] for i in members) / len(members)
        means.append(m)
        sizes.append(len(members))
        zs.append(abs(m) / (sd / len(members) ** 0.5) if sd else 0.0)
    out.update(bin_residual_means=means, bin_sizes=sizes,
               max_abs_z=max(zs) if zs else None)

    if not zs or max(zs) <= MAX_RESIDUAL_Z:
        out.update(linear_fit_adequate=True,
                   reason="binned residuals are consistent with zero, so the "
                          "linear identity control is usable")
        return out
    out["reason"] = (
        f"the line does not fit: r2 = {out['r2']:.3f} and the binned residual "
        f"means are systematic (non-monotonic, worst bin {max(zs):.1f} standard "
        "errors from zero). Residuals from a misspecified fit measure the "
        "misspecification, not the structures"
    )
    return out


UNCONTROLLED_WARNING = (
    "UNCONTROLLED: this ratio is not adjusted for sequence identity. Foldseek "
    "scores rise with sequence similarity, so within > between can mean no "
    "more than 'closer relatives look more alike'. Re-run foldseek with "
    "fident in --format-output and pass identities= to compare "
    "identity-adjusted residuals instead."
)


def _mean(values: List[float]) -> Optional[float]:
    return sum(values) / len(values) if values else None


def structure_coherence(
    groups: Dict[str, List[str]],
    pair_scores: Dict[frozenset, float],
    identities: Optional[Dict[frozenset, float]] = None,
    metric: str = "bits",
) -> List[dict]:
    """Mean within-subfamily vs between-subfamily structural score.

    Only OBSERVED pairs count — foldseek omits pairs below its reporting
    threshold, and inventing zeros for them would fake incoherence.
    mean_between is None (and coherent is None) when no cross-subfamily
    pair was observed for the group.

    **The raw comparison is confounded** (issue #39). Foldseek scores track
    sequence similarity closely, so "within > between" is also what you get
    when subfamilies are simply groups of close relatives — which they are by
    construction. Measured on the PEPC clan the raw rule called 5 of 6
    subfamilies coherent (ratio median 1.24, only OG4 below at 0.915), and
    that number cannot be cited as structural evidence on its own.

    Pass per-pair `fident` as `identities` to control for it: the score is
    regressed on identity across every observed pair and the within/between
    comparison is repeated on the RESIDUALS.

    **The residual comparison is not automatically citable either, and on the
    PEPC data it is not citable at all.** Two things have to be checked first,
    and both fail there:

    * *Does the line fit?* r2 = 0.430 for bits and 0.497 for alntmscore over
      5,343 pairs, with binned residual means running +0.006 / -0.061 /
      -0.068 / +0.074 / ... rather than sitting at zero. alntmscore saturates
      at 1.0 while fident keeps climbing, so high-identity pairs fall below
      the line whatever their structure — and within-pairs are almost all
      high-identity. A naive reading of the residuals then calls every
      subfamily incoherent, which is the misspecification talking.
      (`residual_diagnostics`)
    * *Do within- and between-pairs overlap in identity at all?* Binned at
      0.05: OG15 2 shared bins (3 within-pairs), OG2 2 (3), OG3 3 (3),
      OG4 11 (55), OG5 1 (1), **OG8 0**. Outside OG4 there is essentially
      nothing to compare, so the honest answer is not "no coherence" but
      "cannot tell".

    So each row carries a `verdict` of `coherent`, `not_coherent` or
    `no_interpretation_available` with a `reason`, following the same
    discipline as `sdp_core_relationship`. `coherent_controlled` is None
    whenever the verdict is `no_interpretation_available`; the diagnostics
    (`fit_r2`, `linear_fit_adequate`, `n_shared_identity_bins`) travel with
    every row so the table can be judged without rerunning anything. Without
    identities at all, every row is `no_interpretation_available` with an
    explicit UNCONTROLLED warning — emitting the bare ratio silently is what
    this issue exists to stop.
    """
    # Group ids arrive in alignment spelling and the pair table is keyed on
    # structure filenames ('.' -> '_'). Comparing them literally observes no
    # pairs at all and reports it as "none observed" (#42), so resolve through
    # the shared matcher and count what still has no structure.
    structure_ids = sorted({g for pair in pair_scores for g in pair})
    groups, id_report = resolve_group_ids(structure_ids, groups)
    member_of = {}
    for group_id, members in groups.items():
        for m in members:
            if m in member_of and member_of[m] != group_id:
                raise ValueError(
                    f"gene {m!r} belongs to both {member_of[m]!r} and "
                    f"{group_id!r} — overlapping groups would silently "
                    "assign it to whichever was processed last"
                )
            member_of[m] = group_id

    fit, diagnostics = None, None
    if identities:
        shared = sorted(set(pair_scores) & set(identities), key=lambda k: sorted(k))
        fit = _fit_line([identities[k] for k in shared],
                        [pair_scores[k] for k in shared])
        diagnostics = residual_diagnostics(pair_scores, identities)
    residual = {}
    if fit is not None:
        slope, intercept = fit
        for pair, score in pair_scores.items():
            if pair in identities:
                residual[pair] = score - (slope * identities[pair] + intercept)

    rows: List[dict] = []
    for group_id, members in sorted(groups.items()):
        within, between = [], []
        within_res, between_res = [], []
        within_bins: Dict[int, int] = {}
        between_bins: Dict[int, int] = {}
        for pair, score in pair_scores.items():
            a, b = tuple(pair)
            ga, gb = member_of.get(a), member_of.get(b)
            if group_id not in (ga, gb):
                continue
            if ga == gb == group_id:
                bucket, res_bucket, bins = within, within_res, within_bins
            elif ga is not None and gb is not None:
                bucket, res_bucket, bins = between, between_res, between_bins
            else:
                continue
            bucket.append(score)
            if identities and pair in identities:
                b_idx = _identity_bin(identities[pair])
                bins[b_idx] = bins.get(b_idx, 0) + 1
            if pair in residual:
                res_bucket.append(residual[pair])

        mean_w, mean_b = _mean(within), _mean(between)
        res_w, res_b = _mean(within_res), _mean(between_res)
        # Residuals existing on both sides means the control RAN, not that it
        # HELD. `sequence_controlled` below reports the latter.
        has_residuals = res_w is not None and res_b is not None
        shared_bins = sorted(set(within_bins) & set(between_bins))
        support = {
            "n_shared_identity_bins": len(shared_bins),
            "n_within_pairs_in_shared_bins": sum(within_bins[b] for b in shared_bins),
            "n_between_pairs_in_shared_bins": sum(between_bins[b] for b in shared_bins),
        }
        verdict, reason = _controlled_verdict(
            has_residuals, res_w, res_b, diagnostics, support, identities, fit)
        # One meaning, one column: the control held iff a verdict came out of
        # it. Reporting True next to an UNCONTROLLED reason let anyone
        # filtering on this column tally these rows as controlled results.
        decided = verdict in ("coherent", "not_coherent")
        rows.append({
            "subfamily": group_id,
            "metric": metric,
            "n_within_pairs": len(within),
            "n_between_pairs": len(between),
            "mean_within": mean_w,
            "mean_between": mean_b,
            "coherent": (mean_w > mean_b) if (mean_w is not None
                                              and mean_b is not None) else None,
            "sequence_controlled": decided,
            "mean_within_residual": res_w,
            "mean_between_residual": res_b,
            "coherent_controlled": (res_w > res_b) if decided else None,
            "verdict": verdict,
            "reason": reason,
            "fit_r2": (diagnostics or {}).get("r2"),
            "linear_fit_adequate": (diagnostics or {}).get("linear_fit_adequate"),
            **support,
            "n_members_without_structure": len(
                id_report["groups"][group_id]["unmatched"]),
            "warning": "" if decided else reason,
        })
    return rows


def _controlled_verdict(has_residuals: bool, res_w: Optional[float],
                        res_b: Optional[float], diagnostics: Optional[dict],
                        support: dict, identities, fit) -> Tuple[str, str]:
    """`coherent` / `not_coherent` / `no_interpretation_available` + why.

    Three things have to hold before the residual comparison may conclude
    anything, and each of them failed somewhere on the real PEPC data:

    1. the control ran at all (identities present, a line fittable, residuals
       on both sides of this subfamily's comparison);
    2. the line FITS — otherwise the residual difference is the
       misspecification (r2 = 0.43/0.50 measured, binned residuals systematic);
    3. within- and between-pairs actually OVERLAP in identity — otherwise the
       comparison extrapolates. Measured per subfamily with 0.05-wide bins:
       OG15 2 shared bins (3 within pairs), OG2 2 (3), OG3 3 (3), OG4 11 (55),
       OG5 1 (1), OG8 **0**. Only OG4 has anything worth calling support.

    The verdict follows the same discipline as `sdp_core_relationship`: an
    unsupported comparison is `no_interpretation_available` with a reason, not
    a quiet negative.
    """
    if not has_residuals:
        return "no_interpretation_available", _uncontrolled_reason(identities, fit)
    if diagnostics and not diagnostics["linear_fit_adequate"]:
        return "no_interpretation_available", (
            "UNCONTROLLED: the identity control cannot be trusted — "
            + diagnostics["reason"]
            + ". Within-pairs cluster at high identity, so the residual "
              "difference reproduces the curvature rather than the structures"
        )
    if (support["n_shared_identity_bins"] < MIN_SHARED_BINS
            or support["n_within_pairs_in_shared_bins"] < MIN_WITHIN_IN_SHARED_BINS):
        return "no_interpretation_available", (
            "within and between pairs do not share an identity range: "
            f"{support['n_shared_identity_bins']} shared "
            f"{IDENTITY_BIN_WIDTH:g}-wide identity bin(s) holding "
            f"{support['n_within_pairs_in_shared_bins']} within-pair(s) "
            f"(need {MIN_SHARED_BINS} bins and {MIN_WITHIN_IN_SHARED_BINS} "
            "pairs). Comparing residuals across disjoint identity ranges "
            "extrapolates the fit instead of controlling with it"
        )
    if res_w > res_b:
        return "coherent", (
            "identity-adjusted: within-pair residuals exceed between-pair "
            "residuals over a shared identity range, with an adequate fit"
        )
    return "not_coherent", (
        "identity-adjusted: within-pair residuals do NOT exceed between-pair "
        "residuals over a shared identity range, with an adequate fit"
    )


def _uncontrolled_reason(identities: Optional[Dict[frozenset, float]],
                         fit) -> str:
    """Why the identity control could not be applied — never just silence."""
    if not identities:
        return UNCONTROLLED_WARNING
    if fit is None:
        return (
            "UNCONTROLLED: sequence identities were supplied but could not be "
            "regressed out — fewer than three shared pairs, or every pair has "
            "the same identity, so there is no spread to adjust for. "
            + UNCONTROLLED_WARNING.split(". ", 1)[1]
        )
    return (
        "UNCONTROLLED for this subfamily: the regression was fitted, but this "
        "group has no residual on one side of the comparison (no observed "
        "within- or between-pair carrying an identity). "
        + UNCONTROLLED_WARNING.split(". ", 1)[1]
    )


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


# --- how confident is one subfamily assignment? ----------------------------
#
# Conventions ADOPTED HERE, not derived from these data. SH-aLRT >= 80 and
# UFboot >= 95 are the bars IQ-TREE's own manual recommends for calling a
# clade well supported; UFboot >= 70 is the usual "worth reporting, not
# established" band. Two independent datasets recovering the same membership
# is our own minimum, taken from the #40 rebuild where amino acid, codon 1+2
# and codon 1+2+3 were run separately. Nothing below is a measured property
# of the PEPC clan; change the numbers here and say so when you do.
GRADE_MIN_SH_ALRT = 80
GRADE_HIGH_MIN_UFBOOT = 95
GRADE_PROVISIONAL_MIN_UFBOOT = 70
GRADE_MIN_CONSISTENT_DATASETS = 2

GRADE_HIGH = "HIGH"
GRADE_PROVISIONAL = "PROVISIONAL"
GRADE_UNRESOLVED = "UNRESOLVED"


def _support_pair(
        raw: Optional[str]) -> Tuple[Optional[float], Optional[float]]:
    """Node label -> (SH-aLRT, UFboot).

    IQ-TREE writes "SH-aLRT/UFBoot" ("96.2/94") with `-alrt -bb`, and
    "SH-aLRT/aBayes/UFBoot" when `-abayes` is on too, so the first and last
    numbers are the two we grade on. A single number is not self-describing:
    it is held to BOTH bars rather than assumed to be the lenient one.
    """
    if not raw:
        return None, None
    vals = []
    for part in str(raw).split("/"):
        try:
            vals.append(float(part))
        except ValueError:
            continue
    if not vals:
        return None, None
    return vals[0], vals[-1]


def grade_assignment(sh_alrt: Optional[float], ufboot: Optional[float],
                     n_consistent: Optional[int] = None,
                     consistent: Optional[List[str]] = None,
                     inconsistent: Optional[List[str]] = None
                     ) -> Tuple[str, str]:
    """(grade, reason) for one subfamily assignment. See the constants above.

    The WEAKER support decides. Each of the two supports is held to its own
    bar and the one that fails is named in the reason — a strong SH-aLRT
    never buys a weak UFboot through, and the two are never averaged
    (mean(100, 90) = 95 would clear the HIGH bar that neither 90 nor the pair
    deserves).

    Cross-tree consistency is EVIDENCE. `n_consistent=None` means no
    independent dataset was supplied, which is not agreement and not
    disagreement: it caps the grade at PROVISIONAL and says in the reason
    that the condition was never evaluated. Only a consistency that was
    measured and came out short reads as a failure.
    """
    named = ""
    if inconsistent:
        named = " — " + ", ".join(inconsistent) + " recovered other members"
    if sh_alrt is None or ufboot is None:
        return GRADE_UNRESOLVED, ("no support value on the clade — nothing to "
                                  "grade the assignment with")
    if sh_alrt < GRADE_MIN_SH_ALRT:
        return GRADE_UNRESOLVED, (
            f"SH-aLRT {sh_alrt:g} is below {GRADE_MIN_SH_ALRT:g} "
            f"(UFboot {ufboot:g}) — the weaker support decides")
    if ufboot < GRADE_PROVISIONAL_MIN_UFBOOT:
        return GRADE_UNRESOLVED, (
            f"UFboot {ufboot:g} is below {GRADE_PROVISIONAL_MIN_UFBOOT:g} "
            f"(SH-aLRT {sh_alrt:g}) — the weaker support decides")
    if (n_consistent is not None
            and n_consistent < GRADE_MIN_CONSISTENT_DATASETS):
        return GRADE_UNRESOLVED, (
            f"cross-tree consistency was evaluated and FAILED: {n_consistent} "
            f"dataset(s) recovered this membership, need "
            f"{GRADE_MIN_CONSISTENT_DATASETS}" + named)
    if ufboot >= GRADE_HIGH_MIN_UFBOOT and n_consistent is not None:
        return GRADE_HIGH, (
            f"SH-aLRT {sh_alrt:g} and UFboot {ufboot:g} both clear "
            f"{GRADE_MIN_SH_ALRT:g}/{GRADE_HIGH_MIN_UFBOOT:g}, and "
            f"{n_consistent} datasets recover the same membership"
            + (" (" + ", ".join(consistent) + ")" if consistent else ""))
    caveats = []
    if ufboot < GRADE_HIGH_MIN_UFBOOT:
        caveats.append(
            f"UFboot {ufboot:g} is below {GRADE_HIGH_MIN_UFBOOT:g} "
            f"(SH-aLRT {sh_alrt:g}) — the weaker support decides")
    if n_consistent is None:
        caveats.append(
            "cross-tree consistency was NOT evaluated — no independent "
            "dataset was supplied, so the same membership has never been "
            "recovered twice. HIGH is unreachable without that check")
    return GRADE_PROVISIONAL, "; ".join(caveats)


def _cross_tree_agreement(
        label: str, clade: set, tree_name: str,
        cross_tree_members: Optional[Dict[str, Dict[str, object]]]):
    """(n_consistent, agreeing, disagreeing) for one label.

    Agreement is exact membership equality: the same leaves, recovered from a
    dataset that was not this one. The tree being analysed counts as the
    first dataset, so one agreeing partner is enough to reach
    GRADE_MIN_CONSISTENT_DATASETS. A dataset that does not carry the label at
    all did not recover the membership either, and is listed as disagreeing.
    """
    if not cross_tree_members:
        return None, [], []
    agree, disagree = [tree_name], []
    for name in cross_tree_members:      # caller's order, not ours
        other = cross_tree_members[name].get(label)
        if other is not None and set(other) == clade:
            agree.append(name)
        else:
            disagree.append(name)
    return len(agree), agree, disagree


def anchor_transferability(
    newick: str,
    anchor_labels: Dict[str, str],
    query_prefixes: Optional[List[str]] = None,
    min_support: Optional[float] = None,
    cross_tree_members: Optional[Dict[str, Dict[str, object]]] = None,
    tree_name: str = "this tree",
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

    Every row also carries a `grade` — HIGH / PROVISIONAL / UNRESOLVED, see
    `grade_assignment` — because `transferable` alone cannot tell a clade at
    100/98 from one at 100/76 and the manuscript has to say which it has.
    `cross_tree_members` is the consistency evidence: ``{dataset name: {label:
    members}}`` from INDEPENDENT datasets (amino acid, codon 1+2, codon
    1+2+3). Supplying nothing leaves that condition unevaluated rather than
    satisfied, and the grade is capped at PROVISIONAL with that said in
    `grade_reason`.
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

        # The grade is about how well established the assignment is, so a
        # label with no assignment at all — blocked by a rival, or naming an
        # empty query set — is UNRESOLVED whatever the node support says.
        sh_alrt, ufboot = _support_pair(raw_support)
        n_consistent, agree, disagree = _cross_tree_agreement(
            label, clade, tree_name, cross_tree_members)
        if others or dup_with or n_query == 0:
            grade, grade_reason = GRADE_UNRESOLVED, (
                "no assignment to grade — " + verdict)
        else:
            grade, grade_reason = grade_assignment(
                sh_alrt, ufboot, n_consistent, agree, disagree)
        rows.append({
            "label": label,
            "n_anchors": len(own),
            "clade_size": len(clade),
            "n_query_in_clade": n_query,
            "support": raw_support,
            "sh_alrt": sh_alrt,
            "ufboot": ufboot,
            "blocked_by": sorted(set(others) | set(dup_with)),
            "transferable": transferable,
            "grade": grade,
            "grade_reason": grade_reason,
            "cross_tree_evaluated": n_consistent is not None,
            "n_consistent_datasets": n_consistent,
            "consistent_datasets": agree,
            "inconsistent_datasets": disagree,
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
            "sh_alrt": None,
            "ufboot": None,
            "blocked_by": [],
            "transferable": False,
            "grade": GRADE_UNRESOLVED,
            "grade_reason": ("no assignment to grade — no reference name "
                             "applies to this clade"),
            "cross_tree_evaluated": False,
            "n_consistent_datasets": None,
            "consistent_datasets": [],
            "inconsistent_datasets": [],
            "verdict": ("query-only clade — lineage-specific expansion, no "
                        "reference name applies"),
            "members": sorted(clade),
        })
    return rows
