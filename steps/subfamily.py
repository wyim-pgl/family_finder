"""Subfamily diagnostics for ANY family (generalized from the PEPC pilot).

Three independent, pure layers — callers supply plain dicts so everything
tests without ete4 or external tools:

  * sdp_scan: specificity-determining positions — alignment columns where
    one subfamily is conserved on a residue the rest of the family rarely
    uses. This is the "WHY does it split" evidence at residue level.
  * taxonomic_composition: is each subfamily a lineage-specific block
    (taxonomy-driven — species/genus/family-level expansion) or does it
    span orders (an ancient paralog split that exists "by itself")?
  * structure_coherence: given foldseek all-vs-all scores, are subfamilies
    also structurally coherent (mean within-score > mean between-score)?

Foldseek's role: the pipeline already uses foldseek for tier-3 assignment
(gene vs family representative, steps/esm.py) and clan-merge evidence.
Here it is reused a third way — all-vs-all WITHIN a family to ask whether
the sequence-defined subfamilies are visible in structure space at all.
"""

import logging
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

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


def taxonomic_composition(
    groups: Dict[str, List[str]],
    taxonomy: Dict[str, Dict[str, str]],
    delimiter: str = "_",
) -> List[dict]:
    """Per subfamily: does the split follow taxonomy, or cross it?

    For each rank (species < genus < family < order) the purity is the
    largest fraction of members belonging to one taxon. The verdict is
    "lineage-specific (<rank>)" for the LOWEST rank that is pure — i.e.
    the split is explainable as a taxonomic block at that rank — and
    "paralog-split (crosses order)" when even the order rank is mixed:
    the subfamily exists across lineages, so the division is a property
    of the gene family itself, not of taxonomy.

    Species absent from the taxonomy table still count (species-level
    identity comes from the gene-id prefix); higher ranks treat them as
    their own singleton taxa and the row notes them.
    """
    rows: List[dict] = []
    for group_id, members in sorted(groups.items()):
        species = [m.split(delimiter, 1)[0] for m in members]
        unknown = sorted({s for s in species if s not in taxonomy})

        row: dict = {
            "subfamily": group_id,
            "n_members": len(members),
            "n_species": len(set(species)),
        }
        verdict = "paralog-split (crosses order)"
        for rank in TAXONOMY_RANKS:
            if rank == "species":
                taxa = species
            else:
                taxa = [taxonomy.get(s, {}).get(rank, f"unknown:{s}")
                        for s in species]
            counts = Counter(taxa)
            dominant, dom_n = counts.most_common(1)[0]
            purity = dom_n / len(taxa)
            row[f"{rank}_purity"] = round(purity, 3)
            row[f"{rank}_dominant"] = dominant
        for rank in TAXONOMY_RANKS:  # lowest pure rank wins
            if row[f"{rank}_purity"] == 1.0:
                verdict = f"lineage-specific ({rank})"
                break
        row["verdict"] = verdict
        row["notes"] = (
            f"unknown taxonomy for: {','.join(unknown)}" if unknown else ""
        )
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
