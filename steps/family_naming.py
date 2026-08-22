"""Name pipeline families from ONE annotated species (iceplant pattern).

The Mcry/iceplant workflow: a single species has curated functional
annotations (gene_id + description); every confirmed family containing an
annotated gene inherits a name from its members' descriptions, and each
annotated gene gets its family (and subfamily) id back. Pure dict layers;
the CLI lives in name_families.py.

Gene-id matching is deliberately forgiving: annotation tables often lack
the pipeline's SpeciesPrefix_ convention and disagree on isoform suffixes
(.t1 / .1). match_gene tries exact, prefixed, and suffix variants in a
deterministic order and reports what it matched — never a silent guess.

Issue #30 adds a second evidence layer: PLAZA-orthology transfer. Families
without a direct annotation inherit ATH ortholog symbols (map_orthology),
every row carries its provenance (source: direct | orthology), and
name_groups weights the majority vote — a direct annotation outweighs an
orthology transfer (DEFAULT_WEIGHTS).
"""

import logging
from collections import Counter
from typing import Dict, Iterable, List, Optional, Set, Tuple

logger = logging.getLogger("family_finder")

# Isoform suffixes tried when adding/stripping variants, most common first.
_ISOFORM_SUFFIXES = [".t1", ".1", ".t2", ".2"]

# Majority-vote weights per evidence source (issue #30): a curated direct
# annotation counts double an orthology-transferred description.
DEFAULT_WEIGHTS = {"direct": 1.0, "orthology": 0.5}


def index_families(families: Dict[str, Set[str]]) -> Dict[str, str]:
    """gene_id -> family_id over all confirmed families."""
    index: Dict[str, str] = {}
    for family_id, members in families.items():
        for gene in members:
            index[gene] = family_id
    return index


def _candidates(gene: str, species: Optional[str]) -> Iterable[str]:
    """Deterministic id variants: as-is, prefixed, +/- isoform suffixes."""
    bases = [gene]
    if species and not gene.startswith(f"{species}_"):
        bases.append(f"{species}_{gene}")
    for base in bases:
        yield base
        for suffix in _ISOFORM_SUFFIXES:      # annotation locus, pipeline isoform
            yield base + suffix
        for suffix in _ISOFORM_SUFFIXES:      # annotation isoform, pipeline locus
            if base.endswith(suffix):
                yield base[: -len(suffix)]


def match_gene(
    gene: str, index: Dict[str, str], species: Optional[str] = None
) -> Optional[Tuple[str, str]]:
    """Resolve an annotation gene id to (pipeline_gene_id, family_id)."""
    for candidate in _candidates(gene, species):
        family = index.get(candidate)
        if family is not None:
            return candidate, family
    return None


def map_annotations(
    annotations: Dict[str, str],
    families: Dict[str, Set[str]],
    groups: Optional[Dict[str, str]] = None,
    species: Optional[str] = None,
) -> Tuple[List[dict], List[str]]:
    """Per annotated gene: its family (and subfamily) in the pipeline run.

    groups maps pipeline_gene_id -> subfamily label (Possvm / TreeCluster
    output); genes without a label get "". Returns (rows, unmatched_ids).
    """
    index = index_families(families)
    groups = groups or {}
    rows: List[dict] = []
    unmatched: List[str] = []
    for gene, description in sorted(annotations.items()):
        hit = match_gene(gene, index, species)
        if hit is None:
            unmatched.append(gene)
            continue
        pipeline_id, family_id = hit
        rows.append({
            "annotation_gene": gene,
            "pipeline_gene": pipeline_id,
            "family_id": family_id,
            "subfamily": groups.get(pipeline_id, ""),
            "description": description,
            "source": "direct",
            "ortholog": "",
        })
    if unmatched:
        logger.warning(
            f"{len(unmatched)} annotated gene(s) not found in any family "
            f"(first: {unmatched[:3]})"
        )
    return rows, unmatched


def map_orthology(
    orthology: Dict[str, Tuple[str, str]],
    families: Dict[str, Set[str]],
    groups: Optional[Dict[str, str]] = None,
    species: Optional[str] = None,
) -> Tuple[List[dict], List[str]]:
    """Per gene with an ATH ortholog: its family/subfamily, source=orthology.

    orthology maps gene_id -> (ortholog_id, description), e.g. the output
    of extract_plaza_orthology.py. Same row shape as map_annotations so
    both layers feed the same name_groups vote.
    """
    index = index_families(families)
    groups = groups or {}
    rows: List[dict] = []
    unmatched: List[str] = []
    for gene, (ortholog, description) in sorted(orthology.items()):
        hit = match_gene(gene, index, species)
        if hit is None:
            unmatched.append(gene)
            continue
        pipeline_id, family_id = hit
        rows.append({
            "annotation_gene": gene,
            "pipeline_gene": pipeline_id,
            "family_id": family_id,
            "subfamily": groups.get(pipeline_id, ""),
            "description": description,
            "source": "orthology",
            "ortholog": ortholog,
        })
    if unmatched:
        logger.info(
            f"{len(unmatched)} ortholog-mapped gene(s) not found in any "
            f"family (first: {unmatched[:3]})"
        )
    return rows, unmatched


def combine_sources(
    direct_rows: List[dict], ortho_rows: List[dict]
) -> List[dict]:
    """Merge the two evidence layers; per gene, direct suppresses orthology.

    A gene annotated directly already has curated evidence — its orthology
    transfer would only echo (or contradict) it with weaker support, so it
    is dropped rather than double-counted.
    """
    directly_annotated = {row["pipeline_gene"] for row in direct_rows}
    kept = [
        row for row in ortho_rows
        if row["pipeline_gene"] not in directly_annotated
    ]
    return list(direct_rows) + kept


def name_groups(
    rows: List[dict],
    key: str = "family_id",
    weights: Optional[Dict[str, float]] = None,
    group_sizes: Optional[Dict[str, int]] = None,
) -> List[dict]:
    """Name each group by its weighted-majority member description.

    Each row votes with the weight of its evidence source (rows without a
    "source" are direct — the pre-#30 single-species layer). support =
    winning weight / total weight — a 0.5 support name is a coin flip,
    report it as such.

    ⚠️ support is computed over the members that CARRY an annotation, so on
    its own it cannot distinguish "the clade agrees" from "the two members
    who voted agreed". Pass `group_sizes` (group_id -> total members) to get
    `coverage` = annotated / total alongside it. Measured on the PEPC clan:
    naming subfamilies by Arabidopsis symbol gave OG15 (20 members) support
    1.00 from 2 votes — coverage 0.10 is the number that says so. Ties break toward a direct-backed description, then
    lexicographically (deterministic output). Rows with an empty group
    value (e.g. genes without a subfamily label when key="subfamily") are
    skipped. provenance says what carried the winning name: direct when
    any direct row voted for it, orthology otherwise.
    """
    weights = weights if weights is not None else DEFAULT_WEIGHTS
    per_group: Dict[str, List[dict]] = {}
    for row in rows:
        group = row.get(key, "")
        if not group:
            continue
        per_group.setdefault(group, []).append(row)

    named: List[dict] = []
    for group_id, members in sorted(per_group.items()):
        weight_of: Counter = Counter()
        direct_backed: Set[str] = set()
        n_direct = n_orthology = 0
        for row in members:
            source = row.get("source", "direct")
            weight_of[row["description"]] += weights.get(source, 1.0)
            if source == "orthology":
                n_orthology += 1
            else:
                n_direct += 1
                direct_backed.add(row["description"])
        top = min(
            weight_of,
            key=lambda d: (-weight_of[d], d not in direct_backed, d),
        )
        total = sum(weight_of.values())
        size = (group_sizes or {}).get(group_id)
        coverage = min(1.0, len(members) / size) if size else None
        named.append({
            "group_id": group_id,
            "name": top,
            "support": weight_of[top] / total,
            "coverage": coverage,
            "n_annotated": len(members),
            "group_size": size,
            "n_direct": n_direct,
            "n_orthology": n_orthology,
            "provenance": "direct" if top in direct_backed else "orthology",
        })
    return named
