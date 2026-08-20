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
"""

import logging
from collections import Counter
from typing import Dict, Iterable, List, Optional, Set, Tuple

logger = logging.getLogger("family_finder")

# Isoform suffixes tried when adding/stripping variants, most common first.
_ISOFORM_SUFFIXES = [".t1", ".1", ".t2", ".2"]


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
        })
    if unmatched:
        logger.warning(
            f"{len(unmatched)} annotated gene(s) not found in any family "
            f"(first: {unmatched[:3]})"
        )
    return rows, unmatched


def name_groups(rows: List[dict], key: str = "family_id") -> List[dict]:
    """Name each group by its most common member description.

    support = fraction of annotated members carrying the winning
    description — a 0.5 support name is a coin flip, report it as such.
    Rows with an empty group value (e.g. genes without a subfamily label
    when key="subfamily") are skipped.
    """
    per_group: Dict[str, List[str]] = {}
    for row in rows:
        group = row.get(key, "")
        if not group:
            continue
        per_group.setdefault(group, []).append(row["description"])

    named: List[dict] = []
    for group_id, descriptions in sorted(per_group.items()):
        top, count = Counter(descriptions).most_common(1)[0]
        named.append({
            "group_id": group_id,
            "name": top,
            "support": count / len(descriptions),
            "n_annotated": len(descriptions),
        })
    return named
