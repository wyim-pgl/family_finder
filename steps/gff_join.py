"""Joining a family's gene ids onto the GFF models that describe them.

Split out of `steps.gene_structure` so that module stays within the project's
file-size limit; the two are one story. The hand-rolled version of this join
lost 7 of 77 PEPC ids to `.t1` / `.v1.0` / `.EL10_1.0` / `evm.model`-vs-
`evm.TU` differences and said nothing about it, which is indistinguishable
from a family with fewer members (#42). Everything therefore goes through
`utils.gene_ids.match_ids`, the one matcher, and the loss is reported and
can be capped.
"""

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

from utils.gene_ids import canon_gene_id, match_ids, squash_gene_id, strip_isoform

from steps.gene_structure_model import GeneModel

_CANON = canon_gene_id


def _ISOFORM(gene_id: str) -> str:
    return strip_isoform(canon_gene_id(gene_id))


_SQUASH = squash_gene_id


def _aliases(models: Sequence[GeneModel]) -> Dict[str, GeneModel]:
    """Every string a GFF gene can be looked up by.

    A gene is addressable by its own id and by its transcript's — pep files
    carry variously either (`Dcar_evm.model.*` is a transcript id, Ococ's are
    gene ids). The first model to claim a spelling keeps it.
    """
    alias_gene: Dict[str, GeneModel] = {}
    for model in models:
        for alias in (model.gene_id, model.transcript_id):
            if alias and alias not in alias_gene:
                alias_gene[alias] = model
    return alias_gene


def _candidates(alias_gene: Dict[str, GeneModel],
                queries: Sequence[str]) -> Dict[str, GeneModel]:
    """Narrow a whole genome's aliases to the ones a query could reach.

    A GFF holds tens of thousands of genes and a family holds tens. Indexing
    all of them makes every unrelated pair of ids a possible collision — that
    is how `Dcar_evm.model.Hap1Chr1.244` was lost: it squashes to the same
    string as `evm.model.Hap1Chr12.44`, a gene nobody asked about. Only
    aliases a query actually names are considered.
    """
    keys = set()
    for query in queries:
        forms = [query, *_dot_trims(query)]
        for form in forms:
            keys |= {form, _CANON(form), _ISOFORM(form)}
    return {alias: model for alias, model in alias_gene.items()
            if {alias, _CANON(alias), _ISOFORM(alias)} & keys}


def _colliding(aliases: Sequence[str],
               alias_gene: Dict[str, GeneModel]) -> set:
    """Aliases that make one normalisation level ambiguous, deepest first.

    Returned only when `match_ids` has actually refused the reference set:
    an alias is never given up for a collision at a level the join did not
    need to reach.
    """
    for transform in (_SQUASH, _ISOFORM, _CANON):
        by_key: Dict[str, List[str]] = defaultdict(list)
        for alias in sorted(aliases):
            by_key[transform(alias)].append(alias)
        bad = set()
        for group in by_key.values():
            if len(group) < 2:
                continue
            genes = {alias_gene[a].gene_id for a in group}
            keep = min(group, key=len) if len(genes) == 1 else None
            bad |= {a for a in group if a != keep}
        if bad:
            return bad
    return set()


def _match_with_guard(queries: Sequence[str],
                      alias_gene: Dict[str, GeneModel]):
    """`match_ids` against the alias set, giving up ambiguous aliases lazily.

    Returns the match, how many aliases had to be given up, and the alias set
    that survived — the dot-trim retry must use the SAME set, or it walks back
    into the collision that has just been resolved.
    """
    aliases = set(alias_gene)
    dropped = 0
    while True:
        try:
            return match_ids(queries, sorted(aliases)), dropped, sorted(aliases)
        except ValueError:
            bad = _colliding(sorted(aliases), alias_gene)
            if not bad:
                raise
            aliases -= bad
            dropped += len(bad)


def _dot_trims(gene_id: str) -> List[str]:
    """`AH000245.v2.1` -> `AH000245.v2`, `AH000245`.

    pep ids carry an assembly or release suffix the GFF id does not
    (`.v2.1`, `.EL10_1.0`). Trimming trailing dot-fields is tried only after
    `match_ids` has failed on its own levels, so the milder normalisations
    always win first.
    """
    parts = gene_id.split(".")
    return [".".join(parts[:i]) for i in range(len(parts) - 1, 0, -1)]


def strip_species(gene_id: str) -> str:
    """`Bvul_EL10Ac4g08333.EL10_1.0` -> `EL10Ac4g08333.EL10_1.0`.

    Project rule: the species is the prefix before the first `_`. A leading
    underscore survives in reconstructed ids (`Cjam__JBOBQE...`) and is
    dropped, because the GFF spells the contig without it.
    """
    _species, _, rest = gene_id.partition("_")
    return rest.lstrip("_")


def join_gff_models(gene_ids: Sequence[str], models: Sequence[GeneModel],
                    max_unmatched: Optional[float] = None
                    ) -> Tuple[Dict[str, GeneModel], dict]:
    """Map family gene ids onto GFF models through `utils.gene_ids.match_ids`.

    The hand-rolled version of this join lost 7 of 77 PEPC ids to `.t1` /
    `.v1.0` / `.EL10_1.0` / `evm.model`-vs-`evm.TU` differences, and said
    nothing about it — a shrinking gene set looks exactly like a family that
    has fewer members. Everything goes through the one matcher, and
    `max_unmatched` turns the loss into a refusal.
    """
    queries = {gene: strip_species(gene) for gene in gene_ids}
    alias_gene = _candidates(_aliases(models), list(queries.values()))
    matched, n_ambiguous, references = _match_with_guard(
        list(queries.values()), alias_gene)
    mapping = {gene: alias_gene[matched.mapping[key]]
               for gene, key in queries.items() if key in matched.mapping}
    level = matched.level

    unmatched = [gene for gene in gene_ids if gene not in mapping]
    if unmatched and references:
        recovered = {}
        for gene in unmatched:
            for trimmed in _dot_trims(queries[gene]):
                hit = match_ids([trimmed], references).mapping.get(trimmed)
                if hit is not None:
                    recovered[gene] = alias_gene[hit]
                    break
        if recovered:
            mapping.update(recovered)
            level = "dot-trimmed"
            unmatched = [gene for gene in gene_ids if gene not in mapping]

    fraction = len(unmatched) / len(gene_ids) if gene_ids else 0.0
    report = {
        "level": level,
        "n_requested": len(gene_ids),
        "n_matched": len(mapping),
        "n_unmatched": len(unmatched),
        "unmatched": unmatched,
        "unmatched_fraction": fraction,
        "n_ambiguous_aliases": n_ambiguous,
    }
    if max_unmatched is not None and fraction > max_unmatched:
        raise ValueError(
            f"{len(unmatched)} of {len(gene_ids)} gene ids stayed unmatched "
            f"against the GFF models ({fraction:.0%} > {max_unmatched:.0%}): "
            f"{unmatched[:5]} — refusing rather than reporting a structure "
            "axis over a silently reduced gene set"
        )
    return mapping, report


def join_by_species(gene_ids: Sequence[str],
                    models_by_species: Dict[str, Sequence[GeneModel]],
                    max_unmatched: Optional[float] = None
                    ) -> Tuple[Dict[str, GeneModel], dict]:
    """Join one species at a time, then report the overall loss.

    Pooling every species' GFF into one reference set would let two genomes
    that happen to spell a locus the same way (`g11025` is an AUGUSTUS id, not
    a unique one) collide, and the ambiguity filter would then drop BOTH — a
    coverage loss caused entirely by pooling. The species prefix is already
    part of every gene id, so the pooling is unnecessary.
    """
    by_species: Dict[str, List[str]] = defaultdict(list)
    for gene in gene_ids:
        by_species[gene.split("_", 1)[0]].append(gene)

    mapping: Dict[str, GeneModel] = {}
    per_species: Dict[str, dict] = {}
    for species in sorted(by_species):
        models = models_by_species.get(species)
        if not models:
            per_species[species] = {
                "level": "no GFF", "n_requested": len(by_species[species]),
                "n_matched": 0, "n_unmatched": len(by_species[species]),
                "unmatched": list(by_species[species]),
                "unmatched_fraction": 1.0, "n_ambiguous_aliases": 0,
            }
            continue
        hits, report = join_gff_models(by_species[species], models)
        mapping.update(hits)
        per_species[species] = report

    unmatched = [gene for gene in gene_ids if gene not in mapping]
    fraction = len(unmatched) / len(gene_ids) if gene_ids else 0.0
    overall = {
        "n_requested": len(gene_ids),
        "n_matched": len(mapping),
        "n_unmatched": len(unmatched),
        "unmatched": unmatched,
        "unmatched_fraction": fraction,
        "species_without_gff": sorted(s for s, r in per_species.items()
                                      if r["level"] == "no GFF"),
        "by_species": per_species,
    }
    if max_unmatched is not None and fraction > max_unmatched:
        raise ValueError(
            f"{len(unmatched)} of {len(gene_ids)} gene ids stayed unmatched "
            f"against the GFF models ({fraction:.0%} > {max_unmatched:.0%}): "
            f"{unmatched[:5]} — refusing rather than reporting a structure "
            "axis over a silently reduced gene set"
        )
    return mapping, overall


