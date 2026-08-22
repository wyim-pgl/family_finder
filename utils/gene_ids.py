"""One canonical form for gene identifiers, and one place that matches them.

Issue #42. Three separate silent coverage losses in a single day traced back
here. Structure filenames replace '.' with '_', so `Obas__JBFLFP010000003.1_000519`
becomes `Obas__JBFLFP010000003_1_000519` and stops matching the alignment.
Transcript suffixes appear on one side and not the other, so
`Ococ_OcoChr10G09070.t1` misses `Ococ_OcoChr10G09070`. And every module carries
its own normalisation, so what one fixes another does not.

None of it errors. A mismatch just drops genes: 111 structures once matched 29
alignment names, and rerunning region_disorder on the PEPC clan lost 16 of 111
the same way. Shrinking coverage then looks exactly like an absent signal.

`match_ids` therefore reports what it did — which normalisation level was
needed, what stayed unmatched, and what fraction that is — and can be given a
ceiling above which it refuses rather than quietly returning less.
"""
import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional

_EXTENSIONS = (".pdb", ".fa", ".fasta", ".afa", ".aln")
# `.t1` / `_t1` / `.1` are isoform or version suffixes. A bare `_000001` is
# NOT: gene ids like Cgig_..._SGP5p_1338_000001 carry the gene number after an
# underscore, and stripping it merges distinct loci. The collision guard in
# _index caught exactly that on the real PEPC structure set.
_ISOFORM = re.compile(r"(?:[._]t\d+|\.\d+)$")


def canon_gene_id(gene_id: str) -> str:
    """Normalise separators and drop a file extension. Idempotent."""
    out = gene_id.strip()
    for ext in _EXTENSIONS:
        if out.endswith(ext):
            out = out[: -len(ext)]
            break
    return out.replace(".", "_")


def strip_isoform(gene_id: str) -> str:
    """Remove one trailing transcript suffix (`.t1`, `_t1`, `.1`).

    Only one, and only at the end: a locus name may legitimately end in digits
    (`OcoChr10G09070`), and stripping those would merge distinct genes.
    """
    return _ISOFORM.sub("", gene_id.strip(), count=1)


@dataclass(frozen=True)
class IdMatch:
    mapping: Dict[str, str]
    level: str
    unmatched: List[str] = field(default_factory=list)
    unmatched_fraction: float = 0.0


def _index(reference: Iterable[str], transform) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for ref in reference:
        key = transform(ref)
        if key in out and out[key] != ref:
            raise ValueError(
                f"Reference ids {out[key]!r} and {ref!r} collide under "
                f"normalisation (both become {key!r}) — the mapping would be "
                "ambiguous"
            )
        out[key] = ref
    return out


def match_ids(queries: Iterable[str], reference: Iterable[str],
              max_unmatched: Optional[float] = None) -> IdMatch:
    """Map query ids onto reference ids, loosening only as far as needed.

    Tries exact identity first, then canonical form, then canonical form with
    the transcript suffix removed, and reports which level it settled on. Pass
    `max_unmatched` (a fraction) to make excessive loss raise instead of
    silently reducing coverage.
    """
    queries = list(queries)
    reference = list(reference)
    levels = (
        ("exact", lambda s: s),
        ("canonical", canon_gene_id),
        ("isoform", lambda s: strip_isoform(canon_gene_id(s))),
    )

    best: Optional[IdMatch] = None
    for level, transform in levels:
        index = _index(reference, transform)
        mapping, unmatched = {}, []
        for q in queries:
            hit = index.get(transform(q))
            if hit is None:
                unmatched.append(q)
            else:
                mapping[q] = hit
        frac = len(unmatched) / len(queries) if queries else 0.0
        best = IdMatch(mapping, level, unmatched, frac)
        if not unmatched:
            break

    assert best is not None
    if max_unmatched is not None and best.unmatched_fraction > max_unmatched:
        raise ValueError(
            f"{len(best.unmatched)} of {len(queries)} ids stayed unmatched "
            f"({best.unmatched_fraction:.0%} > {max_unmatched:.0%}) even at the "
            f"{best.level!r} level: {best.unmatched[:5]} — refusing rather than "
            "returning reduced coverage that reads as an absent signal"
        )
    return best
