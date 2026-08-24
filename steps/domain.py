"""Pfam domain architecture as a scalable source of NEGATIVES (issue #47).

The tree can only say "these are one family" when it is given something known
to lie OUTSIDE the family. PEPC was the lucky case: the literature had already
split plant-type from bacterial-type PEPC, and `anchors.tsv` already labelled
ATH_AT1G68750.1 bacterial. Scored on 15-species data the sign flipped exactly
once across 128 genes, so PPC4 was a defensible negative and the five PTPC
fragments merged into one 111-gene family. The other 23,743 families have no
such gift, and drawing negatives from below-cut vote edges does not work -
that is precisely where membership is ambiguous.

Pfam supplies negatives for every family at once, and supplies them from
curation rather than from the graph whose boundaries are in question: two
proteins built from different domain architectures are different families.

Everything here is written around one hazard. This repository's signature
failure is the truncated or mis-annotated model, and a lost domain looks
exactly like a different architecture. So no function compares two genes. A
family's signature contains only the domains a MAJORITY of its members carry,
and two families are incompatible only when a domain is in one signature and
absent from the other - a short model cannot move a majority by itself.

The rule is deliberately one-directional, like the structural channel:

    different architecture  =>  different family        (usable)
    same architecture       =>  nothing                 (many families share
                                                         one architecture)

So this yields outgroups, never merges.
"""
from collections import Counter
from pathlib import Path
from typing import Dict, List, Sequence, Set

# domtblout columns. The Pfam identifier is the ACCESSION, not the name:
# names are editable, PF00311 is not, and a renamed domain must not read as a
# different architecture.
#   hmmscan   gene = query name 3,   domain = target accession 1
#   hmmsearch gene = target name 0,  domain = query accession 4
# ali from/to (17, 18) are on the sequence in both layouts.
_MIN_FIELDS = 19
_IEVALUE, _ALI_FROM, _ALI_TO = 12, 17, 18
_SCAN_GENE, _SCAN_ACC = 3, 1
_SEARCH_GENE, _SEARCH_ACC = 0, 4

DEFAULT_MAX_EVALUE = 1e-5
MAJORITY = 0.5


def _strip_version(acc: str) -> str:
    """PF00311.12 -> PF00311. A version bump is not a new domain."""
    return acc.split(".", 1)[0]


def parse_pfam_domtblout(path, max_evalue: float = DEFAULT_MAX_EVALUE,
                         swap: bool = False) -> Dict[str, List[str]]:
    """gene_id -> domain accessions, ordered along the protein.

    hmmscan emits rows grouped by profile, so they arrive in arbitrary
    positional order; an architecture read in that order would be noise.

    Overlapping hits are collapsed to the strongest: two profiles matching the
    same span describe one region, and counting both would invent a domain the
    protein does not have. Genuine tandem repeats survive, because they do not
    overlap.

    `swap` reads hmmsearch layout, where query and target are reversed.
    """
    q_col, t_col = ((_SEARCH_GENE, _SEARCH_ACC) if swap
                    else (_SCAN_GENE, _SCAN_ACC))
    hits: Dict[str, list] = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < _MIN_FIELDS:
                continue
            try:
                evalue = float(f[_IEVALUE])
                start, end = int(f[_ALI_FROM]), int(f[_ALI_TO])
            except ValueError:
                continue
            if evalue > max_evalue:
                continue
            hits.setdefault(f[q_col], []).append(
                (start, end, evalue, _strip_version(f[t_col])))

    out: Dict[str, List[str]] = {}
    for gene, rows in hits.items():
        kept: List[tuple] = []
        for start, end, evalue, acc in sorted(rows, key=lambda r: r[2]):
            if any(start <= k_end and k_start <= end
                   for k_start, k_end, _, _ in kept):
                continue
            kept.append((start, end, evalue, acc))
        out[gene] = [acc for _, _, _, acc in sorted(kept)]
    return out


def architecture_signature(
    genes: Sequence[str], architecture: Dict[str, List[str]],
) -> Set[str]:
    """The domains a majority of a family's members carry.

    A majority, not a union: one truncated model that lost its domains, or one
    mis-annotated model that gained a spurious one, must not define what the
    family is made of. Genes with no domain hits count in the denominator -
    a family that is mostly fragments should not get a confident signature.
    """
    if not genes:
        return set()
    counts = Counter()
    for gene in genes:
        counts.update(set(architecture.get(gene, ())))
    threshold = len(genes) * MAJORITY
    return {acc for acc, n in counts.items() if n > threshold}


def incompatible_families(family: str,
                          signatures: Dict[str, Set[str]]) -> List[str]:
    """Families whose signature is disjoint from this one.

    Disjoint, not merely different. Sharing a domain says related, which is
    not a claim this makes; only complete architectural disjunction is taken
    as "outside the family". An empty signature - a family with no majority
    domain, usually one built largely of fragments - never makes anything
    incompatible in either direction: absence of domain evidence must not
    become evidence of absence.
    """
    own = signatures.get(family) or set()
    if not own:
        return []
    return sorted(
        other for other, sig in signatures.items()
        if other != family and sig and not (sig & own)
    )


def pick_domain_outgroup(cluster: Sequence[str],
                         signatures: Dict[str, Set[str]],
                         proximity: Dict[str, float],
                         n: int = 3) -> List[str]:
    """The nearest architecturally incompatible families outside the cluster.

    Nearest, because an outgroup should sit just outside. A distant panel is
    monophyletic with respect to everything it brackets, and a monophyletic
    outgroup always returns "one family" - correct but tautological. Proximity
    is whatever the caller measured (vote fraction, bit score); this only
    orders by it.

    Cluster members are excluded even when incompatible: an outgroup drawn
    from inside the group under test is circular. If that empties the result,
    the honest answer is no outgroup rather than a convenient one.
    """
    members = set(cluster)
    candidates: Set[str] = set()
    for family in cluster:
        candidates.update(incompatible_families(family, signatures))
    candidates -= members
    return sorted(candidates,
                  key=lambda f: (-proximity.get(f, 0.0), f))[:n]
