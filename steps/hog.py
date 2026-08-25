"""OrthoFinder HOG levels, exposed as the SUBFAMILY axis (issue #47).

OrthoFinder writes one hierarchical-orthogroup table per internal node of the
species tree - N0 (root) through N13 in the 15-species round-1 run. It is
tempting to read that as the family-rank hierarchy this pipeline lacks, so
that "family" could be a declared cut rather than a per-cluster judgement
needing an outgroup.

Measured on that run, it is not:

    OG 28,150  ->  N0 HOG 28,161      (more, not fewer)
    OGs split into >1 HOG at N0:  3,153
    HOGs holding genes from >1 OG:    0   <- decisive

A HOG never crosses an orthogroup boundary, so the hierarchy SUBDIVIDES and
can never merge - and merging is what fragmentation needs. It is also written
per ROUND: the Ococ PEPC flagship sits in round 1's N0 while the Mcry flagship
is a round-2 gene absent from that table entirely, so no level can hold both.
Fragmentation and the hierarchy are split along the same axis.

What the levels are exactly right for is the rank BELOW family. Once a family
is settled - by an outgroup, as PTPC was - they give nested subfamily
partitions at every species-tree node, already computed, with depth as a free
parameter and no new compute at all.

The levels are NOT a uniform depth knob. Every table carries all 15 species
columns, but a gene appears at node N only if its lineage traces to an
ancestor present there, so coverage collapses with depth. Measured on the
outgroup-confirmed PTPC family (111 genes):

    N0   6 HOGs, 105 covered,   6 uncovered
    N1   3 HOGs,  16 covered,  95 uncovered
    N5   6 HOGs,  88 covered,  23 uncovered
    N13  7 HOGs,  20 covered,  91 uncovered

So N0 is the operative subfamily level and the deeper nodes are partial
refinements of the subset that reaches them - useful detail, not a deeper cut
of the whole family. Nesting does hold where coverage overlaps (`refines`
found no violation across all 14 levels), which is what lets the deeper views
be read as refinements at all.

Two consequences the callers must keep visible:
  * a family merged across rounds has members no single table covers, so
    `subfamily_partition` reports them under None rather than dropping them;
  * a merged family's partition can never be coarser than its OG boundaries,
    because HOGs do not cross them. The subfamily axis therefore inherits the
    family axis's fragmentation and cannot be used to check it.
"""
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# Species columns start after HOG, OG and the gene-tree parent clade.
_FIXED_COLUMNS = 3
_LEVEL_RE = re.compile(r"^N(\d+)\.tsv$")


def hog_level_files(directory) -> List[Path]:
    """The per-node HOG tables, shallowest node first.

    Sorted numerically: a lexicographic sort puts N10 before N2 and would
    silently reorder the hierarchy. The `N0.ids.tsv` variant is skipped - it
    carries OrthoFinder's internal sequence ids, not gene ids.
    """
    directory = Path(directory)
    numbered = []
    for path in directory.iterdir():
        m = _LEVEL_RE.match(path.name)
        if m:
            numbered.append((int(m.group(1)), path))
    return [p for _, p in sorted(numbered)]


def parse_hog_table(path) -> Dict[str, str]:
    """gene_id -> HOG id, from one N<node>.tsv.

    Raises when a gene appears in two HOGs at the same level: the partition
    would be ill-defined and the gene would be counted twice downstream.
    """
    path = Path(path)
    hog_of: Dict[str, str] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if len(header) <= _FIXED_COLUMNS:
            raise ValueError(f"{path}: no species columns")
        for lineno, line in enumerate(fh, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            # A truncated row (killed writer, format shift) silently drops
            # its trailing species columns — those genes would vanish from
            # the partition instead of the file being refused.
            if len(parts) != len(header):
                raise ValueError(
                    f"{path}:{lineno}: row has {len(parts)} columns, header "
                    f"has {len(header)}"
                )
            hog = parts[0]
            for cell in parts[_FIXED_COLUMNS:]:
                for gene in cell.split(","):
                    gene = gene.strip()
                    if not gene:
                        continue
                    previous = hog_of.get(gene)
                    if previous is not None and previous != hog:
                        raise ValueError(
                            f"{path}:{lineno}: gene {gene!r} is in both "
                            f"{previous} and {hog}")
                    hog_of[gene] = hog
    return hog_of


def subfamily_partition(
    genes: Sequence[str], hog_of: Dict[str, str],
) -> Dict[Optional[str], List[str]]:
    """Split a family's genes by HOG at one level.

    Genes the table does not cover land under None instead of vanishing.
    That bucket is not an error to hide: HOG tables are per round, so a family
    merged across rounds always has some, and its size is how much of the
    family this level can actually speak about.
    """
    out: Dict[Optional[str], List[str]] = {}
    for gene in genes:
        out.setdefault(hog_of.get(gene), []).append(gene)
    return {k: sorted(v) for k, v in out.items()}


def refines(deeper: Dict[str, str], shallow: Dict[str, str]) -> bool:
    """Whether `deeper` only splits `shallow` - never joins across it.

    Nesting is what makes "subfamily depth" a single number. If a deeper level
    ever united two shallow HOGs, choosing a depth would stop being a cut and
    the levels could not be compared.
    """
    owner: Dict[str, str] = {}
    for gene, deep in deeper.items():
        shallow_hog = shallow.get(gene)
        if shallow_hog is None:
            continue
        previous = owner.setdefault(deep, shallow_hog)
        if previous != shallow_hog:
            return False
    return True


def nested_partitions(
    genes: Sequence[str], levels: Sequence[Tuple[str, Dict[str, str]]],
) -> List[Tuple[str, Dict[Optional[str], List[str]]]]:
    """One subfamily partition per level, in the order given."""
    return [(name, subfamily_partition(genes, hog_of))
            for name, hog_of in levels]
