"""Intron structure as a GENE-MODEL QUALITY axis — not a subfamily one.

The plan was the one the comparative-genomics literature uses: exon/intron
structure is an axis independent of sequence, so when the tree cannot resolve
a subfamily boundary the structure ought to. Measured on the PEPC clan (issue
#38) it does not, and the measurement is worth stating because it is what
redefines this module:

  * 70 of the 77 clan members reach a GFF (13 species; the 7 that do not are
    3 ATH/Aco anchors and 4 Sof). Intron positions are mapped from codon
    coordinates onto protein-alignment columns.
  * Nine positions are conserved in >=60% of members: columns 118/145 (one
    position the alignment splits), 298, 406, 481, 515, 559, 783, 1174, 1305.
    ppc-1E1 (n=58) and ppc-1E2 (n=12) share ALL NINE at 91-100%.
  * Positions present in one subfamily and absent from the other: ZERO.

The paper's claim that plant PEPC gene structure has been conserved since the
monocot-dicot split is correct, and that is precisely the problem: the
conservation is OLDER than the duplication we are trying to split, so its
diagnostic power for subfamily membership is exactly nil.

What the same measurement DOES support is the opposite reading. Twenty-four
of the sixty-one genes with a named annotation programme deviate from the
conserved set, and the deviation rate tracks the programme, not the clade:

    Helixer        29 genes   0.52 off-set introns/gene   0.10 missing
    EVM            18         0.39                        0.11
    phytozome v13   7         0.29                        0.00
    phytozome v12   3         1.33                        1.33
    AUGUSTUS/UNR    4         0.00                        0.00
    (unknown)       9         not scored -- see below

So a gene whose introns fall outside the family's conserved set is evidence
about the MODEL, and the annotation programme is a confounder that has to
travel with the number. Every row this module emits therefore carries GFF
column 2, and a gene whose programme is unknown gets no intron count at all —
an uncontrolled count reads as biology and is not citable (same discipline as
`steps.subfamily.structure_coherence`, issue #39-3).

Genomic intron coordinates do not survive divergence; the position of an
intron WITHIN the coding sequence does, which is why the codon coordinate is
the one carried across species and then mapped onto alignment columns.
"""

import logging
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from steps.gene_structure_model import GeneModel
from steps.gff_join import (join_by_species, join_gff_models,  # noqa: F401
                            strip_species)

logger = logging.getLogger("family_finder")

# A conserved position is one at least this fraction of the family carries.
CONSERVED_MIN_FRACTION = 0.6
# Alignment columns drift by an indel or two between species; positions within
# this many columns of each other are the same intron.
COLUMN_WINDOW = 2
# A position is diagnostic only if some group nearly always has it and another
# nearly never does. On the PEPC clan nothing comes close.
DIAGNOSTIC_PRESENT_MIN = 0.6
DIAGNOSTIC_ABSENT_MAX = 0.2
# Spread in off-set introns per gene across annotation programmes, above which
# the per-gene deviation must not be read as biology.
# An alignment indel can displace one intron by more columns than `window`
# tolerates (the first PEPC intron moves 27). Beyond this many columns the two
# halves are not called the same position -- a distant pair of rare introns
# would otherwise be merged into a conserved one that never existed.
MAX_SPLIT_GAP = 60
# Each half of a split position must still be carried by this fraction of the
# family. A column one gene has and sixty-nine do not is a model defect, and
# absorbing it into the conserved set would hide exactly what this axis
# reports. Both PEPC halves are near a half of the family (34 and 36 of 70).
MIN_SPLIT_HALF = 0.1
# One intron cannot be in two places, so the halves must be essentially
# complementary -- but one gene of seventy carrying both is a model with a
# genuine extra intron there, and it must not veto a position the other
# sixty-nine share.
MAX_SPLIT_OVERLAP = 0.05
PROGRAMME_SPREAD_MAX = 0.2
MIN_GENES_PER_PROGRAMME = 3

UNKNOWN_PROGRAMME = "unknown"
# GFF3 spells an undefined column-2 value as ".", which is not the name of
# an annotation programme; treating it as one would let a whole genome
# (Sund's Hun.hic.gff3, Ococ's Oco.filter.tidy) form its own "programme"
# and hide the very confounder this axis exists to expose.
_NO_SOURCE = ("", ".")
UNCONTROLLED_STATUS = (
    "uncontrolled: the annotation programme is unknown, so an intron count "
    "cannot be told apart from an annotator's habit"
)

_ID = re.compile(r"ID=([^;]+)")
_PARENT = re.compile(r"Parent=([^;]+)")
_TRANSCRIPT_TYPES = ("mRNA", "transcript")


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------

def cds_structure(blocks: Sequence[Tuple[int, int]],
                  strand: str) -> Tuple[int, int, List[float]]:
    """(n_cds, total nt, intron positions in codon coordinates).

    Blocks are ordered by TRANSCRIPTION, so a minus-strand gene is read from
    the highest coordinate down; ordering by start would place every intron at
    the wrong codon. A position of 56.0 means the intron follows codon 56
    exactly (phase 0); 56.67 means it interrupts codon 57.
    """
    if not blocks:
        return 0, 0, []
    ordered = sorted(blocks, reverse=(strand == "-"))
    lengths = [end - start + 1 for start, end in ordered]
    positions, cumulative = [], 0
    for length in lengths[:-1]:
        cumulative += length
        positions.append(cumulative / 3.0)
    return len(lengths), sum(lengths), positions


def codon_to_column(aligned_seq: str, codon_position: float) -> Optional[int]:
    """Alignment column (1-based) holding the codon an intron sits at.

    Returns None when the position falls past the end of this sequence rather
    than inventing a column — a truncated model must not silently borrow one.
    """
    residue = _residue_of(codon_position)
    if residue < 1:
        return None
    seen = 0
    for column, char in enumerate(aligned_seq, start=1):
        if char == "-":
            continue
        seen += 1
        if seen == residue:
            return column
    return None


def _residue_of(codon_position: float) -> int:
    """Codon coordinate -> the residue the intron is attached to."""
    whole = math.floor(codon_position + 1e-9)
    fraction = codon_position - whole
    return int(whole) if fraction < 1e-6 else int(whole) + 1


def intron_columns(aligned_seq: str,
                   codon_positions: Iterable[float]) -> List[int]:
    """Every intron position of one gene, as alignment columns."""
    out = []
    for position in codon_positions:
        column = codon_to_column(aligned_seq, position)
        if column is not None:
            out.append(column)
    return sorted(out)


# ---------------------------------------------------------------------------
# the family's conserved set, and departures from it
# ---------------------------------------------------------------------------

def _cluster(per_gene: Dict[str, List[int]],
             window: int) -> List[Dict[str, int]]:
    """Greedily group observed columns into positions, gene -> its column."""
    remaining = {gene: list(columns) for gene, columns in per_gene.items()}
    clusters: List[Dict[str, int]] = []
    while True:
        candidates = sorted({c for cols in remaining.values() for c in cols})
        if not candidates:
            break
        best, best_support = None, 0
        for candidate in candidates:
            support = sum(
                1 for cols in remaining.values()
                if any(abs(c - candidate) <= window for c in cols)
            )
            if support > best_support:
                best, best_support = candidate, support
        members: Dict[str, int] = {}
        for gene, cols in remaining.items():
            hits = [c for c in cols if abs(c - best) <= window]
            if hits:
                members[gene] = min(hits, key=lambda c: abs(c - best))
                remaining[gene] = [c for c in cols if c not in hits]
        clusters.append(members)
    return clusters


def _merge_split_positions(clusters: List[Dict[str, int]],
                           n_genes: int,
                           min_fraction: float,
                           max_gap: int = MAX_SPLIT_GAP) -> List[Dict[str, int]]:
    """Rejoin one intron that the ALIGNMENT, not evolution, put in two places.

    Measured case: the first PEPC intron sits in alignment column 118 in 32
    genes and column 145 in 28 — 27 columns apart, far outside any sane
    jitter window — because an N-terminal indel is not aligned across the
    clan. Neither column reaches the 60% floor on its own, so the position
    disappears from the conserved set and every one of the 60 genes is then
    charged with an extra intron it does not have.

    A merge is allowed only when the evidence says one intron:

      * the two columns are within `max_gap` of each other,
      * each half is carried by at least `MIN_SPLIT_HALF` of the family,
      * almost no gene carries both columns (one intron cannot be in two
        places; `MAX_SPLIT_OVERLAP` allows for the odd model that does),
      * no OTHER position carried by more than one gene lies between them
        (that would make them two positions with something in the middle; a
        single gene's stray intron does not veto a family-level position --
        that stray IS the model defect this axis reports), and
      * together they clear the conservation floor.

    Two positions that each clear the floor cannot both be missing from the
    other's genes, so this can only ever rejoin halves, never fuse two real
    positions.
    """
    def fraction(members):
        return len(members) / n_genes if n_genes else 0.0

    merged = True
    while merged:
        merged = False
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                left, right = clusters[i], clusters[j]
                if min(fraction(left), fraction(right)) < MIN_SPLIT_HALF:
                    continue
                union = {**left, **right}
                overlap = len(set(left) & set(right))
                if union and overlap / len(union) > MAX_SPLIT_OVERLAP:
                    continue
                if fraction(union) < min_fraction:
                    continue
                lo, hi = sorted([_median(list(left.values())),
                                 _median(list(right.values()))])
                if hi - lo > max_gap:
                    continue
                blocked = any(
                    len(other) > 1 and lo < _median(list(other.values())) < hi
                    for k, other in enumerate(clusters) if k not in (i, j)
                )
                if blocked:
                    continue
                clusters[i] = union
                clusters.pop(j)
                merged = True
                break
            if merged:
                break
    return clusters


def conserved_columns(per_gene: Dict[str, List[int]],
                      min_fraction: float = CONSERVED_MIN_FRACTION,
                      window: int = COLUMN_WINDOW,
                      max_split_gap: int = MAX_SPLIT_GAP) -> List[dict]:
    """Positions at least `min_fraction` of the family carries.

    Columns within `window` of each other are one position: an indel anywhere
    upstream shifts the whole downstream alignment, so demanding an exact
    column would split one ancestral intron into as many positions as there
    are species. A position the alignment has split MUCH further than that is
    rejoined by `_merge_split_positions`, under conditions that require the
    two halves to behave like one intron.

    Each entry carries `columns`: every anchor the position occupies. A gene
    matches the position if it has an intron within `window` of ANY of them.
    """
    n_genes = len(per_gene)
    if not n_genes:
        return []
    clusters = _merge_split_positions(_cluster(per_gene, window),
                                      n_genes, min_fraction,
                                      max_gap=max_split_gap)

    out = []
    for members in clusters:
        if len(members) < min_fraction * n_genes:
            continue
        columns = sorted(members.values())
        anchors = _anchors(columns, window)
        out.append({
            "column": int(_median(columns)),
            "columns": anchors,
            "n_genes": len(members),
            "fraction": len(members) / n_genes,
            "column_range": [min(columns), max(columns)],
            "split_by_alignment": len(anchors) > 1,
        })
    return sorted(out, key=lambda c: c["column"])


def _anchors(columns: Sequence[int], window: int) -> List[int]:
    """One representative per run of columns that `window` already links."""
    anchors, run = [], []
    for column in sorted(columns):
        if run and column - run[-1] > window:
            anchors.append(int(_median(run)))
            run = []
        run.append(column)
    if run:
        anchors.append(int(_median(run)))
    return anchors


def _median(values: Sequence[int]) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[mid])
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def deviations(columns: Sequence[int], conserved: Sequence[dict],
               window: int = COLUMN_WINDOW) -> dict:
    """This gene's introns outside the conserved set, and the set's gaps."""
    extra = [c for c in columns
             if not any(_hits(c, entry, window) for entry in conserved)]
    missing = [entry["column"] for entry in conserved
               if not any(_hits(c, entry, window) for c in columns)]
    return {"extra": sorted(extra), "missing": sorted(missing)}


def _hits(column: int, entry: dict, window: int) -> bool:
    """Does one intron column fall on this conserved position?"""
    anchors = entry.get("columns") or [entry["column"]]
    return any(abs(column - anchor) <= window for anchor in anchors)


def diagnostic_columns(per_gene: Dict[str, List[int]],
                       groups: Dict[str, List[str]],
                       conserved: Sequence[dict],
                       window: int = COLUMN_WINDOW,
                       present_min: float = DIAGNOSTIC_PRESENT_MIN,
                       absent_max: float = DIAGNOSTIC_ABSENT_MAX) -> List[dict]:
    """Per conserved position, how often each group carries it — and whether
    that separates the groups at all.

    The PEPC answer is no: nine positions, every one of them at 92-100% in
    both subfamilies. The function exists so that claim is measured on each
    new family instead of assumed either way.
    """
    rows = []
    for entry in conserved:
        anchor = entry["column"]
        by_group, present_counts = {}, {}
        for name, members in groups.items():
            seen = [g for g in members if g in per_gene]
            if not seen:
                by_group[name] = None
                continue
            present = sum(
                1 for g in seen
                if any(_hits(c, entry, window) for c in per_gene[g])
            )
            by_group[name] = present / len(seen)
            present_counts[name] = (present, len(seen))
        rates = [v for v in by_group.values() if v is not None]
        diagnostic = bool(rates) and (max(rates) >= present_min
                                      and min(rates) <= absent_max)
        rows.append({
            "column": anchor,
            "by_group": by_group,
            "counts": present_counts,
            "diagnostic": diagnostic,
        })
    return rows


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def parse_gff_cds(path: Path) -> List[GeneModel]:
    """One GeneModel per gene: its longest coding transcript, with the source.

    GFF column 2 is kept because it is the covariate the whole axis depends
    on. It is read from the transcript line, which is where the programme that
    produced the model is named.
    """
    cds: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    transcript_gene: Dict[str, str] = {}
    transcript_source: Dict[str, str] = {}

    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature, attributes = fields[2], fields[8]
            if feature in _TRANSCRIPT_TYPES:
                found = _ID.search(attributes)
                if not found:
                    continue
                transcript = found.group(1)
                parent = _PARENT.search(attributes)
                transcript_gene[transcript] = (parent.group(1) if parent
                                               else transcript)
                transcript_source[transcript] = fields[1]
            elif feature == "CDS":
                parent = _PARENT.search(attributes)
                if not parent:
                    continue
                for transcript in parent.group(1).split(","):
                    cds[transcript].append((int(fields[3]), int(fields[4]),
                                            fields[6]))

    by_gene: Dict[str, List[Tuple[str, List[Tuple[int, int, str]]]]] = defaultdict(list)
    for transcript, blocks in cds.items():
        gene = transcript_gene.get(transcript, transcript).lstrip("_")
        by_gene[gene].append((transcript, blocks))

    models = []
    for gene, transcripts in by_gene.items():
        transcript, blocks = max(
            transcripts,
            key=lambda item: sum(end - start + 1 for start, end, _ in item[1]),
        )
        models.append(GeneModel(
            gene_id=gene,
            transcript_id=transcript,
            source=transcript_source.get(transcript, ""),
            strand=blocks[0][2],
            blocks=tuple(sorted((start, end) for start, end, _ in blocks)),
        ))
    return sorted(models, key=lambda m: m.gene_id)


# ---------------------------------------------------------------------------
# rows, the programme covariate, and the report
# ---------------------------------------------------------------------------

def build_rows(alignment: Dict[str, str], models: Dict[str, GeneModel],
               conserved: Sequence[dict],
               groups: Optional[Dict[str, List[str]]] = None,
               window: int = COLUMN_WINDOW) -> List[dict]:
    """One row per gene, each carrying its annotation programme.

    A gene whose programme is unknown is emitted with NO intron numbers. The
    count on its own would be read as biology, and the PEPC measurement is
    that the count is largely a property of the annotator.
    """
    group_of = {}
    for name, members in (groups or {}).items():
        for member in members:
            group_of[member] = name

    rows = []
    for gene in sorted(models):
        model = models[gene]
        seq = alignment.get(gene)
        source = model.source.strip()
        source = UNKNOWN_PROGRAMME if source in _NO_SOURCE else source
        n_cds, cds_len, positions = cds_structure(model.blocks, model.strand)
        columns = intron_columns(seq, positions) if seq is not None else []
        row = {
            "gene": gene,
            "species": gene.split("_", 1)[0],
            "subfamily": group_of.get(gene, ""),
            "annotation_source": source,
            "gff_gene_id": model.gene_id,
            "n_cds": n_cds,
            "cds_len_nt": cds_len,
            "intron_columns": columns,
            "n_introns": None,
            "n_extra": None,
            "n_missing": None,
            "extra_columns": [],
            "missing_columns": [],
            "status": "",
        }
        if source == UNKNOWN_PROGRAMME:
            row["status"] = UNCONTROLLED_STATUS
            row["intron_columns"] = []
            rows.append(row)
            continue
        if seq is None:
            row["status"] = "not in alignment"
            rows.append(row)
            continue
        deviation = deviations(columns, conserved, window=window)
        row.update({
            "n_introns": len(columns),
            "n_extra": len(deviation["extra"]),
            "n_missing": len(deviation["missing"]),
            "extra_columns": deviation["extra"],
            "missing_columns": deviation["missing"],
            "status": "ok" if not (deviation["extra"] or deviation["missing"])
                      else "deviates from the family's conserved intron set",
        })
        rows.append(row)
    return rows


def deviation_by_programme(rows: Sequence[dict]) -> List[dict]:
    """Deviation rate per annotation programme — the control, always emitted."""
    buckets: Dict[str, List[dict]] = defaultdict(list)
    for row in rows:
        buckets[row.get("annotation_source") or UNKNOWN_PROGRAMME].append(row)

    out = []
    for source in sorted(buckets):
        group = buckets[source]
        scored = [r for r in group if r.get("n_extra") is not None]
        out.append({
            "annotation_source": source,
            "n_genes": len(group),
            "n_scored": len(scored),
            "extra_per_gene": (sum(r["n_extra"] for r in scored) / len(scored)
                               if scored else None),
            "missing_per_gene": (sum(r["n_missing"] for r in scored) / len(scored)
                                 if scored else None),
            "n_deviating": sum(1 for r in scored
                               if r["n_extra"] or r["n_missing"]),
        })
    return out


def programme_confounding(by_programme: Sequence[dict],
                          spread_max: float = PROGRAMME_SPREAD_MAX,
                          min_genes: int = MIN_GENES_PER_PROGRAMME) -> dict:
    """Is the per-gene deviation readable as biology, or as annotator habit?"""
    usable = [p for p in by_programme
              if p["extra_per_gene"] is not None and p["n_scored"] >= min_genes]
    if len(usable) < 2:
        return {
            "confounded": None,
            "reason": (
                "UNCONTROLLED: fewer than two annotation programmes have "
                f"{min_genes} scored genes, so the programme effect cannot be "
                "estimated on this family. Do not cite a per-gene intron "
                "deviation as biology"
            ),
            "spread": None,
        }
    rates = [p["extra_per_gene"] for p in usable]
    spread = max(rates) - min(rates)
    if spread > spread_max:
        return {
            "confounded": True,
            "spread": spread,
            "reason": (
                f"off-set introns per gene range {min(rates):.2f}-"
                f"{max(rates):.2f} across {len(usable)} annotation programmes "
                f"(spread {spread:.2f} > {spread_max:.2f}). Deviation from the "
                "conserved set flags MODEL quality here; it is not comparable "
                "between species annotated by different programmes"
            ),
        }
    return {
        "confounded": False,
        "spread": spread,
        "reason": (
            f"off-set introns per gene vary by only {spread:.2f} across "
            f"{len(usable)} annotation programmes, so a deviating gene is not "
            "explained by its annotator"
        ),
    }


def structure_report(alignment: Dict[str, str], models: Dict[str, GeneModel],
                     groups: Optional[Dict[str, List[str]]] = None,
                     min_fraction: float = CONSERVED_MIN_FRACTION,
                     window: int = COLUMN_WINDOW) -> dict:
    """The whole axis: conserved set, per-gene deviation, programme control."""
    per_gene: Dict[str, List[int]] = {}
    for gene, model in models.items():
        seq = alignment.get(gene)
        if seq is None:
            continue
        _n, _len, positions = cds_structure(model.blocks, model.strand)
        per_gene[gene] = intron_columns(seq, positions)

    conserved = conserved_columns(per_gene, min_fraction=min_fraction,
                                  window=window)
    rows = build_rows(alignment, models, conserved, groups=groups,
                      window=window)
    by_programme = deviation_by_programme(rows)
    diagnostics = (diagnostic_columns(per_gene, groups, conserved, window=window)
                   if groups else [])
    scored = [r for r in rows if r["n_extra"] is not None]

    return {
        "rows": rows,
        "conserved_columns": conserved,
        "diagnostic_columns": diagnostics,
        "summary": {
            "n_genes": len(rows),
            "n_scored": len(scored),
            "n_conserved_positions": len(conserved),
            "n_diagnostic_positions": sum(1 for d in diagnostics
                                          if d["diagnostic"]),
            "n_deviating_genes": sum(1 for r in scored
                                     if r["n_extra"] or r["n_missing"]),
            "by_programme": by_programme,
            "confounded": programme_confounding(by_programme),
        },
    }


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

COLUMNS = ("gene", "species", "subfamily", "annotation_source", "gff_gene_id",
           "n_cds", "cds_len_nt", "n_introns", "intron_columns", "n_extra",
           "n_missing", "extra_columns", "missing_columns", "status")


def _cell(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, (list, tuple)):
        return ",".join(str(v) for v in value)
    return str(value)


def write_table(rows: Sequence[dict], path: Path) -> None:
    """TSV for `annotation_matrix.py --gene-structure`."""
    with open(path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(_cell(row.get(c)) for c in COLUMNS) + "\n")


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alignment", required=True,
                        help="Family protein alignment (FASTA)")
    parser.add_argument("--gff", action="append", default=[], required=True,
                        metavar="SPECIES=PATH",
                        help="GFF3 for one species prefix; repeatable. A "
                             "species may be given more than once when its "
                             "annotation is split over files")
    parser.add_argument("--groups", default=None,
                        help="TSV of gene<TAB>subfamily (enables the "
                             "diagnostic-power measurement)")
    parser.add_argument("-o", "--out", required=True, help="Output TSV")
    parser.add_argument("--json", default=None, help="Also write the summary")
    parser.add_argument("--min-fraction", type=float,
                        default=CONSERVED_MIN_FRACTION)
    parser.add_argument("--window", type=int, default=COLUMN_WINDOW)
    parser.add_argument("--max-unmatched", type=float, default=0.15,
                        help="Refuse when more than this fraction of family "
                             "gene ids fail to reach a GFF model")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    alignment = read_fasta(Path(args.alignment))
    models_by_species: Dict[str, List[GeneModel]] = defaultdict(list)
    for spec in args.gff:
        species, _, path = spec.partition("=")
        if not path:
            parser.error(f"--gff expects SPECIES=PATH, got {spec!r}")
        models_by_species[species].extend(parse_gff_cds(Path(path)))

    groups: Optional[Dict[str, List[str]]] = None
    if args.groups:
        groups = defaultdict(list)
        for line in open(args.groups):
            if not line.strip():
                continue
            gene, _, label = line.rstrip("\n").partition("\t")
            groups[label].append(gene)
        groups = dict(groups)

    mapping, join = join_by_species(list(alignment), models_by_species,
                                    max_unmatched=args.max_unmatched)
    report = structure_report(alignment, mapping, groups=groups,
                              min_fraction=args.min_fraction,
                              window=args.window)
    report["summary"]["id_join"] = join

    write_table(report["rows"], Path(args.out))
    if args.json:
        Path(args.json).write_text(json.dumps(
            {"summary": report["summary"],
             "conserved_columns": report["conserved_columns"],
             "diagnostic_columns": report["diagnostic_columns"]}, indent=2))

    summary = report["summary"]
    logger.info("%d of %d gene(s) joined to a GFF model, %d unmatched (%s)",
                join["n_matched"], join["n_requested"], join["n_unmatched"],
                ", ".join(join["unmatched"][:8]) or "none")
    if join["species_without_gff"]:
        logger.info("no GFF at all for: %s",
                    ", ".join(join["species_without_gff"]))
    logger.info("%d conserved intron position(s); %d diagnostic; "
                "%d gene(s) deviate",
                summary["n_conserved_positions"],
                summary["n_diagnostic_positions"],
                summary["n_deviating_genes"])
    for entry in summary["by_programme"]:
        logger.info("  %-16s n=%-3d extra/gene=%s missing/gene=%s",
                    entry["annotation_source"], entry["n_genes"],
                    "NA" if entry["extra_per_gene"] is None
                    else f"{entry['extra_per_gene']:.2f}",
                    "NA" if entry["missing_per_gene"] is None
                    else f"{entry['missing_per_gene']:.2f}")
    logger.info("programme control: %s", summary["confounded"]["reason"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
