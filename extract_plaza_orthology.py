#!/usr/bin/env python3
"""Build the --plaza-orthology TSV for name_families.py (issue #30).

Our species are not in PLAZA (Dicots 5.0 carries 131 proteomes; the
closest Caryophyllales are bvu/cqu/spo), so the our-gene -> ATH-ortholog
link is computed by DIAMOND against the PLAZA ath proteome and this
script only PARSES the results — it runs anywhere, no network, no tools.

Cluster-side inputs (pronghorn; PLAZA=/data/gpfs/assoc/pgl/data/
sequence_data/plaza, DIAMOND from the orthofinder conda env):

  zcat $PLAZA/protein/selected/proteome.selected_transcript.ath.fasta.gz \\
      > ath_selected.fa
  diamond makedb --in ath_selected.fa -d ath_selected
  cat data/pep/*.pep.fa > all.pep.fa
  diamond makedb --in all.pep.fa -d all
  # forward: our proteins -> ATH
  diamond blastp -q all.pep.fa -d ath_selected -o fwd.tsv \\
      --max-target-seqs 5 -e 1e-10 --threads 8 --quiet
  # reverse (optional, enables the rbh flag): ATH -> our proteins
  diamond blastp -q ath_selected.fa -d all -o rev.tsv \\
      --max-target-seqs 5 -e 1e-10 --threads 8 --quiet

Then:
  python extract_plaza_orthology.py --forward fwd.tsv --reverse rev.tsv \\
      --ath-annotations ath_symbols.tsv -o plaza_orthology.tsv

--ath-annotations is a 2-col TSV: ATH gene id<TAB>symbol/description
(e.g. from the PLAZA annotation download or TAIR functional descriptions,
gene-level ids like AT1G53310). Genes without an entry keep the ATH gene
id as their description — never blank, a name you can still look up.

Output columns: gene_id, ath_gene, description, pident, bitscore,
evalue, rbh (yes/no, empty without --reverse). name_families.py reads
the first three.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

# PLAZA ath ids are transcript-level (AT1G53310.2); annotations and the
# PPC anchors are gene-level.
_TRANSCRIPT_SUFFIX = re.compile(r"\.\d+$")


def strip_transcript(seq_id: str) -> str:
    """AT1G53310.2 -> AT1G53310 (no-op on gene-level ids)."""
    return _TRANSCRIPT_SUFFIX.sub("", seq_id)


def best_hits(lines: Iterable[str]) -> Dict[str, dict]:
    """DIAMOND outfmt-6 lines -> query -> its best (highest-bitscore) hit.

    Ties keep the first-seen hit (DIAMOND emits hits best-first, so this
    is DIAMOND's own ranking, and deterministic on re-parse).
    """
    hits: Dict[str, dict] = {}
    for line in lines:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        query, subject = fields[0], fields[1]
        bitscore = float(fields[11])
        current = hits.get(query)
        if current is None or bitscore > current["bitscore"]:
            hits[query] = {
                "subject": subject,
                "pident": float(fields[2]),
                "evalue": fields[10],
                "bitscore": bitscore,
            }
    return hits


def best_hits_per_species(lines: Iterable[str]) -> Dict[str, Dict[str, dict]]:
    """Reverse search vs a multi-species DB -> query -> species -> best hit.

    Subject ids follow the pipeline's SpeciesPrefix_GeneID convention;
    RBH is judged within one species (one ATH gene has one ortholog PER
    species — a concatenated-DB overall best would let the closest
    species mask everyone else's reciprocal hit).
    """
    hits: Dict[str, Dict[str, dict]] = {}
    for line in lines:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        query, subject = fields[0], fields[1]
        species = subject.split("_", 1)[0]
        bitscore = float(fields[11])
        per_species = hits.setdefault(query, {})
        current = per_species.get(species)
        if current is None or bitscore > current["bitscore"]:
            per_species[species] = {"subject": subject, "bitscore": bitscore}
    return hits


def reciprocal_best(
    gene: str, ath_transcript: str, reverse: Dict[str, Dict[str, dict]]
) -> bool:
    """True when the ATH transcript's best hit in gene's species is gene."""
    species = gene.split("_", 1)[0]
    back = reverse.get(ath_transcript, {}).get(species)
    return back is not None and back["subject"] == gene


def build_orthology_rows(
    forward: Dict[str, dict],
    reverse: Optional[Dict[str, Dict[str, dict]]],
    annotations: Dict[str, str],
    rbh_only: bool = False,
) -> List[dict]:
    """Join best hits with ATH annotations into name_families input rows."""
    rows: List[dict] = []
    for gene in sorted(forward):
        hit = forward[gene]
        ath_transcript = hit["subject"]
        ath_gene = strip_transcript(ath_transcript)
        if reverse is None:
            rbh = ""
        elif reciprocal_best(gene, ath_transcript, reverse):
            rbh = "yes"
        else:
            rbh = "no"
        if rbh_only and rbh != "yes":
            continue
        rows.append({
            "gene_id": gene,
            "ath_gene": ath_gene,
            "description": annotations.get(ath_gene, ath_gene),
            "pident": f"{hit['pident']:.1f}",
            "bitscore": f"{hit['bitscore']:.0f}",
            "evalue": hit["evalue"],
            "rbh": rbh,
        })
    return rows


def read_annotations(path: Path) -> Dict[str, str]:
    """2-col TSV: ATH gene id -> symbol/description."""
    annotations: Dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        annotations[strip_transcript(fields[0])] = fields[1]
    return annotations


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--forward", required=True,
                    help="DIAMOND outfmt 6: our proteins vs PLAZA ath")
    ap.add_argument("--reverse", default=None,
                    help="DIAMOND outfmt 6: PLAZA ath vs our proteins "
                         "(enables the rbh flag)")
    ap.add_argument("--ath-annotations", default=None,
                    help="TSV: ATH gene id<TAB>symbol/description")
    ap.add_argument("--rbh-only", action="store_true",
                    help="Keep only reciprocal best hits (needs --reverse)")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    if args.rbh_only and not args.reverse:
        ap.error("--rbh-only needs --reverse")

    with open(args.forward) as f:
        forward = best_hits(f)
    reverse = None
    if args.reverse:
        with open(args.reverse) as f:
            reverse = best_hits_per_species(f)
    annotations = (
        read_annotations(Path(args.ath_annotations))
        if args.ath_annotations else {}
    )

    rows = build_orthology_rows(forward, reverse, annotations,
                                rbh_only=args.rbh_only)
    with open(args.output, "w") as f:
        f.write("gene_id\tath_gene\tdescription\tpident\tbitscore\t"
                "evalue\trbh\n")
        for row in rows:
            f.write("\t".join(row.values()) + "\n")
    n_named = sum(1 for r in rows if r["description"] != r["ath_gene"])
    print(f"{len(rows)} orthology rows ({n_named} with a symbol, "
          f"{len(forward) - len(rows)} dropped)")


if __name__ == "__main__":
    main()
