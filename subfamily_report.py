#!/usr/bin/env python3
"""Subfamily diagnostics for one family: WHY does it split?

Standalone (no ete4, no external tools): takes a family's protein
alignment plus a gene->subfamily table (from Possvm / TreeCluster / any
splitter) and writes three TSVs into --outdir:

  sdp_residues.tsv          diagnostic residues per subfamily (with
                            reference numbering when --ref-seq is given)
  taxonomy_attribution.tsv  per subfamily: lineage-specific (species/
                            genus/family/order) vs paralog-split — is the
                            division taxonomy-driven or intrinsic?
  structure_coherence.tsv   (with --pairs) foldseek all-vs-all within- vs
                            between-subfamily coherence

Example (PEPC clan):
  python subfamily_report.py \
    --alignment clan_anchor.aln --groups possvm_groups.tsv \
    --taxonomy taxonomy.tsv --pairs foldseek_allvsall.tsv \
    --ref-seq ATH_AT1G53310.2 --outdir subfam_report
"""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from steps.subfamily import (  # noqa: E402
    load_taxonomy,
    parse_pairwise_scores,
    sdp_scan,
    structure_coherence,
    taxonomic_composition,
)
from utils.seqio import read_fasta  # noqa: E402


def read_groups(path: Path) -> dict:
    groups: dict = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        gene, subfamily = line.rstrip("\n").split("\t")[:2]
        groups.setdefault(subfamily, []).append(gene)
    return groups


def write_tsv(rows, path: Path):
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--alignment", required=True,
                    help="Family protein alignment (FASTA)")
    ap.add_argument("--groups", required=True,
                    help="TSV: gene_id<TAB>subfamily")
    ap.add_argument("--taxonomy", default=None,
                    help="TSV: species<TAB>genus<TAB>family<TAB>order")
    ap.add_argument("--pairs", default=None,
                    help="foldseek easy-search all-vs-all TSV "
                         "(query,target,evalue,bits,alntmscore)")
    ap.add_argument("--ref-seq", default=None,
                    help="Alignment sequence id for reference numbering")
    ap.add_argument("--min-group", type=int, default=5)
    ap.add_argument("--delimiter", default="_")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    alignment = read_fasta(args.alignment)
    groups = read_groups(Path(args.groups))
    print(f"{len(alignment)} sequences, {len(groups)} subfamilies")

    sdp = sdp_scan(alignment, groups, min_group=args.min_group,
                   ref_seq_id=args.ref_seq)
    write_tsv(sdp, outdir / "sdp_residues.tsv")
    print(f"sdp_residues.tsv: {len(sdp)} diagnostic columns")

    taxonomy = load_taxonomy(Path(args.taxonomy)) if args.taxonomy else {}
    tax_rows = taxonomic_composition(groups, taxonomy, args.delimiter)
    write_tsv(tax_rows, outdir / "taxonomy_attribution.tsv")
    for r in tax_rows:
        print(f"  {r['subfamily']}: {r['verdict']}")

    if args.pairs:
        scores = parse_pairwise_scores(Path(args.pairs))
        coh = structure_coherence(groups, scores)
        write_tsv(coh, outdir / "structure_coherence.tsv")
        print(f"structure_coherence.tsv: {len(coh)} subfamilies "
              f"({len(scores)} observed pairs)")


if __name__ == "__main__":
    main()
