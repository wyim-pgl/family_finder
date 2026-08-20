#!/usr/bin/env python3
"""Map one species' annotations to family/subfamily ids and name families.

The iceplant (Mcry) pattern: you have a curated annotation table for ONE
species (gene_id<TAB>description). This tool answers two questions against
a finished pipeline run:

  1. For each annotated gene: which family (and subfamily) did it land in?
  2. For each family/subfamily touched: what should it be CALLED, based on
     its annotated members' descriptions (majority vote with support)?

Outputs in --outdir:
  gene_to_family.tsv     annotation_gene, pipeline_gene, family_id,
                         subfamily, description
  family_names.tsv       family_id, name, support, n_annotated
  subfamily_names.tsv    (when --groups given)
  unmatched_genes.txt    annotation ids found in no family

Example:
  python name_families.py --run-dir output_5sp_v2 \
    --annotations mcry_descriptions.tsv --species Mcry \
    --groups possvm_groups.tsv --outdir naming_mcry
"""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from annotate_families import read_families  # noqa: E402
from steps.family_naming import map_annotations, name_groups  # noqa: E402


def read_two_col_tsv(path: Path) -> dict:
    out = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        out[fields[0]] = fields[1]
    return out


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
    ap.add_argument("--run-dir", required=True,
                    help="Pipeline output dir (has summary.tsv)")
    ap.add_argument("--annotations", required=True,
                    help="TSV: gene_id<TAB>description (one species)")
    ap.add_argument("--species", default=None,
                    help="Species prefix to try when annotation ids lack it")
    ap.add_argument("--groups", default=None,
                    help="TSV: pipeline_gene_id<TAB>subfamily (optional)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    families = read_families(Path(args.run_dir))
    annotations = read_two_col_tsv(Path(args.annotations))
    groups = read_two_col_tsv(Path(args.groups)) if args.groups else None
    print(f"{len(families)} families, {len(annotations)} annotated genes")

    rows, unmatched = map_annotations(
        annotations, families, groups=groups, species=args.species
    )
    write_tsv(rows, outdir / "gene_to_family.tsv")
    (outdir / "unmatched_genes.txt").write_text(
        "\n".join(unmatched) + ("\n" if unmatched else "")
    )
    print(f"mapped {len(rows)}, unmatched {len(unmatched)}")

    write_tsv(name_groups(rows, key="family_id"),
              outdir / "family_names.tsv")
    if groups:
        write_tsv(name_groups(rows, key="subfamily"),
                  outdir / "subfamily_names.tsv")
    named = name_groups(rows, key="family_id")
    print(f"named {len(named)} families "
          f"(low-support <0.5: {sum(1 for n in named if n['support'] < 0.5)})")


if __name__ == "__main__":
    main()
