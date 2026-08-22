#!/usr/bin/env python3
"""Map annotations to family/subfamily ids and name families.

Two evidence layers over a finished pipeline run:

  direct     — the iceplant (Mcry) pattern: a curated annotation table
               for ONE species (gene_id<TAB>description).
  orthology  — PLAZA-orthology transfer (issue #30): ATH ortholog
               symbols/descriptions per gene, from
               extract_plaza_orthology.py (--plaza-orthology).

Each family/subfamily touched is named by a weighted majority vote of
its members' descriptions — direct annotation outweighs orthology
transfer (--orthology-weight, default 0.5) — and every row and name
carries its provenance (source: direct | orthology).

Outputs in --outdir:
  gene_to_family.tsv     annotation_gene, pipeline_gene, family_id,
                         subfamily, description, source, ortholog
  family_names.tsv       family_id, name, support, n_annotated,
                         n_direct, n_orthology, provenance
  subfamily_names.tsv    (when --groups given)
  unmatched_genes.txt    annotation/orthology ids found in no family

Examples:
  python name_families.py --run-dir output_5sp_v2 \
    --annotations mcry_descriptions.tsv --species Mcry \
    --groups possvm_groups.tsv --outdir naming_mcry
  python name_families.py --run-dir output_5sp_v2 \
    --plaza-orthology plaza_orthology.tsv --outdir naming_plaza
"""

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from annotate_families import read_families  # noqa: E402
from steps.family_naming import (  # noqa: E402
    combine_sources,
    map_annotations,
    map_orthology,
    name_groups,
)


def read_orthology_tsv(path: Path) -> dict:
    """extract_plaza_orthology.py output -> gene_id -> (ortholog, description).

    Needs the first three columns (gene_id, ath_gene/ortholog,
    description); extra columns (pident, rbh, ...) are ignored. A header
    line is detected by its first field and skipped.
    """
    out = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3 or fields[0] == "gene_id":
            continue
        out[fields[0]] = (fields[1], fields[2])
    return out


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
    ap.add_argument("--annotations", default=None,
                    help="TSV: gene_id<TAB>description (one species)")
    ap.add_argument("--plaza-orthology", default=None,
                    help="TSV from extract_plaza_orthology.py: gene_id, "
                         "ortholog, description[, ...]")
    ap.add_argument("--orthology-weight", type=float, default=0.5,
                    help="Vote weight of an orthology transfer vs 1.0 "
                         "for a direct annotation (default 0.5)")
    ap.add_argument("--species", default=None,
                    help="Species prefix to try when annotation ids lack it")
    ap.add_argument("--groups", default=None,
                    help="TSV: pipeline_gene_id<TAB>subfamily (optional)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    if not args.annotations and not args.plaza_orthology:
        ap.error("need --annotations and/or --plaza-orthology")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    families = read_families(Path(args.run_dir))
    groups = read_two_col_tsv(Path(args.groups)) if args.groups else None

    rows, unmatched = [], []
    if args.annotations:
        annotations = read_two_col_tsv(Path(args.annotations))
        rows, unmatched = map_annotations(
            annotations, families, groups=groups, species=args.species
        )
        print(f"{len(families)} families, {len(annotations)} annotated "
              f"genes -> {len(rows)} direct rows")
    if args.plaza_orthology:
        orthology = read_orthology_tsv(Path(args.plaza_orthology))
        ortho_rows, ortho_unmatched = map_orthology(
            orthology, families, groups=groups, species=args.species
        )
        n_before = len(rows)
        rows = combine_sources(rows, ortho_rows)
        unmatched = unmatched + ortho_unmatched
        print(f"{len(orthology)} orthology genes -> "
              f"{len(rows) - n_before} rows kept "
              f"(direct suppresses same-gene orthology)")

    write_tsv(rows, outdir / "gene_to_family.tsv")
    (outdir / "unmatched_genes.txt").write_text(
        "\n".join(unmatched) + ("\n" if unmatched else "")
    )
    print(f"mapped {len(rows)}, unmatched {len(unmatched)}")

    weights = {"direct": 1.0, "orthology": args.orthology_weight}
    fam_sizes = {fid: len(members) for fid, members in families.items()}
    named = name_groups(rows, key="family_id", weights=weights,
                        group_sizes=fam_sizes)
    write_tsv(named, outdir / "family_names.tsv")
    if groups:
        sub_sizes: dict = {}
        for sub in groups.values():
            sub_sizes[sub] = sub_sizes.get(sub, 0) + 1
        named_sub = name_groups(rows, key="subfamily", weights=weights,
                                group_sizes=sub_sizes)
        write_tsv(named_sub, outdir / "subfamily_names.tsv")
        thin = [n for n in named_sub
                if n["coverage"] is not None and n["coverage"] < 0.5]
        if thin:
            print(f"  WARNING: {len(thin)} subfamily names rest on <50% of "
                  "their members — support alone will look confident. "
                  "Check the coverage column: "
                  + ", ".join(f"{n['group_id']}({n['n_annotated']}/"
                              f"{n['group_size']})" for n in thin[:5]))
    print(f"named {len(named)} families "
          f"(low-support <0.5: {sum(1 for n in named if n['support'] < 0.5)})")


if __name__ == "__main__":
    main()
