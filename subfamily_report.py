#!/usr/bin/env python3
"""Subfamily diagnostics for one family: WHY does it split?

Standalone (no ete4, no external tools): takes a family's protein
alignment plus a gene->subfamily table (from Possvm / TreeCluster / any
splitter) and writes three TSVs into --outdir:

  sdp_residues.tsv          diagnostic residues per subfamily (with
                            reference numbering when --ref-seq is given)
  taxonomy_attribution.tsv  per subfamily: lineage-specific vs
                            paralog-split. Default evidence (issue #27)
                            is the species tree (--species-tree, the
                            same data/species_tree.nwk the pipeline
                            requires): a subfamily whose species set is
                            a clade strictly inside the family span is
                            lineage-specific; non-monophyletic or
                            family-span sets are paralog splits. A
                            taxonomy TSV (--taxonomy) is an optional
                            label layer (rank names); without a tree it
                            alone drives the legacy rank-purity verdict.
  structure_coherence.tsv   (with --pairs) foldseek all-vs-all within- vs
                            between-subfamily coherence

With --focal-subfamily plus the selection-evidence inputs (--relax-json,
--meme-json/--meme-region, --absrel-json, --expression-share,
--signal-partition), also writes subfunctionalization.md — a narrative
explaining HOW the subfamily diverged (sub- vs neo-functionalization
verdict with the per-axis evidence) — and subfunctionalization.tsv.

Example (PEPC clan):
  python subfamily_report.py \
    --alignment clan_anchor.aln --groups possvm_groups.tsv \
    --species-tree data/species_tree.nwk --taxonomy taxonomy.tsv \
    --pairs foldseek_allvsall.tsv \
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
from steps.subfunctionalization import (  # noqa: E402
    apply_branch_names,
    classify,
    narrative,
    parse_absrel,
    parse_meme,
    parse_relax,
)
from steps.codeml import lrt_pvalue  # noqa: E402
from utils.newick import parse_newick  # noqa: E402
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
    ap.add_argument("--species-tree", default=None,
                    help="Species tree (Newick, leaves = species "
                         "prefixes) — default verdict evidence: clade "
                         "test on each subfamily's species set")
    ap.add_argument("--taxonomy", default=None,
                    help="TSV: species<TAB>genus<TAB>family<TAB>order — "
                         "label layer with --species-tree, sole verdict "
                         "evidence without it (legacy)")
    ap.add_argument("--pairs", default=None,
                    help="foldseek easy-search all-vs-all TSV "
                         "(query,target,evalue,bits,alntmscore)")
    ap.add_argument("--ref-seq", default=None,
                    help="Alignment sequence id for reference numbering")
    ap.add_argument("--min-group", type=int, default=5)
    ap.add_argument("--delimiter", default="_")
    ap.add_argument("--outdir", required=True)
    # --- subfunctionalization narrative (optional evidence axes) ---
    ap.add_argument("--focal-subfamily", default=None,
                    help="Subfamily to explain (writes subfunctionalization.md)")
    ap.add_argument("--family-name", default="family",
                    help="Family label used in the narrative")
    ap.add_argument("--relax-json", default=None, help="HyPhy RELAX json")
    ap.add_argument("--meme-json", default=None, help="HyPhy MEME json")
    ap.add_argument("--absrel-json", default=None, help="HyPhy aBSREL json")
    ap.add_argument("--branch-name-map", default=None,
                    help="TSV short<TAB>real to restore HyPhy-renamed leaves")
    ap.add_argument("--meme-region", default=None,
                    help="Alignment-column range 'LO-HI' of the subfamily-"
                         "specific signal region (MEME sites counted inside)")
    ap.add_argument("--expression-share", type=float, default=None,
                    help="Focal subfamily's share of family expression (0-1)")
    ap.add_argument("--signal-partition", default=None,
                    help="One-line description of the signal-region partition")
    ap.add_argument("--branchsite-mlc", default=None,
                    help="codeml branch-site Model A output (BEB section)")
    ap.add_argument("--branchsite-lnl", default=None,
                    help="'ALT,NULL' lnL pair for the branch-site LRT")
    ap.add_argument("--retargeting-events", type=int, default=0,
                    help="Fitch retargeting events gained by the focal clade")
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
    species_tree = (parse_newick(Path(args.species_tree).read_text())
                    if args.species_tree else None)
    tax_rows = taxonomic_composition(groups, taxonomy, args.delimiter,
                                     species_tree=species_tree)
    write_tsv(tax_rows, outdir / "taxonomy_attribution.tsv")
    source = "species tree" if species_tree else "taxonomy rank purity"
    print(f"taxonomy_attribution.tsv (verdict evidence: {source})")
    for r in tax_rows:
        label = f"  [{r['clade_label']}]" if r.get("clade_label") else ""
        print(f"  {r['subfamily']}: {r['verdict']}{label}")

    if args.pairs:
        scores = parse_pairwise_scores(Path(args.pairs))
        coh = structure_coherence(groups, scores)
        write_tsv(coh, outdir / "structure_coherence.tsv")
        print(f"structure_coherence.tsv: {len(coh)} subfamilies "
              f"({len(scores)} observed pairs)")

    if args.focal_subfamily:
        evidence: dict = {
            "expression_share": args.expression_share,
            "signal_partition": args.signal_partition or "",
            "retargeting_events": args.retargeting_events,
        }
        if args.relax_json:
            evidence["relax"] = parse_relax(Path(args.relax_json))
        if args.meme_json:
            sites = parse_meme(Path(args.meme_json))
            evidence["meme_sites_total"] = len(sites)
            if args.meme_region:
                lo, hi = (int(x) for x in args.meme_region.split("-"))
                evidence["meme_sites_in_region"] = sum(
                    1 for s, _ in sites if lo <= s <= hi
                )
        if args.absrel_json:
            branches = parse_absrel(Path(args.absrel_json))
            if args.branch_name_map:
                name_map = {}
                with open(args.branch_name_map) as f:
                    for line in f:
                        short, real = line.rstrip("\n").split("\t")[:2]
                        name_map[short] = real
                branches = apply_branch_names(branches, name_map)
            # heuristic: HyPhy internal nodes are Node-prefixed; leaves are
            # gene ids (possibly renamed g###). Stem = significant internals.
            evidence["absrel_stem"] = [
                b for b in branches if b[0].startswith("Node")
            ]
            evidence["absrel_terminal"] = [
                b for b in branches if not b[0].startswith("Node")
            ]
        if args.branchsite_mlc:
            from beb_cross import parse_beb  # noqa: E402
            beb = parse_beb(Path(args.branchsite_mlc).read_text())
            sig = [(s_, p_) for s_, _a, p_ in beb if p_ >= 0.95]
            bs: dict = {"beb_sites_total": len(sig)}
            if args.meme_region:
                lo, hi = (int(x) for x in args.meme_region.split("-"))
                bs["beb_sites_in_region"] = sum(1 for s_, _ in sig if lo <= s_ <= hi)
            if args.branchsite_lnl:
                alt, null = (float(x) for x in args.branchsite_lnl.split(","))
                bs["lrt"] = 2 * (alt - null)
                bs["p"] = lrt_pvalue(alt, null)
            evidence["branchsite"] = bs

        verdict = classify(evidence)
        text = narrative(args.family_name, args.focal_subfamily,
                         verdict, evidence)
        (outdir / "subfunctionalization.md").write_text(text + "\n")
        with open(outdir / "subfunctionalization.tsv", "w") as f:
            f.write("subfamily\tverdict\tn_evidence_for\tn_evidence_against\t"
                    "evidence_for\tevidence_against\n")
            f.write("\t".join([
                args.focal_subfamily, verdict["verdict"],
                str(len(verdict["evidence_for"])),
                str(len(verdict["evidence_against"])),
                "; ".join(verdict["evidence_for"]),
                "; ".join(verdict["evidence_against"]),
            ]) + "\n")
        print(f"subfunctionalization.md: {args.focal_subfamily} -> "
              f"{verdict['verdict']}")


if __name__ == "__main__":
    main()
