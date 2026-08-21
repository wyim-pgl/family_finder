#!/usr/bin/env python3
"""Family annotation layers + tier-3 prep (issues #17, #20).

Thin standalone CLI over a finished pipeline run (pattern of
find_neofunctionalization.py) — heavy compute (ESM-2, ESMFold, foldseek,
ESM-ECForest) happens on the cluster; this script runs anywhere.

Modes (combinable):
  --ec              EC annotation layer (issues #20/#28): per-family EC
                    consensus with advisory disagreement flags ->
                    ec_consensus.tsv, and EC-switch events on family gene
                    trees (Fitch, same machinery as retargeting) ->
                    ec_switch_events.tsv. Predictions come from any of
                    --emapper (eggNOG-mapper .emapper.annotations),
                    --clean-csv (CLEAN maxsep CSV) — both gate-validated,
                    merged when both given (emapper EC authoritative,
                    CLEAN confidence attached) — or the legacy
                    --ecforest-cache TSV.
  --prep-fold-list  Tier-3 stage-1 prep (issue #17): list the genes to fold —
                    every unplaced gene (unplaced_proteins.fa) plus one
                    representative per family -> fold_list.tsv. Feed this to
                    the cluster-side ESMFold SLURM array.

Usage:
  python annotate_families.py --run-dir output_5sp --outdir annot \\
      --ec --ecforest-cache ecforest_cache.tsv
  python annotate_families.py --run-dir output_5sp --outdir annot \\
      --prep-fold-list
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Set

from find_neofunctionalization import family_base_dir
from steps.ec_sources import merge_ec_predictions, parse_clean, parse_emapper
from steps.ecforest import (
    consensus_disagreements,
    ec_switch_events,
    family_ec_consensus,
    load_cache,
    write_events_tsv,
)
from utils.seqio import read_fasta

logger = logging.getLogger("family_finder")


def read_families(run_dir: Path) -> Dict[str, Set[str]]:
    """summary.tsv -> family_id -> member gene ids (gene_list column)."""
    families: Dict[str, Set[str]] = {}
    with open(run_dir / "summary.tsv") as f:
        next(f)  # header
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            families[fields[0]] = set(fields[4].split(","))
    return families


def ec_layer(
    run_dir: Path,
    families: Dict[str, Set[str]],
    predictions: Dict[str, dict],
    exclude: Set[str],
    min_agree: float,
    outdir: Path,
) -> None:
    """EC consensus table + EC-switch events across all family trees."""
    events = []
    n_trees = 0
    with open(outdir / "ec_consensus.tsv", "w") as f:
        f.write(
            "family_id\tn_members\tn_ec_predicted\tconsensus_ec\t"
            "n_disagree\tdisagreeing_genes\n"
        )
        for family_id in sorted(families):
            members = families[family_id]
            consensus = family_ec_consensus(members, predictions, min_agree)
            disagree = consensus_disagreements(members, predictions, consensus)
            n_predicted = sum(
                1 for m in members
                if predictions.get(m, {}).get("is_enzyme") and predictions[m].get("ec")
            )
            f.write(
                f"{family_id}\t{len(members)}\t{n_predicted}\t"
                f"{consensus or ''}\t{len(disagree)}\t{','.join(disagree)}\n"
            )

            tree = family_base_dir(run_dir, family_id) / "confirmed_tree.nwk"
            if not tree.exists():
                continue
            n_trees += 1
            try:
                events.extend(
                    ec_switch_events(family_id, tree.read_text(), predictions, exclude)
                )
            except ValueError as e:
                logger.warning(f"{family_id}: tree parse failed ({e})")

    write_events_tsv(events, outdir / "ec_switch_events.tsv")
    logger.info(
        f"EC layer: {len(families)} families -> {outdir / 'ec_consensus.tsv'}; "
        f"{n_trees} trees scanned -> {len(events)} EC-switch events "
        f"({outdir / 'ec_switch_events.tsv'})"
    )


def pick_representative(members: Set[str]) -> str:
    """One deterministic representative per family: first sorted member.

    Deliberately simple prep-stage choice — the cluster-side step may
    substitute the longest/medoid member before folding; the fold list is
    an input manifest, not a scientific decision.
    """
    return sorted(members)[0]


def prep_fold_list(
    run_dir: Path, families: Dict[str, Set[str]], outdir: Path
) -> None:
    """fold_list.tsv: unplaced genes + one representative per family."""
    rows: List[tuple] = []

    unplaced_fa = run_dir / "unplaced_proteins.fa"
    n_unplaced = 0
    if unplaced_fa.exists():
        for gene_id in sorted(read_fasta(str(unplaced_fa))):
            rows.append((gene_id, "unplaced", ""))
            n_unplaced += 1
    else:
        logger.warning(f"No {unplaced_fa} — fold list will contain only representatives")

    for family_id in sorted(families):
        rows.append((pick_representative(families[family_id]), "representative", family_id))

    path = outdir / "fold_list.tsv"
    with open(path, "w") as f:
        f.write("gene_id\trole\tfamily_id\n")
        for gene_id, role, family_id in rows:
            f.write(f"{gene_id}\t{role}\t{family_id}\n")
    logger.info(
        f"Fold list: {n_unplaced} unplaced + {len(families)} representatives "
        f"-> {path}"
    )


def load_ec_predictions(args) -> Dict[str, dict]:
    """Resolve EC predictions from the CLI sources (issue #28 tiering).

    emapper and/or CLEAN are merged (emapper EC authoritative, CLEAN
    confidence attached); the legacy ECForest cache is used only when it is
    the sole source, and warned about otherwise — it failed its
    known-answer gate.
    """
    emapper = parse_emapper(Path(args.emapper)) if args.emapper else {}
    clean = parse_clean(Path(args.clean_csv)) if args.clean_csv else {}
    if emapper or clean:
        if args.ecforest_cache:
            logger.warning("--ecforest-cache ignored: gate-validated sources "
                           "(--emapper/--clean-csv) provided")
        if emapper and clean:
            return merge_ec_predictions(emapper, clean)
        return emapper or clean
    return load_cache(Path(args.ecforest_cache))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--run-dir", required=True, help="Pipeline output dir (has summary.tsv)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ec", action="store_true",
                    help="EC consensus + EC-switch events (needs --emapper, "
                         "--clean-csv and/or --ecforest-cache)")
    ap.add_argument("--emapper", default=None,
                    help="eggNOG-mapper .emapper.annotations (issue #28)")
    ap.add_argument("--clean-csv", default=None,
                    help="CLEAN maxsep CSV (issue #28)")
    ap.add_argument("--ecforest-cache", default=None,
                    help="Legacy ESM-ECForest cache TSV "
                         "(steps.ecforest.save_cache format)")
    ap.add_argument("--exclude-genes", default=None,
                    help="File of gene IDs to exclude from tree mapping "
                         "(truncated models fake EC switches — pass the "
                         "pseudogene-flagged list)")
    ap.add_argument("--min-agree", type=float, default=0.5,
                    help="Min voting-member fraction for a family EC consensus")
    ap.add_argument("--prep-fold-list", action="store_true",
                    help="Write fold_list.tsv (tier-3 stage-1 prep)")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    if not args.ec and not args.prep_fold_list:
        ap.error("Nothing to do: pass --ec and/or --prep-fold-list")
    if args.ec and not (args.emapper or args.clean_csv or args.ecforest_cache):
        ap.error("--ec requires --emapper, --clean-csv and/or --ecforest-cache")

    run_dir, outdir = Path(args.run_dir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    families = read_families(run_dir)
    logger.info(f"Families: {len(families)} from {run_dir / 'summary.tsv'}")

    if args.ec:
        predictions = load_ec_predictions(args)
        logger.info(f"EC predictions: {len(predictions)} genes")
        exclude: Set[str] = set()
        if args.exclude_genes:
            exclude = {l.strip() for l in open(args.exclude_genes) if l.strip()}
            logger.info(f"Excluding {len(exclude)} truncated/flagged genes")
        else:
            logger.warning("No --exclude-genes: truncated models can fake EC "
                           "switches. Pass the pseudogene-flagged list.")
        ec_layer(run_dir, families, predictions, exclude, args.min_agree, outdir)

    if args.prep_fold_list:
        prep_fold_list(run_dir, families, outdir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
