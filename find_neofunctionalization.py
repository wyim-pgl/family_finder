#!/usr/bin/env python3
"""Neofunctionalization detection (issue #16).

Stage 1 (cheap, anywhere): map DeepLoc localizations onto family gene trees,
find retargeting branches -> retargeting_events.tsv
Stage 2 (needs PAML): branch-site Model A vs null on each event branch,
LRT + BH-FDR -> retargeting_selection.tsv

Retargeting + significant positive selection = neofunctionalization candidate.

Usage:
  # stage 1 only (no codeml runs)
  python find_neofunctionalization.py --run-dir output_5sp \\
      --deeploc-cache deeploc_cache.tsv --outdir neofunc --events-only

  # both stages
  python find_neofunctionalization.py --run-dir output_5sp \\
      --deeploc-cache deeploc_cache.tsv --exclude-genes truncated.txt \\
      --outdir neofunc --config config_5sp.json
"""

import argparse
import logging
import sys
from pathlib import Path

from config import Config
from steps.codeml import (
    START_OMEGAS,
    best_lnl_run,
    benjamini_hochberg,
    generate_branch_site_ctl,
    lrt_pvalue,
    parse_lnl,
    run_codeml,
    write_marked_tree,
)
from steps.deeploc import load_cache
from steps.retargeting import find_retargeting_events, write_events_tsv

logger = logging.getLogger("family_finder")


def family_base_dir(run_dir: Path, family_id: str) -> Path:
    """R{round}_{og} -> <run_dir>/round_XX/orthogroups/<og> (same resolution
    as steps.codeml.run_codeml_on_families)."""
    parts = family_id.split("_", 1)
    round_num = int(parts[0][1:])
    return run_dir / f"round_{round_num:02d}" / "orthogroups" / parts[1]


def collect_events(run_dir: Path, predictions, exclude, min_prob, max_families=None):
    """Stage 1: retargeting events across all families in summary.tsv."""
    events = []
    n_scanned = 0
    with open(run_dir / "summary.tsv") as f:
        next(f)
        for line in f:
            family_id = line.split("\t", 1)[0]
            base = family_base_dir(run_dir, family_id)
            tree = base / "confirmed_tree.nwk"
            if not tree.exists():
                continue
            n_scanned += 1
            if max_families and n_scanned > max_families:
                break
            try:
                events.extend(
                    find_retargeting_events(
                        family_id, tree.read_text(), predictions, exclude, min_prob
                    )
                )
            except ValueError as e:
                logger.warning(f"{family_id}: tree parse failed ({e})")
    logger.info(f"Scanned {n_scanned} family trees -> {len(events)} retargeting events")
    return events


def test_events(events, run_dir: Path, outdir: Path, config: Config):
    """Stage 2: branch-site Model A vs null per event. Returns result rows."""
    rows = []
    for i, ev in enumerate(events):
        base = family_base_dir(run_dir, ev.family_id)
        codon_aln = base / "confirmed_codon.afa"
        tree = base / "confirmed_tree.nwk"
        if not codon_aln.exists():
            logger.warning(f"{ev.family_id}: no codon alignment, skipping event")
            continue
        work = outdir / "codeml" / f"{ev.family_id}_ev{i}"
        try:
            marked = write_marked_tree(tree, ev.clade_leaves, work / "marked.nwk")
            # Null: omega fixed at 1, one fit. Alternative: refit from
            # several starting omegas and keep the best lnL — codeml's
            # branch-site optimizer lands in local optima otherwise.
            null_dir = work / "null"
            ctl = generate_branch_site_ctl(ev.family_id, codon_aln, marked,
                                           null_dir, True)
            run_codeml(ctl, null_dir, config)
            lnl_null = parse_lnl(null_dir / "results.txt")

            alt_runs = []
            for w0 in START_OMEGAS:
                sub = work / f"alt_w{w0:g}"
                ctl = generate_branch_site_ctl(ev.family_id, codon_aln, marked,
                                               sub, False, start_omega=w0)
                run_codeml(ctl, sub, config)
                alt_runs.append((f"w{w0:g}", parse_lnl(sub / "results.txt")))
            best_label, lnl_alt = best_lnl_run(alt_runs)
            spread = lnl_alt - min(l for _, l in alt_runs)
            if spread > 0.5:
                logger.warning(
                    f"{ev.family_id}: branch-site alt lnL spread {spread:.2f} "
                    f"across starting omegas (best {best_label}) — local optima "
                    "present, multi-start was necessary"
                )
            lnls = {False: lnl_alt, True: lnl_null}
            p = lrt_pvalue(lnls[False], lnls[True])
            rows.append((ev, lnls[False], lnls[True], p))
            logger.info(f"{ev.family_id} {ev.parent_state}->{ev.child_state}: p={p:.3g}")
        except Exception as e:
            logger.error(f"{ev.family_id} event {i} failed: {e}")
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--run-dir", required=True, help="Pipeline output dir (has summary.tsv)")
    ap.add_argument("--deeploc-cache", required=True,
                    help="DeepLoc cache TSV (steps.deeploc.save_cache format)")
    ap.add_argument("--exclude-genes", default=None,
                    help="File of gene IDs to exclude (truncated models — one per line). "
                         "MANDATORY safeguard for real analyses: truncated models "
                         "produce fake retargeting calls.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--config", default=None, help="Pipeline config JSON (for codeml_bin)")
    ap.add_argument("--min-prob", type=float, default=None,
                    help="Override deeploc_min_prob from config")
    ap.add_argument("--events-only", action="store_true",
                    help="Stage 1 only: write retargeting_events.tsv, skip codeml")
    ap.add_argument("--max-families", type=int, default=None, help="Limit (testing)")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    config = Config.from_json(args.config) if args.config else Config()
    min_prob = args.min_prob if args.min_prob is not None else config.deeploc_min_prob
    run_dir, outdir = Path(args.run_dir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    predictions = load_cache(Path(args.deeploc_cache))
    logger.info(f"DeepLoc predictions: {len(predictions)} genes")

    exclude = set()
    if args.exclude_genes:
        exclude = {l.strip() for l in open(args.exclude_genes) if l.strip()}
        logger.info(f"Excluding {len(exclude)} truncated/flagged genes")
    else:
        logger.warning("No --exclude-genes: truncated models WILL produce fake "
                       "retargeting calls. Pass the pseudogene-flagged list.")

    events = collect_events(run_dir, predictions, exclude, min_prob, args.max_families)
    write_events_tsv(events, outdir / "retargeting_events.tsv")
    logger.info(f"Wrote {outdir / 'retargeting_events.tsv'}")

    if args.events_only or not events:
        return 0

    rows = test_events(events, run_dir, outdir, config)
    qvals = benjamini_hochberg([r[3] for r in rows])
    with open(outdir / "retargeting_selection.tsv", "w") as f:
        f.write("family_id\tparent_state\tchild_state\tn_clade_genes\tclade_genes\t"
                "lnl_alt\tlnl_null\tlrt_p\tfdr_q\n")
        for (ev, lnl_a, lnl_0, p), q in sorted(
            zip(rows, qvals), key=lambda x: x[0][3]
        ):
            f.write(f"{ev.family_id}\t{ev.parent_state}\t{ev.child_state}\t"
                    f"{len(ev.clade_leaves)}\t{ev.clade_id}\t"
                    f"{lnl_a:.4f}\t{lnl_0:.4f}\t{p:.4g}\t{q:.4g}\n")
    n_sig = sum(1 for q in qvals if q < 0.05)
    logger.info(f"Wrote {outdir / 'retargeting_selection.tsv'} "
                f"({len(rows)} tested, {n_sig} at FDR<0.05)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
