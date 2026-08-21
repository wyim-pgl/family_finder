#!/usr/bin/env python3
"""codeml branch-site BEB sites x DeepLoc signal windows (issue #24).

The branch-site run (seltest/) used ``cleandata = 0``, so BEB site numbers
are amino-acid positions in the codeml alignment — the same coordinate
system as the clan protein-alignment columns that
``extract_signal_windows.py`` wrote into ``signal_windows.tsv``
(aln_col_start/aln_col_end). The cross is therefore a direct interval
overlap, no remapping.

Usage (when the pronghorn bs_codeml jobs land):
  python beb_cross.py --mlc seltest/alt/mlc \\
      --windows signal_windows.tsv -o beb_cross.tsv [--min-prob 0.5]

Output: one row per BEB site — site, aa, prob, overlapping window count,
window seq ids, their DeepLoc signal types. Sites with no window overlap
are kept (the negative result is a result: MEME already found zero
episodic sites in the SF3 attention region col 165-209).
"""

import argparse
import csv
import re
from pathlib import Path
from typing import List, Tuple

_BEB_HEADER = "Bayes Empirical Bayes"
_SITE_RE = re.compile(r"^\s*(\d+)\s+([A-Z*-])\s+([01]\.\d+)\**\s*$")


def parse_beb(text: str) -> List[Tuple[int, str, float]]:
    """BEB 'Positive sites' block -> [(site, aa, prob)] (stars stripped).

    The block ends at the first line that is neither a site row nor the
    'Positive sites' sub-header — codeml separates it with a blank line.
    """
    sites: List[Tuple[int, str, float]] = []
    in_beb = in_sites = False
    for line in text.splitlines():
        if _BEB_HEADER in line:
            in_beb = True
            continue
        if in_beb and line.strip().startswith("Positive sites"):
            in_sites = True
            continue
        if in_sites:
            m = _SITE_RE.match(line)
            if m:
                sites.append((int(m.group(1)), m.group(2), float(m.group(3))))
            else:
                break  # blank line / next section: block over
    return sites


def cross_windows(sites, windows_tsv: Path, min_prob: float = 0.0) -> List[dict]:
    """Overlap each BEB site with the signal windows' alignment columns."""
    windows = []
    with open(windows_tsv, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                lo = int(row["aln_col_start"])
                hi = int(row["aln_col_end"])
            except (KeyError, ValueError):
                continue
            windows.append((lo, hi, row["seq_id"], row.get("deeploc_signals", "")))

    rows = []
    for site, aa, prob in sites:
        if prob < min_prob:
            continue
        hit = [(s, sig) for lo, hi, s, sig in windows if lo <= site <= hi]
        seqs = sorted({s for s, _ in hit})
        signals = sorted({sig for _, sig in hit if sig})
        rows.append({
            "site": site, "aa": aa, "prob": prob,
            "n_windows": len(hit),
            "window_seqs": ",".join(seqs),
            "signals": ";".join(signals),
        })
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--mlc", required=True, help="codeml Model-A (alt) main output")
    ap.add_argument("--windows", required=True, help="signal_windows.tsv")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--min-prob", type=float, default=0.0,
                    help="Keep BEB sites with P(w>1) >= this (0 = all)")
    args = ap.parse_args()

    sites = parse_beb(Path(args.mlc).read_text())
    rows = cross_windows(sites, Path(args.windows), args.min_prob)
    with open(args.out, "w") as f:
        f.write("site\taa\tprob\tn_windows\twindow_seqs\tsignals\n")
        for r in rows:
            f.write(f"{r['site']}\t{r['aa']}\t{r['prob']:g}\t"
                    f"{r['n_windows']}\t{r['window_seqs']}\t{r['signals']}\n")
    n_in = sum(1 for r in rows if r["n_windows"])
    print(f"{len(sites)} BEB sites -> {len(rows)} kept, "
          f"{n_in} overlap a signal window -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
