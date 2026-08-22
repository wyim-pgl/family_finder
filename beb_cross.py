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
import sys
from pathlib import Path as _P
sys.path.insert(0, str(_P(__file__).resolve().parent))
from utils.alignment import translate_columns
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


def cross_windows(sites, windows_tsv: Path, min_prob: float = 0.0,
                  site_stamp: str = None, allow_unverified: bool = False,
                  site_alignment: dict = None, window_alignment: dict = None,
                  bridge: str = None) -> List[dict]:
    """Overlap each BEB site with the signal windows' alignment columns.

    The two inputs do not necessarily share a coordinate system. BEB site
    numbers are columns of the family's codon alignment (untrimmed, family
    taxon set); the windows were written in whatever alignment
    extract_signal_windows.py was given, which for the PEPC clan was trimmed
    from 1,468 columns to 876 over a different taxon set. Overlapping them by
    interval alone always returns a number, and zero overlaps is exactly what
    subfunctionalization.classify() reports as evidence against
    neofunctionalization — so a coordinate mismatch manufactures a conclusion.

    Pass `site_stamp` (see utils.alignment.alignment_id) and let the windows
    file carry an `aln_stamp` column. When they cannot be compared, either
    translate first (utils.alignment.translate_columns) or set
    allow_unverified, which marks every row rather than hiding the doubt.
    """
    windows = []
    stamps = set()
    with open(windows_tsv, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            try:
                lo = int(row["aln_col_start"])
                hi = int(row["aln_col_end"])
            except (KeyError, ValueError):
                continue
            stamp = (row.get("aln_stamp") or "").strip()
            if stamp and stamp != "-":
                stamps.add(stamp)
            windows.append((lo, hi, row["seq_id"], row.get("deeploc_signals", "")))

    window_stamp = stamps.pop() if len(stamps) == 1 else None
    verified = bool(site_stamp and window_stamp and site_stamp == window_stamp)

    # The honest resolution when the two coordinate systems really do differ:
    # carry the sites across through a sequence present in both alignments,
    # rather than refusing or pretending they agree. steps/codeml.py runs on
    # the family codon alignment while the windows come from the clan protein
    # alignment, so this is the normal case for the pipeline path, not an edge.
    translate = site_alignment is not None and window_alignment is not None
    if translate:
        if bridge is None:
            raise ValueError(
                "Translating between alignments needs a bridge sequence "
                "present in both (bridge=...)"
            )
        verified = True
    elif not verified and not allow_unverified:
        raise ValueError(
            "Refusing to cross BEB sites with signal windows whose coordinate "
            f"system is unconfirmed (site_stamp={site_stamp!r}, "
            f"window_stamp={window_stamp!r}). Column numbers from different "
            "alignments are not comparable, and a mismatch here silently "
            "returns zero overlaps, which is read downstream as evidence "
            "against neofunctionalization. Translate with "
            "utils.alignment.translate_columns, or pass allow_unverified=True "
            "to record the doubt in the output."
        )

    kept = [(site, aa, prob) for site, aa, prob in sites if prob >= min_prob]
    if translate:
        moved = translate_columns([s for s, _, _ in kept],
                                  site_alignment, window_alignment, via=bridge)
    else:
        moved = [s for s, _, _ in kept]

    rows = []
    for (site, aa, prob), target in zip(kept, moved):
        if target is None:
            rows.append({
                "site": site, "aa": aa, "prob": prob,
                "n_windows": 0, "window_seqs": "", "signals": "",
                "coordinates_verified": verified,
                "translated_site": None, "untranslatable": True,
            })
            continue
        hit = [(s, sig) for lo, hi, s, sig in windows if lo <= target <= hi]
        seqs = sorted({s for s, _ in hit})
        signals = sorted({sig for _, sig in hit if sig})
        rows.append({
            "site": site, "aa": aa, "prob": prob,
            "n_windows": len(hit),
            "coordinates_verified": verified,
            "translated_site": target if translate else None,
            "untranslatable": False,
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
