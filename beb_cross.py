#!/usr/bin/env python3
"""codeml branch-site BEB sites x DeepLoc signal windows (issue #24).

The branch-site run (seltest/) used ``cleandata = 0``, so BEB site numbers
are amino-acid positions in the codeml alignment. Whether those are the same
columns as the clan protein-alignment positions ``extract_signal_windows.py``
wrote into ``signal_windows.tsv`` (aln_col_start/aln_col_end) is a question,
not an assumption (#42): the PEPC pilot happened to run both on one
102-taxon / 1,371-column matrix, while the pipeline path runs codeml on the
family codon alignment (five species) and would not. So the cross either
verifies the two stamps agree, or translates through a bridge sequence, or
refuses.

Usage (when the pronghorn bs_codeml jobs land):
  python beb_cross.py --mlc seltest/alt/mlc \\
      --windows signal_windows.tsv -o beb_cross.tsv [--min-prob 0.5] \\
      --site-alignment clan_anchor.aln --window-alignment clan_anchor.aln
Different matrices on the two sides need ``--bridge <seq id present in both>``.

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
from utils.alignment import alignment_id, translate_columns
from utils.seqio import read_fasta
import re
from pathlib import Path
from typing import Dict, List, Tuple

_BEB_HEADER = "Bayes Empirical Bayes"
# Standard genetic code, built once rather than shipped as 64 literals. Stop
# codons become '*' so a premature stop stays visible instead of ending the
# sequence silently.
_BASES = "TCAG"
_AMINO = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
_CODON_TABLE = {
    a + b + c: _AMINO[i * 16 + j * 4 + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
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


def codon_alignment_as_protein(alignment: Dict[str, str]) -> Dict[str, str]:
    """Collapse a codon alignment to one column per codon.

    `steps/codeml.py` runs on `confirmed_codon.afa`, so the pipeline's site
    alignment is on disk as nucleotides while BEB numbers count amino acids.
    Without this step a family codeml result cannot be lined up against a
    protein alignment at all, and "just use the codon file" would put every
    site three times too far along.

    An all-gap codon becomes a gap; a codon that is part gap or contains an
    ambiguity becomes 'X' rather than a guessed residue — the bridge check in
    translate_columns compares residues, and a wrong guess there would pass a
    comparison that ought to fail.
    """
    width = {len(s) for s in alignment.values()}
    if len(width) > 1:
        raise ValueError(f"Alignment has unequal lengths: {sorted(width)}")
    n = width.pop() if width else 0
    if n % 3:
        raise ValueError(
            f"A codon alignment must be a multiple of three columns wide, "
            f"got {n} — this is not the codon matrix codeml was given"
        )
    out = {}
    for name, seq in alignment.items():
        residues = []
        for i in range(0, n, 3):
            codon = seq[i:i + 3].upper()
            if codon == "---":
                residues.append("-")
            elif "-" in codon or codon not in _CODON_TABLE:
                residues.append("X")
            else:
                residues.append(_CODON_TABLE[codon])
        out[name] = "".join(residues)
    return out


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
    ap.add_argument("--site-alignment", default=None,
                    help="FASTA of the alignment the BEB site numbers index "
                         "(the pipeline path: the family codon alignment as "
                         "protein columns). With --window-alignment and "
                         "--bridge the sites are TRANSLATED into the windows' "
                         "coordinate system instead of assumed to match.")
    ap.add_argument("--site-alignment-is-codon", action="store_true",
                    help="--site-alignment is a codon (nucleotide) matrix, as "
                         "steps/codeml.py writes it — collapse it to one "
                         "column per codon first, which is what BEB numbers")
    ap.add_argument("--window-alignment", default=None,
                    help="FASTA of the alignment signal_windows.tsv columns "
                         "index (the clan protein alignment)")
    ap.add_argument("--bridge", default=None,
                    help="Sequence id present in BOTH alignments, carrying "
                         "the same residues — the translation bridge")
    ap.add_argument("--allow-unverified", action="store_true",
                    help="Cross without establishing a shared coordinate "
                         "system. Every row is then marked "
                         "coordinates_verified=False; a zero overlap from an "
                         "unverified cross is not evidence of anything.")
    args = ap.parse_args()

    site_aln = read_fasta(args.site_alignment) if args.site_alignment else None
    if site_aln and args.site_alignment_is_codon:
        site_aln = codon_alignment_as_protein(site_aln)
    window_aln = (read_fasta(args.window_alignment)
                  if args.window_alignment else None)
    # With only the site alignment named, the stamp is still worth emitting:
    # it lets the windows file be checked against it when it carries one.
    site_stamp = str(alignment_id(site_aln)) if site_aln else None

    sites = parse_beb(Path(args.mlc).read_text())
    rows = cross_windows(sites, Path(args.windows), args.min_prob,
                         site_stamp=site_stamp,
                         allow_unverified=args.allow_unverified,
                         site_alignment=site_aln if window_aln else None,
                         window_alignment=window_aln if site_aln else None,
                         bridge=args.bridge)
    header = ["site", "aa", "prob", "translated_site", "n_windows",
              "coordinates_verified", "untranslatable", "window_seqs", "signals"]
    with open(args.out, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write(f"{r['site']}\t{r['aa']}\t{r['prob']:g}\t"
                    f"{'' if r['translated_site'] is None else r['translated_site']}\t"
                    f"{r['n_windows']}\t{r['coordinates_verified']}\t"
                    f"{r['untranslatable']}\t{r['window_seqs']}\t{r['signals']}\n")
    n_in = sum(1 for r in rows if r["n_windows"])
    n_lost = sum(1 for r in rows if r["untranslatable"])
    print(f"{len(sites)} BEB sites -> {len(rows)} kept, "
          f"{n_in} overlap a signal window -> {args.out}")
    if n_lost:
        print(f"{n_lost} site(s) had no position in the window alignment (the "
              "bridge is gapped there) and were counted as untranslatable, "
              "not as outside every window")
    if not rows or not rows[0]["coordinates_verified"]:
        print("WARNING: the two coordinate systems were never shown to match — "
              "these overlaps, including any zero, are not evidence")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
