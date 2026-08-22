#!/usr/bin/env python3
"""DeepLoc attention windows, mapped to alignment columns (issue #24 step 2).

DeepLoc run with per-residue attribution (`-p`) writes one `alpha_<id>.csv` per
sequence: the attention the localization head paid to each residue. A run of
unusually high attention is a candidate localization-signal region. This script
turns those profiles into windows, expresses each window in ALIGNMENT columns,
and writes `signal_windows.tsv`.

Window rule (unchanged from the PEPC pilot that produced the SF3 result):
per-sequence threshold = mean + 2SD, runs above it, gaps of <= 2 aa merged,
windows shorter than 5 aa dropped. Coordinates are 1-based throughout.

**Why this file is in the repo** (issue #42). It used to live only on the GPU
box (`~/pepc_pilot/`), so `aln_col_start` had a provenance nothing could check.
`beb_cross.py` asserted in prose that those columns were the same coordinate
system as codeml's BEB sites; for the PEPC pilot that happened to be true
(one 102-taxon / 1,371-column matrix ran every axis), and for the pipeline
path — codeml on the five-species family codon alignment — it is not. So every
row now carries `aln_stamp`, the `utils.alignment.alignment_id()` of the
alignment its columns refer to, and a consumer can verify instead of assume.

Two other things the original did quietly and this one does not: identifiers
are matched through `utils.gene_ids.match_ids` rather than a local
lowercase-and-strip rule (DeepLoc mangles `Aco_Aco010025.1` into
`aco_aco0100251`), and a sequence whose windows cannot be placed on the
alignment is reported rather than written out with an empty column.

Usage:
  python extract_signal_windows.py --alignment clan_anchor.aln \\
      --attention-dir deeploc_p_out --results deeploc_p_out/results_*.csv \\
      -o signal_windows.tsv
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils.alignment import alignment_id  # noqa: E402
from utils.gene_ids import match_ids  # noqa: E402
from utils.seqio import read_fasta  # noqa: E402

MIN_LEN = 5
MERGE_GAP = 2
N_SD = 2.0

COLUMNS = ["seq_id", "start", "end", "peak_pos", "peak_attention",
           "mean_attention", "aln_col_start", "aln_col_end",
           "deeploc_signals", "aln_stamp"]


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values)


def _pstdev(values: Sequence[float]) -> float:
    m = _mean(values)
    return (sum((v - m) ** 2 for v in values) / len(values)) ** 0.5


def read_attention(path: Path) -> List[float]:
    """One DeepLoc `alpha_<id>.csv` -> the per-residue attention profile.

    The layout is exactly two columns (AA, Alpha). Anything else is a
    different file than this parser understands, and guessing at it would put
    silently wrong numbers into every window below.
    """
    values: List[float] = []
    with open(path) as f:
        next(f, None)  # header
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split(",")
            if len(fields) != 2:
                raise ValueError(
                    f"{path}: expected two columns (AA, Alpha), got "
                    f"{len(fields)}: {line[:60]!r}"
                )
            values.append(float(fields[1]))
    return values


def extract_windows(values: Sequence[float], *, n_sd: float = N_SD,
                    merge_gap: int = MERGE_GAP,
                    min_len: int = MIN_LEN) -> List[Tuple[int, int]]:
    """Runs of above-threshold attention as 1-based `(start, end)` residues.

    The threshold is per-sequence (mean + n_sd * SD) rather than global,
    because attention profiles are normalised per protein — a family-wide
    cutoff would just rank proteins by how peaked their profile is.
    """
    if len(values) < min_len:
        return []
    sd = _pstdev(values)
    if sd == 0.0:
        return []  # a flat profile has no peaks, whatever the mean is
    threshold = _mean(values) + n_sd * sd
    above = [i + 1 for i, v in enumerate(values) if v > threshold]

    merged: List[List[int]] = []
    for pos in above:
        if merged and pos - merged[-1][1] <= merge_gap + 1:
            merged[-1][1] = pos
        else:
            merged.append([pos, pos])
    return [(s, e) for s, e in merged if e - s + 1 >= min_len]


def _residue_columns(aligned_seq: str) -> Dict[int, int]:
    """This sequence's ungapped residue number -> alignment column (1-based)."""
    out, n = {}, 0
    for col, ch in enumerate(aligned_seq, start=1):
        if ch != "-":
            n += 1
            out[n] = col
    return out


def map_windows_to_columns(windows: Dict[str, List[Tuple[int, int]]],
                           alignment: Dict[str, str],
                           attention: Optional[Dict[str, Sequence[float]]] = None,
                           signals: Optional[Dict[str, str]] = None,
                           report_skipped: bool = False):
    """Express each sequence's windows in alignment columns, stamped.

    A sequence the alignment does not have, or a window running past its
    residues, is SKIPPED AND NAMED — the original wrote an empty
    `aln_col_start` and carried on, which is how a shrinking window set became
    indistinguishable from an absent signal.

    Set `report_skipped` to get `(rows, skipped)` instead of just the rows.
    """
    stamp = str(alignment_id(alignment))
    matched = match_ids(sorted(windows), list(alignment)).mapping

    rows: List[dict] = []
    skipped: Dict[str, str] = {}
    for seq_id in sorted(windows):
        ref = matched.get(seq_id)
        if ref is None:
            skipped[seq_id] = "not in the alignment"
            continue
        cols = _residue_columns(alignment[ref])
        profile = (attention or {}).get(seq_id)
        for start, end in windows[seq_id]:
            if start not in cols or end not in cols:
                skipped[seq_id] = (
                    f"window {start}-{end} runs past the {len(cols)} aligned "
                    "residues of this sequence"
                )
                continue
            row = {
                "seq_id": ref,
                "start": start,
                "end": end,
                "aln_col_start": cols[start],
                "aln_col_end": cols[end],
                "deeploc_signals": (signals or {}).get(seq_id, ""),
                "aln_stamp": stamp,
            }
            if profile:
                segment = list(profile[start - 1:end])
                peak = max(range(start, end + 1), key=lambda p: profile[p - 1])
                row.update(peak_pos=peak,
                           peak_attention=f"{profile[peak - 1]:.6f}",
                           mean_attention=f"{_mean(segment):.6f}")
            else:
                row.update(peak_pos="", peak_attention="", mean_attention="")
            rows.append(row)

    return (rows, skipped) if report_skipped else rows


def write_windows(rows: List[dict], path: Path) -> None:
    """signal_windows.tsv, header always present so an empty result is visible."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in COLUMNS})


def read_signals(path: Path) -> Dict[str, str]:
    """DeepLoc `results_*.csv` -> gene id -> its predicted Signals string."""
    out: Dict[str, str] = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            out[row["Protein_ID"]] = row.get("Signals", "")
    return out


def collect_attention(attention_dir: Path,
                      seq_ids: Sequence[str]) -> Tuple[Dict[str, List[float]],
                                                       List[str]]:
    """Load one attention profile per sequence, naming the ones with none.

    DeepLoc derives the filename from the identifier by lowercasing it and
    deleting everything outside [a-z0-9_], so the match goes through
    `utils.gene_ids.match_ids` — the one normaliser — rather than a copy of
    that rule kept here.
    """
    files = {p.name[len("alpha_"):-len(".csv")]: p
             for p in sorted(Path(attention_dir).glob("alpha_*.csv"))}
    matched = match_ids(seq_ids, list(files))
    profiles = {seq_id: read_attention(files[key])
                for seq_id, key in matched.mapping.items()}
    return profiles, list(matched.unmatched)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--alignment", required=True,
                    help="Protein alignment (FASTA) whose columns the output "
                         "reports — its content stamp goes in every row")
    ap.add_argument("--attention-dir", required=True,
                    help="Directory of DeepLoc -p alpha_<id>.csv profiles")
    ap.add_argument("--results", default=None,
                    help="DeepLoc results_*.csv, for the Signals column")
    ap.add_argument("-o", "--out", required=True, help="signal_windows.tsv")
    ap.add_argument("--n-sd", type=float, default=N_SD,
                    help=f"Threshold = mean + this many SD (default {N_SD})")
    ap.add_argument("--min-len", type=int, default=MIN_LEN,
                    help=f"Shortest window kept, in residues (default {MIN_LEN})")
    ap.add_argument("--merge-gap", type=int, default=MERGE_GAP,
                    help=f"Gap in residues still merged (default {MERGE_GAP})")
    ap.add_argument("--max-unmatched", type=float, default=None,
                    help="Refuse when more than this fraction of sequences "
                         "have no attention profile (default: report and "
                         "continue)")
    args = ap.parse_args()

    alignment = read_fasta(args.alignment)
    profiles, no_profile = collect_attention(Path(args.attention_dir),
                                             list(alignment))
    if args.max_unmatched is not None and alignment:
        lost = len(no_profile) / len(alignment)
        if lost > args.max_unmatched:
            raise SystemExit(
                f"{len(no_profile)} of {len(alignment)} sequences have no "
                f"attention profile ({lost:.0%} > {args.max_unmatched:.0%}): "
                f"{no_profile[:5]} — refusing rather than writing a window set "
                "whose gaps read as absent signals"
            )

    signals = read_signals(Path(args.results)) if args.results else {}
    windows = {seq_id: extract_windows(profile, n_sd=args.n_sd,
                                       merge_gap=args.merge_gap,
                                       min_len=args.min_len)
               for seq_id, profile in profiles.items()}
    windows = {k: v for k, v in windows.items() if v}

    rows, skipped = map_windows_to_columns(windows, alignment,
                                           attention=profiles, signals=signals,
                                           report_skipped=True)
    write_windows(rows, Path(args.out))

    stamp = alignment_id(alignment)
    print(f"alignment {args.alignment}: {stamp}")
    print(f"{len(rows)} windows across {len({r['seq_id'] for r in rows})}/"
          f"{len(alignment)} sequences -> {args.out}")
    if no_profile:
        print(f"{len(no_profile)} sequence(s) had no attention profile: "
              f"{no_profile[:5]}")
    for seq_id, reason in sorted(skipped.items()):
        print(f"  skipped {seq_id}: {reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
