"""Chunked ESMFold structure prediction for the tier-3 screen (issue #34).

ESMFold was chosen over ProstT5 after measuring both. The deciding fact is not
throughput: `tier3_assign` gates entirely on `alntmscore`, and a ProstT5
database has no CA coordinates, so it cannot produce one. Using ProstT5 would
mean a new bits-based decision path whose thresholds have nothing established to
anchor them, while ESMFold produces real structures that the existing TM-score
gate consumes unchanged — and its pLDDT feeds region_disorder.py and
flag_low_plddt as a second use.

Measured 2026-08-22 on the unplaced pool (RTX 4090, fp16, chunk_size 64):

- **1.68 s/sequence at 210 aa**, 20/20 folded, no OOM.
- **Peak VRAM 13.3 GB regardless of length** — the model weights dominate, the
  protein barely contributes, so batching gains little and length is safe until
  it is not.
- Combined with the recorded 75 s at 930 aa, cost goes as **length^2.55**.

That exponent is why MAX_LEN exists. In the 9,786-sequence pool the 45 sequences
above 1,100 aa cost as much as the other 9,741 combined (8.7 h versus 4.1 h);
the longest is 3,956 aa. Capping is not an optimisation, it is most of the bill.
Capped, the whole pool folds in about 4.4 hours — cheaper than ProstT5's 11.7,
which is the opposite of what per-sequence figures at 930 aa suggest.

Structure of the run follows steps/hmm_chunks.py and steps/prostt5_chunks.py:
split, fold per chunk, verify every output, mark the chunk done. ESMFold has no
internal checkpoint, so an unverified chunk must never be treated as finished —
a partially written directory looks exactly like a complete one.
"""
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

MAX_LEN = 1100
DONE_NAME = "chunk.done"

# Measured anchor: 20 sequences of ~210 aa folded in 33.6 s wall clock.
REF_LENGTH = 210.0
REF_SECONDS = 1.68
LENGTH_EXPONENT = 2.55

PLDDT_MIN, PLDDT_MAX = 0.0, 100.0
_B_FACTOR = slice(60, 66)


def partition_by_length(sequences: Dict[str, str],
                        max_len: int = MAX_LEN) -> Tuple[Dict[str, str], Dict[str, int]]:
    """Split into what will be folded and what is too long, keeping both.

    Over-length sequences are returned with their lengths rather than dropped:
    they are a reportable fraction of the pool, not noise.
    """
    foldable = {k: v for k, v in sequences.items() if len(v) <= max_len}
    oversize = {k: len(v) for k, v in sequences.items() if len(v) > max_len}
    return foldable, oversize


def estimate_runtime(lengths: Iterable[int],
                     max_len: Optional[int] = None) -> float:
    """Seconds to fold these sequences, from the measured power law.

    `max_len` excludes over-length sequences from the estimate, which is the
    point of having it — the difference between capped and uncapped is where
    the cost of this screen actually lives.
    """
    total = 0.0
    for length in lengths:
        if max_len is not None and length > max_len:
            continue
        total += REF_SECONDS * (length / REF_LENGTH) ** LENGTH_EXPONENT
    return total


def verify_pdb(path: Union[str, Path],
               expected_residues: Optional[int] = None) -> float:
    """Mean pLDDT, raising when the structure is not usable.

    ESMFold writes per-residue pLDDT into the B-factor column. A structure that
    is missing, empty, has no ATOM records, carries something other than pLDDT
    there, or stops short of the sequence is a failure that must not be counted
    as coverage.
    """
    path = Path(path)
    if not path.exists():
        raise ValueError(f"Structure missing: {path}")
    text = path.read_text()
    if not text.strip():
        raise ValueError(f"Structure file is empty: {path}")

    values: List[float] = []
    seen_residues = set()
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            b = float(line[_B_FACTOR])
        except (ValueError, IndexError):
            raise ValueError(f"Unparseable B-factor column in {path}: {line!r}")
        if not PLDDT_MIN <= b <= PLDDT_MAX:
            raise ValueError(
                f"B-factor {b} outside the pLDDT range in {path} — the column "
                "holds something other than pLDDT"
            )
        values.append(b)
        seen_residues.add(line[22:26].strip())

    if not values:
        raise ValueError(f"Found no ATOM records in {path}")
    if expected_residues is not None and len(seen_residues) != expected_residues:
        raise ValueError(
            f"Structure {path} has {len(seen_residues)} residues, expected "
            f"{expected_residues} — truncated fold"
        )
    return sum(values) / len(values)


def run_esmfold(sequences: Dict[str, str], out_dir: Path, **kwargs) -> None:
    """Fold sequences, writing <out_dir>/<name>.pdb.

    Module-level so tests monkeypatch it; the real implementation lives on the
    GPU host (see gpu:~/pepc_pilot/fold_pepc.py) and is invoked by the caller.
    """
    raise NotImplementedError(
        "run_esmfold must be provided by the caller or monkeypatched; folding "
        "runs on the GPU host, not inside the pipeline process"
    )


def fold_chunk(sequences: Dict[str, str], out_dir: Union[str, Path],
               **kwargs) -> dict:
    """Fold one chunk, verify every structure, then mark it done.

    Resumable: sequences whose PDB already exists are skipped. The sentinel is
    written only after every structure verifies, so a killed run leaves a
    directory that is visibly unfinished rather than one that merely looks
    complete.
    """
    out_dir = Path(out_dir)
    done = out_dir / DONE_NAME
    pending = {k: v for k, v in sequences.items()
               if not (out_dir / f"{k}.pdb").exists()}
    n_skipped = len(sequences) - len(pending)

    run_esmfold(pending, out_dir, **kwargs)

    plddt: Dict[str, float] = {}
    for name, seq in sequences.items():
        plddt[name] = verify_pdb(out_dir / f"{name}.pdb",
                                 expected_residues=len(seq))

    done.write_text("")
    return {"n_folded": len(sequences), "n_skipped": n_skipped, "plddt": plddt}


def completed_chunks(dirs: Sequence[Union[str, Path]]) -> List[Path]:
    """Chunk directories carrying the sentinel, i.e. safe to merge."""
    return [Path(d) for d in dirs if (Path(d) / DONE_NAME).exists()]
