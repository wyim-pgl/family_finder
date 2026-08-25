"""Chunked ProstT5 3Di database construction for the tier-3 screen (issue #34).

ProstT5 turns a protein sequence into a 3Di structural alphabet without folding
it, which is how tier-3 reaches the unplaced pool at genome scale. Three
measured properties of this foldseek build shape the module:

- **Judge GPU use by the CLOCK, not by the log.** `--gpu 1` is parsed (a bad
  value is rejected: `--gpu 999` errors out) but never reaches this path: the
  log prints `Use GPU 0` and the flag is absent from the echoed command line
  whether or not you passed it. That line belongs to a different mmseqs knob.
  What actually drives ProstT5 is ggml, which registers its own `CUDA0`
  backend and says so one line later.

  So `Use GPU 0` tells you nothing, and it cost this project a wrong diagnosis
  five times over. The measurement that matters is elapsed time. On 112
  sequences, same commit (18478813), same `cudaMalloc` symbol count (39), both
  logging `Use GPU 0` and `CUDA0`:

      ~/bin/foldseek-cuda (provenance unknown)   497 s
      rebuilt with -DENABLE_CUDA=1               12 s     41x

  The 37.6x recorded elsewhere for this binary was right; it was simply never
  being obtained. Rebuilding turns the 459,398-gene panel from ~23 days into
  ~14 hours, which is the difference between impossible and overnight.

  Rebuild recipe (foldseek wiki, plus two traps it omits):

      micromamba create -n fsbuild -c conda-forge -c nvidia \
        cmake ninja cuda-nvcc cuda-cudart-dev libcublas-dev \
        libcublas-static cuda-version=12.6 rust
      cmake -DCMAKE_BUILD_TYPE=RELEASE -DENABLE_CUDA=1 \
        -DCMAKE_CUDA_ARCHITECTURES=native \
        -DCUDAToolkit_ROOT=$E -DCUDA_CUDART=$E/lib/libcudart.so \
        -DCMAKE_EXE_LINKER_FLAGS="-L$E/lib -Wl,-rpath,$E/lib" ..

  Trap 1: `rust` is not in the wiki's conda list but mmseqs' corrosion needs
  it, otherwise cmake stops at `Configuring incomplete`.
  Trap 2: this host carries three CUDA runtimes (11.5 in /usr/lib, 12.6 in
  /usr/local, 12.6 in conda) and plain `nvcc` resolves to 11.5. Without
  `CUDA_CUDART` pointed at a 12.x runtime the link dies on `undefined
  reference to cudaGetDeviceProperties_v2`, a symbol CUDA 12 introduced.

  ⚠️ 2026-08-25 correction, measured on the rebuilt binary: `--gpu 1` DOES
  decide the path there. Without it, createdb enters the CPU
  ProstT5ForkRunner even on a CUDA build — one forked worker per thread,
  EACH loading the ~3 GB model (16 workers on a 62 GB box: the children die,
  the parent blocks forever in msgsnd on a full SysV message queue, 0 s of
  CPU time, no error). With `--gpu 1`: ggml CUDA0 engaged, 0.034 s/seq on
  the 15sp panel (10,000-seq chunk in 340 s). The clock rule above stands —
  the msgsnd hang also logs `Use GPU 0`.
- **~4.3 s/sequence at 16 CPU threads** (PEPC, ~950 aa, so near the upper
  bound); ten thousand sequences is roughly twelve hours. On GPU
  (`--gpu 1`): 0.034 s/seq measured — the same ten thousand in ~6 minutes.
- **No internal checkpoint.** A createdb killed at hour eleven leaves nothing.

Hence chunking, in the shape of steps/hmm_chunks.py: split the FASTA, build one
database per chunk, search each, merge the *search results* — never the
databases.

The rule this module exists to enforce:

    Passing a `.gguf` FILE to `--prostt5-model` silently produces no 3Di.

createdb exits in 0.00 s with no error and no warning, prints the path back as
though it were accepted, and writes a database containing amino acids. A screen
run that way is a sequence search believing it is a structure search, and it
finds nothing. `ls <db>_ss` catches it once; it does not survive dozens of
chunks. So the check lives here and fails hard, and the model path is rejected
before createdb is ever reached if it is not a directory.

One further consequence, recorded for callers: a ProstT5 database has no CA
coordinates, so `alntmscore` cannot be computed from it. Requesting it fails
with `No datafile could be found for <db>_ca`. Decisions over ProstT5 hits must
use bit scores, not TM-scores.
"""

import logging
import subprocess
from pathlib import Path
from typing import List, Optional, Union

logger = logging.getLogger("family_finder")

SS_SUFFIX = "_ss"
CHUNK_GLOB = "chunk_*.m8"
DONE_SUFFIX = ".done"
CHUNK_FASTA_FMT = "chunk_{:04d}.fa"


# ------------------------------------------------------------------ split --

def split_fasta(fasta_path: Union[str, Path], out_dir: Union[str, Path],
                chunk_size: int) -> List[Path]:
    """Split a FASTA into chunks of `chunk_size` records.

    A chunk boundary is only ever placed before a `>` line, so a record whose
    sequence is wrapped across lines is never cut in half.
    """
    if chunk_size < 1:
        raise ValueError(f"chunk_size must be >= 1, got {chunk_size}")
    fasta_path = Path(fasta_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    records: List[List[str]] = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                records.append([line])
            elif records:
                records[-1].append(line)
    if not records:
        raise ValueError(f"{fasta_path} contains no sequences")

    chunks: List[Path] = []
    for i in range(0, len(records), chunk_size):
        path = out_dir / CHUNK_FASTA_FMT.format(len(chunks) + 1)
        with open(path, "w") as out:
            for rec in records[i:i + chunk_size]:
                out.writelines(rec)
        chunks.append(path)
    logger.info(f"Split {len(records)} sequences into {len(chunks)} chunks "
                f"of <= {chunk_size} in {out_dir}")
    return chunks


# ------------------------------------------------------------------ guard --

def _entry_count(index_path: Path) -> int:
    return sum(1 for line in index_path.read_text().splitlines() if line.strip())


def verify_3di_db(db_path: Union[str, Path]) -> int:
    """Assert that `db_path` really carries 3Di, and return its entry count.

    Raises RuntimeError describing what is wrong. This is the trap-1 guard;
    every chunk must pass it before its results are trusted.
    """
    db_path = Path(db_path)
    ss = db_path.with_name(db_path.name + SS_SUFFIX)
    ss_index = db_path.with_name(db_path.name + SS_SUFFIX + ".index")
    index = db_path.with_name(db_path.name + ".index")

    if not ss.exists():
        raise RuntimeError(
            f"{ss} is missing — createdb produced no 3Di database for "
            f"{db_path.name}. This is what passing a .gguf FILE to "
            "--prostt5-model does: it exits in 0.00 s without an error and "
            "writes amino acids. Pass the weights DIRECTORY instead"
        )
    if ss.stat().st_size == 0:
        raise RuntimeError(f"{ss} is empty — no 3Di was written for {db_path.name}")
    if not ss_index.exists():
        raise RuntimeError(f"{ss_index} is missing — the 3Di database is incomplete")

    if db_path.exists() and ss.read_bytes() == db_path.read_bytes():
        raise RuntimeError(
            f"{ss} is byte-identical to {db_path} — it holds amino acid "
            "residues, not 3Di states, so this is not a structure search"
        )

    n_ss = _entry_count(ss_index)
    if index.exists():
        n_seq = _entry_count(index)
        if n_ss != n_seq:
            raise RuntimeError(
                f"{ss_index} has {n_ss} entries but {index} has {n_seq} — "
                "the 3Di database does not cover every sequence"
            )
    return n_ss


# ------------------------------------------------------------------ build --

def run_createdb(fasta: Union[str, Path], db: Union[str, Path],
                 model_dir: Union[str, Path], foldseek_bin: str,
                 threads: int, gpu: bool = False) -> None:
    """foldseek createdb with ProstT5. Module-level so tests can replace it.

    `gpu=True` passes `--gpu 1`, which is the ONLY way onto the CUDA path of
    a -DENABLE_CUDA=1 build (37.6x measured, wiki installs.md 2026-08-22).
    Without it createdb takes the CPU ProstT5ForkRunner, which forks one
    worker per thread and loads the ~3 GB model in EACH: on a 62 GB box with
    --threads 16 the children die and the parent blocks forever in msgsnd on
    a full SysV message queue — 0 s of CPU time, no error, indistinguishable
    from a slow run except by the clock. A CPU-only build rejects --gpu 1
    with "No GPU devices found" (exit 1), which is the loud failure we want.
    """
    cmd = [str(foldseek_bin), "createdb", str(fasta), str(db),
           "--prostt5-model", str(model_dir)]
    if gpu:
        cmd += ["--gpu", "1"]
    cmd += ["--threads", str(threads)]
    logger.info(f"Running foldseek createdb: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"foldseek createdb failed for {fasta} (exit {result.returncode}):\n"
            f"{result.stderr[-2000:]}"
        )


def build_chunk_db(fasta: Union[str, Path], db: Union[str, Path],
                   model_dir: Union[str, Path], foldseek_bin: str = "foldseek",
                   threads: int = 16, gpu: bool = False) -> Path:
    """Build one chunk's ProstT5 database and refuse to return without 3Di."""
    model_dir = Path(model_dir)
    if model_dir.exists() and not model_dir.is_dir():
        raise ValueError(
            f"--prostt5-model must be the weights directory, not a file: "
            f"{model_dir}. Passing the .gguf file makes createdb write amino "
            "acids and no 3Di, without reporting an error"
        )
    db = Path(db)
    run_createdb(fasta, db, model_dir, foldseek_bin, threads, gpu=gpu)
    n = verify_3di_db(db)
    logger.info(f"Chunk database {db.name}: 3Di verified, {n} entries")
    return db


def build_chunk_db_from_config(fasta: Union[str, Path], db: Union[str, Path],
                               model_dir: Union[str, Path], config) -> Path:
    """Config-driven entry point for chunk builds.

    Callers (screen drivers, the planned tier-3 step) go through here so the
    device choice comes from `config.prostt5_gpu` and lands in the run
    manifest, instead of each driver hard-coding its own argv and silently
    defaulting onto the CPU fork path.
    """
    return build_chunk_db(fasta, db, model_dir,
                          foldseek_bin=config.foldseek_bin,
                          threads=config.n_workers,
                          gpu=config.prostt5_gpu)


# ------------------------------------------------------------------ merge --

def merge_search_results(results_dir: Union[str, Path],
                         out_path: Union[str, Path],
                         expected: Optional[int] = None) -> int:
    """Concatenate chunk search results, refusing an incomplete set.

    foldseek writes no terminator of its own, unlike HMMER's `# [ok]`, so a
    task killed by a time limit leaves a file that is syntactically fine and
    silently short. Each chunk therefore gets a `<chunk>.done` sentinel written
    only after its search exits zero, and merging requires it. An empty result
    file with a sentinel is a legitimate "found nothing" and merges normally.
    """
    results_dir = Path(results_dir)
    parts = sorted(results_dir.glob(CHUNK_GLOB))

    if expected is not None and len(parts) != expected:
        raise RuntimeError(
            f"expected {expected} chunk results, found {len(parts)} in "
            f"{results_dir} — refusing to merge a partial screen"
        )
    if not parts:
        raise RuntimeError(f"no chunk results ({CHUNK_GLOB}) in {results_dir}")

    incomplete = [p.name for p in parts
                  if not p.with_name(p.name + DONE_SUFFIX).exists()]
    if incomplete:
        raise RuntimeError(
            f"incomplete chunk searches (no '{DONE_SUFFIX}' sentinel): "
            f"{', '.join(incomplete)} — rerun them before merging"
        )

    with open(out_path, "w") as out:
        for p in parts:
            out.write(p.read_text())
    logger.info(f"Merged {len(parts)} chunk results -> {out_path}")
    return len(parts)
