# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Iterative gene family construction pipeline for comparative genomics. Runs repeated rounds of OrthoFinder clustering → MAFFT alignment → pal2nal codon alignment → FastTree gene tree → species-aware pruning, then a post-convergence HMMER profile rescue step for divergent homologs.

## Running the Pipeline

```bash
# Activate the conda environment (required for OrthoFinder, MAFFT, FastTree, etc.)
conda activate orthofinder
# OR: micromamba activate family_finder

# 5-species run
python family_finder.py \
  --protein-dir data/pep --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output_5sp --config config_5sp.json --threads 8 --verbose

# 11-species run
python family_finder.py \
  --protein-dir data_11sp/pep --cds-dir data_11sp/cds \
  --species-tree data_11sp/species_tree.nwk \
  --outdir output_11sp --config config_11sp.json --threads 16 --verbose
```

There are no automated tests. The `tests/` directory is empty. To verify the pipeline works, run on a small dataset and check `output/summary.tsv`.

## Architecture

The pipeline has a single control flow path:

1. **`family_finder.py`** — CLI entry point. Parses args, loads `Config`, calls `pipeline.run()`.
2. **`pipeline.py`** — Orchestrator. Runs the iterative loop and HMMER rescue. The core worker function `process_single_orthogroup()` is defined here (not in a step module) because it sequences align → tree → prune for a single OG and is called via `parallel_map`.
3. **`steps/`** — Each file wraps one external tool via `subprocess.run()`:
   - `orthofinder.py` — runs OrthoFinder, parses `Orthogroups.tsv`
   - `align.py` — MAFFT protein alignment + pal2nal codon alignment (with internal stop codon filter)
   - `tree.py` — FastTree or IQ-TREE wrapper (auto-detects nucleotide vs protein)
   - `prune.py` — two-stage pruning: TreeShrink (if available) then species-aware median-of-ratios
   - `hmmer_rescue.py` — parallel hmmbuild + hmmsearch + result parsing + re-alignment
   - `codeml.py` — PAML/codeml wrapper (optional)
4. **`utils/`** — Pure Python helpers: FASTA I/O (`seqio.py`), species tree distances (`species.py`), `ProcessPoolExecutor` wrapper (`parallel.py`), checkpoint/resume (`checkpoint.py`).
5. **`config.py`** — Single `@dataclass` with all parameters. Loads from JSON with unknown-key warnings.

### Data flow between rounds

```
Round N: current_pool (protein seqs) → split_by_species → OrthoFinder
  → per-OG: extract seqs → align → tree → prune → (confirmed, outliers)
  → all outliers become current_pool for Round N+1
After convergence: unplaced genes → HMMER profiles from families → rescue
```

### Key design choices to know about

- **Gene IDs must be `SpeciesPrefix_GeneID`**. The species is extracted by splitting on the first `_`. This convention is assumed everywhere (seqio, species utils, pruning).
- **Per-OG sequence subsetting**: Workers receive only their OG's sequences, not the full pool. This was a critical performance fix (143K seqs × pickle per worker was the bottleneck).
- **Protein tree fallback**: When pal2nal fails (internal stop codons), the pipeline builds a protein-only tree instead of crashing.
- **HMMER rescue runs once post-convergence**, not per-round. DIAMOND repeatedly fails on the same divergent genes, so re-searching every round adds cost without benefit.
- **`MAFFT_BINARIES` env var** is popped at module level in both `family_finder.py` and `hmmer_rescue.py` to avoid conda conflicts.
- **OrthoFinder conda env PATH** is prepended in `family_finder.py` for `mcl`/`diamond` discovery.

## HPC Environment

This runs on a CentOS 7 HPC cluster at `/data/gpfs/assoc/pgl/bin/family_finder/`. Tool binaries are in conda envs — paths are set in `config_5sp.json` / `config_11sp.json`. The hardcoded PATH in `family_finder.py` points to `/data/gpfs/assoc/pgl/bin/conda/conda_envs/orthofinder/bin`.

## Common Pitfalls

- **OrthoFinder + DIAMOND version mismatch**: OrthoFinder 3.1.3 passes `--ignore-warnings` which DIAMOND 2.1.24 rejects. Must patch OrthoFinder's `config.json`.
- **TreeShrink requires Python ≤3.9** but OrthoFinder needs 3.11+. Pipeline gracefully skips TreeShrink when unavailable.
- **Species tree must use substitution-rate branch lengths** (IQ-TREE), not coalescent units (ASTRAL). The pruning algorithm divides observed by expected distances.
- **Pep/CDS ID mismatch** will silently reduce orthogroup sizes. Always validate IDs match before running.
