# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Iterative gene family construction pipeline for comparative genomics (Caryophyllales:
cacti + *Mesembryanthemum crystallinum* outgroup). Repeated rounds of OrthoFinder
clustering → MAFFT alignment → pal2nal codon alignment → FastTree gene tree →
species-aware pruning, with an optional per-round HMM profile assignment step and a
post-convergence HMMER rescue. Standalone analysis modules build on the confirmed
families: pseudogene detection and DeepLoc-retargeting + branch-site dN/dS
neofunctionalization detection.

## Commands

```bash
# Run the pipeline (needs the conda env for OrthoFinder, MAFFT, FastTree, mcl)
conda activate orthofinder   # or: micromamba activate family_finder
python family_finder.py \
  --protein-dir data/pep --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output_5sp --config config_5sp.json --threads 8 --verbose

# Tests — run anywhere, NO conda env needed (ete4 is stubbed via sys.modules,
# external tools are mocked; only pytest + Biopython required)
python3 -m pytest tests/ -q
# Single test
python3 -m pytest tests/test_size_gate.py::test_emits_small_family_when_two_survive_pruning -v

# Standalone analyses (post-run)
python find_pseudogenes.py --help
python find_neofunctionalization.py --help   # stage 1 (--events-only) runs anywhere; stage 2 needs PAML
```

## Architecture

Single control-flow path: `family_finder.py` (CLI, env setup) → `pipeline.run()`.

**`pipeline.py`** owns the iterative loop. `process_single_orthogroup()` (align → tree
→ prune for one OG, called via `utils/parallel.parallel_map`) lives here, not in a step
module. Each round: split pool by species → OrthoFinder → per-OG processing → collect
confirmed families and outliers → gene-conservation audit (leaked genes are forced back
into the pool) → optional profile assignment → outliers become next round's pool.
Resume reads per-round `round_NN/confirmed_families.tsv` and the newest **completed**
checkpoint — `summary.tsv` is only written at the very end. Because those per-round
outputs carry no record of the configuration that produced them, `utils/manifest.py`
writes a run manifest at startup and **refuses to resume when the config hash
differs**, naming the settings that changed; the hash covers every `Config` field
except resource knobs and tool paths, plus the species tree's contents and the input
FASTA checksums. `--allow-config-change` overrides it and records that it was used.
An output directory with no manifest (anything produced before this landed) warns and
proceeds, flagged `unverified_resume`.

**Assignment tiers** (a gene can be placed by any of):
1. OrthoFinder DIAMOND+MCL clustering, every round (`steps/orthofinder.py`)
2. Per-round HMM profile assignment (`steps/profile_assign.py`), opt-in via
   `profile_assign_per_round: true` — offers pruned/dissolved genes back to existing
   family profiles BEFORE they re-enter clustering. Decisions use bit scores with
   coverage + best-vs-second margin; **never compare E-values** (HMMER prints ties —
   underflow to 0 and 2-sig-fig rounding — and lexicographic order then decides, with
   `R10_*` sorting before `R1_*`)
3. Legacy post-convergence rescue (`steps/hmmer_rescue.py`), `hmmer_rescue: true`
4. Planned: structural tier via ESM/Foldseek (issue #17)

**Two size gates, deliberately different**: `min_orthogroup_size` (4) is the floor to
*start* align/tree work; `min_family_size` (2) is the floor to *emit* a family after
pruning. Never dissolve a whole OG because pruning shrank it — that discards genes that
passed pruning (the historical bug behind outgroup-only families).

**Pruning** (`steps/prune.py`): `prune_criterion: "relative"` (default — terminal
branches subtracted, per-family median normalization, dual threshold
`prune_relative_threshold`/`prune_score_floor`) or `"absolute"` (legacy
`distance_ratio_threshold`). The math lives in pure ete4-free helper functions so it is
unit-testable. `utils/species.validate_species_tree` runs at pipeline start and raises
when tree leaves share zero names with data prefixes (which would silently disable
pruning: every gene scores 0.0).

**`utils/newick.py`** is a self-contained Newick parser + Fitch parsimony + `#1`
clade-marking used by the retargeting/codeml modules and tests — it exists so that
analysis code runs without ete4. `steps/prune.py` still uses ete4 (`from ete4 import
Tree` at module level — tests must stub it, see below).

**Selection-test policy (2026-08-24, HyPhy-first)**: hypothesis tests run on
HyPhy (RELAX / MEME / aBSREL / BUSTED) — codeml full-clan runs proved
impractical (102 taxa: 100 iterations in 23h, non-convergent; 95 taxa: ~2 days)
and the hypothesis verdicts never depended on them. codeml is reserved for two
uses only: regenerating legacy published numbers, and BEB site posteriors when
a result must be compared against BEB-based literature — and then on a
fast-track taxon subset, never the full clan.

The policy has since been tested by the case it anticipated. On the PEPC
foreground, codeml branch-site returned **significant positive selection**
(LRT 27.11, p 9.6e-08, five BEB sites at P>0.95) while MEME found **zero**
sites at FDR, aBSREL left the same stem non-significant, and RELAX reported
**strong relaxation** (p 5.9e-13). The run was valid — `#1` on the right
internal node, its ten leaves exactly the intended foreground — so this is
not a marking error. Model A's ω₂ ≥ 1 class absorbs sites whose constraint
has lifted toward ω = 1, which is what RELAX exists to distinguish, so the
branch-site significance is best read as relaxation misread rather than as a
second finding. **When the two disagree, HyPhy decides**, and codeml's BEB
sites are not cited without stating that MEME confirmed none of them.

**Neofunctionalization** (`find_neofunctionalization.py`, issue #16):
`steps/deeploc.py` (DeepLoc once per proteome, probability vectors cached — dual
targeting stays ambiguous) → `steps/retargeting.py` (Fitch on family trees, switching
branches; truncated genes excluded — they fake retargeting) → `steps/codeml.py`
branch-site Model A vs null, LRT via `erfc` (no scipy), BH-FDR. DeepLoc is an
annotation layer, never a clustering criterion: families are homology units and
localization changes within them.

### Data flow between rounds

```
Round N: current_pool (protein seqs) → split_by_species → OrthoFinder
  → per-OG: extract seqs → align → tree → prune → (confirmed, outliers)
  → conservation audit → [optional profile assignment against existing families]
  → outliers become current_pool for Round N+1
After convergence: unplaced genes → HMMER profiles from families → rescue
Final output includes unplaced_proteins.fa / unplaced_cds.fa and id_agreement.tsv
```

### Key design choices to know about

- **Gene IDs must be `SpeciesPrefix_GeneID`** — species is the prefix before the first
  `_`. Assumed everywhere (seqio, species utils, pruning, all analysis modules).
- **`codon_align` returns `(path_or_None, removed_ids)`** — the removed set (internal
  stop codons) must be recycled by callers or the genes vanish from every output.
- **Per-OG sequence subsetting**: workers receive only their OG's sequences (pickling
  the 143K-seq pool per worker was the original bottleneck).
- **Worker logging**: `parallel_map` relays worker log records to the parent via
  QueueHandler/QueueListener — without it, ProcessPoolExecutor children log into the void.
- **Protein tree fallback** when pal2nal fails; the same pruning threshold then applies
  to protein distances.
- **`MAFFT_BINARIES` env var** is popped at module level in `family_finder.py` and
  `steps/hmmer_rescue.py` to avoid conda conflicts; `family_finder.py` also prepends the
  orthofinder conda env to PATH for `mcl`/`diamond` discovery.

## Test conventions

`tests/` uses pytest. ete4 is not installed locally: stub it before importing
pipeline/prune (`sys.modules["ete4"] = types.ModuleType(...)` — copy the pattern from
`tests/test_size_gate.py`). External tools (MAFFT, FastTree, hmmer, codeml, DeepLoc)
are never invoked in tests — wrappers are module-level functions monkeypatched with
fakes, and parsers are fed synthetic text fixtures. Pure-math helpers (pruning scores,
coverage merging, Fitch) take plain dicts precisely so they test without dependencies.

## Investigation record

**Read `resume.md` before re-diagnosing any clustering problem — start at §1.5.**
The early verdicts stand (pep/CDS integrity ✗, DIAMOND sensitivity ✗, pruning ✗ for
Mcry, taxon sampling ✗), but the blanket conclusion this file used to carry — "the
confirmed cause is MCL/graph fragmentation" — was overturned on 2026-08-22 for the
RESIDUAL unplaced pool: classified exhaustively, Mcry's unplaced genes are mostly
model-quality artifacts and true orphans (graph cuts 15.8%, the LOWEST of the four
species), and among comparable models Mcry's orphan rate equals Cgig's. MCL
fragmentation remains the valid diagnosis for SPLIT FAMILIES (the PEPC clan in six
fragments) — a different set than the unplaced pool; resume.md "서로 다른 집합을
잰 것이다" is the reconciliation. Diagnostic scripts live on the cluster
(`forensics_r8.py`, `score_recluster.py`, `classify_unplaced.py`, `task*.py` under
`~/scratch/bin/family_finder/`).

## HPC Environment

Runs on the Pronghorn cluster (`ssh pronghorn`); run directory
`~/scratch/bin/family_finder/` (= `/data/gpfs/home/wyim/scratch/bin/family_finder`),
5-species artifacts in `output_5sp/`. Tool paths are set in `config_5sp.json`
(conda envs under `/data/gpfs/.../conda_envs/orthofinder`; a separate `iqtree2` env
exists). SLURM: partition `cpu-s2-core-0`, `--account=cpu-s2-pgl-0`.
`config_11sp.json` and the data dirs exist only on the cluster, not in this repo.

## Common Pitfalls

- **OrthoFinder 3.1.3 + DIAMOND 2.1.24**: OrthoFinder passes `--ignore-warnings` which
  this DIAMOND rejects — patch OrthoFinder's `config.json`.
- **`-og` is accepted by OrthoFinder 3.1.3 but hidden from `-h`** — never probe support
  by parsing help output. `-og` also suppresses `Comparative_Genomics_Statistics/`
  (per-species stats); `Orthogroups*.tsv` survive.
- **Standalone `orthofinder -b` runs need the conda env on PATH** (`mcl` discovery) —
  only pipeline runs get it injected by `family_finder.py`. sbatch scripts must
  `export PATH=$HOME/scratch/bin/conda/conda_envs/orthofinder/bin:$PATH`.
- **Hardlink copies of WorkingDirectory** (`cp -al`, for `-b` sweeps): unlink any
  pre-existing `OrthoFinder_graph.txt` / `clusters_*` before re-running, or OrthoFinder
  may reuse a stale graph or open hardlinked originals for write.
- **Species tree must use substitution-rate branch lengths** (IQ-TREE nucleotide GTR+G
  to match FastTree `-nt -gtr -gamma`), NOT ASTRAL coalescent units.
  `validate_species_tree` warns on suspicious scales. The committed
  `data/species_tree.nwk` values were hand-written; a data-estimated replacement is
  tracked in issue #14.
- **Cgig and CgigH are the same *Carnegiea gigantea* genome annotated twice** (MAKER vs
  Helixer). Their species-tree distance (0.02) makes the legacy absolute pruning
  criterion eject ~10% of Cgig leaves — measured, see resume.md.
- **TreeShrink requires Python ≤3.9** (OrthoFinder needs 3.11+); the pipeline warns and
  skips when unavailable.
- **OrthoFinder v3 default MCL inflation is 1.2** (not v2's 1.5) — the config default
  `orthofinder_inflation: 1.2` reproduces it explicitly.
