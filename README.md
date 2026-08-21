# family_finder

Iterative gene family construction pipeline with species-aware pruning.

## Overview

`family_finder` builds gene families by repeatedly running OrthoFinder, aligning sequences, building trees, and pruning outliers. Outlier sequences from each round are re-clustered in the next round, allowing displaced genes to find their true families.

```
Round 1: All seqs → OrthoFinder → per-OG align/tree/prune → confirmed families + outliers
           └─ optional tier-2: offer outliers to existing family HMM profiles first
Round 2: Outliers → OrthoFinder → ... → new families + outliers
Round N: Repeat until convergence
                                          ↓
HMMER Rescue: Unplaced genes → hmmsearch vs family HMM profiles → rescued into families
                                          ↓
Pseudogene Detection (optional): 6 evidence types → high/medium/low candidates
                                          ↓
Post-run analysis (standalone CLIs, see "Annotation & Naming Stack"):
  subfamilies · EC · structure transfer · localization · naming · selection
```

A gene can be placed by any of four tiers: OrthoFinder clustering (every
round), per-round HMM profile assignment (`profile_assign_per_round`),
post-convergence HMMER rescue, and — for ambiguous cases — EPA-ng
phylogenetic placement into an IQ-TREE reference. Everything after that
is annotation, never a membership filter.

### Why iterative?

A single OrthoFinder run can mis-cluster sequences — especially fast-evolving genes or genes from incomplete annotations. By pruning outliers and re-running clustering, these sequences get a second chance to group with their true orthologs.

## Pipeline Steps (per orthogroup)

1. **Protein alignment** — MAFFT
2. **Codon alignment** — pal2nal (protein-guided CDS alignment)
3. **Gene tree** — FastTree or IQ-TREE on CDS alignment (`-nt -gtr -gamma`)
4. **Pruning** — Two-stage outlier detection:
   - **Stage 1: TreeShrink** — statistical branch-length outlier removal
   - **Stage 2: Species-aware ratio** — observed/expected distance ratio using species tree
5. **Re-alignment** — clean alignment and tree from confirmed members only

If codon alignment fails (e.g., CDS annotation errors with internal stop codons), the pipeline falls back to protein-based tree building.

## Species-Aware Pruning Algorithm

For each gene in the tree:

```
outlier_score(gene_i) = median over all gene_j of:
    observed_distance(gene_i, gene_j) / expected_distance(species_i, species_j)
```

- `observed_distance`: pairwise distance in the gene tree
- `expected_distance`: pairwise distance in the species tree
- Same-species comparisons (paralogs) are skipped
- `outlier_score > threshold` (default 5.0) → gene is removed

This approach is rooting-independent, robust to outliers (median), and normalizes for phylogenetic distance so distant species don't trigger false positives.

## Convergence

The pipeline stops when any of:
- `max_rounds` reached (default: 10)
- No new families for `convergence_no_new_families` consecutive rounds (default: 2)
- Outlier pool drops below `convergence_threshold` (default: 5)

## HMMER Profile-Based Rescue

After iterative clustering converges, genes that remain unplaced may still be true members of confirmed families but too divergent for DIAMOND pairwise similarity to detect. The HMMER rescue step addresses this by leveraging profile-based search.

**How it works:**

1. Build HMM profiles from the protein alignments of all confirmed families (using `hmmbuild`)
2. Search unplaced gene sequences against the profile database (using `hmmsearch`)
3. Assign each gene to the best-matching family if E-value < threshold (default: 1e-5)
4. Re-align and rebuild trees for affected families

**Why HMMER succeeds where DIAMOND fails:** DIAMOND compares individual sequence pairs, so a divergent gene may not match any single family member well enough. HMM profiles capture position-specific conservation patterns across the entire family, detecting conserved domain architecture even in highly divergent sequences.

### Results (5-species run)

- After convergence (round 10): 13,963 genes remained unplaced
- HMMER rescue recovered **6,358 genes (45.5%)** into existing families
- Breakdown by species:

| Species | Genes Rescued |
|---|---|
| *O. cochenillifera* (Ococ) | 1,692 |
| *C. gigantea* Helixer (CgigH) | 1,573 |
| *O. basilaris* (Obas) | 1,262 |
| *C. gigantea* (Cgig) | 1,099 |
| *M. crystallinum* (Mcry) | 732 |

### CAM pathway gene rescue

Of 9 initially unplaced CAM pathway genes from *M. crystallinum*, 8 were rescued by HMMER:

| Gene ID | Gene Name | Assigned Family | E-value |
|---|---|---|---|
| Mcr2G22880 | CKB1 | R1_OG0000560 | 1.7e-147 |
| Mcr4G10250 | Pfk5 | R1_OG0004219 | 2e-158 |
| Mcr1G03910 | Nst | R1_OG0015572 | 5.2e-76 |
| Mcr1G01450 | VhaE | R1_OG0000561 | 8.5e-73 |
| Mcr5G24020 | VhaB | R1_OG0001531 | 1.4e-59 |
| Mcr1G10640 | Cbl3 | R3_OG0000142 | 3.7e-34 |
| Mcr9G21470 | Lda1 | R1_OG0012544 | 1.2e-9 |
| Mcr5G21920 | FAR1 | R1_OG0000605 | 2.6e-8 |

The single unrescued gene, Mcr4G19500 (Gln1/Gln2, 78 aa), is a truncated gene model too short for reliable profile matching.

## Pseudogene Detection (Optional)

After iterative clustering and HMMER rescue, the pipeline can optionally identify pseudogene candidates across the genome. Pseudogenes -- genes that have lost their protein-coding function -- are common in plant genomes but are often mis-annotated as functional genes. Detecting them improves downstream analyses (e.g., gene family size comparisons, selection tests) by separating functional genes from non-functional copies.

Pseudogene detection is **enabled by default** but can be disabled with `--no-pseudogene-detection` or by setting `"pseudogene_detection": false` in the config. It can also be run standalone via `find_pseudogenes.py` on a completed pipeline output.

### Evidence types

The detector examines six independent lines of evidence. Each gene accumulates an evidence vector, and only genes with at least one positive signal are reported as candidates.

| Evidence | What it detects | Threshold |
|---|---|---|
| Internal stop codons | Premature termination codons in the protein sequence — direct loss-of-function evidence | Any internal `*` |
| CDS/protein length discrepancy | CDS length differs from protein×3 by >10% — frameshifts or assembly errors | \|ratio - 1.0\| > 0.1 |
| Truncated gene | Gene shorter than 50% of its family's median length — partial gene model | < 50% of family median |
| Orphan gene | Unplaced after all clustering rounds + HMMER rescue — likely pseudogene or foreign (HGT) | Not in any family |
| GC3 composition outlier | 3rd-codon-position GC% deviates >3 SD from species mean — compositional drift from pseudogenization or HGT | \|z-score\| > 3.0 |
| Long branch length | Branch length in gene tree >3x median — accelerated evolution from relaxed selection | distance ratio > 3.0 |

### Confidence scoring

Each evidence type carries a weight reflecting its biological informativeness:

| Evidence | Weight | Rationale |
|---|---|---|
| Internal stop codons | 0.40 | Direct evidence of loss of function |
| CDS/protein length mismatch | 0.35 | Frameshift or assembly error |
| Truncated gene | 0.30 | Partial gene model |
| Orphan gene | 0.25 | Unplaced genes are usually pseudogenes or foreign (HGT) |
| GC3 composition outlier | 0.25 | Compositional drift indicates pseudogenization or HGT |
| Long branch length | 0.20 | Accelerated evolution (could also be positive selection) |

The **confidence score** is the sum of applicable weights, capped at 1.0. Classification uses score-based thresholds:

| Classification | Criteria | Interpretation |
|---|---|---|
| `pseudogene_high` | score ≥ 0.50 | Multiple evidence lines — highly likely pseudogene |
| `pseudogene_medium` | score ≥ 0.25 | Single meaningful evidence (orphan, truncated, stop codon, GC3 outlier) |
| `pseudogene_low` | score > 0 | Only long branch (0.20) — weak signal, may be under positive selection |
| `functional` | score = 0 | No pseudogene evidence |

### Usage

**Integrated mode (optional, enabled by default):**

```bash
# Standard run — pseudogene detection runs automatically
python family_finder.py \
  --protein-dir data/pep --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output_5sp --threads 8 --verbose

# Restrict pseudogene analysis to one species
python family_finder.py ... --pseudogene-species Ococ

# Disable pseudogene detection entirely
python family_finder.py ... --no-pseudogene-detection
```

**Standalone mode:** Run on an already-completed pipeline output directory (must contain `summary.tsv`):

```bash
# Single species
python find_pseudogenes.py \
  --protein-dir data/pep --cds-dir data/cds \
  --outdir output_5sp --species Ococ

# All species with custom truncation threshold
python find_pseudogenes.py \
  --protein-dir data/pep --cds-dir data/cds \
  --outdir output_5sp --truncation-threshold 0.4
```

### Output files

All pseudogene output is written to `<outdir>/pseudogene_analysis/`:

| File | Description |
|---|---|
| `pseudogene_candidates.tsv` | Full candidate list with all evidence columns and confidence score |
| `pseudogene_summary.txt` | Human-readable statistics: classification breakdown, evidence type counts |
| `pseudogene_candidates.pep.fa` | Protein sequences of all candidates |
| `pseudogene_candidates.cds.fa` | CDS sequences of all candidates |
| `pseudogene_candidates.bed` | BED file for genome browser (red=high, orange=medium, yellow=low) |
| `family_pseudogene_enrichment.tsv` | Per-family pseudogene concentration |
| `chromosomal_distribution.tsv` | Per-chromosome pseudogene density |
| `species_comparison.tsv` | Cross-species pseudogene rates (all-species mode only) |

### GFF3 filtering

Use the pseudogene results to create a clean GFF3 without pseudogenes:

```python
# Example: filter pseudogenes from GFF3
pseudo_ids = set()
with open("output_5sp/pseudogene_analysis/pseudogene_candidates_Ococ.tsv") as f:
    f.readline()  # skip header
    for line in f:
        parts = line.strip().split("\t")
        if parts[2].startswith("pseudogene"):  # all confidence levels
            pseudo_ids.add(parts[0].split("_", 1)[1])  # strip species prefix

# Then filter GFF3 lines where gene ID is in pseudo_ids
```

### Config parameters

| Parameter | Type | Default | CLI flag | Description |
|---|---|---|---|---|
| `pseudogene_detection` | bool | `true` | `--no-pseudogene-detection` | Enable/disable pseudogene detection |
| `pseudogene_truncation_threshold` | float | `0.5` | — | Flag genes shorter than this fraction of family median |
| `pseudogene_species_filter` | string | `""` | `--pseudogene-species` | Restrict to one species (e.g., `"Ococ"`) |

### Results (5-species run, *O. cochenillifera*)

| Classification | Count | % of 33,745 genes |
|---|---|---|
| `pseudogene_high` | 181 | 0.5% |
| `pseudogene_medium` | 4,020 | 11.9% |
| `pseudogene_low` | 142 | 0.4% |
| **Total candidates** | **4,343** | **12.9%** |

Top evidence types: orphan/unplaced (2,639), truncated (1,200), GC3 outlier (365), long branch (166), internal stops (161).

## Installation

All environments here use **micromamba** — a single static binary, no base
environment, and the fastest solver of the three. (`conda` and `mamba` accept
the same subcommands if you already have one; substitute freely.)

### 1. Install micromamba

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# one-time shell setup, so `micromamba activate` works in new shells
micromamba shell init -s bash -r ~/micromamba
exec $SHELL
```

### 2. Create the pipeline environment

```bash
micromamba create -y -f environment.yml
micromamba activate family_finder
```

Or build it explicitly, without `environment.yml`:

```bash
micromamba create -y -n family_finder -c conda-forge -c bioconda \
  python=3.11 orthofinder mafft fasttree iqtree paml pal2nal ete4 biopython rich
micromamba activate family_finder

# Optional: TreeShrink needs Python <=3.9, so give it its own environment
micromamba create -y -n treeshrink -c conda-forge -c bioconda treeshrink
```

**Put `conda-forge` before `bioconda`.** With strict channel priority the
reverse order can pin a years-old build — measured: HyPhy resolved to 2.5.29
(2021) that way, and pinning a newer version made the solve unsatisfiable.

### 3. Verify

```bash
python -c "from ete4 import Tree; from Bio import SeqIO; print('OK')"
orthofinder -h | head -1
mafft --version
FastTree 2>&1 | head -1
pal2nal.pl 2>&1 | head -1
```

Running a command without activating first (handy in scripts and SLURM jobs):

```bash
micromamba run -n family_finder python family_finder.py --help
```

### Annotation stack (optional, per axis)

The post-run annotation CLIs each need their own tool. They are independent —
install only the axes you want; `annotation_matrix.py` merges whatever is
present and reports what is missing. Every gotcha below was hit in practice;
full install logs live in the lab wiki (`guide/installs.md`).

Each axis gets its **own environment** — DeepLoc, SignalP and CLEAN pin
mutually incompatible torch/numpy versions, so sharing one environment means
one tool silently overwrites another's dependencies.

| Axis | Tool | Device | Install |
|---|---|---|---|
| EC (orthology) | eggNOG-mapper 2.1.15 | CPU | `micromamba create -y -n emapper -c conda-forge -c bioconda python=3.11 eggnog-mapper`, then download the DB (~48 GB extracted) |
| EC (embedding) | CLEAN | GPU | `git clone https://github.com/tttianhao/CLEAN`, follow its README; pulls ESM-1b weights (7.8 GB) |
| Structure transfer | Foldseek + AFDB-SwissProt | CPU | static binary from the Foldseek releases page, then `foldseek databases Alphafold/Swiss-Prot afdb_swissprot tmp` (3.9 GB) |
| Structures (input to the above) | ESMFold *or* ProstT5 3Di | GPU | ESMFold via `transformers.EsmForProteinFolding` in a CUDA-torch environment; ProstT5 is the fold-free alternative when no GPU is available |
| Localization | DeepLoc 2.1 | GPU | academic licence from DTU HealthTech (registration-gated); `pip install --no-deps <package>`, then add its dependencies individually. Run with `-d cuda` |
| Signal peptide | SignalP 6.0 fast | CPU | academic licence from DTU HealthTech; `pip install <package>/`, then copy the 1.6 GB model weights into the installed package. The CPU torch build is the right choice here |
| Subfamilies | Possvm, TreeCluster | CPU | `pip install possvm treecluster` (or from source) |
| Selection | HyPhy ≥2.5.60, PAML | CPU | `micromamba create -y -n hyphy -c conda-forge -c bioconda hyphy paml` |

#### CPU or GPU?

**The pipeline itself is CPU-only.** OrthoFinder, MAFFT, FastTree, IQ-TREE,
pal2nal, HMMER, codeml and HyPhy have no GPU path — a CPU cluster runs
everything in the sections above. GPUs only matter for three annotation axes
that run protein language models.

**Reference machine** (what the measurements below were taken on):

| | |
|---|---|
| GPU | NVIDIA RTX 4090, 24 GB VRAM |
| Driver / CUDA | 560.35.03 / CUDA 12.6 |
| Host | 16 cores, 62 GB RAM |

**GPU requirements, by axis:**

| Axis | VRAM | Host RAM | Notes |
|---|---|---|---|
| CLEAN (ESM-1b, 650M) | modest — fits well inside 24 GB | — | not separately profiled; the model is an order of magnitude smaller than the two below |
| DeepLoc `Accurate` (ProtT5-XL, ~3B) | **5.7 GB measured** | **16.3 GB** (model load, not accumulation) | an 8 GB card is enough |
| ESMFold | **the constraint** — 24 GB covers proteins up to ~1000 aa in fp16 with `chunk_size=64` | — | longer proteins need more VRAM or smaller chunks |

So **8 GB of VRAM runs everything except ESMFold**, and ESMFold is itself
optional — Foldseek can take ProstT5 3Di sequences instead of predicted
structures when no large GPU is available (at the cost of losing TM-score,
since ProstT5-built databases carry no CA coordinates).

**Driver vs torch build.** CUDA minor-version compatibility means a driver
new enough for the CUDA *major* version runs any build of it: this machine's
CUDA 12.6 driver runs cu121, cu124 and cu126 wheels side by side (verified —
three environments here hold torch 1.13.1+cpu, 2.5.1/cu124 and 2.13.0/cu126).
Match the major version to your driver and do not chase the minor one.

**Disk is the requirement people forget.** Model weights and databases
dominate — measured on the reference machine:

| Item | Size |
|---|---|
| eggNOG DB (extracted) | 48 GB |
| Foldseek AFDB-SwissProt | 3.9 GB |
| ESM weights (torch hub cache) | 7.4 GB |
| SignalP 6 model weights | 1.6 GB |
| CLEAN checkout + checkpoints | 0.9 GB |
| ProtT5-XL (HuggingFace cache, DeepLoc `Accurate`) | ~11 GB |

Budget **~75 GB** for the full stack, and put the caches somewhere with room
(`HF_HOME`, `TORCH_HOME`) rather than letting them land in a small `$HOME`.

| Axis | Device | Measured |
|---|---|---|
| eggNOG-mapper | **CPU only** (DIAMOND) | — |
| SignalP 6 (fast model) | **CPU is fine** | 0.23 s/seq; 1.6 s/seq on an older cluster node |
| Foldseek search | **CPU** | seconds for ~100 structures against AFDB-SwissProt |
| CLEAN | **GPU strongly preferred** (ESM-1b) | 102 seqs ≈ 2 min on an RTX 4090 |
| DeepLoc 2.x `Accurate` | **GPU strongly preferred** (ProtT5-XL, ~3B params) | 29,405 proteins in 21 min on a 4090; 5.7 GB VRAM, 16.3 GB RSS |
| ESMFold | **GPU effectively required** | ~75 s per ~930 aa protein on a 4090; 24 GB VRAM covers ~1000 aa |

Install the matching PyTorch build — this is the one choice that decides
whether an axis can use the GPU at all:

```bash
# CPU-only host (smaller download, no driver needed)
pip install "torch==1.13.1" --index-url https://download.pytorch.org/whl/cpu

# CUDA host — pick the cuXXX that matches your driver (nvidia-smi shows it)
pip install "torch==2.1.2" --index-url https://download.pytorch.org/whl/cu121
```

`numpy<2` applies to **both** builds whenever torch is <2.x (see the traps
below). Confirm what you actually got before running anything long:

```bash
python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
```

**⚠️ Silent CPU fallback is the recurring failure here — a run that is 50×
slower looks exactly like a run that is working.** Three measured cases:

- **DeepLoc 2.0 never touches the GPU as shipped.** There is no `.cuda()` or `.to(device)` anywhere in its model code, so the 3B-parameter encoder runs on CPU at batch size 1. Moving just the frozen encoder to CUDA fixes it and produces **byte-identical** output (verified against a `CUDA_VISIBLE_DEVICES=` run) — the patch changes speed only. DeepLoc **2.1** accepts `-d cuda` directly and needs no patch.
- **A stray copy of the package in the working directory wins over the installed one.** After unpacking the tarball, `import DeepLoc2` resolves to `./DeepLoc2` and silently reverts to the unpatched CPU code. Delete the extracted directory and check with `python -c "import DeepLoc2; print(DeepLoc2.__file__)"`.
- **Foldseek's ProstT5 GPU flag fell back to CPU without an error** (measured: 95 sequences took 13 min on CPU). Read the log, do not assume.

While a long job runs, `nvidia-smi` should show real utilisation — for the
DeepLoc/ESMFold axes we measured 88–100%. Near-zero means you are on the CPU
path regardless of what you passed.

**Channel order matters.** On a machine with strict channel priority,
`-c bioconda -c conda-forge` can pin an ancient build — HyPhy resolved to a
2021 release (2.5.29) that way. Put **conda-forge first**: `-c conda-forge -c bioconda`.

**Known install traps** (all measured, not hypothetical):

- **SignalP 6 / DeepLoc: pin `numpy<2`.** Their `torch<2` requirement is ABI-incompatible with NumPy 2.x, which pip installs by default. Symptom: `RuntimeError: Numpy is not available` at prediction time, not install time.
- **DeepLoc on Python 3.12: pin `setuptools<81`** (`pkg_resources` removed) and `pip install sentencepiece` (`T5Tokenizer requires the SentencePiece library`). Install with `--no-deps` or it will overwrite the env's torch.
- **eggNOG-mapper: the hardcoded download host is dead.** `download_eggnog_data.py` points at `eggnogdb.embl.de`, which no longer resolves ([#575](https://github.com/eggnogdb/eggnog-mapper/issues/575)). Fetch the same files from the `eggnog5.embl.de` mirror manually. Also pin **Python 3.11** — 3.12+ fails on the removed `distutils` ([#516](https://github.com/eggnogdb/eggnog-mapper/issues/516)).
- **Old glibc (RHEL7 / CentOS 7) breaks pip.** With glibc 2.17 there is no matching manylinux wheel, so pip builds from source and the system gcc (gnu89) rejects C99 code — `pillow` fails with `'for' loop initial declarations are only allowed in C99 mode`. Install C-extension packages (pillow, matplotlib, numpy, scipy) **from conda-forge** and add only pure-Python packages with `pip install --no-deps`.
- **Old libstdc++ shadows the env's.** Even then, conda-forge matplotlib may load the system library and die with `GLIBCXX_3.4.20 not found`. `micromamba install libstdcxx-ng` alone does not fix it — wrap the entry point with `LD_LIBRARY_PATH=$ENV/lib`.
- **ProstT5-built Foldseek databases have no CA coordinates.** Asking for TM-score output (`qtmscore`) against one dies with `Cannot open db_ca`. Use bits/E-value with `--alignment-type 2` there; TM-score works against real structure databases such as AFDB.

### Verify dependencies

```bash
# All of these must succeed:
orthofinder -h          # OrthoFinder (includes diamond, mcl, famsa)
mafft --version         # MAFFT aligner
FastTree 2>&1 | head    # FastTree (note: capital F and T)
pal2nal.pl 2>&1 | head  # pal2nal codon aligner
python -c "from ete4 import Tree"    # ete4 tree library
python -c "from Bio import SeqIO"    # Biopython

# Optional:
iqtree --version                     # IQ-TREE (alternative tree builder)
codeml 2>&1 | head                   # PAML/codeml (selection analysis)
run_treeshrink.py -h                 # TreeShrink (outlier detection)
```

### Troubleshooting

| Problem | Solution |
|---|---|
| `diamond: Invalid option: ignore-warnings` | OrthoFinder 3.1.3 + diamond <2.1.25 incompatibility. Remove `--ignore-warnings` from `orthofinder/run/config.json` |
| `Cannot run MCL` | Ensure `mcl` is in PATH: `which mcl` or install via `micromamba install -c conda-forge -c bioconda mcl` |
| `MAFFT_BINARIES` conflict | Pipeline auto-clears this env var. If issues persist: `unset MAFFT_BINARIES` |
| `No module named 'Bio'` | Run with the correct Python: the one inside your micromamba environment (`micromamba run -n family_finder python ...` avoids the mistake) |
| TreeShrink won't install | Requires Python <=3.9. Use a separate env or skip (Stage 2 pruning still works) |
| `ModuleNotFoundError: rich` | `pip install rich` in the OrthoFinder environment |
| `RuntimeError: Numpy is not available` (SignalP/DeepLoc) | `torch<2` is ABI-incompatible with NumPy 2.x: `pip install "numpy<2"` |
| `GLIBCXX_3.4.20 not found` | System libstdc++ is older than the environment's. Wrap the entry point with `LD_LIBRARY_PATH=$ENV/lib` |
| pip builds `pillow` from source and fails on C99 | glibc too old for manylinux wheels. Install C-extension packages from conda-forge, add pure-Python ones with `pip install --no-deps` |
| eggNOG-mapper DB download hangs/fails | The hardcoded host is dead — download from the `eggnog5.embl.de` mirror manually |
| A codeml job shows `FAILED` but wrote results | codeml exits 1 on `error: end of tree file` **after** writing everything. Judge by the output: an `lnL` line and (if enabled) a complete `Bayes Empirical Bayes` block mean the run is usable. A missing `Time used:` line is only the trace of that exit |
| BEB probabilities disagree with what you parsed | codeml prints a **NEB block before the BEB block**, with different values for the same site (measured: 0.931 vs 0.970). Parse only after the `Bayes Empirical Bayes` header — codeml itself says "please use the BEB results" |

## Quick Start

```bash
# 1. Install
micromamba create -y -f environment.yml && micromamba activate family_finder

# 2. Prepare input (see Input Preparation below)

# 3. Run
python family_finder.py \
  --protein-dir data/pep \
  --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output \
  --threads 8 \
  --verbose

# 4. Check results
cat output/summary.tsv | head
ls output/final_families/
```

## Input Preparation

### 1. Protein and CDS FASTA files

Each species needs one protein file and one CDS file in separate directories:

```
data/pep/
  Mcry.pep.fa      # Mesembryanthemum crystallinum protein sequences
  Ococ.pep.fa      # Opuntia cochenillifera proteins
  Cgig.pep.fa      # Carnegiea gigantea proteins

data/cds/
  Mcry.cds.fa      # Mesembryanthemum crystallinum CDS sequences
  Ococ.cds.fa      # Opuntia cochenillifera CDS
  Cgig.cds.fa      # Carnegiea gigantea CDS
```

### 2. Gene ID format

Gene IDs **must** follow `SpeciesPrefix_GeneID` format. The species is extracted from the prefix before the first `_`:

```
>Mcry_Mcr1G24690          → species "Mcry"
>CgigH_Cgig_v2_SGP5p_1.1 → species "CgigH"
```

If your original gene IDs don't have species prefixes, add them:

```bash
# Add species prefix to all gene IDs in a FASTA file
sed 's/^>/>Mcry_/' original.pep.fa > data/pep/Mcry.pep.fa
sed 's/^>/>Mcry_/' original.cds.fa > data/cds/Mcry.cds.fa
```

**Important**: Protein and CDS files must use **identical gene IDs**. The pipeline matches them by ID.

### 3. Extracting CDS/protein from GFF3 + genome

If you have a GFF3 annotation and genome FASTA:

```bash
# Using gffread (recommended)
gffread annotation.gff3 -g genome.fa -x cds.fa -y pep.fa

# Or using AGAT (filters longest isoform first)
agat_sp_keep_longest_isoform.pl -gff annotation.gff3 -o longest.gff3
agat_sp_extract_sequences.pl -g longest.gff3 -f genome.fa -t cds -o cds.fa
agat_sp_extract_sequences.pl -g longest.gff3 -f genome.fa -p -o pep.fa
```

### 4. Species tree

Provide a Newick-format species tree with branch lengths. Leaf names must match species prefixes in gene IDs:

```
((Mcry:0.5,Cgig:0.5):0.3,(Ococ:0.2,Obas:0.2):0.3);
```

Recommended approach: estimate branch lengths from a concatenated single-copy ortholog (e.g. BUSCO) codon alignment with IQ-TREE. Do NOT use ASTRAL output directly — its coalescent-unit branch lengths break the observed/expected distance ratio used for pruning.

To include two annotations of the same genome (e.g., MAKER vs Helixer), add them as sister taxa with near-zero distance:

```
((Mcry:0.5,(Cgig:0.01,CgigH:0.01):0.49):0.3,(Ococ:0.2,Obas:0.2):0.3);
```

### 5. Validate inputs before running

```bash
# Check gene ID consistency between protein and CDS
for sp in Mcry Ococ Obas Cgig; do
  pep_ids=$(grep "^>" data/pep/${sp}.pep.fa | sort)
  cds_ids=$(grep "^>" data/cds/${sp}.cds.fa | sort)
  diff <(echo "$pep_ids") <(echo "$cds_ids") | head -5
  echo "${sp}: $(grep -c '^>' data/pep/${sp}.pep.fa) pep, $(grep -c '^>' data/cds/${sp}.cds.fa) cds"
done

# Check for internal stop codons (annotation quality)
python -c "
from Bio import SeqIO
for f in ['Mcry','Ococ','Obas','Cgig']:
    bad = sum(1 for r in SeqIO.parse(f'data/pep/{f}.pep.fa','fasta') if '*' in str(r.seq).rstrip('*'))
    total = sum(1 for _ in SeqIO.parse(f'data/pep/{f}.pep.fa','fasta'))
    print(f'{f}: {bad}/{total} genes with internal stop codons')
"
```

## Usage

```bash
python family_finder.py \
  --protein-dir data/pep \
  --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output \
  --threads 8 \
  --verbose
```

### Resume after interruption

```bash
# Pipeline saves checkpoints after each round
# Resume from the last completed round:
python family_finder.py \
  --protein-dir data/pep \
  --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output \
  --resume \
  --threads 8
```

### Required arguments

| Argument | Description |
|---|---|
| `--protein-dir` | Directory of per-species protein FASTA files (e.g., `Mcry.pep.fa`) |
| `--cds-dir` | Directory of per-species CDS FASTA files (e.g., `Mcry.cds.fa`) |
| `--species-tree` | Newick species tree with substitution-rate branch lengths (e.g. IQ-TREE concatenation; NOT ASTRAL coalescent units) |
| `--outdir` | Output directory |

### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `--config` | — | JSON config file (overrides defaults) |
| `--resume` | — | Resume from last checkpoint |
| `--max-rounds` | 10 | Maximum iterative rounds |
| `--threshold` | 5.0 | Distance ratio threshold for pruning |
| `--threads` | 8 | Parallel workers / OrthoFinder threads |
| `--tree-builder` | fasttree | `fasttree` or `iqtree` |
| `--run-codeml` | — | Run PAML/codeml on confirmed families |
| `--no-hmmer-rescue` | — | Disable HMMER rescue step after convergence |
| `--hmmer-evalue` | 1e-5 | E-value threshold for HMMER rescue |
| `--no-pseudogene-detection` | — | Disable post-convergence pseudogene detection |
| `--pseudogene-species` | — | Restrict pseudogene analysis to one species (e.g., `Ococ`) |
| `--verbose` | — | Debug logging |

### JSON config

All parameters can be set via a JSON config file:

```json
{
  "orthofinder_bin": "orthofinder",
  "mafft_bin": "mafft",
  "fasttree_bin": "FastTree",
  "pal2nal_bin": "pal2nal.pl",
  "hmmbuild_bin": "hmmbuild",
  "hmmsearch_bin": "hmmsearch",
  "hmmpress_bin": "hmmpress",
  "hmmer_rescue": true,
  "hmmer_evalue": 1e-5,
  "hmmer_chunk_size": 0,
  "orthofinder_threads": 8,
  "n_workers": 8,
  "max_rounds": 10,
  "min_orthogroup_size": 4,
  "min_family_size": 2,
  "distance_ratio_threshold": 5.0,
  "treeshrink_quantile": 0.05,
  "tree_builder": "fasttree",
  "orthofinder_inflation": 1.2,
  "clustering_species_exclude": [],
  "profile_assign_per_round": false,
  "prune_criterion": "relative"
}
```

### Key v2 parameters

| Parameter | Default | Description |
|---|---|---|
| `orthofinder_inflation` | `1.2` | MCL inflation (`-I`); reproduces the OrthoFinder v3 default explicitly. Sweep {1.05–1.3} measured flat — treat as a dead knob (issue #10/#25) |
| `clustering_species_exclude` | `[]` | Species prefixes kept OUT of tier-1 clustering every round (e.g. `["CgigH"]` — duplicate annotation distorting MCL, +9.0% co-clustering when excluded). Excluded genes rejoin the unplaced pool after convergence for profile mapping / HMMER rescue (issue #12) |
| `profile_assign_per_round` | `false` | Tier-2: offer each round's outliers to existing family HMM profiles before re-clustering (issue #13). Bit-score based — never compares E-values |
| `prune_criterion` | `"relative"` | Terminal-branch-subtracted, family-median-normalized pruning with dual threshold (`prune_relative_threshold`=3.0, `prune_score_floor`=2.0); `"absolute"` = legacy `distance_ratio_threshold` (issue #14) |
| `min_family_size` | `2` | Floor to EMIT a family after pruning — deliberately separate from `min_orthogroup_size` (issue #11) |
| `epa_min_lwr` | `0.2` | EPA-ng adjudicator: min like-weight ratio to accept a placement. Calibrated by PEWO prune-and-place on the PEPC clan (issue #23): thresholds 0.1–0.3 give recall 1.0 / precision(nd≤2) 0.97; the old 0.8 lost ~9% recall for nothing |
| `epa_lwr_margin` | `0.0` | Best-vs-second LWR margin — measured as a dead knob in calibration; exact ties remain ambiguous regardless. Catastrophic misplacements are short fragments (<150 aa), so gate on length, not margin |
| `epa_lwr_aggregate` | `"edge"` | `"family"` (opt-in) accepts on the SUM of LWR over each family's edges — recovers correct placements whose mass splits across adjacent edges of the same family (the dominant failure mode in calibration) |
| `epa_min_query_len` | `150` | Reject fragment queries (ungapped aa) before placement — every catastrophic misplacement in calibration was an 80–129 aa fragment. `0` disables |
| `hmmer_chunk_size` | `0` | `0` = one hmmsearch call (unchanged behaviour). `>0` splits the profile database into chunks of this many profiles and runs them through a generated **SLURM-optional** script (`sbatch --wait --array` when available, bounded local pool otherwise), then merges. Measured need: rescue cost is profiles × sequences, so the 5sp run's 5.8 h extrapolates to 2.7–4.1 days at 15sp — past a 3-day limit, with no checkpoint inside hmmsearch (issue #31). Companion knobs: `hmmer_chunk_concurrent` (10), `hmmer_chunk_sbatch_extra` (partition/account flags) |
| `clustering_method` | `"orthofinder"` | Tier-1 backend; `"sonicparanoid"` adapter exists but was rejected on measured assignment rates (issue #22) |

## Output

```
outdir/
  round_01/
    input/                     # Per-species FASTA for OrthoFinder
    orthofinder/               # OrthoFinder results
    orthogroups/OG0000000/     # Per-OG outputs:
      proteins.afa             #   protein alignment
      codon.afa                #   codon alignment
      tree.nwk                 #   gene tree
      confirmed_proteins.fa    #   confirmed member sequences
      confirmed_cds.fa         #   confirmed CDS sequences
      confirmed_proteins.afa   #   re-aligned (clean)
      confirmed_codon.afa      #   re-aligned codon (clean)
      confirmed_tree.nwk       #   re-built tree (clean)
    outlier_pool.fa            # Sequences for next round
    round_stats.json           # Round statistics
  round_02/ ...
  hmmer_rescue/
    family_profiles/           # HMM profiles per family
    rescue_summary.tsv         # gene_id, best_family, evalue
  pseudogene_analysis/
    pseudogene_candidates.tsv  # All candidates with evidence columns
    pseudogene_summary.txt     # Human-readable statistics
    pseudogene_candidates.pep.fa  # Protein sequences of candidates
    pseudogene_candidates.cds.fa  # CDS sequences of candidates
    pseudogene_candidates.bed  # BED file for genome browser
    family_pseudogene_enrichment.tsv  # Per-family pseudogene rates
    chromosomal_distribution.tsv      # Per-chromosome density
    species_comparison.tsv            # Cross-species rates
  final_families/              # All confirmed families (incl. HMMER-rescued genes)
  summary.tsv                  # family_id, round, n_genes, n_species, gene_list
  pipeline.log
```

### summary.tsv

Tab-separated file with one row per confirmed family:

```
family_id    round    n_genes    n_species    gene_list
R1_OG0002940    1    9    5    CgigH_...,Cgig_...,Mcry_...,Obas_...,Ococ_...
```

## Examples

### Example 1: Standard 5-species run

Build gene families across 5 species with all defaults (HMMER rescue + pseudogene detection enabled):

```bash
python family_finder.py \
  --protein-dir data/pep \
  --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output_5sp \
  --config config_5sp.json \
  --threads 8 \
  --verbose
```

### Example 2: Resume after interruption

```bash
python family_finder.py \
  --protein-dir data/pep --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output_5sp --resume --threads 8
```

### Example 3: Pseudogene detection only (standalone)

Run on a completed pipeline output without re-running the full pipeline:

```bash
# All species
python find_pseudogenes.py \
  --protein-dir data/pep --cds-dir data/cds \
  --outdir output_5sp

# Single species (Opuntia cochenillifera)
python find_pseudogenes.py \
  --protein-dir data/pep --cds-dir data/cds \
  --outdir output_5sp --species Ococ
```

### Example 4: Filter pseudogenes from GFF3

Create a clean GFF3 by removing pseudogene candidates detected by the pipeline:

```python
# filter_pseudogenes.py
pseudo_ids = set()
with open("output_5sp/pseudogene_analysis/pseudogene_candidates_Ococ.tsv") as f:
    f.readline()  # skip header
    for line in f:
        parts = line.strip().split("\t")
        cls = parts[2]
        if cls.startswith("pseudogene"):  # all confidence levels
            # Strip species prefix: "Ococ_OcoChr01G00010" -> "OcoChr01G00010"
            pseudo_ids.add(parts[0].split("_", 1)[1])

skip_gene = False
with open("Ococ/Oco_clean.gff3") as fin, open("Ococ/Oco_nopseudo.gff3", "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line); continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            fout.write(line); continue
        if parts[2] == "gene":
            gene_id = next((a[3:] for a in parts[8].split(";") if a.startswith("ID=")), None)
            skip_gene = gene_id in pseudo_ids
            if not skip_gene:
                fout.write(line)
        elif not skip_gene:
            fout.write(line)
```

To extract a pseudogene-only GFF3 (with classification annotations):

```python
# Same setup, but invert the logic:
if parts[2] == "gene" and gene_id in pseudo_ids:
    cls, score = pseudo_meta[gene_id]
    parts[8] = f"{parts[8]};pseudogene_class={cls};confidence={score}"
    fout.write("\t".join(parts) + "\n")
```

### Example 5: Disable pseudogene detection

For runs where you only need gene families (no pseudogene analysis):

```bash
python family_finder.py \
  --protein-dir data/pep --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output_5sp --no-pseudogene-detection
```

### Example 6: Adjust pruning strictness

For divergent species (e.g., paleopolyploids), loosen the distance ratio threshold:

```bash
python family_finder.py \
  --protein-dir data/pep --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output_loose --threshold 8.0 --threads 16
```

### Example 7: IQ-TREE for higher-quality trees

Slower but more accurate gene trees with bootstrap support:

```bash
python family_finder.py \
  --protein-dir data/pep --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output_iqtree --tree-builder iqtree --threads 16
```

### Example 8: Run codeml selection analysis

Compute dN/dS for each confirmed family using PAML:

```bash
python family_finder.py \
  --protein-dir data/pep --cds-dir data/cds \
  --species-tree data/species_tree.nwk \
  --outdir output_5sp --run-codeml
# Results in: output_5sp/codeml/<family_id>/<model>/results.txt
```

---

## Example Results: 5-species cactus CAM gene analysis

Input: 4 cactus species + 1 alternative annotation

| Species | Prefix | Genes | Source |
|---|---|---|---|
| *Mesembryanthemum crystallinum* | Mcry | 25,226 | MAKER |
| *Opuntia cochenillifera* | Ococ | 33,745 | MAKER |
| *Opuntia basilaris* | Obas | 28,244 | MAKER |
| *Carnegiea gigantea* | Cgig | 29,163 | MAKER |
| *Carnegiea gigantea* | CgigH | 27,583 | Helixer |

Results (10 rounds, 143,961 input sequences):

| Round | New Families | Outlier Pool |
|---|---|---|
| 1 | 16,250 | 23,707 |
| 2 | 964 | 16,741 |
| 3 | 188 | 15,466 |
| 4 | 73 | 14,883 |
| 5 | 28 | 14,735 |
| 6 | 59 | 14,180 |
| 7 | 16 | 14,081 |
| 8 | 11 | 14,000 |
| 9 | 6 | 13,968 |
| 10 | 1 | 13,963 |

**Total: 17,596 gene families**

### CAM gene families confirmed

All key CAM (Crassulacean Acid Metabolism) genes clustered correctly across 5 species:

| Gene | Family | Genes | Species | Notes |
|---|---|---|---|---|
| Ppck1/2 | R1_OG0002940 | 9 | 5 | Phosphoenolpyruvate carboxylase kinase |
| Ppc4 | R1_OG0010756 | 5 | 5 | PEP carboxylase 4 |
| Ppc1/2 | R1_OG0000093 | 43 | 5 | PEP carboxylase 1/2 |
| Ppcrk1 | R1_OG0005451 | 7 | 5 | PPCK-related kinase 1 |
| Ppcrk2/3 | R1_OG0008534 | 6 | 5 | PPCK-related kinase 2/3 |

CgigH (Helixer annotation) recovered PPC4 and PPCK genes that were present but unannotated in the MAKER annotation, confirming their presence in the *Carnegiea gigantea* genome.

### Value of iteration: CAM genes found only through later rounds

33 CAM pathway genes were placed into correct families only through iterative rounds (R2--R10), having been pruned or misassigned in R1. Examples include:

- **R2--R4:** NADP-Me3, Snrk2.3, CIPK24, Ers1/Etr1, Myb61
- **R9:** Hkl1

These genes would have been lost in a single-pass OrthoFinder analysis.

### Pseudogene detection results

Cross-species comparison from `species_comparison.tsv`:

| Species | Total Genes | Placed | Orphan | Pseudo High | Pseudo Medium | Pseudo Low | Pseudo Total | Pseudo Rate |
|---|---|---|---|---|---|---|---|---|
| *C. gigantea* (MAKER) | 29,163 | 27,465 | 1,698 | 0 | 916 | 2,110 | 3,026 | 10.4% |
| *C. gigantea* (Helixer) | 27,583 | 24,086 | 3,497 | 0 | 1,451 | 3,792 | 5,243 | 19.0% |
| *M. crystallinum* | 25,226 | 21,461 | 3,765 | 0 | 1,522 | 4,027 | 5,549 | 22.0% |
| *O. basilaris* | 28,244 | 25,880 | 2,364 | 0 | 1,414 | 2,597 | 4,011 | 14.2% |
| *O. cochenillifera* | 33,745 | 31,106 | 2,639 | 4 | 2,174 | 2,930 | 5,108 | 15.1% |

(values from initial 6-evidence run with strong/weak counting; rerun with the score-based classification produces 181 high / 4,020 medium / 142 low for Ococ — see updated table below).

### Top high-confidence pseudogenes (Ococ)

After re-running with score-based classification (`score >= 0.50`), 181 Ococ genes were classified as `pseudogene_high`. Examples:

| Gene | Confidence | Evidence |
|---|---|---|
| Ococ_OcoScaf00040937_210000_219999G00010 | 0.90 | 6 internal stops + truncated (32% of median) + long branch (16.7x) |
| Ococ_OcoScaf00009120G00010 | 0.90 | 1 internal stop + truncated (38%) + extreme long branch |
| Ococ_OcoScaf00040926_306000_316999G00010 | 0.90 | 1 internal stop + truncated (21%) + long branch (7.2x) |
| Ococ_OcoScaf00009117G00010 | 0.85 | 1 internal stop + truncated (30%) + GC3 outlier (z=3.36) |

### Pseudogene-enriched gene families (Ococ)

Top families with the highest pseudogene fractions, from `family_pseudogene_enrichment_Ococ.tsv`:

| Family | Members | Pseudo | Fraction | Avg Confidence |
|---|---|---|---|---|
| R2_OG0000680 | 5 | 5 | 100.0% | 0.40 |
| R2_OG0000040 | 24 | 21 | 87.5% | 0.42 |
| R2_OG0001015 | 4 | 3 | 75.0% | 0.15 |
| R1_OG0000467 | 18 | 11 | 61.1% | 0.15 |
| R2_OG0000678 | 5 | 3 | 60.0% | 0.27 |

These families likely represent dying gene families undergoing recent pseudogenization.

### GFF3 filtering results (Ococ)

Using the pseudogene calls to filter the original annotation:

| File | Genes | Size | Description |
|---|---|---|---|
| `Oco_clean.gff3` | 33,745 | 39 MB | Original annotation |
| `Oco_nopseudo.gff3` | 29,402 | 36 MB | Functional only (all 4,343 pseudogenes removed) |
| `Oco_pseudogenes.gff3` | 4,343 | 3.3 MB | Pseudogene candidates (with class + confidence attrs) |

This produces a cleaner gene set for downstream analyses (BUSCO, synteny, ortholog inference) by removing 12.9% of mis-annotated or non-functional gene models.

## Development History & Iterative Improvements

This pipeline was developed through multiple iterations, each addressing issues discovered during real data testing on cactus genomes.

### v1: Basic pipeline (`27a747b`)

- Initial implementation: OrthoFinder → MAFFT → FastTree → species-aware ratio pruning
- Config via YAML (later switched to JSON due to PyYAML dependency issue)
- Single-pass pruning using observed/expected distance ratio

**Problem discovered**: Ratio-only pruning caught only 18 out of 14,712 orthogroups (0.12%). The threshold of 5.0 was too permissive, and the placeholder species tree with uniform distances didn't reflect true evolutionary distances.

### v2: CDS-based trees + ete4 fixes (`f2b5fcb`)

- Switched from protein trees to **CDS-based trees** (`FastTree -nt -gtr -gamma`) for better resolution
- Added pal2nal for protein-guided codon alignment
- Fixed ete4 API (v4 vs v3): `tree.leaves()` returns generator, `get_distance()` requires tree method
- Fixed MAFFT_BINARIES environment conflict
- Fixed OrthoFinder output directory handling

**Lesson**: CDS trees provide better branch length estimates for pruning because synonymous substitutions add signal that protein trees lose.

### v3: Micromamba environment (`22df265`)

- Created `environment.yml` with all dependencies pinned
- Resolved OrthoFinder `rich` dependency issue
- Fixed OrthoFinder shebang pointing to wrong Python

### v4: TreeShrink integration (`e3e9ad8`)

- Added **two-stage pruning**:
  - Stage 1: TreeShrink for statistical branch-length outlier detection
  - Stage 2: Species-aware ratio filter on TreeShrink survivors
- TreeShrink uses quantile-based detection (default q=0.05) — catches outliers that ratio-based methods miss

**Problem discovered**: TreeShrink requires Python <=3.9, conflicts with OrthoFinder env (Python 3.11+). Pipeline gracefully skips TreeShrink when unavailable.

### v5: Robustness + performance (`b3dfd99`)

Addressed issues found during the 5-species (143K sequences) production run:

- **pal2nal failures diagnosed**: All 17 failures traced to Ococ (Opuntia cochenillifera) — 161 genes with internal stop codons in CDS, indicating annotation errors (frameshifts, pseudogenes)
- **Internal stop codon filter**: `_filter_internal_stops()` removes problematic genes before pal2nal, logs which genes and why
- **Protein tree fallback**: When codon alignment fails entirely, builds tree from protein alignment instead of crashing
- **Auto-detect sequence type**: FastTree automatically uses `-nt -gtr -gamma` for nucleotide or protein mode based on input
- **Pickle overhead fix**: Workers now receive only their OG's sequences (~5-50 seqs) instead of the entire pool (143K seqs). Reduced inter-process data transfer by ~1000x
- **Lightweight returns**: Workers return gene ID sets instead of sequence dictionaries, further reducing memory
- **OrthoFinder diamond compatibility**: Removed `--ignore-warnings` flag incompatible with diamond 2.1.24
- **PATH management**: Auto-adds the OrthoFinder environment to PATH for mcl/diamond

### Performance impact of optimizations

| Version | Round 1 time (est.) | Memory per worker |
|---|---|---|
| v1-v4 | ~7 hours | ~700MB (full pool copied) |
| v5 | ~5 hours (projected) | ~10MB (per-OG seqs only) |

### Lessons learned

1. **Test with real data early** — synthetic tests pass but production data reveals annotation quality issues, tool incompatibilities, and performance bottlenecks
2. **Annotation quality varies widely** — Ococ had 161 genes with internal stops; always validate input data
3. **Pickle serialization matters** — Python multiprocessing copies all arguments to each worker; keep payloads small
4. **Graceful degradation** — fall back to protein trees when CDS fails, skip TreeShrink when unavailable, log everything for diagnosis
5. **Two annotations of the same genome** can be treated as separate "species" with near-zero distance to validate gene recovery

## Known Issues

- **Ococ annotation quality**: 161 genes have internal stop codons in CDS, causing pal2nal failures. The pipeline auto-filters these and falls back to protein trees.
- **TreeShrink**: Requires Python <=3.9; may not install in newer environments. Pipeline works without it (Stage 2 pruning only).
- **Large outlier pools** *(v1 run; the v2 figure is pending)*: After HMMER rescue, 7,605 genes remain unplaced. Most are from orthogroups with <4 members (below `min_orthogroup_size`) or are species-specific orphans.

## Annotation & Naming Stack (post-run)

Standalone CLIs over a finished run; every axis is gate-validated on the PEPC
clan (known answer: EC 4.1.1.31). All are annotation layers — never
clustering criteria or membership filters. See methods.md §2.X.6–2.X.7.

| CLI | What it does |
|---|---|
| `annotate_families.py --ec` | Family EC consensus + Fitch EC-switch events. Sources: `--emapper` (eggNOG-mapper annotations) and/or `--clean-csv` (CLEAN maxsep), merged with the emapper EC authoritative and the CLEAN confidence attached (`steps/ec_sources.py`). ESM-ECForest failed the gate and remains only as a legacy `--ecforest-cache` fallback. |
| `annotation_matrix.py` | Merges every axis (emapper, CLEAN, foldseek→AFDB-SwissProt transfer TSV, DeepLoc, SignalP 6) into one per-gene TSV; with `--expected-ec` adds a membership verdict — member (≥2 axes support), intruder (0 support + confident contradiction), review (abstentions are not evidence). PEPC clan: 106 member / 11 review / 1 intruder. |
| `name_families.py` | Family/subfamily naming: direct annotation table + `--plaza-orthology` layer (DIAMOND best-hit/RBH vs PLAZA ath, built by `extract_plaza_orthology.py`), weighted majority (direct 1.0 > orthology 0.5) with provenance columns. |
| `subfamily_report.py` | Subfamily diagnostics: SDP scan, `--species-tree` monophyly/MRCA attribution (root-span rule; taxonomy TSV = optional labels), foldseek structural coherence. With `--focal-subfamily` + selection evidence (RELAX/MEME/aBSREL JSONs, expression share, signal region) also writes `subfunctionalization.md` — a narrative verdict (sub- vs neo-functionalization) explaining HOW the subfamily diverged. |
| `steps/hmm_chunks.py` | Optional chunked hmmsearch for the rescue (`hmmer_chunk_size`, default 0 = single run). Splits the profile DB, writes a **SLURM-optional** runner (`sbatch --wait --array` when available, bounded local pool otherwise), and merges — refusing any chunk missing HMMER's `# [ok]` terminator. Needed at 15sp scale: cost is profiles x sequences, so one run extrapolates past the 3-day limit. |
| `beb_cross.py` | codeml branch-site BEB sites × DeepLoc signal windows (alignment-column overlap; requires the branch-site run to use `cleandata = 0`). |

Caveat that motivated the design: no sequence-, orthology-, or
structure-based tool sees residue-level catalytic loss (the SF2 tandem-array
lesson) — EC and structure calls are read alongside domain/catalytic-residue
evidence.

### Running the stack on one family

`annotate_stack.py` declares every axis as one replayable plan instead of
five hand-assembled invocations. `--dry-run` prints the commands and names
whatever is blocked, without touching anything:

```bash
# what would run, and what is missing
python annotate_stack.py --family-fasta clan.fa --workdir '~/annot/clan' \
    --structures '~/pdb' --expected-ec 4.1.1.31 --dry-run

# run the heavy axes on a GPU host, then merge locally
python annotate_stack.py --family-fasta clan.fa --workdir '~/annot/clan' \
    --structures '~/pdb' --expected-ec 4.1.1.31 --host gpu --outdir annot_clan
```

The merge produces `annotation_matrix.tsv` — one row per gene, one column
block per axis, plus a membership verdict when `--expected-ec` is given:

| verdict | meaning |
|---|---|
| `member` | ≥2 axes support the expected EC |
| `intruder` | no axis supports it AND at least one **confidently** contradicts it (a wrong-EC structural hit at qTM ≥ 0.6, or a wrong EC at confidence ≥ 0.5) |
| `review` | everything else — an abstention (CLEAN confidence ≈ 0, no structural hit, gene absent from an axis) is **not** evidence of intrusion |

That last row is the point of the design: on the PEPC clan, the remote
members abstain on every axis while the one true intruder is contradicted by
two independent axes (106 member / 11 review / 1 intruder).

### Explaining a subfamily split

With the selection evidence in hand, `subfamily_report.py` writes a narrative
verdict rather than only diagnostics:

```bash
python subfamily_report.py \
    --alignment clan.aln --groups subfamilies.tsv --ref-seq ATH_AT1G53310.2 \
    --species-tree data/species_tree.nwk --outdir subfam_report \
    --focal-subfamily SF3 --family-name "PEPC clan" \
    --relax-json relax.json --absrel-json absrel.json \
    --meme-json meme.json --meme-region 165-209 \
    --branchsite-mlc alt/results.txt --branchsite-lnl=-49044.15,-49051.83 \
    --expression-share 0.74 --signal-partition "N-terminal NLS region"
```

`subfunctionalization.md` then states the verdict and the reasoning per axis:
relaxed selection **plus** a partitioned ancestral role **plus** no positive
selection on the stem ⇒ subfunctionalization; positive selection on the stem
or a gained retargeting event ⇒ neofunctionalization; missing tests ⇒
`insufficient evidence` rather than a guess.

## Project Structure

```
family_finder/
  family_finder.py          # CLI entry point
  find_pseudogenes.py       # Standalone pseudogene detection script
  find_neofunctionalization.py  # DeepLoc retargeting + branch-site dN/dS
  annotate_families.py      # EC layer + tier-3 fold-list prep
  annotation_matrix.py      # All-axis merge + membership verdicts
  name_families.py          # Family/subfamily naming (direct + PLAZA orthology)
  subfamily_report.py       # Subfamily diagnostics report
  beb_cross.py              # BEB sites x signal windows
  extract_plaza_orthology.py  # DIAMOND fwd/rev -> orthology TSV (RBH per species)
  config.py                 # Config dataclass + JSON loader
  pipeline.py               # Iterative loop orchestrator
  steps/
    orthofinder.py           # OrthoFinder wrapper + parser
    align.py                 # MAFFT protein + pal2nal codon alignment
    tree.py                  # FastTree / IQ-TREE wrapper
    prune.py                 # TreeShrink + species-aware pruning
    codeml.py                # PAML/codeml wrapper
    hmmer_rescue.py          # HMMER profile-based rescue of unplaced genes
    profile_assign.py        # Per-round HMM profile assignment tier
    pseudogene.py            # Pseudogene detection (evidence collection + reporting)
    ec_sources.py            # eggNOG-mapper + CLEAN parsers/merge
    family_naming.py         # Naming vote logic (direct + orthology weights)
    subfamily.py             # HOG subfamilies, SDP scan, monophyly attribution
    deeploc.py / retargeting.py / esm.py / epa.py  # annotation & placement tiers
  utils/
    seqio.py                 # FASTA I/O, species splitting
    species.py               # Species tree loading, pairwise distances
    parallel.py              # ProcessPoolExecutor wrapper
    checkpoint.py            # Resume/checkpoint logic
    logging_setup.py         # Logging configuration
```

## License

MIT
