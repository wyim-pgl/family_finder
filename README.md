# family_finder

Iterative gene family construction pipeline with species-aware pruning.

## Overview

`family_finder` builds gene families by repeatedly running OrthoFinder, aligning sequences, building trees, and pruning outliers. Outlier sequences from each round are re-clustered in the next round, allowing displaced genes to find their true families.

```
Round 1: All seqs → OrthoFinder → per-OG align/tree/prune → confirmed families + outliers
Round 2: Outliers → OrthoFinder → per-OG align/tree/prune → new families + outliers
Round N: Repeat until convergence
                                          ↓
HMMER Rescue: Unplaced genes → hmmsearch vs family HMM profiles → rescued into existing families
                                          ↓
Pseudogene Detection (optional): All genes scored across 6 evidence types → candidates classified high/medium/low
```

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

### Option 1: micromamba (recommended)

```bash
# Install micromamba if not available
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# Create environment from environment.yml
micromamba create -f environment.yml
micromamba activate family_finder

# Verify installation
python -c "from ete4 import Tree; from Bio import SeqIO; print('OK')"
orthofinder -h | head -1
mafft --version
FastTree 2>&1 | head -1
pal2nal.pl 2>&1 | head -1
```

### Option 2: conda

```bash
conda env create -f environment.yml
conda activate family_finder
```

### Option 3: Manual installation

```bash
# Create conda env with bioconda tools
conda create -n family_finder -c conda-forge -c bioconda \
  python=3.11 orthofinder mafft fasttree iqtree paml pal2nal ete4 biopython rich

conda activate family_finder

# Optional: TreeShrink (requires Python <=3.9, separate env recommended)
conda create -n treeshrink -c bioconda treeshrink
```

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
| `Cannot run MCL` | Ensure `mcl` is in PATH: `which mcl` or install via `conda install -c bioconda mcl` |
| `MAFFT_BINARIES` conflict | Pipeline auto-clears this env var. If issues persist: `unset MAFFT_BINARIES` |
| `No module named 'Bio'` | Run with the correct Python: the one inside your conda/micromamba env |
| TreeShrink won't install | Requires Python <=3.9. Use a separate env or skip (Stage 2 pruning still works) |
| `ModuleNotFoundError: rich` | `pip install rich` in the OrthoFinder environment |

## Quick Start

```bash
# 1. Install
micromamba create -f environment.yml && micromamba activate family_finder

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
  Mcry.pep.fa      # Mammillaria protein sequences
  Ococ.pep.fa      # Opuntia cochenillifera proteins
  Cgig.pep.fa      # Carnegiea gigantea proteins

data/cds/
  Mcry.cds.fa      # Mammillaria CDS sequences
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

Recommended approach: build from low-copy-number orthologs using ASTRAL.

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
| `--species-tree` | Newick species tree (e.g., from ASTRAL) |
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
  "orthofinder_threads": 8,
  "n_workers": 8,
  "max_rounds": 10,
  "min_orthogroup_size": 4,
  "distance_ratio_threshold": 5.0,
  "treeshrink_quantile": 0.05,
  "tree_builder": "fasttree"
}
```

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

## Example: 5-species cactus CAM gene analysis

Input: 4 cactus species + 1 alternative annotation

| Species | Prefix | Genes | Source |
|---|---|---|---|
| *Mammillaria cristata* | Mcry | 25,226 | MAKER |
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
- **PATH management**: Auto-adds orthofinder conda env to PATH for mcl/diamond

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
- **TreeShrink**: Requires Python <=3.9; may not install in newer conda environments. Pipeline works without it (Stage 2 pruning only).
- **Large outlier pools**: After HMMER rescue, 7,605 genes remain unplaced. Most are from orthogroups with <4 members (below `min_orthogroup_size`) or are species-specific orphans.

## Project Structure

```
family_finder/
  family_finder.py          # CLI entry point
  find_pseudogenes.py       # Standalone pseudogene detection script
  config.py                 # Config dataclass + JSON loader
  pipeline.py               # Iterative loop orchestrator
  steps/
    orthofinder.py           # OrthoFinder wrapper + parser
    align.py                 # MAFFT protein + pal2nal codon alignment
    tree.py                  # FastTree / IQ-TREE wrapper
    prune.py                 # TreeShrink + species-aware pruning
    codeml.py                # PAML/codeml wrapper
    hmmer_rescue.py          # HMMER profile-based rescue of unplaced genes
    pseudogene.py            # Pseudogene detection (evidence collection + reporting)
  utils/
    seqio.py                 # FASTA I/O, species splitting
    species.py               # Species tree loading, pairwise distances
    parallel.py              # ProcessPoolExecutor wrapper
    checkpoint.py            # Resume/checkpoint logic
    logging_setup.py         # Logging configuration
```

## License

MIT
