# family_finder

Iterative gene family construction pipeline with species-aware pruning.

## Overview

`family_finder` builds gene families by repeatedly running OrthoFinder, aligning sequences, building trees, and pruning outliers. Outlier sequences from each round are re-clustered in the next round, allowing displaced genes to find their true families.

```
Round 1: All seqs → OrthoFinder → per-OG align/tree/prune → confirmed families + outliers
Round 2: Outliers → OrthoFinder → per-OG align/tree/prune → new families + outliers
Round N: Repeat until convergence
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

## Installation

### Dependencies (micromamba)

```bash
micromamba create -f environment.yml
micromamba activate family_finder
```

Required tools:
- OrthoFinder (+ diamond, mcl, famsa, FastTree)
- MAFFT
- FastTree or IQ-TREE
- pal2nal
- TreeShrink (optional, for Stage 1 pruning)
- PAML/codeml (optional, for selection analysis)

Python packages: ete4, Biopython

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

### Required inputs

| Argument | Description |
|---|---|
| `--protein-dir` | Directory of per-species protein FASTA files (e.g., `Mcry.pep.fa`) |
| `--cds-dir` | Directory of per-species CDS FASTA files (e.g., `Mcry.cds.fa`) |
| `--species-tree` | Newick species tree (e.g., from ASTRAL) |
| `--outdir` | Output directory |

### Gene ID format

Gene IDs must follow `SpeciesPrefix_GeneID` format. The species is extracted from the prefix before the first `_`:
- `Mcry_Mcr1G24690` → species `Mcry`
- `CgigH_Cgig_v2_SGP5p_31_000132.1` → species `CgigH`

Protein and CDS files must use matching gene IDs.

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
| `--verbose` | — | Debug logging |

### JSON config

All parameters can be set via a JSON config file:

```json
{
  "orthofinder_bin": "orthofinder",
  "mafft_bin": "mafft",
  "fasttree_bin": "FastTree",
  "pal2nal_bin": "pal2nal.pl",
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
  final_families/              # All confirmed families (copied)
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

## Known Issues

- **Ococ annotation quality**: 161 genes have internal stop codons in CDS, causing pal2nal failures. The pipeline auto-filters these and falls back to protein trees.
- **TreeShrink**: Requires Python <=3.9; may not install in newer conda environments. Pipeline works without it (Stage 2 pruning only).
- **Large outlier pools**: Most remaining outliers after convergence are from orthogroups with <4 members (below `min_orthogroup_size`), not from pruning failures.

## Project Structure

```
family_finder/
  family_finder.py          # CLI entry point
  config.py                 # Config dataclass + JSON loader
  pipeline.py               # Iterative loop orchestrator
  steps/
    orthofinder.py           # OrthoFinder wrapper + parser
    align.py                 # MAFFT protein + pal2nal codon alignment
    tree.py                  # FastTree / IQ-TREE wrapper
    prune.py                 # TreeShrink + species-aware pruning
    codeml.py                # PAML/codeml wrapper
  utils/
    seqio.py                 # FASTA I/O, species splitting
    species.py               # Species tree loading, pairwise distances
    parallel.py              # ProcessPoolExecutor wrapper
    checkpoint.py            # Resume/checkpoint logic
    logging_setup.py         # Logging configuration
```

## License

MIT
