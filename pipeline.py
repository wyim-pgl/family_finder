"""Iterative pipeline orchestrator.

Runs repeated rounds of OrthoFinder → align → tree → prune.
Outliers from each round feed into the next round.
"""

import json
import logging
import shutil
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

from config import Config
from steps.orthofinder import run_orthofinder, parse_orthogroups
from steps.align import align_protein, codon_align
from steps.tree import build_tree
from steps.prune import prune_orthogroup
from utils.seqio import write_fasta, split_by_species
from utils.species import (
    load_species_tree,
    compute_pairwise_distances,
    get_species,
)
from utils.checkpoint import save_checkpoint, find_latest_checkpoint
from utils.parallel import parallel_map

logger = logging.getLogger("family_finder")


def process_single_orthogroup(args: tuple):
    """Process a single orthogroup: align → tree → prune.

    Called by parallel_map. Takes a single tuple of all arguments.
    """
    (og_id, gene_ids, protein_pool, cds_pool,
     expected_distances, config, round_dir) = args

    og_dir = round_dir / "orthogroups" / og_id
    og_dir.mkdir(parents=True, exist_ok=True)

    # Skip if too small
    if len(gene_ids) < config.min_orthogroup_size:
        outlier_seqs = {gid: protein_pool[gid] for gid in gene_ids if gid in protein_pool}
        return (og_id, None, outlier_seqs)

    try:
        # 1. Extract protein sequences
        prot_seqs = {gid: protein_pool[gid] for gid in gene_ids if gid in protein_pool}
        if len(prot_seqs) < config.min_orthogroup_size:
            return (og_id, None, prot_seqs)

        # 2. Align
        prot_aln = align_protein(prot_seqs, og_dir / "proteins.afa", config)

        # 3. Build tree
        tree_path = build_tree(prot_aln, og_dir / "tree.nwk", config)

        # 4. Species-aware pruning
        gene_to_species = {
            gid: get_species(gid, config.species_delimiter) for gid in prot_seqs
        }
        confirmed, outliers = prune_orthogroup(
            str(tree_path), gene_to_species, expected_distances, config
        )

        # 5. If confirmed set is large enough, produce final outputs
        if len(confirmed) >= config.min_orthogroup_size:
            confirmed_prots = {gid: protein_pool[gid] for gid in confirmed}
            write_fasta(confirmed_prots, str(og_dir / "confirmed_proteins.fa"))

            # Re-align without outliers
            confirmed_aln = align_protein(
                confirmed_prots, og_dir / "confirmed_proteins.afa", config
            )

            # Re-build tree
            build_tree(confirmed_aln, og_dir / "confirmed_tree.nwk", config)

            # Codon alignment if CDS available
            confirmed_cds = {gid: cds_pool[gid] for gid in confirmed if gid in cds_pool}
            if len(confirmed_cds) >= config.min_orthogroup_size:
                codon_align(
                    og_dir / "confirmed_proteins.afa",
                    confirmed_cds,
                    og_dir / "confirmed_codon.afa",
                    config,
                )

            outlier_seqs = {gid: protein_pool[gid] for gid in outliers if gid in protein_pool}
            return (og_id, confirmed, outlier_seqs)
        else:
            # Too small after pruning → all go to outlier pool
            all_seqs = {gid: protein_pool[gid] for gid in gene_ids if gid in protein_pool}
            return (og_id, None, all_seqs)

    except Exception as e:
        logger.error(f"Failed to process {og_id}: {e}")
        all_seqs = {gid: protein_pool[gid] for gid in gene_ids if gid in protein_pool}
        return (og_id, None, all_seqs)


def run(
    protein_dir: str,
    cds_dir: str,
    species_tree_path: str,
    outdir: str,
    config: Config,
    resume: bool = False,
):
    """Run the iterative gene family construction pipeline.

    Args:
        protein_dir: Directory of per-species protein FASTA files.
        cds_dir: Directory of per-species CDS FASTA files.
        species_tree_path: Path to Newick species tree.
        outdir: Output directory.
        config: Pipeline configuration.
        resume: Whether to resume from checkpoint.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load species tree and compute expected distances
    logger.info(f"Loading species tree from {species_tree_path}")
    species_tree = load_species_tree(species_tree_path)
    expected_distances = compute_pairwise_distances(species_tree)
    logger.info(f"Computed pairwise distances for {len(set(k[0] for k in expected_distances))} species")

    # Load all sequences
    from utils.seqio import build_seq_map
    logger.info(f"Loading protein sequences from {protein_dir}")
    current_pool = build_seq_map(protein_dir)
    logger.info(f"Loaded {len(current_pool)} protein sequences")

    logger.info(f"Loading CDS sequences from {cds_dir}")
    cds_pool = build_seq_map(cds_dir)
    logger.info(f"Loaded {len(cds_pool)} CDS sequences")

    # Resume handling
    start_round = 1
    all_confirmed_families: Dict[str, Set[str]] = {}

    if resume:
        cp = find_latest_checkpoint(outdir)
        if cp and cp["status"] == "completed":
            start_round = cp["round_number"] + 1
            # Reload the outlier pool from the completed round
            pool_fasta = outdir / f"round_{cp['round_number']:02d}" / "outlier_pool.fa"
            if pool_fasta.exists():
                from utils.seqio import read_fasta
                current_pool = read_fasta(str(pool_fasta))
                logger.info(f"Resuming from round {start_round} with {len(current_pool)} sequences")

    round_num = start_round - 1
    rounds_with_no_new = 0

    while round_num < config.max_rounds:
        round_num += 1
        round_dir = outdir / f"round_{round_num:02d}"
        round_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"=== Round {round_num} === ({len(current_pool)} sequences)")
        save_checkpoint(round_dir, round_num, len(current_pool), "started")

        # Step 1: Split pool by species
        input_dir = round_dir / "input"
        split_by_species(current_pool, str(input_dir), config.species_delimiter)

        # Step 2: Run OrthoFinder
        of_dir = round_dir / "orthofinder"
        try:
            results_dir = run_orthofinder(input_dir, of_dir, config)
        except Exception as e:
            logger.error(f"OrthoFinder failed in round {round_num}: {e}")
            break

        # Step 3: Parse orthogroups
        orthogroups = parse_orthogroups(results_dir)
        logger.info(f"Round {round_num}: {len(orthogroups)} orthogroups found")

        save_checkpoint(round_dir, round_num, len(current_pool), "processing")

        # Step 4: Process each orthogroup
        work_items = [
            (og_id, gene_ids, current_pool, cds_pool,
             expected_distances, config, round_dir)
            for og_id, gene_ids in orthogroups.items()
        ]

        results = parallel_map(
            process_single_orthogroup,
            work_items,
            n_workers=config.n_workers,
        )

        # Step 5: Collect confirmed families and outliers
        new_families = {}
        new_outlier_pool = {}

        assigned_genes = set()
        for r in results:
            if r is None:
                continue
            og_id, confirmed, outlier_seqs = r
            assigned_genes.update(outlier_seqs.keys())

            if confirmed is not None:
                family_id = f"R{round_num}_{og_id}"
                new_families[family_id] = confirmed
                all_confirmed_families[family_id] = confirmed
                assigned_genes.update(confirmed)

            new_outlier_pool.update(outlier_seqs)

        # Add unassigned genes (not in any orthogroup) to outlier pool
        all_og_genes = set()
        for gene_ids in orthogroups.values():
            all_og_genes.update(gene_ids)
        for gid in current_pool:
            if gid not in all_og_genes:
                new_outlier_pool[gid] = current_pool[gid]

        # Save outlier pool for resume
        write_fasta(new_outlier_pool, str(round_dir / "outlier_pool.fa"))

        # Step 6: Log round statistics
        stats = {
            "round": round_num,
            "input_sequences": len(current_pool),
            "orthogroups": len(orthogroups),
            "new_families": len(new_families),
            "outlier_pool_size": len(new_outlier_pool),
        }
        with open(round_dir / "round_stats.json", "w") as f:
            json.dump(stats, f, indent=2)

        logger.info(
            f"Round {round_num} summary: "
            f"{len(new_families)} new families, "
            f"{len(new_outlier_pool)} outliers remaining"
        )

        save_checkpoint(round_dir, round_num, len(new_outlier_pool), "completed")

        # Step 7: Check convergence
        if len(new_families) == 0:
            rounds_with_no_new += 1
        else:
            rounds_with_no_new = 0

        if rounds_with_no_new >= config.convergence_no_new_families:
            logger.info(f"Converged: no new families for {rounds_with_no_new} consecutive rounds")
            break

        if len(new_outlier_pool) < config.convergence_threshold:
            logger.info(f"Converged: outlier pool ({len(new_outlier_pool)}) below threshold ({config.convergence_threshold})")
            break

        current_pool = new_outlier_pool

    # Final: Assemble all confirmed families
    _write_final_output(all_confirmed_families, current_pool, cds_pool, outdir, config)

    logger.info(f"Pipeline complete: {len(all_confirmed_families)} total families across {round_num} rounds")


def _write_final_output(
    families: Dict[str, Set[str]],
    remaining_pool: Dict[str, str],
    cds_pool: Dict[str, str],
    outdir: Path,
    config: Config,
):
    """Write final summary and copy confirmed family files to final_families/."""
    final_dir = outdir / "final_families"
    final_dir.mkdir(parents=True, exist_ok=True)

    # Copy confirmed family outputs
    for family_id, gene_ids in families.items():
        # Parse round and OG from family_id (e.g., "R1_OG0000000")
        parts = family_id.split("_", 1)
        round_num = int(parts[0][1:])
        og_id = parts[1]
        src_dir = outdir / f"round_{round_num:02d}" / "orthogroups" / og_id
        dst_dir = final_dir / family_id
        if src_dir.exists():
            shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)

    # Write summary TSV
    with open(outdir / "summary.tsv", "w") as f:
        f.write("family_id\tround\tn_genes\tn_species\tgene_list\n")
        for family_id, gene_ids in sorted(families.items()):
            parts = family_id.split("_", 1)
            round_num = parts[0][1:]
            species = {
                get_species(gid, config.species_delimiter) for gid in gene_ids
            }
            gene_list = ",".join(sorted(gene_ids))
            f.write(f"{family_id}\t{round_num}\t{len(gene_ids)}\t{len(species)}\t{gene_list}\n")

    logger.info(f"Final output written to {final_dir}")
    logger.info(f"Summary: {outdir / 'summary.tsv'}")
