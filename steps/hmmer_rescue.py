"""HMMER-based rescue step for unplaced genes after iterative pipeline convergence.

After OrthoFinder rounds converge, some genes remain unplaced because DIAMOND
cannot detect distant homologs. This module:
  1. Builds HMM profiles from each confirmed family's protein alignment.
  2. Searches unplaced genes against these profiles with hmmsearch.
  3. Assigns hits (E-value < threshold) to the best-matching family.
  4. Re-aligns and re-builds trees for families that gained new members.
"""

import logging
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

from config import Config
from utils.seqio import read_fasta, write_fasta

# Prevent MAFFT_BINARIES version conflict in conda/micromamba environments
os.environ.pop("MAFFT_BINARIES", None)

logger = logging.getLogger("family_finder")


def _build_hmm_worker(args: tuple) -> Optional[str]:
    """Build a single HMM profile (for parallel execution).

    Returns family_id on success, None on failure.
    """
    family_id, alignment_path, hmm_path, hmmbuild_bin = args
    cmd = [hmmbuild_bin, "--amino", "-n", family_id, str(hmm_path), str(alignment_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return None
    return family_id


def _concat_hmms(hmm_dir: Path, combined_path: Path) -> Optional[Path]:
    """Concatenate individual HMM files into a single HMM database and press it."""
    hmm_files = sorted(hmm_dir.glob("*.hmm"))
    if not hmm_files:
        return None

    with open(combined_path, "wb") as out_f:
        for hf in hmm_files:
            with open(hf, "rb") as inf:
                while True:
                    chunk = inf.read(8192)
                    if not chunk:
                        break
                    out_f.write(chunk)

    return combined_path


def _run_hmmsearch(
    hmm_db: Path, query_fasta: Path, outpath: Path, config: Config
) -> Path:
    """Run hmmsearch of unplaced proteins against family HMM profiles."""
    cmd = [
        config.hmmsearch_bin,
        "--tblout", str(outpath),
        "--noali",
        "-E", str(config.hmmer_evalue),
        "--cpu", str(min(config.n_workers, 8)),
        str(hmm_db),
        str(query_fasta),
    ]
    logger.debug(f"Running hmmsearch: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"hmmsearch failed: {result.stderr[:300]}")
        raise RuntimeError(f"hmmsearch failed with return code {result.returncode}")
    return outpath


def _parse_hmmsearch_tblout(
    tblout_path: Path, evalue_cutoff: float
) -> Dict[str, Tuple[str, float]]:
    """Parse hmmsearch --tblout output.

    Returns dict of gene_id -> (best_family_id, evalue) for genes passing cutoff.
    Only keeps the best hit (lowest E-value) per gene.
    """
    best_hits: Dict[str, Tuple[str, float]] = {}

    # HMMER3 --tblout columns (whitespace-delimited):
    # 0: target name, 1: target accession, 2: query name (HMM),
    # 3: query accession, 4: full-seq E-value, 5: full-seq score, ...
    with open(tblout_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            gene_id = fields[0]       # target name (unplaced gene)
            family_id = fields[2]     # query name (HMM profile = family_id)
            evalue = float(fields[4]) # full-sequence E-value

            if evalue > evalue_cutoff:
                continue

            if gene_id not in best_hits or evalue < best_hits[gene_id][1]:
                best_hits[gene_id] = (family_id, evalue)

    return best_hits


def rescue_unplaced(
    families: Dict[str, Set[str]],
    unplaced_pool: Dict[str, str],
    protein_pool: Dict[str, str],
    cds_pool: Dict[str, str],
    outdir: Path,
    config: Config,
) -> Dict[str, Set[str]]:
    """HMMER rescue: assign unplaced genes to existing families via HMM profile search.

    Args:
        families: Dict of family_id -> set of confirmed gene_ids.
        unplaced_pool: Dict of gene_id -> protein sequence for unplaced genes.
        protein_pool: Full protein pool (all genes, for re-alignment).
        cds_pool: Full CDS pool (all genes).
        outdir: Pipeline output directory.
        config: Pipeline configuration.

    Returns:
        Updated families dict with rescued genes added.
    """
    rescue_dir = outdir / "hmmer_rescue"
    rescue_dir.mkdir(parents=True, exist_ok=True)
    hmm_dir = rescue_dir / "hmm_profiles"
    hmm_dir.mkdir(parents=True, exist_ok=True)

    if not unplaced_pool:
        logger.info("HMMER rescue: no unplaced genes to rescue")
        return families

    logger.info(f"HMMER rescue: {len(unplaced_pool)} unplaced genes, "
                f"{len(families)} family profiles to search against")

    # Step 1: Build HMM profiles in parallel
    work_items = []
    for family_id in families:
        aln_path = _find_family_alignment(family_id, outdir)
        if aln_path is None:
            continue
        hmm_path = hmm_dir / f"{family_id}.hmm"
        if not hmm_path.exists():
            work_items.append((family_id, str(aln_path), str(hmm_path), config.hmmbuild_bin))

    # Already-built HMMs
    n_existing = sum(1 for f in hmm_dir.glob("*.hmm") if f.stat().st_size > 0)

    logger.info(f"HMMER rescue: building {len(work_items)} HMM profiles "
                f"({n_existing} already exist), using {config.n_workers} workers")

    n_built = n_existing
    n_failed = 0
    if work_items:
        n_workers = min(config.n_workers, len(work_items))
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_build_hmm_worker, item): item[0]
                       for item in work_items}
            done = 0
            for future in as_completed(futures):
                try:
                    result = future.result()
                    if result is not None:
                        n_built += 1
                    else:
                        n_failed += 1
                except Exception as e:
                    n_failed += 1
                    fam_id = futures[future]
                    logger.debug(f"  hmmbuild worker exception for {fam_id}: {e}")
                done += 1
                if done % 2000 == 0:
                    logger.info(f"  hmmbuild progress: {done}/{len(work_items)}")
    if n_failed:
        logger.info(f"  hmmbuild: {n_failed} profiles failed")

    logger.info(f"HMMER rescue: {n_built} HMM profiles ready")

    if n_built == 0:
        logger.warning("HMMER rescue: no HMM profiles built, skipping")
        return families

    # Step 2: Concatenate HMMs into single database
    combined_hmm = rescue_dir / "all_families.hmm"
    if _concat_hmms(hmm_dir, combined_hmm) is None:
        logger.warning("HMMER rescue: no HMM files to concatenate")
        return families

    # hmmpress the database (-f to force overwrite existing index)
    press_cmd = [config.hmmpress_bin, "-f", str(combined_hmm)]
    press_result = subprocess.run(press_cmd, capture_output=True, text=True)
    if press_result.returncode != 0:
        logger.error(f"hmmpress failed: {press_result.stderr[:300]}")
        return families

    # Step 3: Write unplaced genes to FASTA
    unplaced_fasta = rescue_dir / "unplaced_proteins.fa"
    write_fasta(unplaced_pool, str(unplaced_fasta))

    # Step 4: Run hmmsearch
    tblout = rescue_dir / "hmmsearch_results.tblout"
    _run_hmmsearch(combined_hmm, unplaced_fasta, tblout, config)

    # Step 5: Parse results and assign genes to families
    hits = _parse_hmmsearch_tblout(tblout, config.hmmer_evalue)

    if not hits:
        logger.info("HMMER rescue: no significant hits found")
        return families

    # Group rescued genes by family
    rescued_by_family: Dict[str, Set[str]] = {}
    for gene_id, (family_id, evalue) in hits.items():
        rescued_by_family.setdefault(family_id, set()).add(gene_id)
        logger.debug(f"  Rescued {gene_id} → {family_id} (E={evalue:.1e})")

    logger.info(f"HMMER rescue: {len(hits)} genes rescued into "
                f"{len(rescued_by_family)} families")

    # Step 6: Add rescued genes to families and re-align affected families
    updated_families = dict(families)
    for family_id, new_genes in rescued_by_family.items():
        if family_id not in updated_families:
            continue

        old_members = updated_families[family_id]
        updated_families[family_id] = old_members | new_genes

        # Re-align the expanded family
        _realign_family(
            family_id, updated_families[family_id],
            protein_pool, cds_pool, outdir, config,
        )

    # Write rescue summary
    _write_rescue_summary(hits, rescue_dir)

    return updated_families


def _find_family_alignment(family_id: str, outdir: Path) -> Optional[Path]:
    """Find the confirmed protein alignment for a family."""
    # Check final_families first
    aln = outdir / "final_families" / family_id / "confirmed_proteins.afa"
    if aln.exists():
        return aln

    # Fall back to round directory
    try:
        parts = family_id.split("_", 1)
        round_num = int(parts[0][1:])
        og_id = parts[1]
    except (ValueError, IndexError):
        return None
    aln = outdir / f"round_{round_num:02d}" / "orthogroups" / og_id / "confirmed_proteins.afa"
    if aln.exists():
        return aln

    return None


def _realign_family(
    family_id: str,
    gene_ids: Set[str],
    protein_pool: Dict[str, str],
    cds_pool: Dict[str, str],
    outdir: Path,
    config: Config,
):
    """Re-align a family after adding rescued genes."""
    from steps.align import align_protein, codon_align
    from steps.tree import build_tree

    rescue_fam_dir = outdir / "hmmer_rescue" / "families" / family_id
    rescue_fam_dir.mkdir(parents=True, exist_ok=True)

    # Collect sequences
    prot_seqs = {gid: protein_pool[gid] for gid in gene_ids if gid in protein_pool}
    cds_seqs = {gid: cds_pool[gid] for gid in gene_ids if gid in cds_pool}

    if len(prot_seqs) < 2:
        return

    try:
        # Write input
        write_fasta(prot_seqs, str(rescue_fam_dir / "proteins.fa"))

        # Protein alignment
        prot_aln = align_protein(prot_seqs, rescue_fam_dir / "proteins.afa", config)

        # Codon alignment
        codon_aln = codon_align(
            rescue_fam_dir / "proteins.afa", cds_seqs,
            rescue_fam_dir / "codon.afa", config,
        )

        # Build tree
        if codon_aln is not None:
            build_tree(codon_aln, rescue_fam_dir / "tree.nwk", config)
        else:
            build_tree(prot_aln, rescue_fam_dir / "tree.nwk", config)

        logger.debug(f"  Re-aligned {family_id} with {len(prot_seqs)} members")
    except Exception as e:
        logger.warning(f"  Re-alignment failed for {family_id}: {e}")


def _write_rescue_summary(
    hits: Dict[str, Tuple[str, float]], rescue_dir: Path
):
    """Write TSV summary of rescued genes."""
    summary_path = rescue_dir / "rescue_summary.tsv"
    with open(summary_path, "w") as f:
        f.write("gene_id\tfamily_id\tevalue\n")
        for gene_id, (family_id, evalue) in sorted(hits.items()):
            f.write(f"{gene_id}\t{family_id}\t{evalue:.2e}\n")
    logger.info(f"Rescue summary: {summary_path}")
