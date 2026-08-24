"""HMMER-based rescue step for unplaced genes after iterative pipeline convergence.

After OrthoFinder rounds converge, some genes remain unplaced because DIAMOND
cannot detect distant homologs. This module:
  1. Builds HMM profiles from each confirmed family's protein alignment.
  2. Searches unplaced genes against these profiles with hmmsearch.
  3. Assigns hits (E-value < threshold) to the best-matching family.
  4. Re-aligns and re-builds trees for families that gained new members.
"""

from dataclasses import dataclass
import logging
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from config import Config
from utils.seqio import write_fasta

# Prevent MAFFT_BINARIES version conflict in conda/micromamba environments
os.environ.pop("MAFFT_BINARIES", None)

from steps.hmm_chunks import (DOM_CHUNK_GLOB, merge_tblouts, split_hmm_file,
                              write_runner_script)

logger = logging.getLogger("family_finder")

_COPY_CHUNK_BYTES = 8192   # streaming unit for the profile-database concatenation


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
    """Concatenate the individual HMM profiles into a single HMM database.

    Streamed rather than read whole: the 5sp panel's database is ~1.6 GB.
    Returns None — leaving any existing combined_path untouched — when the
    directory holds no profile to concatenate.
    """
    hmm_files = sorted(hmm_dir.glob("*.hmm"))
    if not hmm_files:
        return None

    with open(combined_path, "wb") as out_f:
        for hf in hmm_files:
            with open(hf, "rb") as inf:
                shutil.copyfileobj(inf, out_f, _COPY_CHUNK_BYTES)

    return combined_path


def _run_hmmsearch(
    hmm_db: Path, query_fasta: Path, outpath: Path, config: Config
) -> Path:
    """Run hmmsearch of unplaced proteins against family HMM profiles."""
    cmd = [
        config.hmmsearch_bin,
        "--tblout", str(outpath),
        "--domtblout", str(outpath.with_suffix(".domtblout")),
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


def _run_hmmsearch_chunked(
    combined_hmm: Path, query_fasta: Path, tblout: Path,
    rescue_dir: Path, config
) -> None:
    """Split the profile DB, run the chunks (SLURM-optional), merge.

    Identical output to a single run — profiles are independent, so only
    the order of tblout blocks differs and the parser keeps the best hit
    per gene regardless. Fails loudly on any incomplete chunk.
    """
    chunk_dir = rescue_dir / "chunks"
    out_dir = rescue_dir / "chunk_out"
    chunks = split_hmm_file(combined_hmm, chunk_dir, config.hmmer_chunk_size)

    runner = write_runner_script(
        rescue_dir / "run_hmmsearch_chunks.sh",
        n_chunks=len(chunks), chunk_dir=chunk_dir, query_fasta=query_fasta,
        out_dir=out_dir, hmmsearch_bin=config.hmmsearch_bin,
        threads=min(config.n_workers, 8), evalue=config.hmmer_evalue,
        max_concurrent=config.hmmer_chunk_concurrent,
        sbatch_extra=config.hmmer_chunk_sbatch_extra,
    )
    logger.info(f"HMMER rescue: {len(chunks)} chunks of "
                f"{config.hmmer_chunk_size} profiles -> {runner}")
    result = subprocess.run(["bash", str(runner)], capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"chunk runner failed: {result.stderr[-500:]}")
        raise RuntimeError(f"chunked hmmsearch failed (exit {result.returncode})")
    merge_tblouts(out_dir, tblout, expected=len(chunks))
    merge_tblouts(
        out_dir, tblout.with_suffix(".domtblout"), expected=len(chunks),
        glob=DOM_CHUNK_GLOB,
    )


@dataclass(frozen=True)
class RescueHit:
    """One gene's rescue assignment, with the evidence that decided it."""
    family: str
    bits: float
    evalue: float
    margin: Optional[float]   # best - second best, None when only one family hit
    grade: str                # HIGH / PROVISIONAL / UNRESOLVED
    reason: str


def _grade_hit(best_bits: float, second_bits: Optional[float],
               config: Config) -> Tuple[str, str]:
    """Grade an assignment from the best-vs-second bit-score margin.

    Same rule the per-round tier applies (steps/profile_assign.py): the margin
    must clear max(profile_margin_bits, profile_margin_frac x best). A fixed
    bit margin alone is meaningless across scales - ten bits separates two weak
    profiles and is noise between two that score ten thousand.
    """
    if second_bits is None:
        return "HIGH", "only one family profile hit this gene"
    margin = best_bits - second_bits
    if margin <= 0:
        return "UNRESOLVED", (
            f"tie: the two best families score the same ({best_bits:.1f} bits). "
            f"Nothing in the evidence prefers either, and emission order must "
            f"not decide."
        )
    required = max(config.profile_margin_bits,
                   config.profile_margin_frac * best_bits)
    if margin >= required:
        return "HIGH", (f"margin {margin:.1f} bits clears the required "
                        f"{required:.1f}")
    return "PROVISIONAL", (
        f"margin {margin:.1f} bits is below the required {required:.1f} — the "
        f"per-round tier would refuse this assignment"
    )


def _parse_rescue_domtblout(
    dom_path: Path, evalue_cutoff: float, config: Config
) -> Dict[str, RescueHit]:
    """Grade rescue assignments from --domtblout, with a coverage gate.

    Issue #45. The rescue read --tblout, which carries full-sequence scores and
    no domain envelopes, so it could not measure coverage at all: a forty-residue
    match to a thousand-residue profile scored, passed the E-value cutoff and
    became family membership. The per-round tier refuses exactly that case on
    `profile_min_coverage` / `profile_min_query_coverage`, and the two tiers
    disagreeing about what counts as a member is the inconsistency this issue
    is about.

    Domain envelopes are merged before measuring, so a profile matched in two
    halves is not penalised for being split. Both directions are checked - a
    short conserved domain inside a long protein covers the profile while
    covering almost none of the query, and that is not membership either.

    Parsing is `steps.profile_assign.parse_domtblout`, deliberately: one parser,
    one definition of coverage across both tiers.
    """
    from steps.profile_assign import parse_domtblout

    per_gene = parse_domtblout(dom_path)
    hits: Dict[str, RescueHit] = {}
    for gene_id, profile_hits in per_gene.items():
        candidates = [h for h in profile_hits if h.full_evalue <= evalue_cutoff]
        if not candidates:
            continue
        candidates.sort(key=lambda h: -h.full_bits)
        best = candidates[0]
        second_bits = candidates[1].full_bits if len(candidates) > 1 else None
        grade, reason = _grade_hit(best.full_bits, second_bits, config)

        if best.profile_cov < config.rescue_min_profile_coverage:
            grade, reason = "UNRESOLVED", (
                f"profile coverage {best.profile_cov:.2f} is below "
                f"{config.rescue_min_profile_coverage:.2f} — the match is too "
                f"partial to call membership"
            )
        elif best.query_cov < config.profile_min_query_coverage:
            grade, reason = "UNRESOLVED", (
                f"query coverage {best.query_cov:.2f} is below "
                f"{config.profile_min_query_coverage:.2f} — the profile matches "
                f"a small part of this protein"
            )

        hits[gene_id] = RescueHit(
            family=best.family_id,
            bits=best.full_bits,
            evalue=best.full_evalue,
            margin=None if second_bits is None else best.full_bits - second_bits,
            grade=grade,
            reason=reason,
        )
    return hits


def _parse_hmmsearch_tblout(
    tblout_path: Path, evalue_cutoff: float, config: Config
) -> Dict[str, RescueHit]:
    """Parse hmmsearch --tblout into one graded assignment per gene.

    **The decision is the bit score, never the E-value** (issue #45). HMMER
    prints E-values at two significant figures and underflows small ones to
    zero, so equal E-values are common: 1,036 of 29,943 rescued genes in the
    five-species run tied at the best E-value, 1,030 of them at zero. The old
    code kept the first such hit, and first meant tblout order, which is
    `sorted(glob("*.hmm"))` order - where `R10_*` precedes `R1_*` because
    '0' < '_'. Bit scores separated 1,035 of those 1,036, and disagreed with
    the E-value winner outright for 24 genes.

    The E-value is still the pre-filter (hmmsearch's own -E) and is still
    reported; it is simply not what picks the family.
    """
    per_gene: Dict[str, List[Tuple[float, str, float]]] = {}

    # HMMER3 --tblout columns (whitespace-delimited):
    # 0: target name, 1: target accession, 2: query name (HMM),
    # 3: query accession, 4: full-seq E-value, 5: full-seq score, ...
    with open(tblout_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 6:
                continue
            gene_id = fields[0]        # target name (unplaced gene)
            family_id = fields[2]      # query name (HMM profile = family_id)
            evalue = float(fields[4])  # full-sequence E-value
            bits = float(fields[5])    # full-sequence score

            if evalue > evalue_cutoff:
                continue
            per_gene.setdefault(gene_id, []).append((bits, family_id, evalue))

    hits: Dict[str, RescueHit] = {}
    for gene_id, rows in per_gene.items():
        # Sort by bits only. Family id is NOT a tiebreaker - letting it decide
        # is the defect this function exists to remove.
        rows.sort(key=lambda r: -r[0])
        best_bits, best_family, best_evalue = rows[0]
        second_bits = rows[1][0] if len(rows) > 1 else None
        grade, reason = _grade_hit(best_bits, second_bits, config)
        hits[gene_id] = RescueHit(
            family=best_family,
            bits=best_bits,
            evalue=best_evalue,
            margin=None if second_bits is None else best_bits - second_bits,
            grade=grade,
            reason=reason,
        )
    return hits


def _collect_hmm_work_items(
    families: Dict[str, Set[str]], outdir: Path, hmm_dir: Path, config: Config
) -> list:
    """Which families still need an hmmbuild run, as _build_hmm_worker args."""
    work_items = []
    for family_id in families:
        aln_path = _find_family_alignment(family_id, outdir)
        if aln_path is None:
            continue
        hmm_path = hmm_dir / f"{family_id}.hmm"
        if not hmm_path.exists():
            work_items.append((family_id, str(aln_path), str(hmm_path), config.hmmbuild_bin))
    return work_items


def _run_hmmbuild_pool(work_items: list, n_workers: int) -> Tuple[int, int]:
    """Run the hmmbuild work items in a process pool -> (n_ok, n_failed)."""
    n_ok = 0
    n_failed = 0
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(_build_hmm_worker, item): item[0]
                   for item in work_items}
        done = 0
        for future in as_completed(futures):
            try:
                if future.result() is not None:
                    n_ok += 1
                else:
                    n_failed += 1
            except Exception as e:
                n_failed += 1
                logger.debug(f"  hmmbuild worker exception for "
                             f"{futures[future]}: {e}")
            done += 1
            if done % 2000 == 0:
                logger.info(f"  hmmbuild progress: {done}/{len(work_items)}")
    return n_ok, n_failed


def _build_profile_db(
    families: Dict[str, Set[str]], outdir: Path, rescue_dir: Path,
    hmm_dir: Path, config: Config,
) -> Optional[Path]:
    """Build every missing HMM profile and concatenate them into one database.

    Returns the database path, or None when the rescue cannot proceed (no
    profile built, nothing to concatenate, or hmmpress failed).
    """
    work_items = _collect_hmm_work_items(families, outdir, hmm_dir, config)
    n_built = sum(1 for f in hmm_dir.glob("*.hmm") if f.stat().st_size > 0)

    logger.info(f"HMMER rescue: building {len(work_items)} HMM profiles "
                f"({n_built} already exist), using {config.n_workers} workers")

    if work_items:
        n_ok, n_failed = _run_hmmbuild_pool(
            work_items, min(config.n_workers, len(work_items)))
        n_built += n_ok
        if n_failed:
            logger.info(f"  hmmbuild: {n_failed} profiles failed")

    logger.info(f"HMMER rescue: {n_built} HMM profiles ready")

    if n_built == 0:
        logger.warning("HMMER rescue: no HMM profiles built, skipping")
        return None

    combined_hmm = rescue_dir / "all_families.hmm"
    if _concat_hmms(hmm_dir, combined_hmm) is None:
        logger.warning("HMMER rescue: no HMM files to concatenate")
        return None

    # hmmpress is only needed by hmmscan, never by hmmsearch; skip it when
    # chunking (it costs minutes and 1.6 GB on the 5sp panel for nothing).
    if not config.hmmer_chunk_size:
        press_cmd = [config.hmmpress_bin, "-f", str(combined_hmm)]
        press_result = subprocess.run(press_cmd, capture_output=True, text=True)
        if press_result.returncode != 0:
            logger.error(f"hmmpress failed: {press_result.stderr[:300]}")
            return None

    return combined_hmm


def _search_unplaced(
    combined_hmm: Path, unplaced_pool: Dict[str, str], rescue_dir: Path,
    config: Config,
) -> Dict[str, Tuple[str, float]]:
    """Search the unplaced pool against the profile database -> best hits."""
    unplaced_fasta = rescue_dir / "unplaced_proteins.fa"
    write_fasta(unplaced_pool, str(unplaced_fasta))

    # One hmmsearch call, or chunked (issue #31)
    tblout = rescue_dir / "hmmsearch_results.tblout"
    if config.hmmer_chunk_size:
        _run_hmmsearch_chunked(combined_hmm, unplaced_fasta, tblout,
                               rescue_dir, config)
    else:
        _run_hmmsearch(combined_hmm, unplaced_fasta, tblout, config)

    # Prefer --domtblout: it is the only output carrying domain envelopes, and
    # without them coverage cannot be checked at all (issue #45). Falling back
    # is allowed for a rescue directory written before this existed, but it is
    # said out loud rather than degrading quietly.
    dom = tblout.with_suffix(".domtblout")
    if dom.exists():
        return _parse_rescue_domtblout(dom, config.hmmer_evalue, config)
    logger.warning(
        "no --domtblout beside %s: grading on bit scores alone, WITHOUT the "
        "coverage gate the per-round tier applies", tblout)
    return _parse_hmmsearch_tblout(tblout, config.hmmer_evalue, config)


def _group_by_family(
    hits: Dict[str, RescueHit]
) -> Dict[str, Set[str]]:
    """Invert gene -> RescueHit into family -> {genes}."""
    rescued_by_family: Dict[str, Set[str]] = {}
    for gene_id, hit in hits.items():
        # UNRESOLVED means the evidence does not prefer this family: a tie in
        # bits, or coverage too thin to call membership. Placing it anyway
        # would let emission order decide after all, which is the whole defect.
        if hit.grade == "UNRESOLVED":
            logger.debug(f"  Not placed {gene_id}: {hit.reason}")
            continue
        rescued_by_family.setdefault(hit.family, set()).add(gene_id)
        logger.debug(f"  Rescued {gene_id} → {hit.family} "
                     f"({hit.bits:.1f} bits, {hit.grade})")
    return rescued_by_family


def _absorb_rescued(
    families: Dict[str, Set[str]],
    rescued_by_family: Dict[str, Set[str]],
    protein_pool: Dict[str, str],
    cds_pool: Dict[str, str],
    outdir: Path,
    config: Config,
) -> Dict[str, Set[str]]:
    """Add rescued genes to their families and re-align the ones that grew.

    Hits against a family id that is not in `families` are dropped here —
    they still appear in the rescue summary.
    """
    updated_families = dict(families)
    for family_id, new_genes in rescued_by_family.items():
        if family_id not in updated_families:
            continue

        updated_families[family_id] = updated_families[family_id] | new_genes
        _realign_family(
            family_id, updated_families[family_id],
            protein_pool, cds_pool, outdir, config,
        )
    return updated_families


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
        families: family_id -> set of confirmed gene_ids.
        unplaced_pool: gene_id -> protein sequence, for unplaced genes only.
        protein_pool, cds_pool: the full pools, used to re-align a grown family.
        outdir, config: pipeline output directory and configuration.

    Returns:
        `families` with rescued genes added — or `families` itself, unchanged,
        whenever the rescue cannot place anything.
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

    combined_hmm = _build_profile_db(families, outdir, rescue_dir, hmm_dir, config)
    if combined_hmm is None:
        return families

    hits = _search_unplaced(combined_hmm, unplaced_pool, rescue_dir, config)
    if not hits:
        logger.info("HMMER rescue: no significant hits found")
        return families

    rescued_by_family = _group_by_family(hits)
    logger.info(f"HMMER rescue: {len(hits)} genes rescued into "
                f"{len(rescued_by_family)} families")

    updated_families = _absorb_rescued(
        families, rescued_by_family, protein_pool, cds_pool, outdir, config,
    )
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

        # Codon alignment. Genes removed by the internal-stop filter remain
        # family members (they are already placed); they are only excluded
        # from the codon alignment.
        codon_aln, _stop_removed = codon_align(
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
        f.write("gene_id\tfamily_id\tbits\tmargin\tgrade\tevalue\treason\n")
        for gene_id, hit in sorted(hits.items()):
            margin = "" if hit.margin is None else f"{hit.margin:.1f}"
            f.write(f"{gene_id}\t{hit.family}\t{hit.bits:.1f}\t{margin}\t"
                    f"{hit.grade}\t{hit.evalue:.2e}\t{hit.reason}\n")
    logger.info(f"Rescue summary: {summary_path}")
