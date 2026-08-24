"""Iterative pipeline orchestrator.

Runs repeated rounds of OrthoFinder → align → tree → prune.
Outliers from each round feed into the next round.
"""

import json
import logging
import shutil
from pathlib import Path
from typing import Dict, Optional, Set, Tuple, List

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
    validate_species_tree,
)
from utils.checkpoint import save_checkpoint, find_latest_checkpoint
from utils.manifest import start_manifest, finish_manifest
from utils.parallel import parallel_map

logger = logging.getLogger("family_finder")


class ClusteringFailed(RuntimeError):
    """Clustering died mid-run. Terminal: the round is named, so is the cause.

    Raised rather than logged because the alternative was a zero exit code, a
    truncated summary.tsv and a manifest reading "completed" - indistinguishable
    from a finished run in a batch log, and trusted as the baseline by the next
    resume.
    """

    def __init__(self, method: str, round_num: int, cause):
        super().__init__(
            f"clustering ({method}) failed in round {round_num}: {cause}")


def process_single_orthogroup(args: tuple):
    """Process a single orthogroup: align → tree → prune.

    Called by parallel_map. Takes a single tuple of all arguments.
    """
    (og_id, gene_ids, protein_pool, cds_pool,
     expected_distances, config, round_dir) = args

    og_dir = round_dir / "orthogroups" / og_id
    og_dir.mkdir(parents=True, exist_ok=True)

    # Skip if too small — return gene IDs only (not sequences) to avoid pickle overhead
    if len(gene_ids) < config.min_orthogroup_size:
        return (og_id, None, set(gene_ids))

    try:
        # 1. Extract sequences
        prot_seqs = {gid: protein_pool[gid] for gid in gene_ids if gid in protein_pool}
        cds_seqs = {gid: cds_pool[gid] for gid in gene_ids if gid in cds_pool}
        missing = set(gene_ids) - set(prot_seqs)
        if missing:
            logger.debug(f"  {og_id}: {len(missing)} gene(s) missing from protein pool")
        # Need both pep and cds for each gene
        common_ids = set(prot_seqs) & set(cds_seqs)
        if len(common_ids) < config.min_orthogroup_size:
            return (og_id, None, set(gene_ids))
        prot_seqs = {gid: prot_seqs[gid] for gid in common_ids}
        cds_seqs = {gid: cds_seqs[gid] for gid in common_ids}

        # 2. Protein alignment (guide) → codon alignment via pal2nal
        prot_aln = align_protein(prot_seqs, og_dir / "proteins.afa", config)
        codon_aln, stop_removed = codon_align(
            og_dir / "proteins.afa", cds_seqs, og_dir / "codon.afa", config
        )

        # 3. Build tree from CDS (codon) alignment; fall back to protein if pal2nal failed
        if codon_aln is not None:
            tree_path = build_tree(codon_aln, og_dir / "tree.nwk", config)
        else:
            logger.debug(f"  {og_id}: pal2nal failed, using protein alignment for tree")
            tree_path = build_tree(prot_aln, og_dir / "tree.nwk", config)

        # 4. Species-aware pruning (using CDS tree distances)
        gene_to_species = {
            gid: get_species(gid, config.species_delimiter) for gid in common_ids
        }
        confirmed, outliers = prune_orthogroup(
            str(tree_path), gene_to_species, expected_distances, config
        )

        # Genes removed by the internal-stop filter are not leaves of the
        # codon tree: pruning never sees them, so they would land in neither
        # confirmed nor outliers. Recycle them explicitly (issue #15).
        outliers = set(outliers) | (set(stop_removed) - set(confirmed))

        # 5. If enough confirmed members survive pruning, emit the family.
        # Uses min_family_size (not min_orthogroup_size): dissolving a whole OG
        # because pruning left it below the processing floor discards genes that
        # PASSED pruning, and makes small-but-real families (e.g. one outgroup
        # gene + 2-3 ingroup orthologs) impossible to form. See issue #11.
        if len(confirmed) >= config.min_family_size:
            confirmed_prots = {gid: protein_pool[gid] for gid in confirmed if gid in protein_pool}
            confirmed_cds = {gid: cds_pool[gid] for gid in confirmed if gid in cds_pool}
            write_fasta(confirmed_prots, str(og_dir / "confirmed_proteins.fa"))
            write_fasta(confirmed_cds, str(og_dir / "confirmed_cds.fa"))

            # Skip re-alignment if nothing was pruned
            if not outliers:
                shutil.copy(og_dir / "proteins.afa", og_dir / "confirmed_proteins.afa")
                if codon_aln is not None:
                    shutil.copy(codon_aln, og_dir / "confirmed_codon.afa")
                shutil.copy(og_dir / "tree.nwk", og_dir / "confirmed_tree.nwk")
            else:
                # Re-align without outliers (protein guide → codon)
                confirmed_prot_aln = align_protein(
                    confirmed_prots, og_dir / "confirmed_proteins.afa", config
                )
                # Genes removed here stay members of the family (they passed
                # pruning); they are only excluded from the codon alignment.
                confirmed_codon, _ = codon_align(
                    og_dir / "confirmed_proteins.afa",
                    confirmed_cds,
                    og_dir / "confirmed_codon.afa",
                    config,
                )

                # Re-build tree (CDS if available, else protein)
                if confirmed_codon is not None:
                    build_tree(
                        og_dir / "confirmed_codon.afa",
                        og_dir / "confirmed_tree.nwk",
                        config,
                    )
                else:
                    build_tree(
                        og_dir / "confirmed_proteins.afa",
                        og_dir / "confirmed_tree.nwk",
                        config,
                    )

            return (og_id, confirmed, outliers)
        else:
            return (og_id, None, set(gene_ids))

    except Exception as e:
        logger.error(f"Failed to process {og_id}: {e}")
        return (og_id, None, set(gene_ids))


def _audit_gene_conservation(
    pool_ids,
    new_families: Dict[str, Set[str]],
    outlier_gene_ids: Set[str],
    delimiter: str = "_",
) -> Set[str]:
    """Round-level gene-conservation audit (issue #15).

    Every gene entering a round must either be placed in a new family or
    recycled as an outlier. Returns the set of leaked genes (neither placed
    nor recycled), logging any leak at ERROR with a per-species breakdown.
    """
    placed_genes: Set[str] = set()
    for genes in new_families.values():
        placed_genes.update(genes)

    leaked = set(pool_ids) - placed_genes - set(outlier_gene_ids)
    if leaked:
        by_species: Dict[str, int] = {}
        for gid in leaked:
            sp = gid.split(delimiter, 1)[0]
            by_species[sp] = by_species.get(sp, 0) + 1
        breakdown = ", ".join(
            f"{sp}={n}" for sp, n in sorted(by_species.items())
        )
        examples = ", ".join(sorted(leaked)[:10])
        logger.error(
            f"Gene-conservation audit: {len(leaked)} gene(s) neither placed "
            f"nor recycled this round ({breakdown}); forcing back into the "
            f"outlier pool. Examples: {examples}"
        )
    return leaked


def _check_id_agreement(
    protein_pool: Dict[str, str],
    cds_pool: Dict[str, str],
    outdir,
    delimiter: str = "_",
    warn_pct: float = 95.0,
) -> Dict[str, tuple]:
    """Per-species pep/CDS ID agreement report (issue #2).

    A pep/CDS ID mismatch silently shrinks orthogroups (genes lacking a CDS
    partner are dropped from every OG). This makes the dropout visible:
    writes outdir/id_agreement.tsv and warns per species when the shared
    fraction falls below warn_pct of the protein IDs.

    Returns dict of species -> (species, n_pep, n_cds, n_shared, pct_shared).
    """
    pep_by_sp: Dict[str, Set[str]] = {}
    cds_by_sp: Dict[str, Set[str]] = {}
    for gid in protein_pool:
        pep_by_sp.setdefault(gid.split(delimiter, 1)[0], set()).add(gid)
    for gid in cds_pool:
        cds_by_sp.setdefault(gid.split(delimiter, 1)[0], set()).add(gid)

    stats: Dict[str, tuple] = {}
    for sp in sorted(set(pep_by_sp) | set(cds_by_sp)):
        pep = pep_by_sp.get(sp, set())
        cds = cds_by_sp.get(sp, set())
        shared = pep & cds
        pct = 100.0 * len(shared) / len(pep) if pep else 0.0
        stats[sp] = (sp, len(pep), len(cds), len(shared), pct)
        if pct < warn_pct:
            logger.warning(
                f"ID agreement: species {sp} has only {len(shared)}/{len(pep)} "
                f"({pct:.1f}%) protein IDs with a matching CDS ID — genes "
                f"without a CDS partner are silently dropped from orthogroups. "
                f"Check pep/CDS ID consistency."
            )

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "id_agreement.tsv", "w") as f:
        f.write("species\tn_pep\tn_cds\tn_shared\tpct_shared\n")
        for sp in sorted(stats):
            _, n_pep, n_cds, n_shared, pct = stats[sp]
            f.write(f"{sp}\t{n_pep}\t{n_cds}\t{n_shared}\t{pct:.1f}\n")
    return stats


def round_yield(n_genes_placed: int, pool_size: int) -> float:
    """Fraction of the round's input pool that the round actually placed.

    Genes, not families: one family can hold forty genes or two, so a family
    count says nothing about how much of the pool moved.
    """
    if pool_size <= 0:
        return 0.0
    return n_genes_placed / pool_size


def _round_placed_gene_ids(
    new_families: Dict[str, Set[str]],
    profile_assigned_gene_ids: Set[str],
) -> Set[str]:
    """Return every gene removed from this round's input pool exactly once.

    Profile assignment may target a family created in an earlier round, so
    those genes are absent from ``new_families``.  Omitting them makes a
    productive round look barren to the convergence gate even though the
    conservation equation correctly removed them from the outlier pool.
    """
    placed = set().union(*new_families.values()) if new_families else set()
    n_family_members = sum(len(genes) for genes in new_families.values())
    if n_family_members != len(placed):
        raise RuntimeError(
            "round placed the same gene in more than one new family"
        )
    placed.update(profile_assigned_gene_ids)
    return placed


def should_stop_iterating(recent_yields: List[float], config: Config) -> bool:
    """Whether the round loop has stopped paying for itself (issue #46).

    Two rules, and the second contains the first. With
    `convergence_min_yield` at its default 0.0 this is exactly the old
    behaviour - stop after `convergence_no_new_families` consecutive rounds
    that placed nothing. Raise the threshold and the loop also ends when a
    round's yield falls below it, which is what the measured tail asks for:
    on the 5-species run rounds 5 and 6 placed 113 genes of 134,175 between
    them, at roughly two hours each on the 15-species panel.
    """
    if not recent_yields:
        return False
    if config.convergence_min_yield > 0:
        return recent_yields[-1] < config.convergence_min_yield
    barren = 0
    for y in reversed(recent_yields):
        if y > 0:
            break
        barren += 1
    return barren >= config.convergence_no_new_families


def should_skip_profile_tier(assignments_per_round: List[int],
                             config: Config) -> bool:
    """Whether per-round profile assignment has earned another round.

    It rebuilds every family profile from scratch each round - the log says
    `0 cached` every time - so a round that assigns nothing is hours of
    hmmbuild and hmmsearch for no placement. Measured on the 15-species run:
    round 4 placed 4 genes, round 5 placed 0 in 2h11m.
    """
    limit = config.profile_assign_off_after_barren
    if limit <= 0:
        return False
    if len(assignments_per_round) < limit:
        return False
    return all(n == 0 for n in assignments_per_round[-limit:])


def _write_round_families(round_dir, new_families: Dict[str, Set[str]]):
    """Write this round's confirmed families incrementally (resume support).

    Format: family_id<TAB>comma-joined gene ids, one family per line.
    """
    round_dir = Path(round_dir)
    round_dir.mkdir(parents=True, exist_ok=True)
    with open(round_dir / "confirmed_families.tsv", "w") as f:
        for family_id in sorted(new_families):
            f.write(f"{family_id}\t{','.join(sorted(new_families[family_id]))}\n")


def _load_round_families(outdir, max_round: int) -> Dict[str, Set[str]]:
    """Rebuild confirmed families from per-round confirmed_families.tsv files.

    Only rounds <= max_round (the newest COMPLETED round) are loaded, so a
    partially processed later round cannot contaminate the resumed state.
    Profile assignments are replayed from each completed round as well. They
    may target a family created in an earlier round, so recording only the
    current round's tier-1 families loses those genes on resume after they have
    already been removed from outlier_pool.fa.
    """
    families: Dict[str, Set[str]] = {}
    gene_owner: Dict[str, str] = {}
    for round_dir in sorted(Path(outdir).glob("round_*")):
        try:
            round_num = int(round_dir.name.split("_", 1)[1])
        except (ValueError, IndexError):
            continue
        if round_num > max_round:
            continue
        tsv = round_dir / "confirmed_families.tsv"
        if not tsv.exists():
            raise ValueError(
                f"completed round {round_num} has no {tsv.name}; resume would "
                "silently omit its families"
            )
        with open(tsv) as f:
            for lineno, line in enumerate(f, start=1):
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 2 or not parts[0] or not parts[1]:
                    raise ValueError(
                        f"{tsv}:{lineno}: malformed confirmed-family row"
                    )
                family_id = parts[0]
                try:
                    family_round = int(family_id.split("_", 1)[0][1:])
                except (ValueError, IndexError) as exc:
                    raise ValueError(
                        f"{tsv}:{lineno}: invalid family id {family_id!r}"
                    ) from exc
                if family_round != round_num:
                    raise ValueError(
                        f"{tsv}:{lineno}: {family_id} belongs to round "
                        f"{family_round}, not round {round_num}"
                    )
                if family_id in families:
                    raise ValueError(
                        f"{tsv}:{lineno}: duplicate family id {family_id!r}"
                    )
                members = parts[1].split(",")
                if len(members) != len(set(members)):
                    raise ValueError(
                        f"{tsv}:{lineno}: duplicate gene within {family_id}"
                    )
                for gene in members:
                    previous = gene_owner.get(gene)
                    if previous is not None:
                        raise ValueError(
                            f"{tsv}:{lineno}: gene {gene!r} occurs in both "
                            f"{previous} and {family_id}"
                        )
                    gene_owner[gene] = family_id
                families[family_id] = set(members)

        assignments = round_dir / "profile_assign" / "profile_assignments.tsv"
        if not assignments.exists():
            continue
        with open(assignments) as f:
            header = f.readline().rstrip("\n").split("\t")
            try:
                gi = header.index("gene_id")
                fi = header.index("family_id")
            except ValueError as exc:
                raise ValueError(
                    f"{assignments}: missing gene_id/family_id header"
                ) from exc
            assignment_seen = set()
            for lineno, line in enumerate(f, start=2):
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(gi, fi):
                    raise ValueError(
                        f"{assignments}:{lineno}: truncated assignment row"
                    )
                gene, family_id = parts[gi], parts[fi]
                if gene in assignment_seen:
                    raise ValueError(
                        f"{assignments}:{lineno}: duplicate assignment for "
                        f"gene {gene!r}"
                    )
                assignment_seen.add(gene)
                if not gene or family_id not in families:
                    raise ValueError(
                        f"{assignments}:{lineno}: assignment names an empty "
                        f"gene or unknown/future family {family_id!r}"
                    )
                previous = gene_owner.get(gene)
                if previous is not None and previous != family_id:
                    raise ValueError(
                        f"{assignments}:{lineno}: gene {gene!r} is already in "
                        f"{previous}, cannot also assign it to {family_id}"
                    )
                families[family_id].add(gene)
                gene_owner[gene] = family_id
    return families


def _load_completed_round_pool(outdir, round_num: int) -> Dict[str, str]:
    """Load the outlier side of a completed round or refuse the resume.

    Advancing to the next round with the original full input when this file is
    absent reprocesses already placed genes and can overwrite the meaning of a
    resumed run.  An empty FASTA is valid; a missing artifact is not.
    """
    pool_fasta = Path(outdir) / f"round_{round_num:02d}" / "outlier_pool.fa"
    if not pool_fasta.is_file():
        raise ValueError(
            f"completed round {round_num} has no {pool_fasta.name}; resume "
            "cannot reconstruct its unplaced gene set"
        )
    from utils.seqio import read_fasta
    return read_fasta(str(pool_fasta))


def partition_excluded_species(
    pool: Dict[str, str], exclude: list, delimiter: str = "_"
) -> tuple:
    """Split a sequence pool into (clustering_pool, excluded_pool) by the
    species prefixes in `exclude` (issue #12: clustering_species_exclude).

    Prefix matching is exact per get_species (so excluding "Cgig" leaves
    "CgigH_*" alone). A listed species with no genes in the pool is a
    likely config typo and is logged as a warning. Inputs are not mutated.
    """
    if not exclude:
        return dict(pool), {}

    exclude_set = set(exclude)
    kept: Dict[str, str] = {}
    excluded: Dict[str, str] = {}
    for gene_id, seq in pool.items():
        if get_species(gene_id, delimiter) in exclude_set:
            excluded[gene_id] = seq
        else:
            kept[gene_id] = seq

    seen = {get_species(g, delimiter) for g in excluded}
    for sp in sorted(exclude_set - seen):
        logger.warning(
            f"clustering_species_exclude lists '{sp}' but no gene in the "
            f"pool has that species prefix — check the config"
        )
    return kept, excluded


def _prepare_clustering_pools(
    current_pool: Dict[str, str],
    input_pool: Dict[str, str],
    exclude: list,
    delimiter: str = "_",
) -> tuple:
    """Partition a fresh or resumed pool without losing excluded species.

    A completed-round ``outlier_pool.fa`` deliberately contains only genes
    that participated in clustering.  Therefore a resume must recover the
    excluded side from the original input proteomes, not try to rediscover it
    in that checkpoint file.
    """
    unknown = sorted(set(current_pool) - set(input_pool))
    if unknown:
        raise ValueError(
            "resume outlier pool contains genes absent from the input "
            f"proteomes: {unknown[:10]}"
        )
    changed = sorted(
        gene for gene, sequence in current_pool.items()
        if input_pool[gene] != sequence
    )
    if changed:
        raise ValueError(
            "resume outlier pool contains sequences that differ from the "
            f"input proteomes: {changed[:10]}"
        )

    clustering_input, excluded_pool = partition_excluded_species(
        input_pool, exclude, delimiter,
    )
    clustering_pool = {
        gene: sequence for gene, sequence in current_pool.items()
        if gene in clustering_input
    }
    return clustering_pool, excluded_pool


def run(
    protein_dir: str,
    cds_dir: str,
    species_tree_path: str,
    outdir: str,
    config: Config,
    resume: bool = False,
    allow_config_change: bool = False,
):
    """Run the iterative gene family construction pipeline.

    Args:
        protein_dir: Directory of per-species protein FASTA files.
        cds_dir: Directory of per-species CDS FASTA files.
        species_tree_path: Path to Newick species tree.
        outdir: Output directory.
        config: Pipeline configuration.
        resume: Whether to resume from checkpoint.
        allow_config_change: Resume even though the configuration differs from
            the one on record. Records the change in the run manifest.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Record what is producing this output directory, and refuse a resume
    # under a configuration other than the one that produced the existing
    # rounds. First thing in the run: a mismatch must cost nothing.
    start_manifest(
        outdir, config, species_tree_path, protein_dir, cds_dir,
        resume=resume, allow_config_change=allow_config_change,
    )

    # Load species tree and compute expected distances
    logger.info(f"Loading species tree from {species_tree_path}")
    species_tree = load_species_tree(species_tree_path)
    expected_distances = compute_pairwise_distances(species_tree)
    logger.info(f"Computed pairwise distances for {len(set(k[0] for k in expected_distances))} species")

    # Load all sequences
    from utils.seqio import build_seq_map
    logger.info(f"Loading protein sequences from {protein_dir}")
    current_pool = build_seq_map(protein_dir)
    input_protein_pool = dict(current_pool)
    input_gene_ids = set(current_pool)
    logger.info(f"Loaded {len(current_pool)} protein sequences")

    # Validate species tree against data species prefixes (issue #14)
    data_species = {get_species(g, config.species_delimiter) for g in current_pool}
    for problem in validate_species_tree(species_tree, data_species):
        logger.warning(f"Species tree validation: {problem}")
    tree_species = {leaf.name for leaf in species_tree.leaves()}
    if not (tree_species & data_species):
        raise ValueError(
            "Species tree leaf names have NO overlap with data species prefixes "
            f"(tree: {sorted(tree_species)}; data: {sorted(data_species)}). "
            "Pruning would be silently disabled for every gene — fix the tree "
            "leaf names or the gene ID prefixes."
        )

    logger.info(f"Loading CDS sequences from {cds_dir}")
    cds_pool = build_seq_map(cds_dir)
    logger.info(f"Loaded {len(cds_pool)} CDS sequences")

    # Startup validation: pep/CDS ID agreement per species (issue #2)
    _check_id_agreement(current_pool, cds_pool, outdir, config.species_delimiter)
    logger.info(f"ID agreement report: {outdir / 'id_agreement.tsv'}")

    # Resume handling
    start_round = 1
    all_confirmed_families: Dict[str, Set[str]] = {}

    if resume:
        cp = find_latest_checkpoint(outdir)
        if cp and cp["status"] == "completed":
            start_round = cp["round_number"] + 1
            # Reload the outlier pool from the completed round
            current_pool = _load_completed_round_pool(
                outdir, cp["round_number"],
            )
            logger.info(
                f"Resuming from round {start_round} with "
                f"{len(current_pool)} sequences"
            )

            # Reload previously confirmed families from the incremental
            # per-round confirmed_families.tsv files (summary.tsv is only
            # written at the very end, so it cannot support mid-run resume)
            all_confirmed_families = _load_round_families(
                outdir, cp["round_number"]
            )
            logger.info(f"Reloaded {len(all_confirmed_families)} confirmed families from previous rounds")

    # Keep excluded species out of tier-1 clustering for every round
    # (issue #12: CgigH). Their genes rejoin the unplaced pool after
    # convergence for profile mapping / HMMER rescue.
    current_pool, clustering_excluded_pool = _prepare_clustering_pools(
        current_pool, input_protein_pool,
        config.clustering_species_exclude, config.species_delimiter,
    )
    if clustering_excluded_pool:
        logger.info(
            f"Excluding {len(clustering_excluded_pool)} genes from clustering "
            f"(species: {sorted(config.clustering_species_exclude)}) — they "
            f"will be offered to family profiles after convergence"
        )

    round_num = start_round - 1
    round_yields: List[float] = []
    profile_assignments: List[int] = []
    yield_log: List[dict] = []

    while round_num < config.max_rounds:
        round_num += 1
        round_dir = outdir / f"round_{round_num:02d}"
        round_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"=== Round {round_num} === ({len(current_pool)} sequences)")
        save_checkpoint(round_dir, round_num, len(current_pool), "started")

        # Step 1: Split pool by species
        input_dir = round_dir / "input"
        split_by_species(current_pool, str(input_dir), config.species_delimiter)

        # Step 2: Cluster the pool (tier 1) + parse orthogroups. Backend
        # selected by config.clustering_method (issue #22): OrthoFinder
        # (default) or the SonicParanoid2 adapter — both yield the same
        # og_id -> gene_ids dict, everything downstream is unchanged.
        of_dir = round_dir / "orthofinder"
        try:
            if config.clustering_method == "sonicparanoid":
                from steps.sonicparanoid import (
                    find_groups_file, parse_groups, run_sonicparanoid,
                )
                sp_dir = run_sonicparanoid(input_dir, of_dir, config)
                orthogroups = parse_groups(find_groups_file(sp_dir))
            else:
                results_dir = run_orthofinder(input_dir, of_dir, config)
                orthogroups = parse_orthogroups(results_dir)
        except Exception as e:
            # Not a break. A run whose clustering died must not reach
            # _write_final_output and finish "completed": that partial
            # summary.tsv is exactly what the next resume trusts as its
            # baseline, so swallowing this here re-opens #48 through the
            # front door (#49).
            raise ClusteringFailed(config.clustering_method, round_num, e) from e
        logger.info(f"Round {round_num}: {len(orthogroups)} orthogroups found")

        save_checkpoint(round_dir, round_num, len(current_pool), "processing")

        # Step 4: Process each orthogroup
        # Only pass the sequences each OG actually needs (avoids pickling entire pool)
        work_items = []
        for og_id, gene_ids in orthogroups.items():
            og_prots = {gid: current_pool[gid] for gid in gene_ids if gid in current_pool}
            og_cds = {gid: cds_pool[gid] for gid in gene_ids if gid in cds_pool}
            work_items.append(
                (og_id, gene_ids, og_prots, og_cds,
                 expected_distances, config, round_dir)
            )

        results = parallel_map(
            process_single_orthogroup,
            work_items,
            n_workers=config.n_workers,
        )

        # Step 5: Collect confirmed families and outlier gene IDs
        new_families = {}
        outlier_gene_ids = set()

        for r in results:
            if r is None:
                continue
            og_id, confirmed, outlier_ids = r

            if confirmed is not None:
                family_id = f"R{round_num}_{og_id}"
                # Defensive uniqueness check (issue #4): family ids are
                # R{round}_{og_id}; a duplicate means the same round number
                # ran twice (resume/round-numbering bug) and would silently
                # overwrite a previously confirmed family.
                if family_id in all_confirmed_families:
                    logger.error(
                        f"Duplicate family id {family_id}: overwriting a "
                        f"previously confirmed family — this indicates a "
                        f"round-numbering or resume bug"
                    )
                new_families[family_id] = confirmed
                all_confirmed_families[family_id] = confirmed

            outlier_gene_ids.update(outlier_ids)

        # Add unassigned genes (not in any orthogroup) to outlier pool
        all_og_genes = set()
        for gene_ids in orthogroups.values():
            all_og_genes.update(gene_ids)
        for gid in current_pool:
            if gid not in all_og_genes:
                outlier_gene_ids.add(gid)

        # Gene-conservation audit: every gene in this round's pool must be
        # placed or recycled; force any leak back into the pool (issue #15)
        leaked = _audit_gene_conservation(
            current_pool, new_families, outlier_gene_ids,
            config.species_delimiter,
        )
        outlier_gene_ids.update(leaked)

        # Build outlier pool with sequences from current_pool
        new_outlier_pool = {gid: current_pool[gid] for gid in outlier_gene_ids if gid in current_pool}

        # Per-round profile assignment (issue #13): offer this round's outliers
        # to the confirmed families' HMM profiles BEFORE recycling them into
        # the next outliers-only OrthoFinder run. Errors must never kill a round.
        skip_profiles = should_skip_profile_tier(profile_assignments, config)
        if skip_profiles:
            logger.info(
                f"Round {round_num}: skipping profile assignment — the last "
                f"{config.profile_assign_off_after_barren} rounds assigned "
                f"nothing and it rebuilds every profile each round"
            )
        profile_assigned_gene_ids: Set[str] = set()
        if (config.profile_assign_per_round and all_confirmed_families
                and new_outlier_pool and not skip_profiles):
            try:
                from steps.profile_assign import run_profile_assignment
                from steps.hmmer_rescue import _find_family_alignment

                assigned, _ = run_profile_assignment(
                    families=all_confirmed_families,
                    unplaced_pool=new_outlier_pool,
                    protein_lookup=lambda fam_id: _find_family_alignment(fam_id, outdir),
                    outdir=round_dir / "profile_assign",
                    config=config,
                )
                unknown_targets = sorted(set(assigned) - set(all_confirmed_families))
                if unknown_targets:
                    raise ValueError(
                        "profile assignment returned unknown families: "
                        + ",".join(unknown_targets[:10])
                    )
                assignment_owner = {}
                for fam_id, genes in assigned.items():
                    for gid in genes:
                        if gid not in new_outlier_pool:
                            raise ValueError(
                                f"profile assignment returned {gid!r}, which "
                                "was not in this round's outlier pool"
                            )
                        previous = assignment_owner.get(gid)
                        if previous is not None and previous != fam_id:
                            raise ValueError(
                                f"profile assignment placed {gid!r} in both "
                                f"{previous} and {fam_id}"
                            )
                        assignment_owner[gid] = fam_id
                n_assigned = 0
                for fam_id, genes in assigned.items():
                    merged = all_confirmed_families[fam_id] | genes
                    all_confirmed_families[fam_id] = merged
                    if fam_id in new_families:
                        new_families[fam_id] = merged
                    n_assigned += len(genes)
                    for gid in genes:
                        new_outlier_pool.pop(gid, None)
                        profile_assigned_gene_ids.add(gid)
                profile_assignments.append(n_assigned)
                logger.info(
                    f"Round {round_num}: profile assignment placed {n_assigned} "
                    f"outlier genes into {len(assigned)} families"
                )
            except Exception as e:
                logger.error(
                    f"Round {round_num}: profile assignment failed "
                    f"(continuing without it): {e}"
                )

        try:
            placed_this_round = _round_placed_gene_ids(
                new_families, profile_assigned_gene_ids,
            )
        except RuntimeError as exc:
            raise RuntimeError(f"round {round_num}: {exc}") from exc
        overlap = placed_this_round & set(new_outlier_pool)
        accounted = placed_this_round | set(new_outlier_pool)
        if overlap or accounted != set(current_pool):
            missing = sorted(set(current_pool) - accounted)
            extra = sorted(accounted - set(current_pool))
            raise RuntimeError(
                f"round {round_num} post-profile conservation failed: "
                f"placed={len(placed_this_round)} + "
                f"unplaced={len(new_outlier_pool)} != total={len(current_pool)}; "
                f"overlap={overlap and sorted(overlap)[:5]}, "
                f"missing={missing[:5]}, extra={extra[:5]}"
            )

        # Persist the two sides of the completed-round conservation equation
        # only after profile assignment has moved genes between them. The
        # profile_assignments.tsv is replayed on resume for moves into families
        # created in an earlier round; current-round targets are also present in
        # this self-contained family snapshot.
        _write_round_families(round_dir, new_families)

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

        # Step 7: Check convergence — on genes placed, not families created
        placed = len(placed_this_round)
        this_yield = round_yield(placed, len(current_pool))
        round_yields.append(this_yield)
        yield_log.append({
            "round": round_num,
            "pool_size": len(current_pool),
            "new_families": len(new_families),
            "genes_placed": placed,
            "yield": round(this_yield, 6),
            "profile_assigned": profile_assignments[-1] if profile_assignments else None,
        })
        with open(outdir / "round_yield.tsv", "w") as f:
            f.write("round\tpool_size\tnew_families\tgenes_placed\tyield"
                    "\tprofile_assigned\n")
            for r in yield_log:
                f.write(f"{r['round']}\t{r['pool_size']}\t{r['new_families']}\t"
                        f"{r['genes_placed']}\t{r['yield']}\t"
                        f"{'' if r['profile_assigned'] is None else r['profile_assigned']}\n")
        logger.info(
            f"Round {round_num} yield: {placed} gene(s) of {len(current_pool)} "
            f"= {this_yield:.4%}"
        )

        if should_stop_iterating(round_yields, config):
            # Update before leaving, like the other two convergence branches.
            # Breaking first handed HMMER rescue the PREVIOUS round's pool
            # (#51), so genes placed in the final round were offered again and
            # genes newly orphaned by it were never offered at all.
            current_pool = new_outlier_pool
            logger.info(
                f"Converged after round {round_num}: yield {this_yield:.4%} "
                f"(threshold {config.convergence_min_yield:.4%})"
                if config.convergence_min_yield > 0 else
                f"Converged after round {round_num}"
            )
            break

        if len(new_outlier_pool) < config.convergence_threshold:
            logger.info(f"Converged: outlier pool ({len(new_outlier_pool)}) below threshold ({config.convergence_threshold})")
            current_pool = new_outlier_pool
            break

        current_pool = new_outlier_pool

    # Excluded species rejoin the unplaced pool here so profile mapping /
    # HMMER rescue and the final unplaced outputs see them (issue #12).
    if clustering_excluded_pool:
        logger.info(
            f"Merging {len(clustering_excluded_pool)} clustering-excluded "
            f"genes into the unplaced pool for post-run assignment"
        )
        current_pool = {**current_pool, **clustering_excluded_pool}

    # Rescue and the final fragmentation scan both need sequences that have
    # already left the iterative pool.  Load that full map once for either
    # consumer; merge_scan is valid without enabling rescue.
    full_protein_pool = None
    if all_confirmed_families and (
        (config.hmmer_rescue and current_pool) or config.merge_scan
    ):
        full_protein_pool = build_seq_map(protein_dir)

    # HMMER rescue: assign unplaced genes to existing families via HMM profiles
    if config.hmmer_rescue and all_confirmed_families and current_pool:
        logger.info(f"=== HMMER Rescue === ({len(current_pool)} unplaced genes)")
        from steps.hmmer_rescue import rescue_unplaced

        all_confirmed_families = rescue_unplaced(
            families=all_confirmed_families,
            unplaced_pool=current_pool,
            protein_pool=full_protein_pool,
            cds_pool=cds_pool,
            outdir=outdir,
            config=config,
        )

    # Final: Assemble all confirmed families (must run before pseudogene detection
    # and merge scanning so the final, post-rescue alignments are available).
    # Writing this first preserves the completed base family result if the
    # optional but configured scan later fails; the exception still propagates
    # so the manifest cannot call that partial run completed.
    _write_final_output(
        all_confirmed_families, current_pool, cds_pool, outdir, config,
        expected_gene_ids=input_gene_ids,
    )

    # Post-convergence: nominate families that look like one family split in two.
    # A configured scan is part of the run, not a best-effort annotation: a
    # failure must not be converted into a completed manifest with empty files.
    if config.merge_scan and all_confirmed_families:
        scan_for_fragmented_families(
            all_confirmed_families, full_protein_pool, outdir, config,
        )

    # Pseudogene detection (post-convergence, after final output is written)
    if config.pseudogene_detection:
        logger.info("=== Pseudogene Detection ===")
        from steps.pseudogene import (
            detect_pseudogenes,
            write_pseudogene_report,
            write_pseudogene_summary,
            write_pseudogene_fasta,
            write_family_pseudogene_report,
            write_species_comparison,
            write_pseudogene_bed,
            write_chromosomal_distribution,
        )

        # Rebuild full protein pool for analysis
        from utils.seqio import build_seq_map as _build_seq_map
        full_prot = _build_seq_map(protein_dir)

        sp_filter = config.pseudogene_species_filter or None

        n_analyzed = len(full_prot)
        if sp_filter:
            n_analyzed = sum(
                1 for gid in full_prot
                if gid.split(config.species_delimiter, 1)[0] == sp_filter
            )

        evidence = detect_pseudogenes(
            protein_seqs=full_prot,
            cds_seqs=cds_pool,
            families=all_confirmed_families,
            outdir=outdir,
            config=config,
            species_filter=sp_filter,
        )

        pseudo_dir = outdir / "pseudogene_analysis"
        write_pseudogene_report(evidence, pseudo_dir / "pseudogene_candidates.tsv", sp_filter)
        write_pseudogene_summary(evidence, pseudo_dir / "pseudogene_summary.txt", n_analyzed, sp_filter)
        write_pseudogene_fasta(evidence, full_prot, cds_pool, pseudo_dir)
        write_family_pseudogene_report(evidence, all_confirmed_families, pseudo_dir / "family_pseudogene_enrichment.tsv")
        if not sp_filter:
            write_species_comparison(evidence, full_prot, all_confirmed_families, pseudo_dir / "species_comparison.tsv", config.species_delimiter)
        write_pseudogene_bed(evidence, pseudo_dir / "pseudogene_candidates.bed")
        write_chromosomal_distribution(evidence, full_prot, pseudo_dir / "chromosomal_distribution.tsv", sp_filter, config.species_delimiter)

    finish_manifest(outdir, "completed")
    logger.info(f"Pipeline complete: {len(all_confirmed_families)} total families across {round_num} rounds")


def write_merge_candidates(candidates, families: Dict[str, Set[str]],
                          outdir: Path) -> Path:
    """Record family pairs that look like one family split in two.

    Nominations only. `detect_merge_candidates` says in its own docstring that
    a candidate needs tree validation - align A and B together, build the codon
    tree, score every leaf - before anything is merged, and merging two
    families a codon tree would keep apart is worse than leaving them reported.

    The file is written even when the list is empty, because a missing file
    cannot be told apart from a scan that never ran, and "no fragmentation
    detected" is itself a result.

    Family sizes travel with each pair: 63 members against a 9-member splinter
    is a different claim from two 40-member families, and the reciprocal
    fraction alone hides which one you are looking at.
    """
    path = outdir / "merge_candidates.tsv"
    with open(path, "w") as f:
        f.write("family_a\tfamily_b\tmin_reciprocal\tn_genes_a\tn_genes_b\n")
        for fam_a, fam_b, frac in candidates:
            f.write(f"{fam_a}\t{fam_b}\t{frac:.3f}\t"
                    f"{len(families.get(fam_a, ()))}\t"
                    f"{len(families.get(fam_b, ()))}\n")
    return path


def write_vote_edges(edges, outdir: Path) -> Path:
    """Record every one-way nearest-neighbour edge, uncut.

    The v3 family table was built from a file with exactly this schema and
    nothing in this repository wrote it, so the merge of 3,611 families could
    not be reproduced. The edges are a NOMINATION layer: the tree decides
    (steps/cluster_validate.fragment_verdict), and an edge at frac 0.22 is a
    question for the tree rather than an answer.

    Written even when empty - a missing file cannot be told apart from a scan
    that never ran.
    """
    path = outdir / "vote_edges.tsv"
    with open(path, "w") as f:
        f.write("from\tto\tvotes\tfrom_size\tfrac\n")
        for fam_from, fam_to, votes, from_size, frac in edges:
            f.write(f"{fam_from}\t{fam_to}\t{votes}\t{from_size}\t{frac:.3f}\n")
    return path


def write_fragmentation_clusters(clusters, families: Dict[str, Set[str]],
                                 outdir: Path) -> Path:
    """Name the connected components and count what they hold.

    Gene counts travel with each cluster because a component of two 40-member
    families is a different claim from a 63-member family and a 9-member
    splinter, and the family count alone hides which one you are looking at.

    Ids are positional (C0001, C0002, ...) over the ordering
    `fragmentation_clusters` returns - largest first, members sorted - so the
    same edges always produce the same ids.
    """
    path = outdir / "fragmentation_clusters.tsv"
    with open(path, "w") as f:
        f.write("cluster_id\tn_families\tn_genes\tfamilies\n")
        for i, members in enumerate(clusters, start=1):
            n_genes = sum(len(families.get(fam, ())) for fam in members)
            f.write(f"C{i:04d}\t{len(members)}\t{n_genes}\t"
                    f"{','.join(members)}\n")
    return path


def scan_for_fragmented_families(families: Dict[str, Set[str]],
                                 protein_pool: Dict[str, str],
                                 outdir: Path, config: Config) -> Path:
    """Search the final families' own members against the family profiles.

    Run once, after convergence and rescue. Not per round: the per-round tier
    switches itself off once it stops placing (issue #46), which would take the
    fragmentation scan with it, and the family set it would be scanning is not
    the one that ships.

    The base family table is written before this function runs, but a configured
    scan is still part of the run: any missing sequence/profile or inflated vote
    raises.  Empty output files are reserved for a genuinely complete scan with
    no cross-family hits.
    """
    from steps.profile_assign import (build_profiles, detect_merge_candidates,
                                      fragmentation_clusters, parse_domtblout,
                                      vote_edges, _hmmsearch_domtblout)
    from steps.hmmer_rescue import _find_family_alignment

    scan_dir = outdir / "merge_scan"
    scan_dir.mkdir(parents=True, exist_ok=True)

    # A retry in the same output directory must not leave an older successful
    # scan looking current when this attempt fails before rewriting it.
    for name in ("vote_edges.tsv", "fragmentation_clusters.tsv",
                 "merge_candidates.tsv"):
        stale = outdir / name
        if stale.exists():
            stale.unlink()

    member_owner = {}
    for family_id, genes in families.items():
        for gene_id in genes:
            previous = member_owner.get(gene_id)
            if previous is not None and previous != family_id:
                raise ValueError(
                    f"merge scan: gene {gene_id!r} belongs to both "
                    f"{previous!r} and {family_id!r}"
                )
            member_owner[gene_id] = family_id
    missing_sequences = sorted(set(member_owner) - set(protein_pool))
    if missing_sequences:
        raise ValueError(
            "merge scan: family members absent from the protein pool: "
            + ",".join(missing_sequences[:10])
        )
    if not member_owner:
        raise ValueError("merge scan: no family members to search")

    members = {gene_id: protein_pool[gene_id] for gene_id in member_owner}

    query = scan_dir / "family_members.fa"
    write_fasta(members, str(query))

    # Build the database here.  The previous wiring merely checked whether
    # this path happened to exist, so every clean merge_scan run returned
    # header-only "no fragmentation" files without ever running hmmsearch.
    alignments = {}
    for family_id in families:
        # Rescued families have a current alignment under hmmer_rescue; the
        # round's confirmed alignment predates those added members.
        rescued = (outdir / "hmmer_rescue" / "families" / family_id /
                   "proteins.afa")
        alignment = rescued if rescued.exists() else _find_family_alignment(
            family_id, outdir,
        )
        if alignment is not None:
            alignments[family_id] = alignment
    missing_alignments = sorted(set(families) - set(alignments))
    if missing_alignments:
        raise RuntimeError(
            "merge scan: no protein alignment for families: "
            + ",".join(missing_alignments[:10])
        )

    hmm_db = build_profiles(
        families=families,
        alignment_lookup=alignments.get,
        hmm_dir=scan_dir / "hmm_profiles",
        config=config,
    )
    if hmm_db is None:
        raise RuntimeError("merge scan: failed to build the family profile database")

    # build_profiles can continue after individual hmmbuild failures.  That is
    # useful for assignment, but a whole-table audit may not silently omit a
    # family profile.  Verify the exact NAME set in the combined ASCII HMM DB.
    profile_names = []
    with open(hmm_db) as fh:
        for line in fh:
            if line.startswith("NAME"):
                parts = line.split(maxsplit=1)
                if len(parts) == 2:
                    profile_names.append(parts[1].strip())
    if (len(profile_names) != len(set(profile_names))
            or set(profile_names) != set(families)):
        missing = sorted(set(families) - set(profile_names))
        extra = sorted(set(profile_names) - set(families))
        raise RuntimeError(
            "merge scan: profile database does not exactly cover the family "
            f"table (profiles={len(profile_names)}, families={len(families)}, "
            f"missing={missing[:10]}, extra={extra[:10]})"
        )

    dom = scan_dir / "members_vs_families.domtblout"
    _hmmsearch_domtblout(hmm_db, query, dom, config)
    hits = parse_domtblout(dom)

    # The uncut edge layer first: it is what the tree validation consumes, and
    # the reciprocity cut below throws away the cases most worth judging - in
    # the 15-species run the cut left 8,902 of 23,744 families with no edge at
    # all, the flagship PEPC splinter among them, at frac 0.22.
    edges = vote_edges(hits, families, config)
    # One gene casts one vote, so total votes cannot exceed the members that
    # were searched. That is the invariant - NOT "edges <= families", which is
    # false: a family's members need not agree on a neighbour, and the 15sp
    # sample produced 353 edges from 200 families. domtblout is ordered by
    # profile, not by gene, so anything that assumes gene-contiguous rows
    # inflates the counts into frac values above 1.0.
    n_members = sum(len(m) for m in families.values())
    total_votes = sum(v for _, _, v, _, _ in edges)
    inflated = [e for e in edges if e[2] > e[3]]
    if total_votes > n_members or inflated:
        raise RuntimeError(
            f"merge scan: {total_votes} votes from {n_members} members, "
            f"{len(inflated)} edge(s) with votes above the source family's "
            "size - vote inflation; refusing to write"
        )
    clusters = fragmentation_clusters(edges)
    candidates = detect_merge_candidates(hits, families, config)
    write_vote_edges(edges, outdir)
    write_fragmentation_clusters(clusters, families, outdir)
    result = write_merge_candidates(candidates, families, outdir)
    logger.info(
        f"Merge scan: {len(edges)} vote edge(s) over {len(families)} families "
        f"(profile_cov >= {config.merge_min_profile_cov}) forming "
        f"{len(clusters)} cluster(s); {len(candidates)} pair(s) also clear "
        f"reciprocal >= {config.merge_min_reciprocal}"
    )
    return result


def _write_final_output(
    families: Dict[str, Set[str]],
    remaining_pool: Dict[str, str],
    cds_pool: Dict[str, str],
    outdir: Path,
    config: Config,
    expected_gene_ids: Optional[Set[str]] = None,
):
    """Write final summary and copy confirmed family files to final_families/."""
    placed_list = [gene for genes in families.values() for gene in genes]
    placed_genes = set(placed_list)
    if len(placed_list) != len(placed_genes):
        raise RuntimeError(
            "final conservation failed: at least one gene occurs in multiple "
            "families"
        )
    unplaced_ids = set(remaining_pool) - placed_genes
    if expected_gene_ids is not None:
        accounted = placed_genes | unplaced_ids
        overlap = placed_genes & unplaced_ids
        if overlap or accounted != expected_gene_ids:
            missing = sorted(expected_gene_ids - accounted)
            extra = sorted(accounted - expected_gene_ids)
            raise RuntimeError(
                "final conservation failed: "
                f"placed={len(placed_genes)} + unplaced={len(unplaced_ids)} "
                f"!= total={len(expected_gene_ids)}; overlap={len(overlap)}, "
                f"missing={missing[:10]}, extra={extra[:10]}"
            )

    final_dir = outdir / "final_families"
    final_dir.mkdir(parents=True, exist_ok=True)

    # Copy confirmed family outputs
    for family_id, gene_ids in families.items():
        # Parse round and OG from family_id (e.g., "R1_OG0000000")
        try:
            parts = family_id.split("_", 1)
            round_num = int(parts[0][1:])
            og_id = parts[1]
        except (ValueError, IndexError):
            continue
        src_dir = outdir / f"round_{round_num:02d}" / "orthogroups" / og_id
        dst_dir = final_dir / family_id
        if src_dir.exists():
            shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)

        # Overlay HMMER rescue re-alignments if they exist
        rescue_dir = outdir / "hmmer_rescue" / "families" / family_id
        if rescue_dir.exists():
            shutil.copytree(rescue_dir, dst_dir, dirs_exist_ok=True)

    # Write summary TSV
    with open(outdir / "summary.tsv", "w") as f:
        f.write("family_id\tround\tn_genes\tn_species\tgene_list\n")
        for family_id, gene_ids in sorted(families.items()):
            try:
                parts = family_id.split("_", 1)
                round_num = parts[0][1:]
            except (ValueError, IndexError):
                round_num = "?"
            species = {
                get_species(gid, config.species_delimiter) for gid in gene_ids
            }
            gene_list = ",".join(sorted(gene_ids))
            f.write(f"{family_id}\t{round_num}\t{len(gene_ids)}\t{len(species)}\t{gene_list}\n")

    # Write on-disk record of what the pipeline failed to place (issue #15).
    # remaining_pool may still contain genes later placed by HMMER rescue,
    # so subtract everything that ended up in a family.
    unplaced_prot = {
        gid: seq for gid, seq in remaining_pool.items() if gid not in placed_genes
    }
    unplaced_cds = {
        gid: cds_pool[gid] for gid in unplaced_prot if gid in cds_pool
    }
    write_fasta(unplaced_prot, str(outdir / "unplaced_proteins.fa"))
    write_fasta(unplaced_cds, str(outdir / "unplaced_cds.fa"))
    logger.info(
        f"Unplaced genes: {len(unplaced_prot)} proteins "
        f"({len(unplaced_cds)} with CDS) written to "
        f"{outdir / 'unplaced_proteins.fa'}"
    )

    logger.info(f"Final output written to {final_dir}")
    logger.info(f"Summary: {outdir / 'summary.tsv'}")
