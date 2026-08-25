"""Nominate structurally co-clustered merge candidates from a ProstT5 3Di DB.

This is a standalone structural analysis step for issue #17's planned tier:
it reads a run directory's CUMULATIVE confirmed families (the union of every
`round_*/confirmed_families.tsv` — each file holds only that round's NEWLY
confirmed families, so any single round under-reports the table), runs
`foldseek cluster` against an already-built ProstT5 3Di database, and emits
CANDIDATE family-merge sets.

Design constraints, enforced here:
  * Reuse steps.prostt5_chunks.verify_3di_db() before any foldseek call.
  * Never merge families automatically. Structure nominates only.
  * Keep foldseek invocations in thin module-level wrappers so tests replace
    them and pure cluster-mapping logic runs without foldseek installed.
  * Guard greedy set-cover runaway clusters with a hard size cap: oversized
    clusters are written to a separate overflow file and excluded from the
    candidate table.
"""

import argparse
import logging
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple

from config import Config, resolve_tool
from steps.prostt5_chunks import verify_3di_db

logger = logging.getLogger("family_finder")


def load_confirmed_families(path: Path) -> Dict[str, Set[str]]:
    """Parse one round's confirmed_families.tsv into family -> member ids."""
    families: Dict[str, Set[str]] = {}
    gene_owner: Dict[str, str] = {}
    with open(path) as f:
        for lineno, line in enumerate(f, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 2 or not parts[0] or not parts[1]:
                raise ValueError(f"{path}:{lineno}: malformed confirmed-family row")
            family_id = parts[0]
            members = parts[1].split(",")
            if len(members) != len(set(members)):
                raise ValueError(f"{path}:{lineno}: duplicate gene within {family_id}")
            if family_id in families:
                raise ValueError(f"{path}:{lineno}: duplicate family id {family_id!r}")
            for gene_id in members:
                previous = gene_owner.get(gene_id)
                if previous is not None:
                    raise ValueError(
                        f"{path}:{lineno}: gene {gene_id!r} occurs in both "
                        f"{previous} and {family_id}"
                    )
                gene_owner[gene_id] = family_id
            families[family_id] = set(members)
    return families


def load_all_round_families(run_dir: Path) -> Dict[str, Set[str]]:
    """Union every round's confirmed_families.tsv in a pipeline run directory.

    pipeline._write_round_families writes each round INCREMENTALLY: a round's
    file holds only the families newly confirmed in that round. Reading a
    single round therefore counts every other round's members as unplaced —
    on a multi-round run a foldseek cluster spanning a round-1 and a round-2
    family would score n_families == 1 and a real merge candidate would be
    silently dropped. The cumulative table is the union, with family ids and
    gene ownership checked across rounds.
    """
    run_dir = Path(run_dir)
    tsvs = sorted(run_dir.glob("round_*/confirmed_families.tsv"))
    if not tsvs:
        raise FileNotFoundError(
            f"no round_*/confirmed_families.tsv under {run_dir}; expected a "
            "pipeline output directory"
        )
    families: Dict[str, Set[str]] = {}
    gene_owner: Dict[str, str] = {}
    for tsv in tsvs:
        for family_id, members in load_confirmed_families(tsv).items():
            if family_id in families:
                raise ValueError(
                    f"{tsv}: family id {family_id!r} already confirmed in an "
                    "earlier round"
                )
            for gene_id in members:
                previous = gene_owner.get(gene_id)
                if previous is not None:
                    raise ValueError(
                        f"{tsv}: gene {gene_id!r} occurs in both {previous} "
                        f"and {family_id}"
                    )
                gene_owner[gene_id] = family_id
            families[family_id] = members
    return families


def gene_to_family_index(families: Dict[str, Set[str]]) -> Dict[str, str]:
    """Invert family membership to gene -> family."""
    index: Dict[str, str] = {}
    for family_id, members in families.items():
        for gene_id in members:
            previous = index.get(gene_id)
            if previous is not None and previous != family_id:
                raise ValueError(
                    f"gene {gene_id!r} appears in both {previous!r} and {family_id!r}"
                )
            index[gene_id] = family_id
    return index


def parse_foldseek_cluster_tsv(path: Path) -> List[Tuple[str, str]]:
    """Parse foldseek createtsv cluster output into (representative, member).

    Malformed rows raise instead of being skipped: a truncated createtsv
    (disk full, node preemption mid-write) would otherwise yield silently
    incomplete clusters, and statuses computed from partial membership look
    exactly like a finished run. Extra columns beyond the first two are
    createtsv's own and are ignored.
    """
    pairs: List[Tuple[str, str]] = []
    with open(path) as f:
        for lineno, line in enumerate(f, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0] or not fields[1]:
                raise ValueError(
                    f"{path}:{lineno}: malformed foldseek cluster row: {line!r}"
                )
            pairs.append((fields[0], fields[1]))
    return pairs


def cluster_pairs_to_clusters(pairs: Sequence[Tuple[str, str]]) -> List[Dict[str, object]]:
    """Group rep/member pairs into stable cluster records.

    Each sequence must belong to exactly one representative. Duplicate rows for
    the same rep/member pair are ignored; a member assigned to two different
    representatives is rejected as malformed input.
    """
    by_rep: Dict[str, List[str]] = {}
    member_owner: Dict[str, str] = {}
    rep_order: List[str] = []

    for representative, member in pairs:
        if representative not in by_rep:
            by_rep[representative] = []
            rep_order.append(representative)
        owner = member_owner.get(member)
        if owner is not None and owner != representative:
            raise ValueError(
                f"member {member!r} appears in both {owner!r} and {representative!r}"
            )
        member_owner[member] = representative
        if member not in by_rep[representative]:
            by_rep[representative].append(member)
        if representative not in member_owner:
            member_owner[representative] = representative
            if representative not in by_rep[representative]:
                by_rep[representative].append(representative)

    clusters: List[Dict[str, object]] = []
    for representative in rep_order:
        members = sorted(set(by_rep[representative]))
        clusters.append({"representative": representative, "members": members})
    return clusters


def _species_of(gene_id: str, delimiter: str) -> str:
    return gene_id.split(delimiter, 1)[0]


def _format_join(items: Iterable[str]) -> str:
    values = sorted(set(items))
    return ",".join(values) if values else "-"


def _format_counts(counts: Dict[str, int]) -> str:
    if not counts:
        return "-"
    return ",".join(f"{key}:{counts[key]}" for key in sorted(counts))


def summarize_clusters(
    clusters: Sequence[Dict[str, object]],
    gene_to_family: Dict[str, str],
    max_cluster_size: int,
    species_delimiter: str = "_",
) -> Tuple[List[Dict[str, object]], Dict[str, int]]:
    """Annotate foldseek clusters against family assignments.

    Returns per-cluster records and an aggregate status counter. Status values:
      * candidate: spans >=2 confirmed families and is within the size cap
      * overflow: exceeds the size cap; never nominated
      * single_family_with_unplaced: one family plus unplaced genes
      * single_family: one family only
      * unplaced_only: no confirmed-family members
    """
    records: List[Dict[str, object]] = []
    summary = {
        "total_clusters": 0,
        "candidate_clusters": 0,
        "overflow_clusters": 0,
        "single_family_clusters": 0,
        "single_family_with_unplaced_clusters": 0,
        "unplaced_only_clusters": 0,
        "candidate_genes": 0,
        "candidate_unplaced_genes": 0,
        "max_cluster_size_seen": 0,
    }

    for idx, cluster in enumerate(clusters, start=1):
        representative = str(cluster["representative"])
        members = sorted(set(cluster["members"]))  # type: ignore[arg-type]
        family_counts: Dict[str, int] = {}
        species_counts: Dict[str, int] = {}
        unplaced: List[str] = []

        for gene_id in members:
            species = _species_of(gene_id, species_delimiter)
            species_counts[species] = species_counts.get(species, 0) + 1
            family_id = gene_to_family.get(gene_id)
            if family_id is None:
                unplaced.append(gene_id)
            else:
                family_counts[family_id] = family_counts.get(family_id, 0) + 1

        n_members = len(members)
        n_families = len(family_counts)
        if n_members > max_cluster_size:
            status = "overflow"
        elif n_families >= 2:
            status = "candidate"
        elif n_families == 1 and unplaced:
            status = "single_family_with_unplaced"
        elif n_families == 1:
            status = "single_family"
        else:
            status = "unplaced_only"

        record = {
            "cluster_id": f"FSCLUST_{idx:06d}",
            "representative_gene_id": representative,
            "n_genes": n_members,
            "n_species": len(species_counts),
            "n_families": n_families,
            "family_ids": sorted(family_counts),
            "family_gene_counts": family_counts,
            "unplaced_gene_ids": sorted(unplaced),
            "n_unplaced": len(unplaced),
            "species_counts": species_counts,
            "status": status,
        }
        records.append(record)

        summary["total_clusters"] += 1
        summary["max_cluster_size_seen"] = max(summary["max_cluster_size_seen"], n_members)
        if status == "candidate":
            summary["candidate_clusters"] += 1
            summary["candidate_genes"] += n_members
            summary["candidate_unplaced_genes"] += len(unplaced)
        elif status == "overflow":
            summary["overflow_clusters"] += 1
        elif status == "single_family_with_unplaced":
            summary["single_family_with_unplaced_clusters"] += 1
        elif status == "single_family":
            summary["single_family_clusters"] += 1
        else:
            summary["unplaced_only_clusters"] += 1

    return records, summary


def write_merge_candidates_tsv(records: Sequence[Dict[str, object]], path: Path) -> None:
    """Write candidate clusters only."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(
            "cluster_id\trepresentative_gene_id\tn_genes\tn_species\tn_families\t"
            "family_ids\tfamily_gene_counts\tn_unplaced\tunplaced_gene_ids\t"
            "species_counts\n"
        )
        for rec in records:
            if rec["status"] != "candidate":
                continue
            f.write(
                f"{rec['cluster_id']}\t{rec['representative_gene_id']}\t"
                f"{rec['n_genes']}\t{rec['n_species']}\t{rec['n_families']}\t"
                f"{_format_join(rec['family_ids'])}\t"
                f"{_format_counts(rec['family_gene_counts'])}\t"
                f"{rec['n_unplaced']}\t"
                f"{_format_join(rec['unplaced_gene_ids'])}\t"
                f"{_format_counts(rec['species_counts'])}\n"
            )


def write_overflow_tsv(records: Sequence[Dict[str, object]], path: Path) -> None:
    """Write oversized clusters that were intentionally skipped."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(
            "cluster_id\trepresentative_gene_id\tn_genes\tn_species\tn_families\t"
            "family_ids\tn_unplaced\tunplaced_gene_ids\tspecies_counts\n"
        )
        for rec in records:
            if rec["status"] != "overflow":
                continue
            f.write(
                f"{rec['cluster_id']}\t{rec['representative_gene_id']}\t"
                f"{rec['n_genes']}\t{rec['n_species']}\t{rec['n_families']}\t"
                f"{_format_join(rec['family_ids'])}\t{rec['n_unplaced']}\t"
                f"{_format_join(rec['unplaced_gene_ids'])}\t"
                f"{_format_counts(rec['species_counts'])}\n"
            )


def write_summary_tsv(summary: Dict[str, int], verified_entries: int, path: Path) -> None:
    """Write aggregate counts as a simple key/value TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"verified_3di_entries\t{verified_entries}\n")
        for key in sorted(summary):
            f.write(f"{key}\t{summary[key]}\n")


def run_foldseek_cluster(
    db_path: Path,
    cluster_db: Path,
    tmp_dir: Path,
    config: Config,
) -> Path:
    """Run foldseek cluster on the prebuilt ProstT5 DB."""
    foldseek_bin = resolve_tool(config.foldseek_bin, "foldseek_bin")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        foldseek_bin, "cluster",
        str(db_path), str(cluster_db), str(tmp_dir),
        "--threads", str(max(1, config.n_workers)),
        "-e", str(config.merge_candidate_evalue),
        "-c", str(config.merge_candidate_min_coverage),
        "--min-seq-id", str(config.merge_candidate_min_seq_id),
    ]
    logger.info(f"Running foldseek cluster: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"foldseek cluster failed for {db_path} (exit {result.returncode}):\n"
            f"{result.stderr[-2000:]}"
        )
    return cluster_db


def run_foldseek_createtsv(
    db_path: Path,
    cluster_db: Path,
    out_tsv: Path,
    config: Config,
) -> Path:
    """Convert foldseek cluster DB output to rep/member TSV."""
    foldseek_bin = resolve_tool(config.foldseek_bin, "foldseek_bin")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        foldseek_bin, "createtsv",
        str(db_path), str(db_path), str(cluster_db), str(out_tsv),
    ]
    logger.info(f"Running foldseek createtsv: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"foldseek createtsv failed for {cluster_db} (exit {result.returncode}):\n"
            f"{result.stderr[-2000:]}"
        )
    return out_tsv


def run_merge_candidates(
    run_dir: Path,
    db_path: Path,
    out_dir: Path,
    config: Config,
) -> Dict[str, Path]:
    """Validate the DB, cluster it, cross-reference families, and write outputs.

    ``run_dir`` is the pipeline OUTPUT directory (the parent of round_*), not
    a single round: the confirmed-family table is cumulative across rounds.
    """
    run_dir = Path(run_dir)
    db_path = Path(db_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    verified_entries = verify_3di_db(db_path)
    families = load_all_round_families(run_dir)
    gene_to_family = gene_to_family_index(families)

    cluster_db = out_dir / "foldseek_clusters"
    tmp_dir = out_dir / "foldseek_tmp"
    cluster_tsv = out_dir / "foldseek_clusters.tsv"
    run_foldseek_cluster(db_path, cluster_db, tmp_dir, config)
    run_foldseek_createtsv(db_path, cluster_db, cluster_tsv, config)

    pairs = parse_foldseek_cluster_tsv(cluster_tsv)
    clusters = cluster_pairs_to_clusters(pairs)
    records, summary = summarize_clusters(
        clusters,
        gene_to_family,
        max_cluster_size=config.merge_candidate_max_cluster_size,
        species_delimiter=config.species_delimiter,
    )

    # NOT merge_candidates.tsv: the pipeline's vote-based merge scan writes
    # (and deletes as stale, pipeline.py run_merge_scan) a file of that exact
    # basename with an incompatible schema. Pointing --outdir at the run
    # directory must not clobber it, and by-filename consumers must never
    # get the wrong schema silently.
    candidates_tsv = out_dir / "structural_merge_candidates.tsv"
    overflow_tsv = out_dir / "structural_merge_overflow.tsv"
    summary_tsv = out_dir / "structural_merge_summary.tsv"
    write_merge_candidates_tsv(records, candidates_tsv)
    write_overflow_tsv(records, overflow_tsv)
    write_summary_tsv(summary, verified_entries, summary_tsv)

    logger.info(
        "Structural merge candidates: %d candidate clusters, %d overflow clusters",
        summary["candidate_clusters"],
        summary["overflow_clusters"],
    )
    return {
        "candidates": candidates_tsv,
        "overflow": overflow_tsv,
        "summary": summary_tsv,
        "cluster_pairs": cluster_tsv,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Nominate merge candidates from a ProstT5 3Di foldseek cluster run.",
    )
    parser.add_argument(
        "--run-dir", required=True,
        help="Pipeline output directory containing round_*/confirmed_families.tsv "
             "(the cumulative family table is the union across rounds)",
    )
    parser.add_argument(
        "--db", required=True,
        help="Existing ProstT5 3Di database path (validated with verify_3di_db)",
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Output directory for merge-candidate tables",
    )
    parser.add_argument(
        "--config", default="",
        help="Optional JSON config; defaults to Config() when omitted",
    )
    parser.add_argument(
        "--threads", type=int, default=0,
        help="Override config.n_workers / foldseek --threads for this run",
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable debug logging",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    config = Config.from_json(args.config) if args.config else Config()
    if args.threads:
        config.n_workers = args.threads

    outputs = run_merge_candidates(
        Path(args.run_dir),
        Path(args.db),
        Path(args.outdir),
        config,
    )
    print(outputs["candidates"])
    print(outputs["overflow"])
    print(outputs["summary"])


if __name__ == "__main__":
    main()
