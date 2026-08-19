"""First-class HMM profile assignment for the iterative pipeline (issue #13).

Promotes profile search from a one-shot post-convergence rescue
(steps/hmmer_rescue.py, kept as the legacy path) to a per-round step:
after each round's pruning, unplaced/outlier genes are offered to every
confirmed family's HMM profile BEFORE they are recycled into the next
outliers-only OrthoFinder run.

Design rules (from the issue #8 forensics):
  * Decisions use --domtblout, never --tblout.
  * Assignment requires profile coverage, query coverage, and a
    best-vs-second-best BIT SCORE margin.
  * E-values are used only as a pass/fail significance gate and are
    NEVER compared between hits. tblout E-values are printed at 2
    significant figures and underflow to 0, so ties are routine, and the
    legacy tie-break (first line of a lexicographically sorted HMM db,
    where R10_* sorts before R1_*) silently favored later rounds.
"""

import hashlib
import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

from config import Config
from utils.seqio import write_fasta

logger = logging.getLogger("family_finder")

# Minimum number of whitespace-delimited columns needed from a domtblout
# row (through column 18, "ali to"). Full rows have 23 columns.
_DOMTBL_MIN_FIELDS = 19

# Profile-coverage floors for the reassignment screen and merge detection.
REASSIGN_MIN_PROFILE_COV = 0.6
MERGE_MIN_PROFILE_COV = 0.5

# Maximum candidates reported per ambiguous gene.
MAX_AMBIGUOUS_CANDIDATES = 5


@dataclass(frozen=True)
class ProfileHit:
    """One gene-vs-family profile hit aggregated over its domains."""

    gene_id: str
    family_id: str
    full_evalue: float
    full_bits: float
    # Per-domain rows: (hmm_from, hmm_to, ali_from, ali_to, dom_bits)
    domains: List[Tuple[int, int, int, int, float]] = field(default_factory=list)
    hmm_len: int = 0
    seq_len: int = 0

    @property
    def profile_cov(self) -> float:
        """Fraction of the HMM covered by merged, non-overlapping domains."""
        if self.hmm_len <= 0:
            return 0.0
        covered = _merged_length([(d[0], d[1]) for d in self.domains])
        return covered / self.hmm_len

    @property
    def query_cov(self) -> float:
        """Fraction of the query sequence covered by merged domains."""
        if self.seq_len <= 0:
            return 0.0
        covered = _merged_length([(d[2], d[3]) for d in self.domains])
        return covered / self.seq_len

    @property
    def nbits(self) -> float:
        """Full-sequence bit score normalized by profile length."""
        if self.hmm_len <= 0:
            return 0.0
        return self.full_bits / self.hmm_len


def _merged_length(intervals: List[Tuple[int, int]]) -> int:
    """Total length of the union of 1-based inclusive intervals."""
    if not intervals:
        return 0
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    total = 0
    cur_start, cur_end = ordered[0]
    for start, end in ordered[1:]:
        if start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            total += cur_end - cur_start + 1
            cur_start, cur_end = start, end
    total += cur_end - cur_start + 1
    return total


def parse_domtblout(path: Path) -> Dict[str, List[ProfileHit]]:
    """Parse hmmsearch --domtblout output into per-gene ProfileHit lists.

    Profiles are the QUERIES, sequences are the TARGETS (same convention
    as hmmer_rescue). Per-domain rows are aggregated per (gene, family).
    Columns used: target name 0, tlen 2, query name 3, qlen 5,
    full E-value 6, full bits 7, dom bits 13, hmm from 15, hmm to 16,
    ali from 17, ali to 18.

    Returns dict of gene_id -> hits sorted by full_bits DESCENDING.
    E-values are stored for the significance gate but never used to rank.
    """
    accum: Dict[Tuple[str, str], dict] = {}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < _DOMTBL_MIN_FIELDS:
                continue
            try:
                gene_id = fields[0]
                seq_len = int(fields[2])
                family_id = fields[3]
                hmm_len = int(fields[5])
                full_evalue = float(fields[6])
                full_bits = float(fields[7])
                dom_bits = float(fields[13])
                hmm_from = int(fields[15])
                hmm_to = int(fields[16])
                ali_from = int(fields[17])
                ali_to = int(fields[18])
            except ValueError:
                logger.debug(f"Skipping malformed domtblout line: {line[:80]!r}")
                continue

            key = (gene_id, family_id)
            rec = accum.setdefault(key, {
                "full_evalue": full_evalue,
                "full_bits": full_bits,
                "hmm_len": hmm_len,
                "seq_len": seq_len,
                "domains": [],
            })
            rec["domains"].append((hmm_from, hmm_to, ali_from, ali_to, dom_bits))

    hits_by_gene: Dict[str, List[ProfileHit]] = {}
    for (gene_id, family_id), rec in accum.items():
        hit = ProfileHit(
            gene_id=gene_id,
            family_id=family_id,
            full_evalue=rec["full_evalue"],
            full_bits=rec["full_bits"],
            domains=rec["domains"],
            hmm_len=rec["hmm_len"],
            seq_len=rec["seq_len"],
        )
        hits_by_gene.setdefault(gene_id, []).append(hit)

    for hits in hits_by_gene.values():
        hits.sort(key=lambda h: -h.full_bits)
    return hits_by_gene


def _member_digest(members: Set[str]) -> str:
    """sha1 of the sorted member ids — cache key for a family's profile."""
    joined = "\n".join(sorted(members))
    return hashlib.sha1(joined.encode("utf-8")).hexdigest()


def _hmmbuild(
    family_id: str, alignment_path: Path, hmm_path: Path, config: Config
) -> bool:
    """Build one HMM profile. Module-level so tests can mock it."""
    cmd = [
        config.hmmbuild_bin, "--amino", "-n", family_id,
        str(hmm_path), str(alignment_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.debug(f"hmmbuild failed for {family_id}: {result.stderr[:200]}")
        return False
    return hmm_path.exists() and hmm_path.stat().st_size > 0


def _hmmpress(hmm_db: Path, config: Config) -> bool:
    """Press the combined HMM database. Module-level so tests can mock it."""
    cmd = [config.hmmpress_bin, "-f", str(hmm_db)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"hmmpress failed: {result.stderr[:300]}")
        return False
    return True


def _hmmsearch_domtblout(
    hmm_db: Path, query_fasta: Path, outpath: Path, config: Config
) -> Path:
    """Run hmmsearch with --domtblout. Module-level so tests can mock it."""
    cmd = [
        config.hmmsearch_bin,
        "--domtblout", str(outpath),
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


def build_profiles(
    families: Dict[str, Set[str]],
    alignment_lookup: Callable[[str], Optional[Path]],
    hmm_dir: Path,
    config: Config,
) -> Optional[Path]:
    """Build (or reuse cached) per-family HMMs and press a combined database.

    A family's profile is reused only when its .hmm file exists AND is
    non-empty AND its `.members` sidecar matches the sha1 of the current
    sorted member ids. This fixes two hmmer_rescue defects: zero-byte
    files from failed builds were treated as cached, and profiles were
    never rebuilt when family membership changed.

    Returns the combined, pressed HMM database path, or None when no
    profiles could be built.
    """
    hmm_dir = Path(hmm_dir)
    hmm_dir.mkdir(parents=True, exist_ok=True)

    ready: List[Path] = []
    work: List[Tuple[str, Path, Path, Path, str]] = []
    for family_id in sorted(families):
        aln_path = alignment_lookup(family_id)
        if aln_path is None:
            continue
        hmm_path = hmm_dir / f"{family_id}.hmm"
        members_path = hmm_dir / f"{family_id}.members"
        digest = _member_digest(families[family_id])
        is_cached = (
            hmm_path.exists()
            and hmm_path.stat().st_size > 0
            and members_path.exists()
            and members_path.read_text().strip() == digest
        )
        if is_cached:
            ready.append(hmm_path)
        else:
            work.append((family_id, Path(aln_path), hmm_path, members_path, digest))

    n_failed = 0
    if work:
        logger.info(
            f"Profile assignment: building {len(work)} HMM profiles "
            f"({len(ready)} cached)"
        )
        n_threads = max(1, min(config.n_workers, len(work)))
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            futures = {
                executor.submit(_hmmbuild, fam, aln, hmm, config): (hmm, members, digest)
                for fam, aln, hmm, members, digest in work
            }
            for future in as_completed(futures):
                hmm_path, members_path, digest = futures[future]
                try:
                    ok = future.result()
                except Exception as e:
                    logger.debug(f"hmmbuild worker exception for {hmm_path.stem}: {e}")
                    ok = False
                if ok:
                    members_path.write_text(digest + "\n")
                    ready.append(hmm_path)
                else:
                    n_failed += 1
    if n_failed:
        logger.info(f"Profile assignment: {n_failed} hmmbuild failures")

    if not ready:
        return None

    # Concatenate the explicit list of ready profiles (never a directory
    # glob — stale files and lexicographic order must not matter).
    combined = hmm_dir.parent / "all_families.hmm"
    with open(combined, "wb") as out_f:
        for hmm_path in sorted(ready):
            with open(hmm_path, "rb") as in_f:
                while True:
                    chunk = in_f.read(8192)
                    if not chunk:
                        break
                    out_f.write(chunk)

    if not _hmmpress(combined, config):
        return None
    return combined


def _margin_required(best_bits: float, config: Config) -> float:
    return max(config.profile_margin_bits, config.profile_margin_frac * best_bits)


def assign(
    hits_by_gene: Dict[str, List[ProfileHit]],
    current_assignment: Dict[str, Optional[str]],
    config: Config,
) -> Tuple[Dict[str, str], List[dict], Dict[str, Tuple[str, str]]]:
    """Decide profile assignments from parsed domtblout hits.

    NEW genes (current_assignment[gene] is None) are accepted into the
    top-bit-score family iff ALL of:
      * full_evalue <= config.hmmer_evalue (significance gate only)
      * profile_cov >= config.profile_min_coverage
      * query_cov  >= config.profile_min_query_coverage
      * no second hit, or best.full_bits - second.full_bits >=
        max(config.profile_margin_bits,
            config.profile_margin_frac * best.full_bits)
    The margin uses BIT SCORES, never E-values: tblout/domtblout E-values
    round to 2 significant figures and underflow to 0, and breaking ties
    by HMM-db order made R10_* families systematically beat R1_*.
    A gene failing ONLY the margin is reported in `ambiguous` with its
    top candidates instead of being force-assigned.

    ALREADY-PLACED genes are screened for reassignment: another family is
    a candidate iff its nbits exceeds the current family's nbits by
    config.profile_reassign_margin_nbits AND its profile_cov >= 0.6.
    These are returned as `moves` and are NOT auto-applied — this screen
    is biased (a gene inflates its own family's profile), so the caller
    must confirm each move with a leave-one-out profile rebuild and tree
    validation (a cluster-side step) before touching family membership.

    Returns:
        (new_assignments: gene -> family,
         ambiguous: [{"gene": ..., "candidates": [...]}, ...],
         moves: gene -> (from_family, to_family))
    """
    new_assignments: Dict[str, str] = {}
    ambiguous: List[dict] = []
    moves: Dict[str, Tuple[str, str]] = {}

    for gene_id, hits in hits_by_gene.items():
        if not hits:
            continue
        current_family = current_assignment.get(gene_id)
        if current_family is None:
            _assign_new_gene(gene_id, hits, config, new_assignments, ambiguous)
        else:
            move = _screen_reassignment(gene_id, hits, current_family, config)
            if move is not None:
                moves[gene_id] = move

    return new_assignments, ambiguous, moves


def _assign_new_gene(
    gene_id: str,
    hits: List[ProfileHit],
    config: Config,
    new_assignments: Dict[str, str],
    ambiguous: List[dict],
) -> None:
    """Apply the acceptance criteria for an unplaced gene (see assign())."""
    best = hits[0]
    passes_gates = (
        best.full_evalue <= config.hmmer_evalue
        and best.profile_cov >= config.profile_min_coverage
        and best.query_cov >= config.profile_min_query_coverage
    )
    if not passes_gates:
        return

    second = hits[1] if len(hits) > 1 else None
    has_margin = (
        second is None
        or best.full_bits - second.full_bits >= _margin_required(best.full_bits, config)
    )
    if has_margin:
        new_assignments[gene_id] = best.family_id
        return

    # Failing ONLY the margin: never pick a winner (ties included) —
    # report the top candidates instead.
    ambiguous.append({
        "gene": gene_id,
        "candidates": [
            {
                "family_id": h.family_id,
                "full_bits": h.full_bits,
                "profile_cov": h.profile_cov,
                "query_cov": h.query_cov,
            }
            for h in hits[:MAX_AMBIGUOUS_CANDIDATES]
        ],
    })


def _screen_reassignment(
    gene_id: str,
    hits: List[ProfileHit],
    current_family: str,
    config: Config,
) -> Optional[Tuple[str, str]]:
    """Biased screen for a better home than the gene's current family."""
    current_nbits = 0.0
    for h in hits:
        if h.family_id == current_family:
            current_nbits = h.nbits
            break

    best_candidate: Optional[ProfileHit] = None
    for h in hits:
        if h.family_id == current_family:
            continue
        if h.profile_cov < REASSIGN_MIN_PROFILE_COV:
            continue
        if h.nbits - current_nbits < config.profile_reassign_margin_nbits:
            continue
        if best_candidate is None or h.nbits > best_candidate.nbits:
            best_candidate = h

    if best_candidate is None:
        return None
    return (current_family, best_candidate.family_id)


def detect_merge_candidates(
    hits_by_gene: Dict[str, List[ProfileHit]],
    families: Dict[str, Set[str]],
    config: Config,
) -> List[Tuple[str, str, float]]:
    """Find family pairs whose members reciprocally hit each other's profile.

    frac(A->B) = share of A's members whose best NON-SELF family hit
    (by full_bits, restricted to hits with profile_cov >= 0.5 and
    full_evalue <= config.hmmer_evalue) is B. A pair is a merge
    candidate iff min(frac(A->B), frac(B->A)) >= config.merge_min_reciprocal.

    Candidates still require tree validation before an actual merge
    (align A∪B, build the codon tree, score every leaf with the pruning
    statistic) — this function only nominates pairs.
    """
    gene_to_family: Dict[str, str] = {}
    for family_id, members in families.items():
        for gene_id in members:
            gene_to_family[gene_id] = family_id

    # Count, per ordered pair (A, B), members of A voting for B.
    votes: Dict[Tuple[str, str], int] = {}
    for gene_id, hits in hits_by_gene.items():
        own_family = gene_to_family.get(gene_id)
        if own_family is None:
            continue
        for h in hits:  # sorted by full_bits desc
            if h.family_id == own_family:
                continue
            if h.profile_cov < MERGE_MIN_PROFILE_COV:
                continue
            if h.full_evalue > config.hmmer_evalue:
                continue
            key = (own_family, h.family_id)
            votes[key] = votes.get(key, 0) + 1
            break  # only the best qualifying non-self hit votes

    candidates: List[Tuple[str, str, float]] = []
    seen: Set[Tuple[str, str]] = set()
    for (fam_a, fam_b), n_ab in votes.items():
        pair = tuple(sorted((fam_a, fam_b)))
        if pair in seen:
            continue
        seen.add(pair)
        n_ba = votes.get((fam_b, fam_a), 0)
        size_a = len(families.get(fam_a, set()))
        size_b = len(families.get(fam_b, set()))
        if size_a == 0 or size_b == 0:
            continue
        frac_ab = n_ab / size_a
        frac_ba = n_ba / size_b
        min_frac = min(frac_ab, frac_ba)
        if min_frac >= config.merge_min_reciprocal:
            candidates.append((pair[0], pair[1], min_frac))

    candidates.sort(key=lambda c: (-c[2], c[0], c[1]))
    return candidates


def run_profile_assignment(
    families: Dict[str, Set[str]],
    unplaced_pool: Dict[str, str],
    protein_lookup: Callable[[str], Optional[Path]],
    outdir: Path,
    config: Config,
) -> Tuple[Dict[str, Set[str]], Set[str]]:
    """Offer unplaced genes to family profiles and accept unambiguous hits.

    Args:
        families: family_id -> confirmed member gene_ids.
        unplaced_pool: gene_id -> protein sequence for unplaced genes.
        protein_lookup: family_id -> path to the family's confirmed
            protein alignment (input to hmmbuild), or None if missing.
        outdir: directory for profiles, search output, and TSV reports.
        config: pipeline configuration.

    Returns:
        (assigned: family_id -> set of newly assigned gene_ids,
         still_unplaced: gene_ids that remain unplaced)
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    still_unplaced: Set[str] = set(unplaced_pool)

    if not families or not unplaced_pool:
        return {}, still_unplaced

    hmm_db = build_profiles(families, protein_lookup, outdir / "hmm_profiles", config)
    if hmm_db is None:
        logger.warning("Profile assignment: no HMM profiles available, skipping")
        return {}, still_unplaced

    query_fasta = outdir / "unplaced_proteins.fa"
    write_fasta(unplaced_pool, str(query_fasta))

    domtblout = outdir / "hmmsearch_results.domtblout"
    _hmmsearch_domtblout(hmm_db, query_fasta, domtblout, config)

    hits_by_gene = parse_domtblout(domtblout)
    current_assignment: Dict[str, Optional[str]] = {g: None for g in unplaced_pool}
    new_assignments, ambiguous, _moves = assign(hits_by_gene, current_assignment, config)

    assigned: Dict[str, Set[str]] = {}
    for gene_id, family_id in new_assignments.items():
        assigned.setdefault(family_id, set()).add(gene_id)
        still_unplaced.discard(gene_id)

    _write_assignments_tsv(new_assignments, hits_by_gene, outdir / "profile_assignments.tsv")
    _write_ambiguous_tsv(ambiguous, outdir / "ambiguous_assignments.tsv")

    logger.info(
        f"Profile assignment: {len(new_assignments)} genes assigned to "
        f"{len(assigned)} families, {len(ambiguous)} ambiguous, "
        f"{len(still_unplaced)} still unplaced"
    )
    return assigned, still_unplaced


def _write_assignments_tsv(
    new_assignments: Dict[str, str],
    hits_by_gene: Dict[str, List[ProfileHit]],
    path: Path,
) -> None:
    with open(path, "w") as f:
        f.write("gene_id\tfamily_id\tfull_evalue\tfull_bits\tprofile_cov\tquery_cov\n")
        for gene_id in sorted(new_assignments):
            family_id = new_assignments[gene_id]
            hit = next(
                (h for h in hits_by_gene.get(gene_id, []) if h.family_id == family_id),
                None,
            )
            if hit is None:
                continue
            f.write(
                f"{gene_id}\t{family_id}\t{hit.full_evalue:.2e}\t"
                f"{hit.full_bits:.1f}\t{hit.profile_cov:.3f}\t{hit.query_cov:.3f}\n"
            )
    logger.debug(f"Profile assignments written: {path}")


def _write_ambiguous_tsv(ambiguous: List[dict], path: Path) -> None:
    with open(path, "w") as f:
        f.write("gene_id\trank\tfamily_id\tfull_bits\tprofile_cov\tquery_cov\n")
        for entry in sorted(ambiguous, key=lambda e: e["gene"]):
            for rank, cand in enumerate(entry["candidates"], start=1):
                f.write(
                    f"{entry['gene']}\t{rank}\t{cand['family_id']}\t"
                    f"{cand['full_bits']:.1f}\t{cand['profile_cov']:.3f}\t"
                    f"{cand['query_cov']:.3f}\n"
                )
    logger.debug(f"Ambiguous assignments written: {path}")
