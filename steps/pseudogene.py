"""Pseudogene detection module for the family_finder pipeline.

Identifies pseudogene candidates using multiple lines of evidence:
  1. Internal stop codons in protein sequences (premature termination)
  2. CDS/protein length discrepancies (potential frameshifts)
  3. Truncated genes (significantly shorter than family homologs)
  4. Orphan genes (unplaced after iterative clustering + HMMER rescue)
  5. Elevated branch lengths in gene trees (relaxed selection pressure)
  6. Intact-but-divergent genes flagged by pruning (ratio outliers)

Each gene receives an evidence vector and a composite classification:
  - "pseudogene_high"   : 2+ strong evidence lines
  - "pseudogene_medium" : 1 strong or 2+ weak evidence lines
  - "pseudogene_low"    : 1 weak evidence line only
  - "functional"        : no pseudogene evidence
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from config import Config

logger = logging.getLogger("family_finder")


def _get_species(gene_id: str, delimiter: str = "_") -> str:
    """Extract species prefix from a gene ID (avoids ete4 import)."""
    return gene_id.split(delimiter, 1)[0]


# ---------------------------------------------------------------------------
# Evidence dataclass
# ---------------------------------------------------------------------------

class PseudogeneEvidence:
    """Accumulates evidence for a single gene being a pseudogene."""

    __slots__ = (
        "gene_id", "species", "internal_stops", "length_ratio",
        "family_id", "is_orphan", "distance_ratio", "is_truncated",
        "gc_content", "gc3_content", "composition_zscore",
        "notes",
    )

    def __init__(self, gene_id: str, species: str):
        self.gene_id = gene_id
        self.species = species
        self.internal_stops: int = 0           # count of internal stop codons
        self.length_ratio: Optional[float] = None  # CDS_len / (prot_len*3)
        self.family_id: Optional[str] = None   # assigned family (if any)
        self.is_orphan: bool = False            # unplaced after all rounds
        self.distance_ratio: Optional[float] = None  # observed/expected distance
        self.is_truncated: bool = False         # shorter than family median
        self.gc_content: Optional[float] = None      # overall GC%
        self.gc3_content: Optional[float] = None     # GC% at 3rd codon position
        self.composition_zscore: Optional[float] = None  # z-score vs species mean
        self.notes: List[str] = []

    @property
    def confidence_score(self) -> float:
        """Compute a weighted confidence score (0-1 range).

        Weights reflect biological informativeness:
          - internal stop codons:  0.40 (direct evidence of loss-of-function)
          - CDS length mismatch:   0.35 (frameshift / assembly error)
          - truncated:             0.30 (partial gene model)
          - orphan:                0.25 (unplaced = likely pseudogene or foreign)
          - GC3 outlier:           0.25 (compositional drift = pseudogene or HGT)
          - long branch:           0.20 (accelerated evolution)
        """
        score = 0.0
        if self.internal_stops > 0:
            score += 0.40
        if self.length_ratio is not None and abs(self.length_ratio - 1.0) > 0.1:
            score += 0.35
        if self.is_truncated:
            score += 0.30
        if self.is_orphan:
            score += 0.25
        if self.composition_zscore is not None and abs(self.composition_zscore) > 3.0:
            score += 0.25
        if self.distance_ratio is not None and self.distance_ratio > 3.0:
            score += 0.20
        return min(score, 1.0)

    @property
    def classification(self) -> str:
        """Classify based on weighted confidence score.

        Thresholds:
          high   : score >= 0.50  (multiple evidence lines or 1 very strong)
          medium : score >= 0.25  (orphan, GC3 outlier, truncated, etc.)
          low    : score >  0     (single weak indicator like long branch)
        """
        s = self.confidence_score
        if s >= 0.50:
            return "pseudogene_high"
        if s >= 0.25:
            return "pseudogene_medium"
        if s > 0:
            return "pseudogene_low"
        return "functional"

    @property
    def evidence_summary(self) -> str:
        """One-line summary of all evidence."""
        parts = []
        if self.internal_stops > 0:
            parts.append(f"internal_stops={self.internal_stops}")
        if self.length_ratio is not None and abs(self.length_ratio - 1.0) > 0.1:
            parts.append(f"length_ratio={self.length_ratio:.3f}")
        if self.is_truncated:
            parts.append("truncated")
        if self.is_orphan:
            parts.append("orphan")
        if self.distance_ratio is not None and self.distance_ratio > 3.0:
            parts.append(f"dist_ratio={self.distance_ratio:.2f}")
        if self.composition_zscore is not None and abs(self.composition_zscore) > 3.0:
            parts.append(f"gc3_zscore={self.composition_zscore:.2f}")
        if self.notes:
            parts.extend(self.notes)
        return ";".join(parts) if parts else "none"


# ---------------------------------------------------------------------------
# Evidence scanners
# ---------------------------------------------------------------------------

def scan_internal_stops(
    protein_seqs: Dict[str, str],
    species_delimiter: str = "_",
) -> Dict[str, PseudogeneEvidence]:
    """Scan protein sequences for internal stop codons.

    Internal '*' characters (excluding terminal) indicate premature
    termination — a hallmark of pseudogenes or annotation errors.
    """
    get_species = _get_species

    results: Dict[str, PseudogeneEvidence] = {}
    n_found = 0

    for gene_id, seq in protein_seqs.items():
        seq_clean = seq.rstrip("*")
        n_stops = seq_clean.count("*")
        if n_stops > 0:
            sp = get_species(gene_id, species_delimiter)
            ev = PseudogeneEvidence(gene_id, sp)
            ev.internal_stops = n_stops
            ev.notes.append(f"premature_stop_codon(n={n_stops})")
            results[gene_id] = ev
            n_found += 1

    logger.info(f"Internal stop scan: {n_found}/{len(protein_seqs)} genes with internal stops")
    return results


def scan_length_discrepancies(
    protein_seqs: Dict[str, str],
    cds_seqs: Dict[str, str],
    species_delimiter: str = "_",
    tolerance: float = 0.1,
) -> Dict[str, PseudogeneEvidence]:
    """Detect CDS/protein length mismatches suggesting frameshifts.

    Expected: len(CDS) == len(protein)*3 or len(protein)*3 + 3 (stop codon).
    Deviations beyond tolerance fraction suggest annotation errors or frameshifts.
    """
    get_species = _get_species

    results: Dict[str, PseudogeneEvidence] = {}

    for gene_id in protein_seqs:
        if gene_id not in cds_seqs:
            continue
        prot_len = len(protein_seqs[gene_id].rstrip("*"))
        cds_len = len(cds_seqs[gene_id].replace("-", ""))
        if prot_len == 0:
            continue

        expected_cds = prot_len * 3
        ratio = cds_len / expected_cds

        if abs(ratio - 1.0) > tolerance:
            sp = get_species(gene_id, species_delimiter)
            ev = PseudogeneEvidence(gene_id, sp)
            ev.length_ratio = ratio
            if ratio < 1.0 - tolerance:
                ev.notes.append(f"cds_too_short(ratio={ratio:.3f})")
            else:
                ev.notes.append(f"cds_too_long(ratio={ratio:.3f})")
            results[gene_id] = ev

    logger.info(f"Length discrepancy scan: {len(results)}/{len(protein_seqs)} genes with mismatches")
    return results


def scan_truncated_genes(
    protein_seqs: Dict[str, str],
    families: Dict[str, Set[str]],
    species_delimiter: str = "_",
    truncation_threshold: float = 0.5,
) -> Dict[str, PseudogeneEvidence]:
    """Identify genes significantly shorter than their family's median length.

    A gene shorter than truncation_threshold * family_median is flagged.
    """
    from statistics import median
    get_species = _get_species

    results: Dict[str, PseudogeneEvidence] = {}

    # Build gene -> family mapping
    gene_to_family: Dict[str, str] = {}
    for fam_id, gene_ids in families.items():
        for gid in gene_ids:
            gene_to_family[gid] = fam_id

    # Compute family median lengths
    family_lengths: Dict[str, List[int]] = defaultdict(list)
    for fam_id, gene_ids in families.items():
        for gid in gene_ids:
            if gid in protein_seqs:
                family_lengths[fam_id].append(len(protein_seqs[gid].rstrip("*")))

    family_medians: Dict[str, float] = {}
    for fam_id, lengths in family_lengths.items():
        if lengths:
            family_medians[fam_id] = median(lengths)

    # Flag truncated genes
    for gene_id, seq in protein_seqs.items():
        fam_id = gene_to_family.get(gene_id)
        if fam_id is None or fam_id not in family_medians:
            continue
        gene_len = len(seq.rstrip("*"))
        med_len = family_medians[fam_id]
        if med_len > 0 and gene_len < truncation_threshold * med_len:
            sp = get_species(gene_id, species_delimiter)
            ev = PseudogeneEvidence(gene_id, sp)
            ev.is_truncated = True
            ev.family_id = fam_id
            ratio = gene_len / med_len
            ev.notes.append(f"truncated({gene_len}aa_vs_median_{med_len:.0f}aa,ratio={ratio:.2f})")
            results[gene_id] = ev

    logger.info(f"Truncation scan: {len(results)} genes flagged as truncated")
    return results


def scan_orphan_genes(
    all_protein_seqs: Dict[str, str],
    families: Dict[str, Set[str]],
    species_delimiter: str = "_",
    species_filter: Optional[str] = None,
) -> Dict[str, PseudogeneEvidence]:
    """Identify genes not assigned to any family (orphan/unplaced genes).

    Orphan genes after exhaustive clustering + HMMER rescue are weak
    pseudogene candidates — they may be too divergent to cluster.
    """
    get_species = _get_species

    # Collect all genes in families
    placed_genes: Set[str] = set()
    for gene_ids in families.values():
        placed_genes.update(gene_ids)

    results: Dict[str, PseudogeneEvidence] = {}

    for gene_id in all_protein_seqs:
        if gene_id in placed_genes:
            continue
        sp = get_species(gene_id, species_delimiter)
        if species_filter and sp != species_filter:
            continue
        ev = PseudogeneEvidence(gene_id, sp)
        ev.is_orphan = True
        ev.notes.append("unplaced_after_rescue")
        results[gene_id] = ev

    logger.info(f"Orphan scan: {len(results)} unplaced genes")
    return results


def scan_gc3_composition(
    cds_seqs: Dict[str, str],
    species_delimiter: str = "_",
    zscore_threshold: float = 3.0,
) -> Dict[str, PseudogeneEvidence]:
    """Detect genes with anomalous GC content at 3rd codon positions (GC3).

    Pseudogenes accumulate mutations freely, causing their GC3 to drift
    from the species-specific mean. Genes with |z-score| > threshold
    relative to the species GC3 distribution are flagged.
    """
    from statistics import mean, stdev
    get_species = _get_species

    results: Dict[str, PseudogeneEvidence] = {}

    # Compute GC3 per gene
    gene_gc3: Dict[str, float] = {}
    gene_gc: Dict[str, float] = {}
    for gene_id, seq in cds_seqs.items():
        seq_upper = seq.upper().replace("-", "").replace("N", "")
        if len(seq_upper) < 9:  # need at least 3 codons
            continue

        # Overall GC
        gc_count = seq_upper.count("G") + seq_upper.count("C")
        gene_gc[gene_id] = gc_count / len(seq_upper) if seq_upper else 0

        # GC at 3rd codon position
        third_bases = [seq_upper[i] for i in range(2, len(seq_upper), 3)
                       if i < len(seq_upper)]
        if not third_bases:
            continue
        gc3 = sum(1 for b in third_bases if b in "GC") / len(third_bases)
        gene_gc3[gene_id] = gc3

    # Compute per-species GC3 stats
    species_gc3: Dict[str, List[float]] = defaultdict(list)
    for gene_id, gc3 in gene_gc3.items():
        sp = get_species(gene_id, species_delimiter)
        species_gc3[sp].append(gc3)

    species_stats: Dict[str, Tuple[float, float]] = {}
    for sp, values in species_gc3.items():
        if len(values) >= 20:
            species_stats[sp] = (mean(values), stdev(values))

    # Flag outliers
    for gene_id, gc3 in gene_gc3.items():
        sp = get_species(gene_id, species_delimiter)
        if sp not in species_stats:
            continue
        sp_mean, sp_std = species_stats[sp]
        if sp_std <= 0:
            continue

        zscore = (gc3 - sp_mean) / sp_std
        if abs(zscore) > zscore_threshold:
            ev = PseudogeneEvidence(gene_id, sp)
            ev.gc_content = gene_gc.get(gene_id)
            ev.gc3_content = gc3
            ev.composition_zscore = zscore
            ev.notes.append(
                f"gc3_outlier(gc3={gc3:.3f},sp_mean={sp_mean:.3f},zscore={zscore:.1f})"
            )
            results[gene_id] = ev

    logger.info(f"GC3 composition scan: {len(results)} outliers (|z|>{zscore_threshold})")
    return results


def scan_distance_outliers(
    outdir: Path,
    species_delimiter: str = "_",
    species_filter: Optional[str] = None,
) -> Dict[str, PseudogeneEvidence]:
    """Identify genes with elevated branch lengths in confirmed family trees.

    Genes with root-to-tip distance > 3x the family median suggest
    accelerated evolution consistent with pseudogenization.
    Falls back to parsing pipeline logs if ete4 is unavailable.
    """
    get_species = _get_species
    import re

    results: Dict[str, PseudogeneEvidence] = {}

    # Strategy 1: Scan confirmed trees in final_families/ using Bio.Phylo
    final_dir = outdir / "final_families"
    if final_dir.exists():
        try:
            from Bio import Phylo
            import io
            has_phylo = True
        except ImportError:
            has_phylo = False

        if has_phylo:
            from statistics import median as _median
            n_trees = 0
            for fam_dir in sorted(final_dir.iterdir()):
                if not fam_dir.is_dir():
                    continue
                tree_path = fam_dir / "confirmed_tree.nwk"
                if not tree_path.exists():
                    continue
                try:
                    tree_text = tree_path.read_text().strip()
                    tree = Phylo.read(io.StringIO(tree_text), "newick")
                except Exception:
                    continue

                n_trees += 1
                terminals = tree.get_terminals()
                if len(terminals) < 3:
                    continue

                # Compute root-to-tip distances
                bl_map = {}
                for term in terminals:
                    if term.name:
                        try:
                            dist = tree.distance(term)
                            bl_map[term.name] = dist
                        except Exception:
                            pass

                if len(bl_map) < 3:
                    continue

                med_bl = _median(list(bl_map.values()))
                if med_bl <= 0:
                    continue

                for gene_id, bl in bl_map.items():
                    ratio = bl / med_bl
                    if ratio <= 3.0:
                        continue
                    sp = get_species(gene_id, species_delimiter)
                    if species_filter and sp != species_filter:
                        continue
                    ev = PseudogeneEvidence(gene_id, sp)
                    ev.distance_ratio = ratio
                    ev.family_id = fam_dir.name
                    ev.notes.append(
                        f"long_branch(bl={bl:.4f},median={med_bl:.4f},ratio={ratio:.1f})"
                    )
                    results[gene_id] = ev

            logger.info(f"Branch length scan: {len(results)} outliers from {n_trees} trees")
            return results

    # Strategy 2: Parse pipeline log for ratio outlier messages
    log_path = outdir / "pipeline.log"
    if log_path.exists():
        pattern = re.compile(
            r"Ratio outlier:\s+(\S+)\s+\(species=(\S+),\s+score=([0-9.]+)\)"
        )
        with open(log_path) as f:
            for line in f:
                m = pattern.search(line)
                if m:
                    gene_id = m.group(1)
                    sp = m.group(2)
                    score = float(m.group(3))
                    if species_filter and sp != species_filter:
                        continue
                    ev = PseudogeneEvidence(gene_id, sp)
                    ev.distance_ratio = score
                    ev.notes.append(f"pruning_outlier(ratio={score:.2f})")
                    results[gene_id] = ev

    logger.info(f"Distance outlier scan: {len(results)} genes")
    return results


# ---------------------------------------------------------------------------
# Composite analysis
# ---------------------------------------------------------------------------

def merge_evidence(
    *evidence_dicts: Dict[str, PseudogeneEvidence],
) -> Dict[str, PseudogeneEvidence]:
    """Merge evidence from multiple scanners into a single dict per gene."""
    merged: Dict[str, PseudogeneEvidence] = {}

    for ev_dict in evidence_dicts:
        for gene_id, ev in ev_dict.items():
            if gene_id not in merged:
                merged[gene_id] = ev
            else:
                existing = merged[gene_id]
                # Merge fields
                if ev.internal_stops > 0:
                    existing.internal_stops = max(existing.internal_stops, ev.internal_stops)
                if ev.length_ratio is not None:
                    existing.length_ratio = ev.length_ratio
                if ev.is_truncated:
                    existing.is_truncated = True
                if ev.is_orphan:
                    existing.is_orphan = True
                if ev.distance_ratio is not None:
                    if existing.distance_ratio is None or ev.distance_ratio > existing.distance_ratio:
                        existing.distance_ratio = ev.distance_ratio
                if ev.gc_content is not None:
                    existing.gc_content = ev.gc_content
                if ev.gc3_content is not None:
                    existing.gc3_content = ev.gc3_content
                if ev.composition_zscore is not None:
                    existing.composition_zscore = ev.composition_zscore
                if ev.family_id is not None:
                    existing.family_id = ev.family_id
                # Merge notes (avoid duplicates)
                existing_notes = set(existing.notes)
                for note in ev.notes:
                    if note not in existing_notes:
                        existing.notes.append(note)

    return merged


def detect_pseudogenes(
    protein_seqs: Dict[str, str],
    cds_seqs: Dict[str, str],
    families: Dict[str, Set[str]],
    outdir: Path,
    config: Config,
    species_filter: Optional[str] = None,
) -> Dict[str, PseudogeneEvidence]:
    """Run all pseudogene evidence scanners and merge results.

    Args:
        protein_seqs: All protein sequences (gene_id -> seq).
        cds_seqs: All CDS sequences (gene_id -> seq).
        families: Confirmed families (family_id -> set of gene_ids).
        outdir: Pipeline output directory (for log parsing).
        config: Pipeline configuration.
        species_filter: Optional species prefix to restrict analysis.

    Returns:
        Dict of gene_id -> PseudogeneEvidence for all candidates.
    """
    delimiter = config.species_delimiter

    # Filter sequences by species if requested
    if species_filter:
        get_species = _get_species
        protein_seqs = {
            gid: seq for gid, seq in protein_seqs.items()
            if get_species(gid, delimiter) == species_filter
        }
        cds_seqs = {
            gid: seq for gid, seq in cds_seqs.items()
            if get_species(gid, delimiter) == species_filter
        }

    logger.info(
        f"Pseudogene detection: analyzing {len(protein_seqs)} genes"
        + (f" (species={species_filter})" if species_filter else "")
    )

    # Run all scanners
    ev_stops = scan_internal_stops(protein_seqs, delimiter)
    ev_length = scan_length_discrepancies(protein_seqs, cds_seqs, delimiter)
    ev_truncated = scan_truncated_genes(
        protein_seqs, families, delimiter,
        truncation_threshold=config.pseudogene_truncation_threshold,
    )
    ev_orphans = scan_orphan_genes(protein_seqs, families, delimiter, species_filter)
    ev_distance = scan_distance_outliers(outdir, delimiter, species_filter)
    ev_gc3 = scan_gc3_composition(cds_seqs, delimiter)

    # Merge all evidence
    merged = merge_evidence(ev_stops, ev_length, ev_truncated, ev_orphans, ev_distance, ev_gc3)

    # Add family assignment info for genes without it
    gene_to_family: Dict[str, str] = {}
    for fam_id, gene_ids in families.items():
        for gid in gene_ids:
            gene_to_family[gid] = fam_id
    for gene_id, ev in merged.items():
        if ev.family_id is None and gene_id in gene_to_family:
            ev.family_id = gene_to_family[gene_id]

    # Summary
    classifications = defaultdict(int)
    for ev in merged.values():
        classifications[ev.classification] += 1

    logger.info(
        f"Pseudogene detection complete: {len(merged)} candidates "
        f"(high={classifications.get('pseudogene_high', 0)}, "
        f"medium={classifications.get('pseudogene_medium', 0)}, "
        f"low={classifications.get('pseudogene_low', 0)})"
    )

    return merged


# ---------------------------------------------------------------------------
# Report writer
# ---------------------------------------------------------------------------

def write_pseudogene_report(
    evidence: Dict[str, PseudogeneEvidence],
    outpath: Path,
    species_filter: Optional[str] = None,
):
    """Write a TSV report of pseudogene candidates.

    Columns:
      gene_id, species, classification, family_id, internal_stops,
      length_ratio, is_truncated, is_orphan, distance_ratio, evidence_summary
    """
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    # Sort: high confidence first, then by gene_id
    priority = {"pseudogene_high": 0, "pseudogene_medium": 1, "pseudogene_low": 2, "functional": 3}
    sorted_evs = sorted(
        evidence.values(),
        key=lambda e: (priority.get(e.classification, 9), e.gene_id),
    )

    with open(outpath, "w") as f:
        f.write(
            "gene_id\tspecies\tclassification\tconfidence_score\tfamily_id\t"
            "internal_stops\tlength_ratio\tis_truncated\tis_orphan\t"
            "distance_ratio\tgc3\tgc3_zscore\tevidence_summary\n"
        )
        for ev in sorted_evs:
            f.write(
                f"{ev.gene_id}\t{ev.species}\t{ev.classification}\t"
                f"{ev.confidence_score:.3f}\t"
                f"{ev.family_id or 'none'}\t"
                f"{ev.internal_stops}\t"
                f"{ev.length_ratio if ev.length_ratio is not None else 'NA'}\t"
                f"{ev.is_truncated}\t{ev.is_orphan}\t"
                f"{ev.distance_ratio if ev.distance_ratio is not None else 'NA'}\t"
                f"{ev.gc3_content if ev.gc3_content is not None else 'NA'}\t"
                f"{ev.composition_zscore if ev.composition_zscore is not None else 'NA'}\t"
                f"{ev.evidence_summary}\n"
            )

    logger.info(f"Pseudogene report written to {outpath}")


def write_pseudogene_summary(
    evidence: Dict[str, PseudogeneEvidence],
    outpath: Path,
    total_genes: int,
    species_filter: Optional[str] = None,
):
    """Write a human-readable summary of pseudogene analysis."""
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    classifications = defaultdict(int)
    species_counts = defaultdict(lambda: defaultdict(int))
    evidence_types = defaultdict(int)

    for ev in evidence.values():
        cls = ev.classification
        classifications[cls] += 1
        species_counts[ev.species][cls] += 1
        if ev.internal_stops > 0:
            evidence_types["internal_stop_codons"] += 1
        if ev.length_ratio is not None and abs(ev.length_ratio - 1.0) > 0.1:
            evidence_types["length_discrepancy"] += 1
        if ev.is_truncated:
            evidence_types["truncated"] += 1
        if ev.is_orphan:
            evidence_types["orphan/unplaced"] += 1
        if ev.distance_ratio is not None and ev.distance_ratio > 3.0:
            evidence_types["distance_outlier"] += 1
        if ev.composition_zscore is not None and abs(ev.composition_zscore) > 3.0:
            evidence_types["gc3_composition_outlier"] += 1

    n_pseudo = sum(
        v for k, v in classifications.items() if k.startswith("pseudogene")
    )

    with open(outpath, "w") as f:
        f.write("=" * 60 + "\n")
        f.write("PSEUDOGENE DETECTION SUMMARY\n")
        f.write("=" * 60 + "\n\n")

        if species_filter:
            f.write(f"Species filter: {species_filter}\n")
        f.write(f"Total genes analyzed: {total_genes}\n")
        f.write(f"Pseudogene candidates: {n_pseudo} ({100*n_pseudo/max(total_genes,1):.1f}%)\n\n")

        f.write("Classification breakdown:\n")
        for cls in ["pseudogene_high", "pseudogene_medium", "pseudogene_low", "functional"]:
            n = classifications.get(cls, 0)
            f.write(f"  {cls:25s}: {n:6d}\n")
        f.write("\n")

        f.write("Evidence type counts:\n")
        for etype, count in sorted(evidence_types.items(), key=lambda x: -x[1]):
            f.write(f"  {etype:25s}: {count:6d}\n")
        f.write("\n")

        if len(species_counts) > 1:
            f.write("Per-species breakdown:\n")
            f.write(f"  {'Species':15s} {'high':>8s} {'medium':>8s} {'low':>8s} {'total':>8s}\n")
            for sp in sorted(species_counts):
                sc = species_counts[sp]
                total = sum(sc.values())
                f.write(
                    f"  {sp:15s} "
                    f"{sc.get('pseudogene_high', 0):8d} "
                    f"{sc.get('pseudogene_medium', 0):8d} "
                    f"{sc.get('pseudogene_low', 0):8d} "
                    f"{total:8d}\n"
                )
            f.write("\n")
        elif species_filter:
            sc = species_counts.get(species_filter, {})
            f.write(f"Species {species_filter}:\n")
            for cls in ["pseudogene_high", "pseudogene_medium", "pseudogene_low"]:
                f.write(f"  {cls:25s}: {sc.get(cls, 0):6d}\n")
            f.write("\n")

    logger.info(f"Pseudogene summary written to {outpath}")


def write_pseudogene_fasta(
    evidence: Dict[str, PseudogeneEvidence],
    protein_seqs: Dict[str, str],
    cds_seqs: Dict[str, str],
    outdir: Path,
    min_confidence: str = "pseudogene_low",
):
    """Write FASTA files of pseudogene candidate sequences.

    Writes separate files for protein and CDS sequences, with
    classification in the FASTA header.
    """
    from utils.seqio import write_fasta

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    priority = {"pseudogene_high": 0, "pseudogene_medium": 1, "pseudogene_low": 2}
    min_priority = priority.get(min_confidence, 2)

    prot_out: Dict[str, str] = {}
    cds_out: Dict[str, str] = {}

    for gene_id, ev in evidence.items():
        cls_priority = priority.get(ev.classification, 99)
        if cls_priority > min_priority:
            continue

        if gene_id in protein_seqs:
            prot_out[gene_id] = protein_seqs[gene_id]
        if gene_id in cds_seqs:
            cds_out[gene_id] = cds_seqs[gene_id]

    if prot_out:
        write_fasta(prot_out, str(outdir / "pseudogene_candidates.pep.fa"))
    if cds_out:
        write_fasta(cds_out, str(outdir / "pseudogene_candidates.cds.fa"))

    logger.info(
        f"Pseudogene FASTA: {len(prot_out)} protein, {len(cds_out)} CDS sequences"
    )


# ---------------------------------------------------------------------------
# Per-family pseudogene enrichment
# ---------------------------------------------------------------------------

def write_family_pseudogene_report(
    evidence: Dict[str, PseudogeneEvidence],
    families: Dict[str, Set[str]],
    outpath: Path,
):
    """Write a per-family report showing pseudogene enrichment.

    For each family, reports:
      - total members, number of pseudogene candidates, pseudogene fraction
      - per-species breakdown of pseudogenes within the family
      - average confidence score of pseudogene members
    Sorted by pseudogene fraction descending.
    """
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for fam_id, gene_ids in families.items():
        n_total = len(gene_ids)
        pseudo_in_fam = []
        species_pseudo = defaultdict(int)
        species_total = defaultdict(int)

        for gid in gene_ids:
            sp = _get_species(gid)
            species_total[sp] += 1
            if gid in evidence and evidence[gid].classification.startswith("pseudogene"):
                pseudo_in_fam.append(evidence[gid])
                species_pseudo[sp] += 1

        n_pseudo = len(pseudo_in_fam)
        if n_pseudo == 0:
            continue

        frac = n_pseudo / n_total if n_total > 0 else 0
        avg_score = sum(e.confidence_score for e in pseudo_in_fam) / n_pseudo

        # Per-species breakdown string
        sp_parts = []
        for sp in sorted(species_pseudo):
            sp_parts.append(f"{sp}:{species_pseudo[sp]}/{species_total[sp]}")

        rows.append({
            "family_id": fam_id,
            "n_total": n_total,
            "n_pseudo": n_pseudo,
            "pseudo_frac": frac,
            "avg_confidence": avg_score,
            "species_breakdown": ",".join(sp_parts),
            "n_species": len(species_total),
        })

    rows.sort(key=lambda r: (-r["pseudo_frac"], -r["n_pseudo"]))

    with open(outpath, "w") as f:
        f.write(
            "family_id\tn_total\tn_pseudo\tpseudo_frac\tavg_confidence\t"
            "n_species\tspecies_breakdown\n"
        )
        for r in rows:
            f.write(
                f"{r['family_id']}\t{r['n_total']}\t{r['n_pseudo']}\t"
                f"{r['pseudo_frac']:.3f}\t{r['avg_confidence']:.3f}\t"
                f"{r['n_species']}\t{r['species_breakdown']}\n"
            )

    logger.info(
        f"Family pseudogene report: {len(rows)} families with pseudogenes "
        f"(out of {len(families)} total)"
    )


# ---------------------------------------------------------------------------
# Cross-species comparison
# ---------------------------------------------------------------------------

def write_species_comparison(
    evidence: Dict[str, PseudogeneEvidence],
    all_protein_seqs: Dict[str, str],
    families: Dict[str, Set[str]],
    outpath: Path,
    species_delimiter: str = "_",
):
    """Write a cross-species pseudogene comparison report.

    For each species, reports total genes, pseudogene counts by confidence,
    orphan count, and pseudogene rate. Highlights species-specific
    pseudogene enrichment.
    """
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    # Count genes per species
    sp_total: Dict[str, int] = defaultdict(int)
    for gid in all_protein_seqs:
        sp = _get_species(gid, species_delimiter)
        sp_total[sp] += 1

    # Count placed genes per species
    placed = set()
    for gene_ids in families.values():
        placed.update(gene_ids)
    sp_placed: Dict[str, int] = defaultdict(int)
    for gid in placed:
        sp = _get_species(gid, species_delimiter)
        sp_placed[sp] += 1

    # Count pseudogenes per species by classification
    sp_pseudo: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    sp_scores: Dict[str, List[float]] = defaultdict(list)
    for gid, ev in evidence.items():
        if ev.classification.startswith("pseudogene"):
            sp_pseudo[ev.species][ev.classification] += 1
            sp_scores[ev.species].append(ev.confidence_score)

    from statistics import mean

    with open(outpath, "w") as f:
        f.write(
            "species\ttotal_genes\tplaced_genes\torphan_genes\t"
            "pseudo_high\tpseudo_medium\tpseudo_low\tpseudo_total\t"
            "pseudo_rate\tavg_confidence\n"
        )
        for sp in sorted(sp_total):
            total = sp_total[sp]
            n_placed = sp_placed.get(sp, 0)
            n_orphan = total - n_placed
            pc = sp_pseudo.get(sp, {})
            n_high = pc.get("pseudogene_high", 0)
            n_med = pc.get("pseudogene_medium", 0)
            n_low = pc.get("pseudogene_low", 0)
            n_total_pseudo = n_high + n_med + n_low
            rate = n_total_pseudo / total if total > 0 else 0
            scores = sp_scores.get(sp, [])
            avg_conf = mean(scores) if scores else 0

            f.write(
                f"{sp}\t{total}\t{n_placed}\t{n_orphan}\t"
                f"{n_high}\t{n_med}\t{n_low}\t{n_total_pseudo}\t"
                f"{rate:.4f}\t{avg_conf:.3f}\n"
            )

    logger.info(f"Species comparison written to {outpath}")


# ---------------------------------------------------------------------------
# BED output for genome browser
# ---------------------------------------------------------------------------

def write_pseudogene_bed(
    evidence: Dict[str, PseudogeneEvidence],
    outpath: Path,
    min_confidence: str = "pseudogene_low",
):
    """Write a BED file of pseudogene candidates for genome browser viewing.

    Extracts chromosome and position from gene IDs when possible.
    Gene IDs like Ococ_OcoChr01G00010 → chrom=OcoChr01, start from gene number.

    BED columns: chrom, start, end, name, score (confidence*1000), strand
    Color: red for high, orange for medium, yellow for low confidence.
    """
    import re
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    priority = {"pseudogene_high": 0, "pseudogene_medium": 1, "pseudogene_low": 2}
    min_priority = priority.get(min_confidence, 2)

    # Color by confidence (RGB for BED)
    colors = {
        "pseudogene_high": "255,0,0",       # red
        "pseudogene_medium": "255,165,0",   # orange
        "pseudogene_low": "255,255,0",      # yellow
    }

    # Pattern to extract chromosome and gene position
    # Handles: Ococ_OcoChr01G00010, Ococ_OcoScaf00000236G00030
    chr_pattern = re.compile(
        r"^[^_]+_([A-Za-z]+(?:Chr\d+|Scaf\d+(?:_\d+_\d+)?))G(\d+)$"
    )

    rows = []
    for gene_id, ev in evidence.items():
        cls_p = priority.get(ev.classification, 99)
        if cls_p > min_priority:
            continue

        m = chr_pattern.match(gene_id)
        if m:
            chrom = m.group(1)
            gene_num = int(m.group(2))
            # Approximate positions (10kb per gene unit)
            start = (gene_num - 1) * 10000
            end = start + 10000
        else:
            chrom = ev.species + "_unknown"
            start = 0
            end = 1000

        score = int(ev.confidence_score * 1000)
        color = colors.get(ev.classification, "128,128,128")
        rows.append((chrom, start, end, gene_id, score, ".", start, end, color))

    rows.sort(key=lambda r: (r[0], r[1]))

    with open(outpath, "w") as f:
        f.write(f'track name="pseudogenes" description="Pseudogene candidates" '
                f'itemRgb="On"\n')
        for chrom, start, end, name, score, strand, tstart, tend, color in rows:
            f.write(
                f"{chrom}\t{start}\t{end}\t{name}\t{score}\t"
                f"{strand}\t{tstart}\t{tend}\t{color}\n"
            )

    logger.info(f"BED file: {len(rows)} entries written to {outpath}")


# ---------------------------------------------------------------------------
# Chromosomal distribution analysis
# ---------------------------------------------------------------------------

def write_chromosomal_distribution(
    evidence: Dict[str, PseudogeneEvidence],
    all_protein_seqs: Dict[str, str],
    outpath: Path,
    species_filter: Optional[str] = None,
    species_delimiter: str = "_",
):
    """Analyze pseudogene distribution across chromosomes.

    Reports pseudogene density per chromosome/scaffold to identify
    hotspots of pseudogenization (e.g., pericentromeric regions).
    """
    import re
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    chr_pattern = re.compile(
        r"^[^_]+_([A-Za-z]+(?:Chr\d+|Scaf[^G]+))G\d+$"
    )

    chrom_total: Dict[str, int] = defaultdict(int)
    chrom_pseudo: Dict[str, int] = defaultdict(int)
    chrom_high: Dict[str, int] = defaultdict(int)

    for gid in all_protein_seqs:
        sp = _get_species(gid, species_delimiter)
        if species_filter and sp != species_filter:
            continue
        m = chr_pattern.match(gid)
        if m:
            chrom = m.group(1)
            chrom_total[chrom] += 1

    for gid, ev in evidence.items():
        if not ev.classification.startswith("pseudogene"):
            continue
        m = chr_pattern.match(gid)
        if m:
            chrom = m.group(1)
            chrom_pseudo[chrom] += 1
            if ev.classification == "pseudogene_high":
                chrom_high[chrom] += 1

    # Build rows and sort by pseudogene rate descending
    rows = []
    for chrom in chrom_total:
        total = chrom_total[chrom]
        n_pseudo = chrom_pseudo.get(chrom, 0)
        n_high = chrom_high.get(chrom, 0)
        rate = n_pseudo / total if total > 0 else 0
        rows.append((chrom, total, n_pseudo, rate, n_high))

    rows.sort(key=lambda r: (-r[3], -r[2], r[0]))

    with open(outpath, "w") as f:
        f.write("chromosome\ttotal_genes\tn_pseudo\tpseudo_rate\tn_high\n")
        for chrom, total, n_pseudo, rate, n_high in rows:
            f.write(f"{chrom}\t{total}\t{n_pseudo}\t{rate:.4f}\t{n_high}\n")

    # Log major chromosomes only (>= 50 genes)
    major_chroms = [(c, t, n, r) for c, t, n, r, _ in rows if t >= 50]
    if major_chroms:
        top = major_chroms[0]
        logger.info(
            f"Chromosomal distribution: {len(rows)} chromosomes/scaffolds. "
            f"Highest rate on major chroms: {top[0]} ({top[3]:.1%}, "
            f"{top[2]}/{top[1]} genes)"
        )
    else:
        logger.info(f"Chromosomal distribution written to {outpath}")
