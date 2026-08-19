"""Plant genomic-LM wrapper: splice/isoform plausibility scoring (issue #18).

Backend-agnostic by design: the cluster picks PlantCaduceus, AgroNT, or Evo 2
after the small benchmark on known plant splice sites (AlphaGenome was
REJECTED — human/mouse-trained, unvalidated on plants). Each backend needs a
thin cluster-side adapter script whose ONLY contract is the interchange TSV:

    seq_id<TAB>score          (header line, one row per scored sequence)

written as {outdir}/glm_scores.tsv. Higher score = more plausible sequence
under the model. This module stays deliberately thin — research scaffolding;
the analysis design (splicing plausibility -> localization change ->
selection signal) lives in issues #18/#3.

Downstream: low-plausibility isoforms are ANNOTATION-ARTIFACT candidates,
filtered before DeepLoc isoform-level localization (issue #16 cache) — an
annotation layer, never a clustering criterion.
"""

import logging
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger("family_finder")

# The interchange filename every backend adapter must produce in {outdir}.
GLM_SCORES_TSV = "glm_scores.tsv"


def run_glm_scores(fasta: Path, outdir: Path, config) -> Path:
    """Run the configured gLM adapter on a FASTA. Returns the scores TSV path.

    Driven by config.glm_score_cmd, e.g.
    "python score_plantcaduceus.py --fasta {fasta} --outdir {outdir}".
    The adapter must write {outdir}/glm_scores.tsv in the interchange format.
    """
    if not config.glm_score_cmd.strip():
        raise ValueError("plant gLM: glm_score_cmd not configured (see config.py)")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = shlex.split(
        config.glm_score_cmd.format(fasta=str(fasta), outdir=str(outdir))
    )
    logger.info(f"Running plant gLM scoring: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"plant gLM scoring failed:\n{result.stderr}")
        raise RuntimeError(
            f"plant gLM scoring failed with return code {result.returncode}"
        )
    scores_tsv = outdir / GLM_SCORES_TSV
    if not scores_tsv.exists():
        raise FileNotFoundError(
            f"gLM adapter did not produce the interchange TSV: {scores_tsv}"
        )
    return scores_tsv


def parse_glm_scores(path: Path) -> Dict[str, float]:
    """Parse the interchange TSV (seq_id<TAB>score, header line) -> scores.

    Tolerates a missing header (a row whose second column parses as a float
    is data). Malformed rows are skipped with a debug log, never guessed.
    """
    scores: Dict[str, float] = {}
    with open(path) as f:
        for i, line in enumerate(f):
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                logger.debug(f"Skipping malformed gLM score line: {line[:80]!r}")
                continue
            try:
                scores[fields[0]] = float(fields[1])
            except ValueError:
                if i == 0:
                    continue  # header line
                logger.debug(f"Skipping malformed gLM score line: {line[:80]!r}")
    return scores


def splice_plausibility_report(
    isoform_scores: Dict[str, float],
    gene_isoform_map: Dict[str, List[str]],
    min_gap: float,
) -> List[dict]:
    """Rank each gene's isoforms by gLM score and flag low-plausibility ones.

    An isoform scoring >= min_gap BELOW its gene's best isoform is flagged
    as a low-plausibility annotation-artifact candidate (the filter feeding
    issue #3's isoform-level DeepLoc analysis). Genes with a single scored
    isoform produce one unflagged row (no within-gene comparison possible).
    Isoforms without a score are skipped (unscored, not artifact).

    Returns rows sorted by (gene, rank):
        {"gene", "isoform", "rank", "score", "best_score", "gap",
         "low_plausibility"}
    """
    rows: List[dict] = []
    for gene in sorted(gene_isoform_map):
        scored = [
            (iso, isoform_scores[iso])
            for iso in gene_isoform_map[gene]
            if iso in isoform_scores
        ]
        if not scored:
            continue
        scored.sort(key=lambda kv: (-kv[1], kv[0]))
        best_score = scored[0][1]
        for rank, (iso, score) in enumerate(scored, start=1):
            gap = best_score - score
            rows.append({
                "gene": gene,
                "isoform": iso,
                "rank": rank,
                "score": score,
                "best_score": best_score,
                "gap": gap,
                "low_plausibility": gap >= min_gap,
            })
    return rows


def write_report_tsv(rows: List[dict], path: Path) -> None:
    with open(path, "w") as f:
        f.write("gene\tisoform\trank\tscore\tbest_score\tgap\tlow_plausibility\n")
        for r in rows:
            f.write(
                f"{r['gene']}\t{r['isoform']}\t{r['rank']}\t{r['score']:.4f}\t"
                f"{r['best_score']:.4f}\t{r['gap']:.4f}\t"
                f"{int(r['low_plausibility'])}\n"
            )
    logger.debug(f"Splice plausibility report written: {path}")
