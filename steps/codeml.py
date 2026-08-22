"""codeml (PAML) integration: .ctl file generation and execution."""

import logging
import subprocess
from pathlib import Path
from typing import Dict, Set

from Bio import AlignIO

from config import Config, resolve_tool

logger = logging.getLogger("family_finder")

# codeml NSsites mapping for common models
MODEL_NSITES = {
    "M0": (0, 0),     # model=0, NSsites=0  (one-ratio)
    "M1a": (0, 1),    # model=0, NSsites=1  (nearly neutral)
    "M2a": (0, 2),    # model=0, NSsites=2  (positive selection)
    "M7": (0, 7),     # model=0, NSsites=7  (beta)
    "M8": (0, 8),     # model=0, NSsites=8  (beta & omega)
}

CTL_TEMPLATE = """\
      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0
       atefun = 0

        model = {model}
      NSsites = {nsites}

        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 1

    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 10

        getSE = 0
 RateAncestor = 0
   Small_Diff = .5e-6
    cleandata = 1
  fix_blength = 0
       method = 0
"""


# Branch-site Model A (issue #16): tests for positive selection on a marked
# foreground branch — the standard test for neofunctionalization after a
# retargeting event. Site models above cannot ask branch-specific questions.
BRANCH_SITE_CTL_TEMPLATE = """\
      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        clock = 0

        model = 2
      NSsites = 2

        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = {fix_omega}
        omega = {omega}

    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 10

        getSE = 0
 RateAncestor = 0
   Small_Diff = .5e-6
    cleandata = 1
  fix_blength = 0
       method = 0
"""


def fasta_to_phylip(fasta_path: Path, phylip_path: Path):
    """Convert FASTA alignment to sequential PHYLIP format for codeml."""
    alignment = AlignIO.read(str(fasta_path), "fasta")
    with open(phylip_path, "w") as f:
        AlignIO.write(alignment, f, "phylip-sequential")


def generate_ctl(
    family_id: str,
    codon_aln: Path,
    tree_path: Path,
    work_dir: Path,
    model_name: str,
) -> Path:
    """Generate a codeml control file for a given model."""
    if model_name not in MODEL_NSITES:
        raise ValueError(f"Unknown model: {model_name}. Supported: {list(MODEL_NSITES.keys())}")

    model, nsites = MODEL_NSITES[model_name]
    work_dir.mkdir(parents=True, exist_ok=True)

    # Convert to PHYLIP
    phylip_path = work_dir / "alignment.phy"
    fasta_to_phylip(codon_aln, phylip_path)

    ctl_content = CTL_TEMPLATE.format(
        seqfile=str(phylip_path),
        treefile=str(tree_path),
        outfile=str(work_dir / "results.txt"),
        model=model,
        nsites=nsites,
    )

    ctl_path = work_dir / "codeml.ctl"
    with open(ctl_path, "w") as f:
        f.write(ctl_content)

    return ctl_path


def write_marked_tree(tree_path: Path, clade_leaves, out_path: Path) -> Path:
    """Write a codeml-ready tree with the given clade tagged '#1' as the
    foreground branch. Branch lengths and internal support labels are dropped
    (PAML does not need the former and chokes on the latter).

    Raises ValueError if the leaves do not form a clade in this tree — events
    produced by steps.retargeting on the same tree always do.
    """
    from utils.newick import mark_clade, parse_newick, write_newick

    root = parse_newick(Path(tree_path).read_text())
    if not mark_clade(root, set(clade_leaves)):
        raise ValueError(
            f"Leaves {sorted(clade_leaves)} do not form a clade in {tree_path}"
        )
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        write_newick(root, with_dist=False, keep_internal_names=False) + "\n"
    )
    return out_path


DEFAULT_START_OMEGA = 1.5
# Standard multi-start set for the branch-site alternative: one value below
# neutrality, one at the historical default, one well above. Measured need:
# the PEPC SF3 fit (LRT 15.36) came from a SINGLE start and could not be
# claimed free of a local optimum until refit from all three.
START_OMEGAS = (0.5, 1.5, 3.0)


def best_lnl_run(runs):
    """Pick the (label, lnL) pair with the highest lnL from multi-start fits."""
    runs = list(runs)
    if not runs:
        raise ValueError("no codeml runs to choose from")
    return max(runs, key=lambda r: r[1])


def generate_branch_site_ctl(
    family_id: str,
    codon_aln: Path,
    marked_tree: Path,
    work_dir: Path,
    null: bool,
    start_omega: float = DEFAULT_START_OMEGA,
) -> Path:
    """Control file for branch-site Model A (null=False) or its null
    (null=True: fix_omega=1, omega=1). LRT df=1 between the pair.

    `start_omega` is the ALTERNATIVE model's starting omega. codeml's
    branch-site optimizer is prone to local optima, so the standard guard
    is to refit from several starting values (see START_OMEGAS) and keep
    the run with the highest lnL (`best_lnl_run`). The null ignores it —
    it fixes omega = 1 by definition.
    """
    work_dir.mkdir(parents=True, exist_ok=True)

    phylip_path = work_dir / "alignment.phy"
    fasta_to_phylip(codon_aln, phylip_path)

    ctl_content = BRANCH_SITE_CTL_TEMPLATE.format(
        seqfile=str(phylip_path),
        treefile=str(marked_tree),
        outfile=str(work_dir / "results.txt"),
        fix_omega=1 if null else 0,
        omega=1 if null else start_omega,
    )
    ctl_path = work_dir / "codeml.ctl"
    ctl_path.write_text(ctl_content)
    return ctl_path


def parse_lnl(results_path: Path) -> float:
    """Extract lnL from a codeml results file (line 'lnL(...): <value> ...')."""
    import re

    for line in open(results_path):
        if line.startswith("lnL"):
            m = re.search(r"lnL.*?:\s*(-?\d+\.\d+)", line)
            if m:
                return float(m.group(1))
    raise ValueError(f"No lnL line found in {results_path}")


def lrt_pvalue(lnl_alt: float, lnl_null: float) -> float:
    """LRT p-value for df=1: 2*(lnL_alt - lnL_null) ~ chi2(1), and
    P(chi2(1) > x) = erfc(sqrt(x/2)) — no scipy needed."""
    import math

    delta = max(0.0, lnl_alt - lnl_null)  # LR statistic = 2*delta
    return math.erfc(math.sqrt(delta))


def benjamini_hochberg(pvalues) -> list:
    """BH FDR q-values, preserving input order."""
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [0.0] * n
    prev = 1.0
    for rank_from_end, idx in enumerate(reversed(order)):
        rank = n - rank_from_end  # 1-based rank of this p-value
        val = min(prev, pvalues[idx] * n / rank)
        q[idx] = val
        prev = val
    return q


def run_codeml(ctl_path: Path, work_dir: Path, config: Config) -> Path:
    """Run codeml with the given control file.

    The binary is resolved through `config.resolve_tool` so a missing or
    unpinned codeml fails naming the setting rather than as a bare OSError
    from subprocess (issue #44).

    **Success is judged from the output, not from the exit status.** codeml
    writes its results and then exits 1 with "error: end of tree file"; this
    was reproduced three times and is recorded in resume.md. Treating a
    non-zero return as failure therefore discards completed analyses. A run
    counts as complete when results.txt carries an lnL line, which is what
    every downstream parser needs; anything less is a real failure and the
    stderr is raised with it.
    """
    binary = resolve_tool(config.codeml_bin, "codeml_bin")
    logger.debug(f"Running codeml: {ctl_path}")

    result = subprocess.run(
        [binary, str(ctl_path.name)],
        cwd=str(work_dir),
        capture_output=True,
        text=True,
        timeout=3600,
    )

    results_file = work_dir / "results.txt"
    if not results_file.exists():
        raise RuntimeError(
            f"codeml wrote no results.txt for {ctl_path} "
            f"(return code {result.returncode}):\n{result.stderr}"
        )

    if "lnL" not in results_file.read_text():
        raise RuntimeError(
            f"codeml results for {ctl_path} carry no lnL line, so the run did "
            f"not finish (return code {result.returncode}):\n{result.stderr}"
        )

    if result.returncode != 0:
        logger.warning(
            f"codeml exited {result.returncode} for {ctl_path} but wrote "
            f"complete results; accepting them. stderr: "
            f"{result.stderr.strip().splitlines()[-1] if result.stderr.strip() else ''}"
        )

    return results_file


def run_codeml_on_families(
    families: Dict[str, Set[str]],
    cds_pool: Dict[str, str],
    outdir: Path,
    config: Config,
):
    """Run codeml on all confirmed families for each specified model."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for family_id in sorted(families.keys()):
        # Locate the codon alignment and tree from the family's round directory
        try:
            parts = family_id.split("_", 1)
            round_num = int(parts[0][1:])
            og_id = parts[1]
        except (ValueError, IndexError):
            logger.warning(f"Skipping {family_id}: cannot parse family_id")
            continue

        base_dir = outdir.parent / f"round_{round_num:02d}" / "orthogroups" / og_id
        codon_aln = base_dir / "confirmed_codon.afa"
        tree_path = base_dir / "confirmed_tree.nwk"

        if not codon_aln.exists() or not tree_path.exists():
            logger.warning(f"Skipping {family_id}: missing codon alignment or tree")
            continue

        for model_name in config.codeml_models:
            work_dir = outdir / family_id / model_name
            try:
                ctl_path = generate_ctl(family_id, codon_aln, tree_path, work_dir, model_name)
                run_codeml(ctl_path, work_dir, config)
                logger.info(f"codeml {model_name} completed for {family_id}")
            except Exception as e:
                logger.error(f"codeml {model_name} failed for {family_id}: {e}")
