"""codeml (PAML) integration: .ctl file generation and execution."""

import logging
import subprocess
from pathlib import Path
from typing import Dict, Set

from Bio import AlignIO

from config import Config
from utils.seqio import write_fasta

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


def run_codeml(ctl_path: Path, work_dir: Path, config: Config) -> Path:
    """Run codeml with the given control file."""
    logger.debug(f"Running codeml: {ctl_path}")

    result = subprocess.run(
        [config.codeml_bin, str(ctl_path.name)],
        cwd=str(work_dir),
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        logger.error(f"codeml failed for {ctl_path}:\n{result.stderr}")
        raise RuntimeError(f"codeml failed with return code {result.returncode}")

    results_file = work_dir / "results.txt"
    if not results_file.exists():
        raise FileNotFoundError(f"codeml output not found: {results_file}")

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
        parts = family_id.split("_", 1)
        round_num = int(parts[0][1:])
        og_id = parts[1]

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
