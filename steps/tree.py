"""Tree building wrappers: FastTree and IQ-TREE."""

import logging
import shutil
import subprocess
from pathlib import Path

from config import Config

logger = logging.getLogger("family_finder")


def build_tree(alignment: Path, outpath: Path, config: Config, nucleotide: bool = True) -> Path:
    """Build a phylogenetic tree from an alignment.

    Args:
        nucleotide: If True, use nucleotide model (-nt -gtr -gamma).
                    If False, use protein model.

    Returns path to the Newick tree file.
    """
    alignment = Path(alignment)
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    # Auto-detect: if alignment has protein extension, use protein mode
    if alignment.name.endswith("proteins.afa"):
        nucleotide = False

    if config.tree_builder == "fasttree":
        return _run_fasttree(alignment, outpath, config, nucleotide)
    elif config.tree_builder == "iqtree":
        return _run_iqtree(alignment, outpath, config, nucleotide)
    else:
        raise ValueError(f"Unknown tree builder: {config.tree_builder}")


def _run_fasttree(alignment: Path, outpath: Path, config: Config, nucleotide: bool = True) -> Path:
    if nucleotide:
        cmd = [config.fasttree_bin, "-nt", "-gtr", "-gamma", str(alignment)]
    else:
        cmd = [config.fasttree_bin, "-gamma", str(alignment)]

    logger.debug(f"Running FastTree: {' '.join(cmd)}")

    tmp_out = outpath.parent / f".{outpath.name}.tmp"
    try:
        with open(tmp_out, "w") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, timeout=1800)

        if result.returncode != 0:
            logger.error(f"FastTree failed:\n{result.stderr}")
            raise RuntimeError(f"FastTree failed with return code {result.returncode}")

        shutil.move(str(tmp_out), str(outpath))
    except BaseException:
        tmp_out.unlink(missing_ok=True)
        raise

    return outpath


def _run_iqtree(alignment: Path, outpath: Path, config: Config, nucleotide: bool = True) -> Path:
    prefix = outpath.parent / outpath.stem
    model = "GTR+G" if nucleotide else "LG+G"
    cmd = [
        config.iqtree_bin,
        "-s", str(alignment),
        "-m", model,
        "-bb", "1000",
        "-nt", "AUTO",
        "--prefix", str(prefix),
    ]

    logger.debug(f"Running IQ-TREE: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)

    if result.returncode != 0:
        logger.error(f"IQ-TREE failed:\n{result.stderr}")
        raise RuntimeError(f"IQ-TREE failed with return code {result.returncode}")

    # IQ-TREE outputs .treefile
    treefile = Path(f"{prefix}.treefile")
    if treefile.exists():
        shutil.move(str(treefile), str(outpath))
    else:
        raise FileNotFoundError(f"IQ-TREE tree file not found at {treefile}")

    return outpath
