"""Tree building wrappers: FastTree and IQ-TREE."""

import logging
import subprocess
from pathlib import Path

from config import Config

logger = logging.getLogger("family_finder")


def build_tree(alignment: Path, outpath: Path, config: Config) -> Path:
    """Build a phylogenetic tree from a protein alignment.

    Returns path to the Newick tree file.
    """
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    if config.tree_builder == "fasttree":
        return _run_fasttree(alignment, outpath, config)
    elif config.tree_builder == "iqtree":
        return _run_iqtree(alignment, outpath, config)
    else:
        raise ValueError(f"Unknown tree builder: {config.tree_builder}")


def _run_fasttree(alignment: Path, outpath: Path, config: Config) -> Path:
    cmd = [config.fasttree_bin, "-gamma", str(alignment)]

    logger.debug(f"Running FastTree: {' '.join(cmd)}")

    with open(outpath, "w") as out_f:
        result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logger.error(f"FastTree failed:\n{result.stderr}")
        raise RuntimeError(f"FastTree failed with return code {result.returncode}")

    return outpath


def _run_iqtree(alignment: Path, outpath: Path, config: Config) -> Path:
    prefix = outpath.parent / outpath.stem
    cmd = [
        config.iqtree_bin,
        "-s", str(alignment),
        "-m", "TEST",
        "-bb", "1000",
        "-nt", "AUTO",
        "--prefix", str(prefix),
    ]

    logger.debug(f"Running IQ-TREE: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"IQ-TREE failed:\n{result.stderr}")
        raise RuntimeError(f"IQ-TREE failed with return code {result.returncode}")

    # IQ-TREE outputs .treefile
    treefile = Path(f"{prefix}.treefile")
    if treefile.exists():
        treefile.rename(outpath)
    else:
        raise FileNotFoundError(f"IQ-TREE tree file not found at {treefile}")

    return outpath
