"""Logging configuration for the family_finder pipeline."""

import logging
import sys
from pathlib import Path


def setup_logging(outdir: Path, level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger("family_finder")
    logger.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    # Console handler
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File handler
    outdir.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(outdir / "pipeline.log")
    fh.setLevel(level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger
