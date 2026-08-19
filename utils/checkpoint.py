"""Checkpoint/resume logic for the iterative pipeline."""

import json
from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger("family_finder")

CHECKPOINT_FILE = "checkpoint.json"


def save_checkpoint(round_dir: Path, round_num: int, pool_size: int,
                    status: str = "started", completed_ogs: Optional[list] = None):
    """Save a checkpoint for the current round."""
    round_dir.mkdir(parents=True, exist_ok=True)
    data = {
        "round_number": round_num,
        "pool_size": pool_size,
        "status": status,
        "completed_orthogroups": completed_ogs or [],
    }
    tmp = round_dir / f"{CHECKPOINT_FILE}.tmp"
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2)
    tmp.rename(round_dir / CHECKPOINT_FILE)


def load_checkpoint(round_dir: Path) -> Optional[dict]:
    """Load checkpoint from a round directory, if it exists."""
    cp = round_dir / CHECKPOINT_FILE
    if cp.exists():
        with open(cp) as f:
            return json.load(f)
    return None


def find_latest_checkpoint(outdir: Path) -> Optional[dict]:
    """Find the newest COMPLETED checkpoint across all rounds.

    Scans backwards and skips checkpoints that are not completed (a run that
    died mid-round leaves a 'started'/'processing' checkpoint; resuming from
    it would silently restart from round 1 with the full pool — issue #15).
    """
    round_dirs = sorted(outdir.glob("round_*"))
    for rd in reversed(round_dirs):
        cp = load_checkpoint(rd)
        if cp is None:
            continue
        if cp.get("status") == "completed":
            logger.info(f"Found checkpoint: round {cp['round_number']}, status={cp['status']}")
            return cp
        logger.warning(
            f"Skipping incomplete checkpoint: round {cp.get('round_number')}, "
            f"status={cp.get('status')}"
        )
    return None
