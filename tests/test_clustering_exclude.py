"""Tests for clustering_species_exclude (issue #26 v2 production, from #12).

Species listed in config.clustering_species_exclude are kept out of the
tier-1 clustering pool every round (CgigH duplicate-annotation case: its
presence distorts MCL cliques, measured co-clustering gain 9.0% when
excluded). Their genes must NOT be lost: they are merged back into the
unplaced pool after convergence so the post-run profile mapping (HMMER
rescue / profile assignment) can place them.

External tools are never invoked; partitioning is a pure function.
"""

import sys
import types
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

# steps.prune imports ete4 at module level; stub it before importing pipeline
if "ete4" not in sys.modules:
    ete4_stub = types.ModuleType("ete4")
    ete4_stub.Tree = object
    sys.modules["ete4"] = ete4_stub

import pipeline  # noqa: E402
from config import Config  # noqa: E402


POOL = {
    "Cgig_g1": "MA", "Cgig_g2": "MB",
    "CgigH_g1": "MC", "CgigH_g2": "MD",
    "Mcry_g1": "ME",
}


def test_config_default_is_empty_list():
    # Arrange / Act
    cfg = Config()
    # Assert
    assert cfg.clustering_species_exclude == []


def test_partition_excludes_listed_species():
    # Arrange / Act
    kept, excluded = pipeline.partition_excluded_species(POOL, ["CgigH"], "_")
    # Assert
    assert set(excluded) == {"CgigH_g1", "CgigH_g2"}
    assert set(kept) == {"Cgig_g1", "Cgig_g2", "Mcry_g1"}
    assert excluded["CgigH_g1"] == "MC"  # sequences travel with the ids


def test_partition_prefix_match_is_exact():
    # Excluding Cgig must NOT drag CgigH along (distinct species prefixes).
    kept, excluded = pipeline.partition_excluded_species(POOL, ["Cgig"], "_")
    assert set(excluded) == {"Cgig_g1", "Cgig_g2"}
    assert "CgigH_g1" in kept and "CgigH_g2" in kept


def test_partition_empty_exclude_is_noop():
    kept, excluded = pipeline.partition_excluded_species(POOL, [], "_")
    assert kept == POOL
    assert excluded == {}


def test_partition_unknown_species_warns_and_keeps_everything(caplog):
    import logging
    with caplog.at_level(logging.WARNING, logger="family_finder"):
        kept, excluded = pipeline.partition_excluded_species(
            POOL, ["Nope"], "_"
        )
    assert kept == POOL
    assert excluded == {}
    assert any("Nope" in r.message for r in caplog.records)


def test_partition_does_not_mutate_input():
    before = dict(POOL)
    pipeline.partition_excluded_species(POOL, ["CgigH"], "_")
    assert POOL == before
