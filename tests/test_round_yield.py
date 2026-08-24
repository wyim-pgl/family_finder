"""Stop iterating when a round stops paying for itself.

Issue #46. Convergence was "N consecutive rounds adding zero new families",
which is a yield threshold of exactly zero: a round that places four genes out
of thirty-five thousand keeps the loop alive. Measured on the five-species v2
run, rounds 4-6 together placed 0.26% of the genes, and on the fifteen-species
run one of those rounds spent 2h11m in profile assignment to place none.

The zero-threshold behaviour is still reachable (yield 0.0), so nothing changes
for a caller that does not ask for it.
"""
import sys
import types
from pathlib import Path

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import Config
from pipeline import (_round_placed_gene_ids, round_yield,
                      should_skip_profile_tier, should_stop_iterating)


def test_yield_is_genes_placed_over_the_pool_it_started_from():
    """Family count is the wrong unit: one family can hold forty genes or two."""
    assert round_yield(n_genes_placed=100, pool_size=1000) == 0.1


def test_an_empty_pool_yields_nothing_rather_than_dividing_by_zero():
    assert round_yield(n_genes_placed=0, pool_size=0) == 0.0


def test_yield_includes_profile_assignment_into_an_older_family():
    # The assigned gene is removed from the round's outlier pool but is not a
    # member of a family created in this round.  It still counts as placement.
    placed = _round_placed_gene_ids(
        {"R2_OG1": {"Sp1_new", "Sp2_new"}}, {"Sp3_assigned_to_R1"},
    )

    assert placed == {"Sp1_new", "Sp2_new", "Sp3_assigned_to_R1"}
    assert round_yield(len(placed), 10) == 0.3
    assert should_stop_iterating([0.0, round_yield(len(placed), 10)], Config()) is False


def test_a_round_below_the_threshold_stops_the_loop():
    cfg = Config()
    cfg.convergence_min_yield = 0.01          # 1% of the pool
    assert should_stop_iterating([0.004], cfg) is True


def test_a_productive_round_keeps_going():
    cfg = Config()
    cfg.convergence_min_yield = 0.01
    assert should_stop_iterating([0.038], cfg) is False


def test_the_default_reproduces_the_old_zero_threshold():
    """Default config must not change any existing run's round count."""
    cfg = Config()
    assert cfg.convergence_min_yield == 0.0
    assert should_stop_iterating([0.0001], cfg) is False   # tiny but non-zero
    assert should_stop_iterating([0.0, 0.0], cfg) is True  # two empty rounds


def test_one_empty_round_is_not_enough_at_the_default():
    assert should_stop_iterating([0.0], Config()) is False


def test_the_five_species_run_would_have_stopped_at_round_four():
    """Measured yields: R1 .953 R2 .038 R3 .007 R4 .002 R5 .001 R6 .00003.

    At 0.1% the loop ends after R4 - the two rounds that follow placed 113 of
    134,175 genes between them.
    """
    cfg = Config()
    cfg.convergence_min_yield = 0.001
    assert should_stop_iterating([0.953], cfg) is False
    assert should_stop_iterating([0.038], cfg) is False
    assert should_stop_iterating([0.007], cfg) is False
    assert should_stop_iterating([0.002], cfg) is False
    assert should_stop_iterating([0.001], cfg) is False     # exactly at the bar
    assert should_stop_iterating([0.0003], cfg) is True


# --- tier-2 auto-off -------------------------------------------------------

def test_the_profile_tier_switches_off_after_two_barren_rounds():
    """R5 spent 2h11m rebuilding 23,744 profiles to place zero genes, having
    placed four the round before."""
    cfg = Config()
    assert should_skip_profile_tier([4, 0], cfg) is False
    assert should_skip_profile_tier([0, 0], cfg) is True


def test_a_single_barren_round_does_not_switch_it_off():
    assert should_skip_profile_tier([0], Config()) is False


def test_one_placement_resets_the_count():
    assert should_skip_profile_tier([0, 0, 1], Config()) is False


def test_the_auto_off_can_be_disabled():
    cfg = Config()
    cfg.profile_assign_off_after_barren = 0
    assert should_skip_profile_tier([0, 0, 0, 0], cfg) is False
