"""Nominate fragmented families once, after convergence (issue #35).

`detect_merge_candidates` has existed in steps/profile_assign.py since #13 and
**nothing in the pipeline ever called it** - the only callers were its own
tests. So the 15-species production run split the PEPC clan into four families
(R1_OG0000440 / R2_OG0000359 / R1_OG0008467 / R1_OG0009826, the flagship CAM
pairs SPLIT 2 of 2) and no output said a word about it.

Placement is not the place to catch this. Rescue only searches genes that were
left unplaced, and a fragmented family's members are all placed, so the
fragmentation is invisible there by construction - measured: none of the four
PEPC families appears among the 1,130 thin-margin family pairs.

The scan runs once, after convergence and rescue, on the final family set, and
only NOMINATES. detect_merge_candidates' own docstring is explicit that a
candidate still needs tree validation before an actual merge, and merging two
families that a codon tree would keep apart is worse than reporting them.
"""
import sys
import types
from pathlib import Path

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import Config
from pipeline import write_merge_candidates


def _rows(tmp_path):
    return (tmp_path / "merge_candidates.tsv").read_text().splitlines()


def test_candidates_are_written_with_their_reciprocal_fraction(tmp_path):
    write_merge_candidates([("R1_OG0000440", "R2_OG0000359", 0.83)],
                           {"R1_OG0000440": {"a"} , "R2_OG0000359": {"b"}},
                           tmp_path)
    rows = _rows(tmp_path)
    assert rows[0].split("\t") == ["family_a", "family_b", "min_reciprocal",
                                   "n_genes_a", "n_genes_b"]
    assert rows[1].startswith("R1_OG0000440\tR2_OG0000359\t0.83")


def test_an_empty_result_still_writes_the_file(tmp_path):
    """A missing file is indistinguishable from a scan that never ran, and
    'no fragmentation detected' is a result worth recording."""
    write_merge_candidates([], {}, tmp_path)
    assert _rows(tmp_path) == ["family_a\tfamily_b\tmin_reciprocal\t"
                               "n_genes_a\tn_genes_b"]


def test_family_sizes_travel_with_the_pair(tmp_path):
    """A 63-member family paired with a 9-member splinter is a different claim
    from two 40-member families, and the reciprocal fraction alone hides it."""
    write_merge_candidates([("A", "B", 0.7)],
                           {"A": set("abcdefghi"), "B": {"x", "y"}}, tmp_path)
    assert _rows(tmp_path)[1].endswith("\t9\t2")


def test_nothing_is_merged_automatically(tmp_path):
    """The scan nominates. Merging without tree validation is how two families
    a codon tree would keep apart become one."""
    families = {"A": {"a"}, "B": {"b"}}
    before = {k: set(v) for k, v in families.items()}
    write_merge_candidates([("A", "B", 0.9)], families, tmp_path)
    assert families == before
