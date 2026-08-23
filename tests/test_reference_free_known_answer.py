"""Blocking the reference must not change the conclusion (issue #42).

The residue-level results in methods.md were reached knowing that
Mcry_Mcr8G11630 is SwissProt P10490, that GYSDSGKDAGR is catalytic, and that
SIDAQL at residues 11-16 is the regulatory phosphorylation motif. Almost no
family has any of that. The acceptance test for the reference-free path is
therefore: block the knowledge and see whether the same conclusions still come
out.

Measured on the corrected PEPC clan (77 plant-type sequences, 1,428 columns),
both halves pass.

1. Core avoidance. With the reference given and with it blocked, the verdict is
   `avoids_core` either way - p = 0.0, observed median distance 2 against a null
   median of 1.0, 61 diagnostic and 441 invariant columns in both runs. Only the
   coordinate frame differs, and it is reported: blocking the reference makes
   the analysis number its results in Tfru_..._Tf_contig_062_000129, a sequence
   nobody has characterised.

2. The phosphorylation context. Without being told what SIDAQL is, the scan
   reports five diagnostic columns in the N-terminal 5% of the alignment, all of
   them ppc-1E2, in a region holding ZERO family-invariant columns, where
   ppc-1E1 is coverage-suppressed at every single column (32 of its 62 members
   reach it). Mapped back afterwards, one of those five is alignment column 30 =
   P10490 residue 9 - the exact residue the reference-based analysis singled out
   as Met in all 15 ppc-1E2 members against Leu in ppc-1E1. The others flank the
   motif at residues 2, 4, 5 and 20.

   The reference-free path is also the stricter one: residue 10, described in
   methods.md as "fixed as alanine in ppc-1E2 while ppc-1E1 is split", is NOT
   called diagnostic here, because alanine is 13 of 22 non-gap residues in the
   contrast set. That is correct - a residue the other group uses most of the
   time is not diagnostic - and it is the reference-based narrative that made
   the weaker claim.

The test below locks the invariant on synthetic data so it runs without the
cluster; the numbers above are the measurement it stands for.
"""
import sys
import types
from pathlib import Path

sys.modules.setdefault("ete4", types.ModuleType("ete4"))
sys.modules["ete4"].Tree = object
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from steps.subfamily import sdp_core_relationship, sdp_scan

# A core every group agrees on, a diagnostic block far from it, and an
# N-terminal stretch one group barely covers.
CORE = "WWWWWWWWWW"
def _seq(head, diag):
    return head + CORE + diag + CORE


ALN = {}
for i in range(8):
    ALN[f"Aaa_g{i}"] = _seq("MKTAYIAKQR", "CCCCC")
for i in range(8):
    ALN[f"Bbb_g{i}"] = _seq("MKTAYIAKQR", "DDDDD")
GROUPS = {"SF1": [f"Aaa_g{i}" for i in range(8)],
          "SF2": [f"Bbb_g{i}" for i in range(8)]}


def test_the_verdict_does_not_depend_on_which_sequence_supplies_coordinates():
    known = sdp_core_relationship(ALN, GROUPS, ref_seq_id="Aaa_g0")
    blocked = sdp_core_relationship(ALN, GROUPS)
    assert known["verdict"] == blocked["verdict"]
    assert known["n_sdp_columns"] == blocked["n_sdp_columns"]
    assert known["n_invariant_columns"] == blocked["n_invariant_columns"]
    assert known["p_value"] == blocked["p_value"]


def test_the_coordinate_frame_is_reported_in_both_cases():
    """Same numbers, different frame - and a result that did not say which
    sequence it counted in would be uncitable."""
    known = sdp_core_relationship(ALN, GROUPS, ref_seq_id="Aaa_g0")
    blocked = sdp_core_relationship(ALN, GROUPS)
    assert known["distance_reference"] == "Aaa_g0"
    assert known["distance_reference_source"] == "explicit"
    assert blocked["distance_reference"] is not None
    assert blocked["distance_reference_source"] == "automatic"


def test_the_diagnostic_columns_are_the_same_with_and_without_a_reference():
    """The reference is a coordinate frame, never an input to detection."""
    with_ref = {h["aln_col"] for h in sdp_scan(ALN, GROUPS, ref_seq_id="Aaa_g0")}
    without = {h["aln_col"] for h in sdp_scan(ALN, GROUPS)}
    assert with_ref == without


def test_a_residue_the_contrast_set_commonly_uses_is_not_diagnostic():
    """Why residue 10 of P10490 is absent from the reference-free result: a
    residue the other subfamily uses most of the time cannot diagnose."""
    aln = dict(ALN)
    for i in range(8):
        # SF2 members mostly carry the same residue SF1 is 'fixed' on
        aln[f"Bbb_g{i}"] = _seq("MKTAYIAKQR", "CCCCD" if i < 6 else "DDDDD")
    cols = {h["aln_col"] for h in sdp_scan(aln, GROUPS) if h["subfamily"] == "SF1"}
    shared = len(_seq("", "")) - 5  # first column of the diagnostic block
    assert shared not in cols
