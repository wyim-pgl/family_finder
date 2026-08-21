"""Tests for the neofunctionalization module (issue #16): newick/Fitch,
DeepLoc parsing, retargeting detection, branch-site helpers."""

import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from steps.codeml import (  # noqa: E402
    benjamini_hochberg,
    generate_branch_site_ctl,
    lrt_pvalue,
    parse_lnl,
    write_marked_tree,
)
from steps.deeploc import parse_deeploc, save_cache, load_cache, top_location  # noqa: E402
from steps.retargeting import find_retargeting_events  # noqa: E402
from utils.newick import fitch, mark_clade, parse_newick, write_newick  # noqa: E402

# FastTree-style tree: internal support labels + branch lengths
NWK = "((Mcry_g1:0.5,Cgig_g1:0.1)0.99:0.2,(Ococ_g1:0.1,Obas_g1:0.1)0.95:0.2);"


# ---------- newick ----------

def test_parse_newick_leaves_and_dists():
    root = parse_newick(NWK)
    assert root.leaf_names() == {"Mcry_g1", "Cgig_g1", "Ococ_g1", "Obas_g1"}
    mcry = next(l for l in root.leaves() if l.name == "Mcry_g1")
    assert mcry.dist == 0.5


def test_write_newick_codeml_mode_drops_support_and_dist():
    root = parse_newick(NWK)
    out = write_newick(root, with_dist=False, keep_internal_names=False)
    assert ":" not in out and "0.99" not in out
    assert out.endswith(";")


def test_mark_clade_emits_hash1():
    root = parse_newick(NWK)
    assert mark_clade(root, {"Ococ_g1", "Obas_g1"})
    out = write_newick(root, with_dist=False, keep_internal_names=False)
    assert "(Ococ_g1,Obas_g1) #1" in out


def test_mark_clade_returns_false_for_non_clade():
    root = parse_newick(NWK)
    assert not mark_clade(root, {"Mcry_g1", "Ococ_g1"})


# ---------- Fitch ----------

def test_fitch_finds_single_switch():
    # Arrange: chloroplast pair vs cytoplasm pair -> one switch branch
    root = parse_newick("((a:1,b:1):1,(c:1,d:1):1);")
    states = {"a": "Plastid", "b": "Plastid", "c": "Cytoplasm", "d": "Cytoplasm"}
    # Act
    fitch(root, states)
    from utils.newick import state_switches
    events = state_switches(root)
    # Assert: exactly one switch, and its clade is one of the two pairs
    assert len(events) == 1
    _, _, clade = events[0]
    assert clade in ({"a", "b"}, {"c", "d"})


def test_fitch_uninformative_leaf_does_not_create_event():
    root = parse_newick("((a:1,b:1):1,(c:1,d:1):1);")
    states = {"a": "Plastid", "b": "Plastid", "c": "Plastid", "d": None}
    fitch(root, states)
    from utils.newick import state_switches
    assert state_switches(root) == []


# ---------- DeepLoc ----------

def test_parse_deeploc_and_top_location(tmp_path):
    csv = tmp_path / "results.csv"
    csv.write_text(
        "Protein_ID,Localizations,Signals,Cytoplasm,Plastid,Mitochondrion\n"
        "g1,Plastid,Transit peptide,0.05,0.91,0.04\n"
        "g2,Cytoplasm|Plastid,,0.45,0.44,0.11\n"
    )
    preds = parse_deeploc(csv)
    assert preds["g1"]["probs"]["Plastid"] == 0.91
    assert top_location(preds["g1"]) == "Plastid"
    # dual-targeting split probabilities -> uninformative, not forced
    assert top_location(preds["g2"], min_prob=0.5) is None
    assert top_location(None) is None


def test_deeploc_cache_roundtrip(tmp_path):
    preds = {"g1": {"label": "Plastid", "probs": {"Plastid": 0.9, "Cytoplasm": 0.1}}}
    p = tmp_path / "cache.tsv"
    save_cache(preds, p)
    assert load_cache(p) == preds


# ---------- retargeting ----------

PREDS = {
    "Mcry_g1": {"label": "Cytoplasm", "probs": {"Cytoplasm": 0.9, "Plastid": 0.1}},
    "Cgig_g1": {"label": "Plastid", "probs": {"Plastid": 0.9, "Cytoplasm": 0.1}},
    "Ococ_g1": {"label": "Plastid", "probs": {"Plastid": 0.9, "Cytoplasm": 0.1}},
    "Obas_g1": {"label": "Plastid", "probs": {"Plastid": 0.9, "Cytoplasm": 0.1}},
}


def test_find_retargeting_events_detects_mcry_switch():
    events = find_retargeting_events("R1_OG1", NWK, PREDS)
    assert len(events) == 1
    ev = events[0]
    assert ev.child_state == "Cytoplasm" and ev.clade_leaves == frozenset({"Mcry_g1"})


def test_truncation_exclude_suppresses_fake_event():
    # The safeguard: excluding the (possibly truncated) switching gene
    # removes the event instead of reporting an artifact
    events = find_retargeting_events("R1_OG1", NWK, PREDS, exclude={"Mcry_g1"})
    assert events == []


# ---------- branch-site helpers ----------

def test_write_marked_tree_and_ctl(tmp_path):
    tree = tmp_path / "t.nwk"
    tree.write_text(NWK)
    marked = write_marked_tree(tree, {"Mcry_g1"}, tmp_path / "marked.nwk")
    assert "Mcry_g1 #1" in marked.read_text()

    aln = tmp_path / "codon.afa"
    aln.write_text(">Mcry_g1\nATGGCT\n>Cgig_g1\nATGGCA\n>Ococ_g1\nATGGCC\n>Obas_g1\nATGGCG\n")
    alt = generate_branch_site_ctl("F", aln, marked, tmp_path / "alt", null=False)
    null = generate_branch_site_ctl("F", aln, marked, tmp_path / "null", null=True)
    alt_txt, null_txt = alt.read_text(), null.read_text()
    for txt in (alt_txt, null_txt):
        assert "model = 2" in txt and "NSsites = 2" in txt
    assert "fix_omega = 0" in alt_txt
    assert "fix_omega = 1" in null_txt and "omega = 1" in null_txt


def test_write_marked_tree_rejects_non_clade(tmp_path):
    tree = tmp_path / "t.nwk"
    tree.write_text(NWK)
    with pytest.raises(ValueError, match="clade"):
        write_marked_tree(tree, {"Mcry_g1", "Ococ_g1"}, tmp_path / "m.nwk")


def test_parse_lnl(tmp_path):
    mlc = tmp_path / "results.txt"
    mlc.write_text("stuff\nlnL(ntime: 25  np: 27):  -2345.678901   +0.000000\n")
    assert parse_lnl(mlc) == -2345.678901


def test_lrt_pvalue_chi2_df1():
    # 2*delta = 3.84 is the classic chi2(1) 5% critical value
    p = lrt_pvalue(lnl_alt=-100.0, lnl_null=-101.92)
    assert abs(p - 0.05) < 0.001
    assert lrt_pvalue(-100.0, -100.0) == 1.0
    assert lrt_pvalue(-100.0, -99.0) == 1.0  # alt worse than null -> LR clamped


def test_benjamini_hochberg():
    q = benjamini_hochberg([0.01, 0.04, 0.03, 0.005])
    assert q[3] == pytest.approx(0.02)   # 0.005*4/1
    assert q[0] == pytest.approx(0.02)   # 0.01*4/2, monotone
    assert q[2] == pytest.approx(0.04)   # 0.03*4/3
    assert q[1] == pytest.approx(0.04)   # 0.04*4/4
    assert benjamini_hochberg([]) == []


# --- multiple starting omega (branch-site local-optima guard) -------------

def test_branch_site_ctl_accepts_start_omega(tmp_path):
    """codeml branch-site is prone to local optima; the standard guard is to
    refit from several starting omegas and keep the best lnL."""
    aln = tmp_path / "codon.fa"
    aln.write_text(">a\nATGATG\n>b\nATGATG\n")
    marked = tmp_path / "t.nwk"
    marked.write_text("(a #1,b);")
    ctl = generate_branch_site_ctl("F", aln, marked, tmp_path / "w3",
                                   null=False, start_omega=3.0)
    text = ctl.read_text()
    assert "omega = 3" in text.replace("omega = 3.0", "omega = 3")
    assert "fix_omega = 0" in text


def test_branch_site_null_ignores_start_omega(tmp_path):
    """The null fixes omega=1 by definition — a start value must not leak in."""
    aln = tmp_path / "codon.fa"
    aln.write_text(">a\nATGATG\n>b\nATGATG\n")
    marked = tmp_path / "t.nwk"
    marked.write_text("(a #1,b);")
    ctl = generate_branch_site_ctl("F", aln, marked, tmp_path / "wn",
                                   null=True, start_omega=3.0)
    text = ctl.read_text()
    assert "fix_omega = 1" in text
    assert "omega = 1" in text


def test_branch_site_default_start_omega_unchanged(tmp_path):
    """Backward compatibility: omitting start_omega keeps the historical 1.5."""
    aln = tmp_path / "codon.fa"
    aln.write_text(">a\nATGATG\n>b\nATGATG\n")
    marked = tmp_path / "t.nwk"
    marked.write_text("(a #1,b);")
    ctl = generate_branch_site_ctl("F", aln, marked, tmp_path / "wd", null=False)
    assert "omega = 1.5" in ctl.read_text()


def test_best_lnl_run_picks_maximum():
    from steps.codeml import best_lnl_run
    runs = [("w0.5", -100.5), ("w1.5", -99.2), ("w3", -101.0)]
    assert best_lnl_run(runs) == ("w1.5", -99.2)


def test_best_lnl_run_rejects_empty():
    from steps.codeml import best_lnl_run
    try:
        best_lnl_run([])
    except ValueError:
        return
    raise AssertionError("empty run list must raise")
