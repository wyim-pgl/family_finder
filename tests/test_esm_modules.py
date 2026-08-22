"""Tests for steps/esm.py (issue #17): embedding cache + cosine tie-break,
synthetic-PDB pLDDT parsing, foldseek parsing, tier-3 assignment.

No ESM/torch/foldseek required: subprocess wrappers are monkeypatched and
parsers are fed synthetic text fixtures (house test convention).
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from config import Config  # noqa: E402
from steps import esm  # noqa: E402
from steps.esm import (  # noqa: E402
    FOLDSEEK_FORMAT_OUTPUT,
    StructHit,
    cosine,
    embedding_tiebreak,
    flag_low_plddt,
    load_embeddings,
    parse_foldseek,
    parse_plddt,
    run_esm_embed,
    run_foldseek,
    save_embeddings,
    tier3_assign,
    write_tier3_tsv,
)


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

def _atom_line(serial, name, resseq, bfactor, record="ATOM"):
    """One fixed-column PDB atom line (B-factor = pLDDT for ESMFold)."""
    return (
        f"{record:<6}{serial:>5} {name:<4} MET A{resseq:>4}    "
        f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.00:6.2f}{bfactor:6.2f}"
        f"           C\n"
    )


def _shit(gene, fam, *, tm, bits=100.0, evalue=1e-10):
    return StructHit(gene_id=gene, family_id=fam, evalue=evalue,
                     bits=bits, alntmscore=tm)


class _FakeCompleted:
    def __init__(self, returncode=0, stderr=""):
        self.returncode = returncode
        self.stderr = stderr
        self.stdout = ""


# ---------------------------------------------------------------------------
# cosine + embedding cache
# ---------------------------------------------------------------------------

def test_cosine_orthogonal_parallel_and_zero():
    assert cosine([1.0, 0.0], [0.0, 1.0]) == pytest.approx(0.0)
    assert cosine([1.0, 2.0], [2.0, 4.0]) == pytest.approx(1.0)
    assert cosine([0.0, 0.0], [1.0, 1.0]) == 0.0  # zero-norm -> uninformative


def test_cosine_rejects_dimension_mismatch():
    with pytest.raises(ValueError, match="dimension"):
        cosine([1.0], [1.0, 2.0])


def test_embeddings_cache_roundtrip(tmp_path):
    embs = {"Sp1_g1": [0.25, -1.5, 3.0], "Sp2_g2": [1e-7, 0.0, 2.5]}
    path = tmp_path / "embeddings.tsv"
    save_embeddings(embs, path)
    assert path.read_text().startswith("gene_id\tembedding\n")
    assert load_embeddings(path) == embs


# ---------------------------------------------------------------------------
# embedding_tiebreak (the #13 ambiguous-case hook)
# ---------------------------------------------------------------------------

TIE_EMBS = {
    "gene": [1.0, 0.0],
    "famA_rep1": [1.0, 0.0],   # cosine 1.0
    "famA_rep2": [1.0, 0.0],   # cosine 1.0 -> famA mean 1.0
    "famB_rep1": [1.0, 1.0],   # cosine ~0.707 -> famB mean ~0.707
}
CANDIDATES = {"famA": ["famA_rep1", "famA_rep2"], "famB": ["famB_rep1"]}


def test_tiebreak_picks_family_clearing_margin():
    assert embedding_tiebreak("gene", CANDIDATES, TIE_EMBS, min_delta=0.02) == "famA"


def test_tiebreak_returns_none_below_margin():
    # famA mean 1.0 vs famB mean ~0.707: delta ~0.293 < 0.5 -> stays ambiguous
    assert embedding_tiebreak("gene", CANDIDATES, TIE_EMBS, min_delta=0.5) is None


def test_tiebreak_exact_tie_is_never_broken():
    embs = {"gene": [1.0, 0.0], "a1": [2.0, 0.0], "b1": [3.0, 0.0]}
    cands = {"famA": ["a1"], "famB": ["b1"]}  # both cosine exactly 1.0
    assert embedding_tiebreak("gene", cands, embs, min_delta=0.02) is None


def test_tiebreak_missing_gene_or_rep_embeddings():
    # gene without an embedding -> None
    assert embedding_tiebreak("nope", CANDIDATES, TIE_EMBS, 0.02) is None
    # no representative has an embedding -> None
    cands = {"famA": ["missing1"], "famB": ["missing2"]}
    assert embedding_tiebreak("gene", cands, TIE_EMBS, 0.02) is None


def test_tiebreak_single_scorable_family_wins_unopposed():
    cands = {"famA": ["famA_rep1"], "famB": ["missing"]}
    assert embedding_tiebreak("gene", cands, TIE_EMBS, 0.02) == "famA"


# ---------------------------------------------------------------------------
# parse_plddt (synthetic PDB text)
# ---------------------------------------------------------------------------

def test_parse_plddt_means_ca_bfactors_only(tmp_path):
    pdb = tmp_path / "gene.pdb"
    pdb.write_text(
        "HEADER    ESMFOLD PREDICTION\n"
        + _atom_line(1, "N", 1, 99.99)          # non-CA: ignored
        + _atom_line(2, "CA", 1, 90.00)
        + _atom_line(3, "CA", 2, 30.00)
        + _atom_line(4, "CA", 3, 55.55, record="HETATM")  # not ATOM: ignored
        + "TER\nEND\n"
    )
    assert parse_plddt(pdb) == pytest.approx(60.0)


def test_parse_plddt_raises_without_ca(tmp_path):
    pdb = tmp_path / "empty.pdb"
    pdb.write_text("HEADER\n" + _atom_line(1, "N", 1, 80.0) + "END\n")
    with pytest.raises(ValueError, match="No CA atoms"):
        parse_plddt(pdb)


def test_flag_low_plddt_uses_config_threshold():
    config = Config()  # plddt_flag_below = 50.0
    plddts = {"g_low": 32.1, "g_high": 81.0, "g_edge": 50.0}
    assert flag_low_plddt(plddts, config) == ["g_low"]  # strict <, sorted


# ---------------------------------------------------------------------------
# parse_foldseek
# ---------------------------------------------------------------------------

def test_parse_foldseek_declared_columns_and_suffix_strip(tmp_path):
    tsv = tmp_path / "hits.tsv"
    tsv.write_text(
        # query, target, evalue, bits, alntmscore (FOLDSEEK_FORMAT_OUTPUT)
        "Sp1_g1\tR1_OG0000001.pdb\t1e-10\t250.0\t0.82\n"
        "Sp1_g1\tR2_OG0000002.pdb\t1e-12\t300.0\t0.60\n"   # better evalue/bits, worse TM
        "Sp1_g1\tR1_OG0000001.pdb\t1e-08\t200.0\t0.79\n"   # same family: keep best TM
        "malformed line without tabs\n"
        "Sp2_g9\tR3_OG0000003\tnot_a_float\t1.0\t0.5\n"    # unparseable: skipped
    )
    hits = parse_foldseek(tsv)
    assert set(hits) == {"Sp1_g1"}
    assert [h.family_id for h in hits["Sp1_g1"]] == ["R1_OG0000001", "R2_OG0000002"]
    assert hits["Sp1_g1"][0].alntmscore == 0.82  # best-TM dedupe within family
    assert hits["Sp1_g1"][0].evalue == 1e-10


# ---------------------------------------------------------------------------
# tier3_assign
# ---------------------------------------------------------------------------

def test_tier3_accepts_clear_winner():
    hits = {"g1": [_shit("g1", "famA", tm=0.80), _shit("g1", "famB", tm=0.55)]}
    assignments, ambiguous = tier3_assign(hits, Config())
    assert assignments == {"g1": "famA"}
    assert ambiguous == []


def test_tier3_single_hit_above_floor_accepted():
    hits = {"g1": [_shit("g1", "famA", tm=0.51)]}
    assignments, ambiguous = tier3_assign(hits, Config())
    assert assignments == {"g1": "famA"}


def test_tier3_rejects_below_min_tmscore():
    hits = {"g1": [_shit("g1", "famA", tm=0.49)]}  # < 0.5 default
    assignments, ambiguous = tier3_assign(hits, Config())
    assert assignments == {} and ambiguous == []  # unplaced, not ambiguous


def test_tier3_margin_failure_goes_ambiguous_with_candidates():
    hits = {"g1": [_shit("g1", "famA", tm=0.80), _shit("g1", "famB", tm=0.78)]}
    assignments, ambiguous = tier3_assign(hits, Config())  # margin 0.02 < 0.05
    assert assignments == {}
    assert len(ambiguous) == 1
    cands = ambiguous[0]["candidates"]
    assert [c["family_id"] for c in cands] == ["famA", "famB"]


def test_tier3_equal_tm_tie_is_ambiguous_even_with_zero_margin():
    config = Config(tier3_margin_tm=0.0)
    hits = {"g1": [_shit("g1", "famA", tm=0.7, bits=300.0),
                   _shit("g1", "famB", tm=0.7, bits=100.0)]}
    assignments, ambiguous = tier3_assign(hits, config)
    assert assignments == {}  # bits/e-value/order never break a TM tie
    assert len(ambiguous) == 1


def test_tier3_evalue_never_decides():
    # famB has an astronomically better e-value but a worse TM: famA wins
    hits = {"g1": [_shit("g1", "famA", tm=0.80, evalue=1e-3),
                   _shit("g1", "famB", tm=0.60, evalue=1e-100)]}
    assignments, _ = tier3_assign(hits, Config())
    assert assignments == {"g1": "famA"}


def test_write_tier3_tsv(tmp_path):
    hits = {"g1": [_shit("g1", "famA", tm=0.80)]}
    path = tmp_path / "tier3.tsv"
    write_tier3_tsv({"g1": "famA"}, hits, path)
    lines = path.read_text().splitlines()
    assert lines[0] == "gene_id\tfamily_id\talntmscore\tbits\tevalue"
    assert lines[1].startswith("g1\tfamA\t0.800\t")


# ---------------------------------------------------------------------------
# Subprocess wrappers (mocked — no ESM/foldseek installed locally)
# ---------------------------------------------------------------------------

def test_run_esm_embed_requires_configured_template(tmp_path):
    with pytest.raises(ValueError, match="not configured"):
        run_esm_embed(tmp_path / "p.fa", tmp_path / "out", Config())


def test_run_esm_embed_formats_template(tmp_path, monkeypatch):
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return _FakeCompleted(0)

    monkeypatch.setattr(esm.subprocess, "run", fake_run)
    config = Config(esm_embed_cmd="esm-extract esm2_t33_650M_UR50D {fasta} {outdir} --include mean")
    out = run_esm_embed(tmp_path / "p.fa", tmp_path / "embed", config)
    assert out == tmp_path / "embed" and out.is_dir()
    assert calls[0][0] == "esm-extract"
    assert str(tmp_path / "p.fa") in calls[0]
    assert calls[0][-2:] == ["--include", "mean"]


def test_run_foldseek_pins_format_output(tmp_path, monkeypatch):
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return _FakeCompleted(0)

    monkeypatch.setattr(esm.subprocess, "run", fake_run)
    out = run_foldseek(tmp_path / "pdbs", tmp_path / "famdb",
                       tmp_path / "out" / "hits.tsv", Config())
    assert out == tmp_path / "out" / "hits.tsv"
    cmd = calls[0]
    assert cmd[0] == "foldseek" and cmd[1] == "easy-search"
    fmt_idx = cmd.index("--format-output")
    assert cmd[fmt_idx + 1] == FOLDSEEK_FORMAT_OUTPUT


def test_wrapper_raises_on_nonzero_exit(tmp_path, monkeypatch):
    monkeypatch.setattr(
        esm.subprocess, "run",
        lambda cmd, **kw: _FakeCompleted(1, stderr="CUDA out of memory"),
    )
    config = Config(esmfold_cmd="esm-fold -i {fasta} -o {outdir}")
    with pytest.raises(RuntimeError, match="return code 1"):
        esm.run_esmfold(tmp_path / "p.fa", tmp_path / "out", config)


# ---------------------------------------------------------------------------
# Config defaults (issues #17/#18/#20 block)
# ---------------------------------------------------------------------------

def test_config_defaults_for_structural_tier():
    config = Config()
    assert config.esm_embed_cmd == "" and config.esmfold_cmd == ""
    assert config.foldseek_bin == "foldseek"
    assert config.plddt_flag_below == 50.0
    assert config.tier3_min_tmscore == 0.5
    assert config.tier3_margin_tm == 0.05
    assert config.embed_tiebreak_min_delta == 0.02
    assert config.glm_score_cmd == "" and config.ecforest_cmd == ""


# ---------------------------------------------------------------------------
# foldseek target -> family id (issue #42: one normaliser, not one per module)
# ---------------------------------------------------------------------------

def test_target_to_family_strips_every_structure_extension():
    from steps.esm import _target_to_family

    for target in ("R1_OG0000134.pdb", "R1_OG0000134.cif",
                   "R1_OG0000134.pdb.gz", "R1_OG0000134.cif.gz"):
        assert _target_to_family(target) == "R1_OG0000134"


def test_target_to_family_leaves_a_bare_family_id_alone():
    from steps.esm import _target_to_family

    assert _target_to_family("R1_OG0000134") == "R1_OG0000134"


def test_target_to_family_uses_the_shared_canonical_form():
    # A representative named from a dotted gene id reaches foldseek with the
    # dot intact; the family index is keyed on the canonical form.
    from steps.esm import _target_to_family

    assert _target_to_family("Obas__JBFLFP010000003.1_000519.pdb") == \
        "Obas__JBFLFP010000003_1_000519"
