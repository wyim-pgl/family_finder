"""steps/hmmer_rescue.rescue_unplaced() — the post-convergence rescue tier.

These lock the behaviour that the refactor into _build_profile_db /
_search_unplaced / _group_by_family / _absorb_rescued had to preserve
byte-for-byte: which genes end up in which family, what lands in
rescue_summary.tsv, and — just as important — every path on which the
rescue must hand `families` back UNCHANGED rather than a rebuilt copy.

No HMMER binary is ever invoked: `subprocess` is swapped for a fake that
emulates hmmbuild/hmmpress/hmmsearch by writing the files those tools
would write, so _concat_hmms, the tblout parser, merge_tblouts and the
summary writer all run for real.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from config import Config
import steps.hmmer_rescue as hr


# --- fake externals ------------------------------------------------------

class _Result:
    def __init__(self, returncode=0, stderr=""):
        self.returncode = returncode
        self.stdout = ""
        self.stderr = stderr


class FakeSubprocess:
    """Emulates the three HMMER binaries the rescue shells out to."""

    def __init__(self, tblout_body="", chunk_bodies=(), dom_chunk_bodies=(),
                 press_rc=0,
                 build_fail=(), build_silent=False):
        self.tblout_body = tblout_body
        self.chunk_bodies = list(chunk_bodies)
        self.dom_chunk_bodies = list(dom_chunk_bodies)
        self.press_rc = press_rc
        self.build_fail = set(build_fail)
        self.build_silent = build_silent
        self.calls = []

    def run(self, cmd, capture_output=False, text=False, **kw):
        self.calls.append(list(cmd))
        exe = Path(cmd[0]).name

        if exe == "hmmbuild":
            family_id, hmm_path = cmd[3], cmd[4]
            if family_id in self.build_fail:
                return _Result(1, stderr="hmmbuild: bad alignment\n")
            if not self.build_silent:
                Path(hmm_path).write_text(
                    f"HMMER3/f [3.2.1]\nNAME  {family_id}\nLENG  12\nHMM\n//\n")
            return _Result(0)

        if exe == "hmmpress":
            return _Result(self.press_rc, stderr="hmmpress: no such file\n")

        if exe == "hmmsearch":
            Path(cmd[cmd.index("--tblout") + 1]).write_text(self.tblout_body)
            return _Result(0)

        if exe == "bash":                       # the generated chunk runner
            out_dir = Path(cmd[1]).parent / "chunk_out"
            out_dir.mkdir(parents=True, exist_ok=True)
            for i, body in enumerate(self.chunk_bodies):
                (out_dir / f"chunk_{i:04d}.tblout").write_text(body)
            for i, body in enumerate(self.dom_chunk_bodies):
                (out_dir / f"chunk_{i:04d}.domtblout").write_text(body)
            return _Result(0)

        raise AssertionError(f"unexpected command: {cmd}")

    def ran(self, exe):
        return [c for c in self.calls if Path(c[0]).name == exe]


AFA = ">Aaa_g1\nMKTAYIAKQR\n>Bbb_g2\nMKTAYIAKQR\n"

PROT = {"Aaa_u1": "MKTAYIAKQRQ", "Bbb_u2": "MKTAYIAKQRW",
        "Ccc_u3": "MKTAYIAKQRY", "Aaa_g1": "MKTAYIAKQR",
        "Bbb_g2": "MKTAYIAKQR"}
CDS = {k: "ATG" * len(v) for k, v in PROT.items()}


def families():
    return {"R1_OG0000001": {"Aaa_g1", "Bbb_g2"},
            "R1_OG0000002": {"Aaa_g3", "Bbb_g4"},
            "R2_OG0000003": {"Ccc_g5", "Bbb_g6"}}


def tblout(rows):
    """A HMMER3 --tblout body from (gene, family, evalue[, bits]) rows.

    Bit scores matter now: the rescue decides on them, not on the E-value
    (issue #45), and rows left at the same score are a genuine tie the code is
    required to refuse rather than break by emission order.
    """
    head = ("#                              --- full sequence ----\n"
            "# target name  accession  query name  accession  E-value\n")
    body = ""
    for row in rows:
        g, f, e = row[:3]
        bits = row[3] if len(row) > 3 else 120.0
        body += f"{g} - {f} - {e} {bits} 3.4 1 1\n"
    return head + body + "#\n# [ok]\n"


def domtblout(rows):
    """Complete synthetic --domtblout from (gene, family, E, bits[, ali_to])."""
    body = "# target name  ...\n"
    for row in rows:
        gene, family, evalue, bits = row[:4]
        ali_to = row[4] if len(row) > 4 else 100
        body += (
            f"{gene} - 100 {family} - 100 {evalue} {bits} 0.0 "
            f"1 1 {evalue} {evalue} {bits} 0.0 "
            f"1 100 1 {ali_to} 1 {ali_to} 0.99 -\n"
        )
    return body + "#\n# [ok]\n"


def make_outdir(tmp_path, with_alignments=True, round_ogs=()):
    d = tmp_path / "out"
    if with_alignments:
        for fam in families():
            p = d / "final_families" / fam / "confirmed_proteins.afa"
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(AFA)
    for rnd, og in round_ogs:
        p = d / f"round_{rnd:02d}" / "orthogroups" / og / "confirmed_proteins.afa"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(AFA)
    return d


def make_config(**kw):
    cfg = Config()
    cfg.n_workers = 2
    cfg.hmmer_evalue = 1e-5
    cfg.hmmer_chunk_size = 0
    for k, v in kw.items():
        setattr(cfg, k, v)
    return cfg


@pytest.fixture
def realigned(monkeypatch):
    """Swallow the mafft/FastTree re-alignment, recording what it was asked to do."""
    seen = []

    def fake(family_id, gene_ids, protein_pool, cds_pool, outdir, config):
        seen.append((family_id, sorted(gene_ids)))

    monkeypatch.setattr(hr, "_realign_family", fake)
    return seen


def install(monkeypatch, fake):
    monkeypatch.setattr(hr, "subprocess", fake)
    return fake


# --- paths that must return `families` untouched -------------------------

def test_returns_the_same_object_when_there_is_nothing_to_rescue(
        tmp_path, monkeypatch, realigned):
    # Arrange
    fake = install(monkeypatch, FakeSubprocess())
    fams = families()

    # Act
    out = hr.rescue_unplaced(fams, {}, PROT, CDS,
                             make_outdir(tmp_path), make_config())

    # Assert
    assert out is fams
    assert fake.calls == []


def test_returns_families_unchanged_when_no_family_has_an_alignment(
        tmp_path, monkeypatch, realigned):
    fake = install(monkeypatch, FakeSubprocess())
    fams = families()

    out = hr.rescue_unplaced(fams, {"Aaa_u1": PROT["Aaa_u1"]}, PROT, CDS,
                             make_outdir(tmp_path, with_alignments=False),
                             make_config())

    assert out is fams
    assert fake.ran("hmmsearch") == []


def test_returns_families_unchanged_when_hmmpress_fails(
        tmp_path, monkeypatch, realigned):
    fake = install(monkeypatch, FakeSubprocess(press_rc=1))
    fams = families()

    out = hr.rescue_unplaced(fams, {"Aaa_u1": PROT["Aaa_u1"]}, PROT, CDS,
                             make_outdir(tmp_path), make_config())

    assert out is fams
    assert fake.ran("hmmsearch") == []          # never got that far


def test_returns_families_unchanged_when_no_profile_reached_disk(
        tmp_path, monkeypatch, realigned):
    # hmmbuild claims success but writes nothing -> nothing to concatenate
    fake = install(monkeypatch, FakeSubprocess(build_silent=True))
    fams = families()

    out = hr.rescue_unplaced(fams, {"Aaa_u1": PROT["Aaa_u1"]}, PROT, CDS,
                             make_outdir(tmp_path), make_config())

    assert out is fams


def test_returns_families_unchanged_when_no_hit_clears_the_evalue(
        tmp_path, monkeypatch, realigned):
    fake = install(monkeypatch, FakeSubprocess(
        tblout_body=tblout([("Aaa_u1", "R1_OG0000001", "1.0")])))
    fams = families()

    out = hr.rescue_unplaced(fams, {"Aaa_u1": PROT["Aaa_u1"]}, PROT, CDS,
                             make_outdir(tmp_path), make_config())

    assert out is fams
    assert realigned == []


# --- the placing path ----------------------------------------------------

def test_places_genes_into_their_best_scoring_family(
        tmp_path, monkeypatch, realigned):
    # Arrange: Aaa_u1 hits two families, Bbb_u2's only hit is above the cutoff
    install(monkeypatch, FakeSubprocess(tblout_body=tblout([
        ("Aaa_u1", "R1_OG0000001", "1.2e-40", 900.0),
        ("Aaa_u1", "R1_OG0000002", "3.0e-12", 120.0),
        ("Bbb_u2", "R1_OG0000002", "5.0e-3"),
        ("Ccc_u3", "R1_OG0000001", "7.7e-08"),
    ])))
    outdir = make_outdir(tmp_path)
    fams = families()
    unplaced = {k: PROT[k] for k in ("Aaa_u1", "Bbb_u2", "Ccc_u3")}

    # Act
    out = hr.rescue_unplaced(fams, unplaced, PROT, CDS, outdir, make_config())

    # Assert
    assert out is not fams
    assert out["R1_OG0000001"] == {"Aaa_g1", "Bbb_g2", "Aaa_u1", "Ccc_u3"}
    assert out["R1_OG0000002"] == {"Aaa_g3", "Bbb_g4"}    # 120 bits lost to 900
    assert "Bbb_u2" not in {g for m in out.values() for g in m}


def test_does_not_mutate_the_families_it_was_given(
        tmp_path, monkeypatch, realigned):
    install(monkeypatch, FakeSubprocess(tblout_body=tblout([
        ("Aaa_u1", "R1_OG0000001", "1.2e-40")])))
    fams = families()
    before = {k: set(v) for k, v in fams.items()}

    hr.rescue_unplaced(fams, {"Aaa_u1": PROT["Aaa_u1"]}, PROT, CDS,
                       make_outdir(tmp_path), make_config())

    assert fams == before


def test_realigns_exactly_the_families_that_grew(
        tmp_path, monkeypatch, realigned):
    install(monkeypatch, FakeSubprocess(tblout_body=tblout([
        ("Aaa_u1", "R1_OG0000001", "1.0e-30"),
        ("Bbb_u2", "R2_OG0000003", "2.0e-20"),
    ])))

    hr.rescue_unplaced(families(),
                       {k: PROT[k] for k in ("Aaa_u1", "Bbb_u2")},
                       PROT, CDS, make_outdir(tmp_path), make_config())

    assert sorted(f for f, _ in realigned) == ["R1_OG0000001", "R2_OG0000003"]
    assert dict(realigned)["R1_OG0000001"] == ["Aaa_g1", "Aaa_u1", "Bbb_g2"]


def test_hit_against_an_unknown_family_is_reported_but_not_placed(
        tmp_path, monkeypatch, realigned):
    install(monkeypatch, FakeSubprocess(tblout_body=tblout([
        ("Aaa_u1", "R9_OG9999999", "1.0e-50"),
        ("Bbb_u2", "R1_OG0000001", "1.0e-11"),
    ])))
    outdir = make_outdir(tmp_path)

    out = hr.rescue_unplaced(families(),
                             {k: PROT[k] for k in ("Aaa_u1", "Bbb_u2")},
                             PROT, CDS, outdir, make_config())

    assert "R9_OG9999999" not in out
    assert [f for f, _ in realigned] == ["R1_OG0000001"]
    summary = (outdir / "hmmer_rescue" / "rescue_summary.tsv").read_text()
    assert "Aaa_u1\tR9_OG9999999\t" in summary


def test_rescue_summary_is_sorted_by_gene_with_two_digit_exponent_notation(
        tmp_path, monkeypatch, realigned):
    install(monkeypatch, FakeSubprocess(tblout_body=tblout([
        ("Ccc_u3", "R1_OG0000001", "7.7e-08"),
        ("Aaa_u1", "R1_OG0000001", "1.2e-40"),
    ])))
    outdir = make_outdir(tmp_path)

    hr.rescue_unplaced(families(),
                       {k: PROT[k] for k in ("Aaa_u1", "Ccc_u3")},
                       PROT, CDS, outdir, make_config())

    # The summary carries the evidence that decided each call, not just the
    # E-value it used to be picked by (issue #45): bits, the best-vs-second
    # margin, and the grade. This format change is deliberate - the earlier
    # three-column form recorded a number that was not the decision.
    lines = (outdir / "hmmer_rescue" / "rescue_summary.tsv").read_text().splitlines()
    assert lines[0] == "gene_id\tfamily_id\tbits\tmargin\tgrade\tevalue\treason"
    assert [l.split("\t")[0] for l in lines[1:]] == ["Aaa_u1", "Ccc_u3"]
    assert all(l.split("\t")[1] == "R1_OG0000001" for l in lines[1:])
    assert [l.split("\t")[5] for l in lines[1:]] == ["1.20e-40", "7.70e-08"]


def test_unplaced_pool_is_written_out_for_the_search(
        tmp_path, monkeypatch, realigned):
    install(monkeypatch, FakeSubprocess(tblout_body=tblout([])))
    outdir = make_outdir(tmp_path)

    hr.rescue_unplaced(families(), {"Aaa_u1": PROT["Aaa_u1"]},
                       PROT, CDS, outdir, make_config())

    assert (outdir / "hmmer_rescue" / "unplaced_proteins.fa").read_text() == (
        ">Aaa_u1\nMKTAYIAKQRQ\n")


# --- profile database construction ---------------------------------------

def test_existing_profiles_are_reused_instead_of_rebuilt(
        tmp_path, monkeypatch, realigned):
    outdir = make_outdir(tmp_path)
    hmm_dir = outdir / "hmmer_rescue" / "hmm_profiles"
    hmm_dir.mkdir(parents=True)
    for fam in families():
        (hmm_dir / f"{fam}.hmm").write_text(f"PLACEHOLDER {fam}\n//\n")
    install(monkeypatch, FakeSubprocess(tblout_body=tblout([
        ("Aaa_u1", "R1_OG0000001", "1.0e-30")])))

    out = hr.rescue_unplaced(families(), {"Aaa_u1": PROT["Aaa_u1"]},
                             PROT, CDS, outdir, make_config())

    # hmmbuild would have overwritten the placeholders; it was not run
    assert (hmm_dir / "R1_OG0000001.hmm").read_text().startswith("PLACEHOLDER")
    assert "Aaa_u1" in out["R1_OG0000001"]


def test_a_failing_hmmbuild_does_not_abort_the_rescue(
        tmp_path, monkeypatch, realigned):
    install(monkeypatch, FakeSubprocess(
        build_fail=["R1_OG0000002", "R2_OG0000003"],
        tblout_body=tblout([("Aaa_u1", "R1_OG0000001", "1.0e-30")])))

    out = hr.rescue_unplaced(families(), {"Aaa_u1": PROT["Aaa_u1"]},
                             PROT, CDS, make_outdir(tmp_path), make_config())

    assert "Aaa_u1" in out["R1_OG0000001"]


def test_profiles_are_found_in_the_round_directory_when_final_families_is_absent(
        tmp_path, monkeypatch, realigned):
    outdir = make_outdir(tmp_path, with_alignments=False,
                         round_ogs=[(1, "OG0000001")])
    install(monkeypatch, FakeSubprocess(tblout_body=tblout([
        ("Aaa_u1", "R1_OG0000001", "1.0e-30")])))

    out = hr.rescue_unplaced(families(), {"Aaa_u1": PROT["Aaa_u1"]},
                             PROT, CDS, outdir, make_config())

    # only the family whose round-directory alignment exists got a profile
    built = sorted(p.stem for p in
                   (outdir / "hmmer_rescue" / "hmm_profiles").glob("*.hmm"))
    assert built == ["R1_OG0000001"]
    assert "Aaa_u1" in out["R1_OG0000001"]


# --- chunked search (issue #31) ------------------------------------------

def test_chunked_search_skips_hmmpress_and_merges_every_chunk(
        tmp_path, monkeypatch, realigned):
    # Arrange
    outdir = make_outdir(tmp_path)
    hmm_dir = outdir / "hmmer_rescue" / "hmm_profiles"
    hmm_dir.mkdir(parents=True)
    for fam in families():
        (hmm_dir / f"{fam}.hmm").write_text(
            f"HMMER3/f [3.2.1]\nNAME  {fam}\nLENG  12\nHMM\n//\n")
    fake = install(monkeypatch, FakeSubprocess(chunk_bodies=[
        tblout([("Aaa_u1", "R1_OG0000001", "1.0e-30")]),
        tblout([("Bbb_u2", "R2_OG0000003", "2.0e-20")]),
    ], dom_chunk_bodies=[
        domtblout([("Aaa_u1", "R1_OG0000001", "1.0e-30", 120.0)]),
        domtblout([("Bbb_u2", "R2_OG0000003", "2.0e-20", 120.0)]),
    ]))

    # Act
    out = hr.rescue_unplaced(families(),
                             {k: PROT[k] for k in ("Aaa_u1", "Bbb_u2")},
                             PROT, CDS, outdir,
                             make_config(hmmer_chunk_size=2,
                                         hmmer_chunk_concurrent=2))

    # Assert — hits from BOTH chunks landed, and hmmpress was never run
    assert fake.ran("hmmpress") == []
    assert fake.ran("hmmsearch") == []          # the runner script does it
    assert "Aaa_u1" in out["R1_OG0000001"]
    assert "Bbb_u2" in out["R2_OG0000003"]


def test_chunked_search_uses_domain_output_for_the_query_coverage_gate(
        tmp_path, monkeypatch, realigned):
    """A strong two-residue domain must not be rescued via tblout fallback."""
    outdir = make_outdir(tmp_path)
    hmm_dir = outdir / "hmmer_rescue" / "hmm_profiles"
    hmm_dir.mkdir(parents=True)
    for fam in families():
        (hmm_dir / f"{fam}.hmm").write_text(
            f"HMMER3/f [3.2.1]\nNAME  {fam}\nLENG  12\nHMM\n//\n")
    install(monkeypatch, FakeSubprocess(
        chunk_bodies=[
            tblout([("Aaa_u1", "R1_OG0000001", "1.0e-30", 120.0)]),
            tblout([]),
        ],
        dom_chunk_bodies=[
            domtblout([
                ("Aaa_u1", "R1_OG0000001", "1.0e-30", 120.0, 2)
            ]),
            domtblout([]),
        ],
    ))

    out = hr.rescue_unplaced(
        families(), {"Aaa_u1": PROT["Aaa_u1"]}, PROT, CDS, outdir,
        make_config(hmmer_chunk_size=2),
    )

    assert "Aaa_u1" not in out["R1_OG0000001"]
    assert (outdir / "hmmer_rescue" /
            "hmmsearch_results.domtblout").exists()


def test_an_incomplete_chunk_fails_loudly_rather_than_under_reporting(
        tmp_path, monkeypatch, realigned):
    outdir = make_outdir(tmp_path)
    hmm_dir = outdir / "hmmer_rescue" / "hmm_profiles"
    hmm_dir.mkdir(parents=True)
    for fam in families():
        (hmm_dir / f"{fam}.hmm").write_text(
            f"HMMER3/f [3.2.1]\nNAME  {fam}\nLENG  12\nHMM\n//\n")
    truncated = "gene1 - R1_OG0000001 - 1e-40 120.0 3.4 1 1\n"   # no '# [ok]'
    install(monkeypatch, FakeSubprocess(chunk_bodies=[
        tblout([("Aaa_u1", "R1_OG0000001", "1.0e-30")]), truncated]))

    with pytest.raises(RuntimeError, match="incomplete hmmsearch chunks"):
        hr.rescue_unplaced(families(), {"Aaa_u1": PROT["Aaa_u1"]},
                           PROT, CDS, outdir,
                           make_config(hmmer_chunk_size=2))


# --- the extracted helpers in isolation ----------------------------------

def test_group_by_family_inverts_the_hit_table():
    def hit(fam, bits):
        return hr.RescueHit(family=fam, bits=bits, evalue=1e-40, margin=None,
                            grade="HIGH", reason="")

    grouped = hr._group_by_family({
        "Aaa_u1": hit("F1", 900.0), "Bbb_u2": hit("F1", 800.0),
        "Ccc_u3": hit("F2", 700.0)})

    assert grouped == {"F1": {"Aaa_u1", "Bbb_u2"}, "F2": {"Ccc_u3"}}


def test_absorb_rescued_leaves_the_input_sets_alone(tmp_path, monkeypatch):
    monkeypatch.setattr(hr, "_realign_family", lambda *a, **k: None)
    original = {"F1": {"a"}}
    members = original["F1"]

    out = hr._absorb_rescued(original, {"F1": {"b"}}, PROT, CDS,
                             tmp_path, make_config())

    assert out["F1"] == {"a", "b"}
    assert members == {"a"}                     # the original set is untouched
    assert out["F1"] is not members


# --- _concat_hmms: streaming the profile database together ---------------
#
# The profile DB is ~1.6 GB on the 5sp panel, so this is streamed rather
# than read whole. These cover what the stream has to get right: exact
# bytes across the chunk boundary, the glob/sort that fixes concatenation
# order, and leaving the target alone when there is nothing to write.

CHUNK = 8192          # the streaming unit _concat_hmms uses


def binary_payload(seed, n):
    """Deterministic pseudo-binary bytes, NULs and 0xFF included."""
    return bytes((seed * 31 + i * 17 + (i >> 8) * 7) & 0xFF for i in range(n))


def profile_dir(tmp_path, files, name="profiles"):
    d = tmp_path / name
    d.mkdir(parents=True, exist_ok=True)
    for fname, data in files:
        (d / fname).write_bytes(data)
    return d


@pytest.mark.parametrize("size", [
    0, 1, CHUNK - 1, CHUNK, CHUNK + 1, CHUNK * 13 + 77])
def test_a_profile_survives_concatenation_byte_for_byte_at_any_size(
        tmp_path, size):
    # Arrange: sizes that straddle the streaming chunk boundary
    data = binary_payload(3, size)
    d = profile_dir(tmp_path, [("only.hmm", data)])
    out = tmp_path / "all_families.hmm"

    # Act
    result = hr._concat_hmms(d, out)

    # Assert
    assert result == out
    assert out.read_bytes() == data


def test_multiple_profiles_are_joined_in_sorted_order_with_no_separator(
        tmp_path):
    parts = [("b.hmm", binary_payload(1, CHUNK * 2)),
             ("a.hmm", binary_payload(2, 5)),
             ("c.hmm", binary_payload(3, CHUNK + 3))]
    d = profile_dir(tmp_path, parts)

    hr._concat_hmms(d, tmp_path / "all.hmm")

    expected = dict(parts)
    assert (tmp_path / "all.hmm").read_bytes() == (
        expected["a.hmm"] + expected["b.hmm"] + expected["c.hmm"])


def test_round_ten_sorts_before_round_one_in_the_concatenated_database(
        tmp_path):
    # '0' < '_', so R10_ precedes R1_ -- the same lexicographic trap that
    # makes E-value ties resolve wrongly elsewhere in this pipeline.
    d = profile_dir(tmp_path, [("R1_OG1.hmm", b"ONE"),
                               ("R2_OG1.hmm", b"TWO"),
                               ("R10_OG1.hmm", b"TEN")])

    hr._concat_hmms(d, tmp_path / "all.hmm")

    assert (tmp_path / "all.hmm").read_bytes() == b"TENONETWO"


def test_an_empty_profile_contributes_nothing_and_breaks_nothing(tmp_path):
    d = profile_dir(tmp_path, [("a.hmm", b"AAA"), ("b.hmm", b""),
                               ("c.hmm", b"CCC")])

    hr._concat_hmms(d, tmp_path / "all.hmm")

    assert (tmp_path / "all.hmm").read_bytes() == b"AAACCC"


def test_only_dot_hmm_files_are_concatenated(tmp_path):
    # hmmpress leaves .h3i/.h3m/.h3f/.h3p siblings next to the database
    d = profile_dir(tmp_path, [("a.hmm", b"KEEP"), ("a.hmm.h3i", b"INDEX"),
                               ("b.HMM", b"WRONGCASE"), ("notes.txt", b"NO")])

    hr._concat_hmms(d, tmp_path / "all.hmm")

    assert (tmp_path / "all.hmm").read_bytes() == b"KEEP"


def test_no_profiles_returns_none_and_leaves_any_existing_database_alone(
        tmp_path):
    # Arrange: the caller treats None as "cannot proceed" and returns the
    # families untouched -- a truncated database here would be worse than
    # no database at all.
    d = profile_dir(tmp_path, [])
    out = tmp_path / "all.hmm"
    out.write_bytes(b"PREVIOUS DATABASE")

    result = hr._concat_hmms(d, out)

    assert result is None
    assert out.read_bytes() == b"PREVIOUS DATABASE"


def test_a_directory_of_non_profiles_is_the_same_as_an_empty_one(tmp_path):
    d = profile_dir(tmp_path, [("a.txt", b"X"), ("b.hmm.bak", b"Y")])

    assert hr._concat_hmms(d, tmp_path / "all.hmm") is None
    assert not (tmp_path / "all.hmm").exists()


def test_a_stale_database_is_overwritten_not_appended_to(tmp_path):
    d = profile_dir(tmp_path, [("a.hmm", b"FRESH")])
    out = tmp_path / "all.hmm"
    out.write_bytes(b"STALE" * 1000)

    hr._concat_hmms(d, out)

    assert out.read_bytes() == b"FRESH"
