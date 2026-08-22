"""annotate_stack.py — one plan that runs every annotation axis for a family.

The point is repeatability: the PEPC clan was annotated by hand-assembling
five tool invocations across two machines. This module turns that into a
declared plan that can be printed, checked, and replayed for any family.

Tests cover plan construction only — no tool is ever invoked.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from annotate_stack import Axis, build_plan, local_merge_command, missing_inputs


def _plan(**kw):
    kw.setdefault("family_fasta", "/data/clan.fa")
    kw.setdefault("workdir", "~/annot/clan")
    kw.setdefault("expected_ec", "4.1.1.31")
    return build_plan(**kw)


def test_plan_covers_every_axis_by_default():
    names = {a.name for a in _plan()}
    assert names == {"signalp", "deeploc", "emapper", "clean", "foldseek",
                     "structure", "expression"}


def test_axes_can_be_selected():
    names = {a.name for a in _plan(axes=["signalp", "emapper"])}
    assert names == {"signalp", "emapper"}


def test_unknown_axis_is_rejected():
    try:
        _plan(axes=["signalp", "nosuchtool"])
    except ValueError as e:
        assert "nosuchtool" in str(e)
        return
    raise AssertionError("unknown axis must raise")


def test_every_axis_declares_its_output_and_command():
    for a in _plan():
        assert isinstance(a, Axis)
        assert a.command.strip(), a.name
        assert a.output, a.name


def test_stop_codons_are_stripped_before_any_tool():
    """Every tool here chokes on '*' (SignalP/DeepLoc tokenizers, DIAMOND);
    the plan must clean the input rather than leave it to the caller."""
    plan = _plan()
    prep = [a for a in plan if a.name == "signalp"][0]
    assert "s/\\*//g" in prep.command or "clean" in prep.command


def test_foldseek_axis_requires_structures():
    plan = _plan(structures="~/pdb")
    fs = [a for a in plan if a.name == "foldseek"][0]
    assert "~/pdb" in fs.command
    assert "afdb_swissprot" in fs.command


def test_foldseek_axis_is_dropped_without_structures():
    names = {a.name for a in _plan(structures=None, axes=None, drop_unavailable=True)}
    assert "foldseek" not in names


def test_missing_inputs_reports_what_blocks_each_axis():
    blocked = missing_inputs(_plan(structures=None))
    assert "foldseek" in blocked
    assert "signalp" not in blocked


def test_local_merge_command_wires_every_axis_output():
    cmd = local_merge_command(_plan(), outdir="annot", expected_ec="4.1.1.31")
    for flag in ("--emapper", "--clean-csv", "--foldseek-tsv",
                 "--deeploc-csv", "--signalp", "--expected-ec 4.1.1.31"):
        assert flag in cmd, flag


def test_local_merge_skips_axes_not_in_plan():
    cmd = local_merge_command(_plan(axes=["signalp"]), outdir="annot",
                              expected_ec=None)
    assert "--signalp" in cmd
    assert "--emapper" not in cmd
    assert "--expected-ec" not in cmd


def test_commands_are_shell_quoted_against_spaces():
    plan = _plan(family_fasta="/data/my clan.fa")
    sp = [a for a in plan if a.name == "signalp"][0]
    assert "'/data/my clan.fa'" in sp.command


# --- path handling (found by dry-running the real PEPC plan) --------------

def test_tilde_paths_stay_expandable():
    """shlex.quote('~/x') yields '~/x' — a literal directory named '~'.
    Remote paths are given as ~/... constantly, so the leading ~ must
    survive quoting or every axis writes into the wrong place."""
    plan = _plan(family_fasta="~/pepc/clan.fa", workdir="~/annot/pepc")
    sp = [a for a in plan if a.name == "signalp"][0]
    assert "'~/annot/pepc'" not in sp.command
    assert "~/annot/pepc" in sp.command
    assert "'~/pepc/clan.fa'" not in sp.command


def test_tilde_path_with_space_is_still_quoted():
    plan = _plan(workdir="~/my annot/pepc")
    sp = [a for a in plan if a.name == "signalp"][0]
    assert "~/'my annot/pepc'" in sp.command


def test_deeploc_output_is_a_deterministic_file():
    """DeepLoc writes results_<timestamp>.csv; the merge needs one fixed
    path, so the axis must normalise the name itself."""
    plan = _plan()
    dl = [a for a in plan if a.name == "deeploc"][0]
    assert dl.output.endswith(".csv")
    assert "results_" in dl.command and "mv" in dl.command


def test_foldseek_axis_calls_the_repo_transfer_script():
    fs = [a for a in _plan(structures="~/pdb") if a.name == "foldseek"][0]
    assert "fs_transfer.py" in fs.command
    assert "hits.tsv" in fs.command


# --- fetch step (found by actually executing the driver) ------------------

def test_merge_command_uses_local_paths_after_fetch(tmp_path):
    """The axes run on --host, so their outputs are REMOTE paths. A merge
    command that pastes remote paths into a local invocation cannot work —
    the fetched local copies are what the merge must reference."""
    from annotate_stack import fetch_plan, local_merge_command
    plan = _plan(axes=["signalp", "emapper"], workdir="~/annot/clan")
    fetched = fetch_plan(plan, local_dir="annot_clan")
    cmd = local_merge_command(fetched, outdir="annot_clan", expected_ec=None)
    assert "~/annot/clan" not in cmd
    assert "annot_clan/signalp/prediction_results.txt" in cmd


def test_fetch_plan_preserves_axis_names_and_commands(tmp_path):
    from annotate_stack import fetch_plan
    plan = _plan(axes=["signalp"])
    fetched = fetch_plan(plan, local_dir="out")
    assert [a.name for a in fetched] == ["signalp"]
    assert fetched[0].command == plan[0].command      # unchanged
    assert fetched[0].output.startswith("out/")


def test_fetch_plan_keeps_basenames_unique_across_axes():
    """Two axes must not collide on a local filename."""
    from annotate_stack import fetch_plan
    fetched = fetch_plan(_plan(), local_dir="out")
    outs = [a.output for a in fetched]
    assert len(set(outs)) == len(outs)


# --- gene-structure and expression axes (issue #38) ----------------------

def test_the_new_axes_are_blocked_and_named_without_their_inputs():
    plan = build_plan("clan.fa", "~/wd")
    blocked = missing_inputs(plan)
    assert "structure" in blocked and "expression" in blocked
    assert "--alignment" in blocked["structure"]
    assert "--matrix" in blocked["expression"]


def test_the_gene_structure_axis_declares_one_gff_per_species():
    plan = {a.name: a for a in build_plan(
        "clan.fa", "~/wd", alignment="~/clan.aln",
        gffs={"Ococ": "~/Oco.gff3", "Mcry": "~/Mcry.gff3"},
        groups="~/sf.tsv")}
    command = plan["structure"].command
    assert "--gff Mcry=~/Mcry.gff3" in command
    assert "--gff Ococ=~/Oco.gff3" in command
    assert "--groups ~/sf.tsv" in command
    assert "structure" not in missing_inputs(build_plan(
        "clan.fa", "~/wd", alignment="~/clan.aln",
        gffs={"Ococ": "~/Oco.gff3"}))


def test_the_expression_axis_declares_one_matrix_per_species():
    plan = {a.name: a for a in build_plan(
        "clan.fa", "~/wd", members="~/members.txt",
        matrices={"Ococ": "~/RSEM_TPM.average.tsv"})}
    command = plan["expression"].command
    assert "--members ~/members.txt" in command
    assert "--matrix Ococ=~/RSEM_TPM.average.tsv" in command


def test_the_merge_command_picks_up_both_new_axes():
    plan = build_plan("clan.fa", "~/wd", alignment="~/clan.aln",
                      gffs={"Ococ": "~/Oco.gff3"}, members="~/m.txt",
                      matrices={"Ococ": "~/tpm.tsv"})
    merge = local_merge_command(plan, "annot", "4.1.1.31")
    assert "--gene-structure ~/wd/gene_structure.tsv" in merge
    assert "--expression ~/wd/expression.tsv" in merge


def test_drop_unavailable_drops_every_blocked_axis_not_only_foldseek():
    names = {a.name for a in build_plan("clan.fa", "~/wd",
                                        drop_unavailable=True)}
    assert names == {"signalp", "deeploc", "emapper", "clean"}
