#!/usr/bin/env python3
"""Run every annotation axis for one family — as one declared, replayable plan.

The PEPC clan was annotated by hand-assembling five tool invocations across
two machines. That is not repeatable. This module declares the same stack
once so any family can be put through it:

  signalp   SignalP 6         secretory signal peptides
  deeploc   DeepLoc 2.1       localization + sorting signals (NES/NLS/mTP)
  emapper   eggNOG-mapper     EC / gene symbol / GO / KEGG by orthology
  clean     CLEAN             EC prediction with calibrated confidence
  foldseek  Foldseek          structure function transfer vs AFDB-SwissProt
                              (needs predicted structures — ESMFold or ProstT5)
  structure steps/gene_structure  intron positions vs the family's conserved
                              set, with the annotation programme as covariate
                              (needs the family alignment and per-species GFFs)
  expression steps/expression  mean TPM and within-species share
                              (needs a per-species expression matrix; only
                              two of the PEPC clan's species have one)

The heavy axes run on the GPU box; the merge (`annotation_matrix.py`) is
pure Python and runs anywhere. `--dry-run` prints the plan without touching
anything, which is also the fastest way to see what a given family is
missing.

Usage:
  # see the plan (nothing runs)
  python annotate_stack.py --family-fasta clan.fa --workdir '~/annot/clan' \\
      --structures '~/pepc_pilot/pdb' --expected-ec 4.1.1.31 --dry-run

  # execute it on the gpu host, then merge locally
  python annotate_stack.py --family-fasta clan.fa --workdir '~/annot/clan' \\
      --structures '~/pepc_pilot/pdb' --expected-ec 4.1.1.31 \\
      --host gpu --outdir annot_clan

Axis tool paths are the installed locations recorded in the lab wiki
(guide/installs.md, 2026-08-20/21 entries).
"""

import argparse
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

# Installed locations on the GPU box (lab wiki guide/installs.md).
SIGNALP_BIN = "/home/wyim/micromamba/envs/signalp6/bin/signalp6"
DEEPLOC_PY = "/home/wyim/micromamba/envs/transgenic/bin/deeploc2"
EMAPPER_ENV = "/home/wyim/micromamba/envs/emapper/bin"
EGGNOG_DATA = "~/scratch/eggnog_data"
CLEAN_APP = "~/scratch/CLEAN/app"
MICROMAMBA = "/home/wyim/scratch/bin/micromamba"
FOLDSEEK_BIN = "~/bin/foldseek/bin/foldseek"
AFDB_DB = "~/foldseek_dbs/afdb_swissprot"

ALL_AXES = ("signalp", "deeploc", "emapper", "clean", "foldseek",
            "structure", "expression")


def qpath(path: str) -> str:
    """Quote a path for the shell WITHOUT killing a leading `~`.

    `shlex.quote('~/x')` returns `'~/x'`, which bash treats as a literal
    directory named `~` — every axis would then write into the wrong
    place. Remote paths are given as `~/...` constantly here, so keep the
    tilde outside the quotes and quote only the remainder.
    """
    if path == "~":
        return "~"
    if path.startswith("~/"):
        rest = path[2:]
        quoted = shlex.quote(rest)
        return f"~/{rest}" if quoted == rest else f"~/{quoted}"
    return shlex.quote(path)


@dataclass(frozen=True)
class Axis:
    name: str
    command: str
    output: str            # path (remote) holding this axis's result
    needs: str = ""        # human-readable prerequisite, "" when always runnable


def build_plan(
    family_fasta: str,
    workdir: str,
    expected_ec: Optional[str] = None,
    axes: Optional[List[str]] = None,
    structures: Optional[str] = None,
    drop_unavailable: bool = False,
    alignment: Optional[str] = None,
    gffs: Optional[Dict[str, str]] = None,
    matrices: Optional[Dict[str, str]] = None,
    members: Optional[str] = None,
    groups: Optional[str] = None,
) -> List[Axis]:
    """Declare the per-axis commands for one family.

    `structures` is a directory of predicted structures (PDB) — without it
    the foldseek axis cannot run and is either reported as blocked
    (default) or dropped (`drop_unavailable`).

    `alignment` + `gffs` (species -> GFF3) enable the gene-structure axis,
    and `matrices` (species -> expression matrix) the expression axis. Both
    are declared the same way as foldseek's structures: absent means BLOCKED
    and named, never silently skipped, because "no rows" from either axis
    reads exactly like a measured negative.
    """
    wanted = list(ALL_AXES) if axes is None else list(axes)
    unknown = [a for a in wanted if a not in ALL_AXES]
    if unknown:
        raise ValueError(f"unknown axis: {', '.join(unknown)}")

    q = qpath
    fa, wd = q(family_fasta), q(workdir)
    clean_fa = f"{wd}/input.clean.fa"

    # Every tool here chokes on trailing '*' (SignalP/DeepLoc tokenizers,
    # DIAMOND); clean once, up front, and let the axes read the clean copy.
    prep = f"mkdir -p {wd} && sed '/^>/! s/\\*//g' {fa} > {clean_fa}"

    catalog: Dict[str, Axis] = {
        "signalp": Axis(
            "signalp",
            f"{prep} && {SIGNALP_BIN} --fastafile {clean_fa} --organism eukarya "
            f"--output_dir {wd}/signalp --format none --mode fast",
            f"{workdir}/signalp/prediction_results.txt",
        ),
        "deeploc": Axis(
            "deeploc",
            f"{prep} && {DEEPLOC_PY} -f {clean_fa} -o {wd}/deeploc "
            f"-m Accurate -d cuda && mv {wd}/deeploc/results_*.csv "
            f"{wd}/deeploc/deeploc_results.csv",
            f"{workdir}/deeploc/deeploc_results.csv",
        ),
        "emapper": Axis(
            "emapper",
            f"{prep} && PATH={EMAPPER_ENV}:$PATH emapper.py -m diamond "
            f"-i {clean_fa} --data_dir {EGGNOG_DATA} -o family "
            f"--output_dir {wd}/emapper --cpu 14 --override",
            f"{workdir}/emapper/family.emapper.annotations",
        ),
        "clean": Axis(
            "clean",
            f"{prep} && cp {clean_fa} {CLEAN_APP}/data/inputs/family.fasta && "
            f"cd {CLEAN_APP} && {MICROMAMBA} run -n clean python "
            f"CLEAN_infer_fasta.py --fasta_data family",
            f"{CLEAN_APP}/results/inputs/family_maxsep.csv",
        ),
        "foldseek": Axis(
            "foldseek",
            (f"mkdir -p {wd}/foldseek && {FOLDSEEK_BIN} easy-search "
             f"{q(structures) if structures else '<STRUCTURES>'} {AFDB_DB} "
             f"{wd}/foldseek/hits.tsv {wd}/foldseek/tmp --format-output "
             "query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,"
             "alntmscore,qcov,tcov,taxname --max-seqs 30 -e 1e-5 --threads 14"
             f" && python3 fs_transfer.py --hits {wd}/foldseek/hits.tsv "
             f"-o {wd}/foldseek/structure_function_transfer.tsv"),
            f"{workdir}/foldseek/structure_function_transfer.tsv",
            needs="predicted structures (--structures)",
        ),
        "structure": Axis(
            "structure",
            (f"mkdir -p {wd} && python3 -m steps.gene_structure "
             f"--alignment {q(alignment) if alignment else '<ALIGNMENT>'} "
             + ("".join(f"--gff {sp}={q(path)} "
                        for sp, path in sorted((gffs or {}).items()))
                or "--gff <SPECIES=GFF3> ")
             + (f"--groups {q(groups)} " if groups else "")
             + f"-o {wd}/gene_structure.tsv "
             f"--json {wd}/gene_structure.json"),
            f"{workdir}/gene_structure.tsv",
            needs="the family alignment and per-species GFF3s "
                  "(--alignment / --gff SPECIES=PATH)",
        ),
        "expression": Axis(
            "expression",
            (f"mkdir -p {wd} && python3 -m steps.expression "
             f"--members {q(members) if members else '<MEMBERS>'} "
             + ("".join(f"--matrix {sp}={q(path)} "
                        for sp, path in sorted((matrices or {}).items()))
                or "--matrix <SPECIES=MATRIX> ")
             + (f"--groups {q(groups)} " if groups else "")
             + f"-o {wd}/expression.tsv --json {wd}/expression.json"),
            f"{workdir}/expression.tsv",
            needs="a member list and at least one per-species expression "
                  "matrix (--members / --matrix SPECIES=PATH)",
        ),
    }

    plan = [catalog[a] for a in wanted]
    if drop_unavailable:
        blocked = set(missing_inputs(plan))
        plan = [a for a in plan if a.name not in blocked]
    return plan


_PLACEHOLDERS = ("<STRUCTURES>", "<ALIGNMENT>", "<SPECIES=GFF3>",
                 "<MEMBERS>", "<SPECIES=MATRIX>")


def missing_inputs(plan: List[Axis]) -> Dict[str, str]:
    """Axes that cannot run yet -> why."""
    return {a.name: a.needs for a in plan
            if a.needs and any(p in a.command for p in _PLACEHOLDERS)}


_MERGE_FLAG = {
    "emapper": "--emapper",
    "clean": "--clean-csv",
    "foldseek": "--foldseek-tsv",
    "deeploc": "--deeploc-csv",
    "signalp": "--signalp",
    "structure": "--gene-structure",
    "expression": "--expression",
}


def local_merge_command(plan: List[Axis], outdir: str,
                        expected_ec: Optional[str]) -> str:
    """The annotation_matrix.py invocation that consumes this plan's outputs."""
    parts = ["python annotation_matrix.py", f"--outdir {qpath(outdir)}"]
    for a in plan:
        flag = _MERGE_FLAG.get(a.name)
        if flag:
            parts.append(f"{flag} {qpath(a.output)}")
    if expected_ec:
        parts.append(f"--expected-ec {expected_ec}")
    return " \\\n    ".join(parts)


def fetch_plan(plan: List[Axis], local_dir: str) -> List[Axis]:
    """Rewrite each axis's output to where its LOCAL copy will live.

    The axes run on --host, so `Axis.output` is a remote path; the merge
    runs locally and cannot read it. Local paths are namespaced by axis so
    two axes cannot collide on a basename.
    """
    out = []
    for a in plan:
        local = f"{local_dir}/{a.name}/{Path(a.output).name}"
        out.append(Axis(a.name, a.command, local, a.needs))
    return out


def fetch_outputs(plan: List[Axis], host: str, local_dir: str) -> List[Axis]:
    """scp each axis output back, returning the plan with local paths."""
    fetched = fetch_plan(plan, local_dir)
    for remote, local_axis in zip(plan, fetched):
        dest = Path(local_axis.output)
        dest.parent.mkdir(parents=True, exist_ok=True)
        r = subprocess.run(["scp", "-q", f"{host}:{remote.output}", str(dest)])
        if r.returncode != 0:
            print(f"could not fetch {remote.name} from {host}:{remote.output} "
                  "— the merge tolerates a missing axis", file=sys.stderr)
    return fetched


def _pairs(parser, flag: str, values: List[str]) -> Dict[str, str]:
    """`SPECIES=PATH` arguments -> {species: path}."""
    out: Dict[str, str] = {}
    for value in values:
        species, _, path = value.partition("=")
        if not path:
            parser.error(f"{flag} expects SPECIES=PATH, got {value!r}")
        out[species] = path
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--family-fasta", required=True,
                    help="Family protein FASTA (path ON THE HOST that runs the axes)")
    ap.add_argument("--workdir", required=True, help="Remote working directory")
    ap.add_argument("--structures", default=None,
                    help="Directory of predicted structures (enables foldseek axis)")
    ap.add_argument("--alignment", default=None,
                    help="Family alignment (enables the gene-structure axis)")
    ap.add_argument("--gff", action="append", default=[],
                    metavar="SPECIES=PATH",
                    help="GFF3 for one species; repeatable")
    ap.add_argument("--members", default=None,
                    help="Family member list (enables the expression axis)")
    ap.add_argument("--matrix", action="append", default=[],
                    metavar="SPECIES=PATH",
                    help="Expression matrix for one species; repeatable")
    ap.add_argument("--groups", default=None,
                    help="TSV of gene<TAB>subfamily, used by both new axes")
    ap.add_argument("--expected-ec", default=None,
                    help="Family EC for the membership verdict in the merge")
    ap.add_argument("--axes", default=None,
                    help=f"Comma-separated subset of: {','.join(ALL_AXES)}")
    ap.add_argument("--outdir", default="annot",
                    help="Local output directory for the merged matrix")
    ap.add_argument("--host", default=None,
                    help="ssh host to run the axes on (omit for --dry-run only)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print the plan; run nothing")
    args = ap.parse_args()

    axes = args.axes.split(",") if args.axes else None
    plan = build_plan(
        family_fasta=args.family_fasta, workdir=args.workdir,
        expected_ec=args.expected_ec, axes=axes, structures=args.structures,
        drop_unavailable=False,
        alignment=args.alignment, gffs=_pairs(ap, "--gff", args.gff),
        matrices=_pairs(ap, "--matrix", args.matrix), members=args.members,
        groups=args.groups,
    )

    blocked = missing_inputs(plan)
    for name, why in blocked.items():
        print(f"[blocked] {name}: needs {why}", file=sys.stderr)

    print("# annotation stack plan")
    for a in plan:
        mark = "  (BLOCKED)" if a.name in blocked else ""
        print(f"\n## {a.name}{mark}\n{a.command}\n#  -> {a.output}")
    print("\n## merge (local)")
    print(local_merge_command(plan, args.outdir, args.expected_ec))

    if args.dry_run or not args.host:
        if not args.dry_run:
            print("\n(no --host given: nothing executed)", file=sys.stderr)
        return 0

    for a in plan:
        if a.name in blocked:
            print(f"skipping blocked axis: {a.name}", file=sys.stderr)
            continue
        print(f"\n=== running {a.name} on {args.host} ===", file=sys.stderr)
        r = subprocess.run(["ssh", args.host, a.command])
        if r.returncode != 0:
            print(f"axis {a.name} failed (exit {r.returncode}) — continuing; "
                  "the merge tolerates missing axes", file=sys.stderr)
    print(f"\n=== fetching outputs to {args.outdir} ===", file=sys.stderr)
    fetched = fetch_outputs([a for a in plan if a.name not in blocked],
                            args.host, args.outdir)
    print("\n## merge (local, fetched paths)")
    print(local_merge_command(fetched, args.outdir, args.expected_ec))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
