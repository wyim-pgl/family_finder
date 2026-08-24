#!/usr/bin/env python3
"""Tree-validate family boundaries and apply the merges (issue #47).

`steps/profile_assign.vote_edges` nominates: it says which families look like
pieces of one family. `steps/cluster_validate.fragment_verdict` decides, on
each cluster's own codon tree. Until this file existed, only the nominating
half had a caller - the deciding half ran from ad-hoc scripts on one cluster,
so the published family table could not be reproduced, re-run at another
threshold, or applied to a different dataset.

The verdict rule (see steps/cluster_validate for the reasoning):

  INTERLEAVED   the fragment does not form its own clade - it mixes with
                another fragment's. Lineage-axis fragmentation. MERGE.
  MONOPHYLETIC  the fragment is its own clade. Topology alone cannot tell a
                subfamily of the same family from a distinct neighbouring
                family, so this is reported and merged by nobody.

Two modes, because the work is embarrassingly parallel but the bookkeeping
must not be:

  --judge      align/tree/judge clusters (all of them, or one with --cluster,
               which is what a SLURM array task calls)
  --apply      fold the verdicts into the family table and write summary_v3

Splitting them is deliberate. A campaign over thousands of clusters is going
to be interrupted; --judge is resumable per cluster and --apply refuses to run
on an incomplete verdict set rather than quietly shipping a partial merge.
"""
import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

from steps.cluster_validate import fragment_verdict
from utils.seqio import read_fasta, write_fasta

logger = logging.getLogger("family_finder")

VERDICT_HEADER = ("cluster_id\tfragment\tstatus\tn_members\tn_missing\t"
                  "merge_group\treason")


# ---------------------------------------------------------------------------
# reading the run's own outputs
# ---------------------------------------------------------------------------

def load_families(summary_path: Path) -> Dict[str, List[str]]:
    """family_id -> member gene ids, from a run's summary.tsv."""
    families: Dict[str, List[str]] = {}
    with open(summary_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi = header.index("gene_list")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= gi:
                continue
            families[parts[0]] = [g for g in parts[gi].split(",") if g]
    return families


def load_clusters(clusters_path: Path) -> Dict[str, List[str]]:
    """cluster_id -> member family ids, from fragmentation_clusters.tsv."""
    clusters: Dict[str, List[str]] = {}
    with open(clusters_path) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            clusters[parts[0]] = [f for f in parts[3].split(",") if f]
    return clusters


def load_verdict_rows(path: Path) -> List[tuple]:
    """Verdict rows as tuples, header skipped."""
    rows = []
    with open(path) as fh:
        first = fh.readline()
        if not first.startswith("cluster_id\t"):
            fh.seek(0)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6:
                rows.append(tuple(parts))
    return rows


# ---------------------------------------------------------------------------
# verdicts -> merges
# ---------------------------------------------------------------------------

def merge_groups_from_rows(rows: Sequence[tuple]) -> List[List[str]]:
    """Distinct merge groups named by INTERLEAVED rows.

    Every member of a group carries the same `A+B+C` string, so the groups
    are deduplicated rather than counted once per row. A group naming a
    single family merges nothing and is dropped: `fragment_verdict` only
    emits multi-member groups, but a hand-edited verdict file might not.
    """
    seen = set()
    groups: List[List[str]] = []
    for row in rows:
        status, group = row[2], row[5]
        if status != "INTERLEAVED" or not group:
            continue
        members = sorted(f for f in group.split("+") if f)
        if len(members) < 2:
            continue
        key = tuple(members)
        if key not in seen:
            seen.add(key)
            groups.append(members)
    return groups


def apply_merges(
    families: Dict[str, List[str]], groups: Sequence[Sequence[str]],
) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """Fold each group into its first family. Returns (families, provenance).

    Overlapping groups are resolved as connected components before anything
    moves. Taking `[[A,B],[B,C]]` at face value would fold B into A and then
    look for a B that is no longer there - or, worse, copy its genes twice.

    An unknown family id raises rather than being skipped: a verdict file
    from a different run would otherwise silently drop genes, and this
    function's whole job is to preserve them.
    """
    if not groups:
        return dict(families), {}

    for group in groups:
        for fam in group:
            if fam not in families:
                raise KeyError(
                    f"merge group names {fam!r}, which is not in the family "
                    f"table - the verdicts and the summary are from "
                    f"different runs")

    parent: Dict[str, str] = {}

    def find(x: str) -> str:
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for group in groups:
        for other in group[1:]:
            a, b = find(group[0]), find(other)
            if a != b:
                parent[a] = b

    comps: Dict[str, List[str]] = {}
    for fam in parent:
        comps.setdefault(find(fam), []).append(fam)

    merged = dict(families)
    provenance: Dict[str, str] = {}
    for members in comps.values():
        members = sorted(members)
        if len(members) < 2:
            continue
        keep = members[0]
        genes = list(families[keep])
        for other in members[1:]:
            genes.extend(families[other])
            del merged[other]
        merged[keep] = genes
        provenance[keep] = "+".join(members)
    return merged, provenance


# ---------------------------------------------------------------------------
# judging one cluster (thin wrappers so tests never invoke the tools)
# ---------------------------------------------------------------------------

def _run(cmd: List[str], stdout_path: Path = None) -> bool:
    """Run one external tool. Module-level so tests can monkeypatch it."""
    if stdout_path:
        with open(stdout_path, "w") as out:
            r = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE)
    else:
        r = subprocess.run(cmd, capture_output=True)
    if r.returncode != 0:
        logger.debug("failed: %s", " ".join(cmd[:3]))
    return r.returncode == 0


def build_cluster_tree(pep: Path, cds: Path, workdir: Path, cfg) -> Path:
    """Align, project to codons, and build the cluster's tree.

    Falls back to the protein tree when pal2nal produces nothing - the same
    fallback the pipeline uses, and the verdict only needs topology.
    """
    aln = workdir / "aln.fa"
    codon = workdir / "codon.fa"
    tree = workdir / "tree.nwk"
    if not _run([cfg.mafft_bin, "--auto", "--thread", str(cfg.threads),
                 str(pep)], aln):
        raise RuntimeError(f"MAFFT failed on {pep}")
    ok = (cds and cds.exists()
          and _run(["perl", cfg.pal2nal, str(aln), str(cds),
                    "-output", "fasta"], codon)
          and codon.exists() and codon.stat().st_size > 0)
    if ok:
        _run([cfg.fasttree_bin, "-nt", "-gtr", "-gamma", str(codon)], tree)
    else:
        _run([cfg.fasttree_bin, str(aln)], tree)
    if not tree.exists() or tree.stat().st_size == 0:
        raise RuntimeError(f"FastTree produced no tree for {pep}")
    return tree


def judge_cluster(cluster_id: str, fam_ids: Sequence[str],
                  families: Dict[str, List[str]], pep_pool, cds_pool,
                  outdir: Path, cfg) -> List[tuple]:
    """Align/tree/judge one cluster and return its verdict rows."""
    workdir = outdir / cluster_id
    workdir.mkdir(parents=True, exist_ok=True)
    fragments = {f: families[f] for f in fam_ids}
    wanted = [g for m in fragments.values() for g in m]

    pep = workdir / "cluster.pep.fa"
    write_fasta({g: pep_pool[g] for g in wanted if g in pep_pool}, str(pep))
    cds = workdir / "cluster.cds.fa"
    have_cds = {g: cds_pool[g] for g in wanted if g in cds_pool}
    if have_cds:
        write_fasta(have_cds, str(cds))

    tree = build_cluster_tree(pep, cds if have_cds else None, workdir, cfg)
    verdict = fragment_verdict(tree.read_text(), fragments)

    group_of = {}
    for grp in verdict["merge_groups"]:
        for fam in grp:
            group_of[fam] = "+".join(grp)
    rows = []
    for fam in sorted(verdict["fragments"]):
        d = verdict["fragments"][fam]
        rows.append((cluster_id, fam, d["status"], str(d["n_members"]),
                     str(d["n_missing"]), group_of.get(fam, ""), d["reason"]))
    return rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def write_summary_v3(families: Dict[str, List[str]],
                     provenance: Dict[str, str],
                     original: Path, out: Path,
                     cluster_of: Dict[str, str]) -> None:
    """Write the merged family table, carrying merge provenance."""
    meta = {}
    with open(original) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ri = header.index("round") if "round" in header else None
        for line in fh:
            p = line.rstrip("\n").split("\t")
            meta[p[0]] = p[ri] if ri is not None else ""
    with open(out, "w") as f:
        f.write("family_id\tround\tn_genes\tn_species\tgene_list\t"
                "merged_from\tcluster_id\n")
        for fam in sorted(families):
            genes = families[fam]
            species = {g.split("_", 1)[0] for g in genes}
            f.write(f"{fam}\t{meta.get(fam, '')}\t{len(genes)}\t"
                    f"{len(species)}\t{','.join(genes)}\t"
                    f"{provenance.get(fam, '')}\t{cluster_of.get(fam, '')}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True,
                    help="pipeline output dir (holds summary.tsv and "
                         "fragmentation_clusters.tsv)")
    ap.add_argument("--judge", action="store_true",
                    help="align/tree/judge clusters")
    ap.add_argument("--apply", action="store_true",
                    help="fold verdicts into the family table")
    ap.add_argument("--cluster", help="judge only this cluster (array task)")
    ap.add_argument("--pep-dir", help="proteome dir, default <run-dir>/../data/pep")
    ap.add_argument("--cds-dir", help="CDS dir")
    ap.add_argument("--workdir", help="default <run-dir>/validate")
    ap.add_argument("--out", help="apply: default <run-dir>/summary_v3.tsv")
    ap.add_argument("--mafft-bin", default="mafft")
    ap.add_argument("--fasttree-bin", default="FastTree")
    ap.add_argument("--pal2nal", default="pal2nal.pl")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--expect", type=int,
                    help="apply: required cluster count; refuses a partial set")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    run_dir = Path(args.run_dir)
    workdir = Path(args.workdir) if args.workdir else run_dir / "validate"
    workdir.mkdir(parents=True, exist_ok=True)
    verdict_dir = workdir / "verdicts"
    verdict_dir.mkdir(exist_ok=True)

    if not (args.judge or args.apply):
        ap.error("nothing to do: pass --judge and/or --apply")

    families = load_families(run_dir / "summary.tsv")
    clusters = load_clusters(run_dir / "fragmentation_clusters.tsv")
    logger.info("%d families, %d clusters", len(families), len(clusters))

    if args.judge:
        args.mafft_bin = args.mafft_bin
        cfg = argparse.Namespace(
            mafft_bin=args.mafft_bin, fasttree_bin=args.fasttree_bin,
            pal2nal=args.pal2nal, threads=args.threads)
        pep_pool, cds_pool = {}, {}
        pep_dir = Path(args.pep_dir) if args.pep_dir else None
        if pep_dir:
            for fa in sorted(pep_dir.glob("*.fa*")):
                pep_pool.update(read_fasta(str(fa)))
        if args.cds_dir:
            for fa in sorted(Path(args.cds_dir).glob("*.fa*")):
                cds_pool.update(read_fasta(str(fa)))

        todo = [args.cluster] if args.cluster else sorted(clusters)
        for cid in todo:
            out = verdict_dir / f"{cid}.rows.tsv"
            if out.exists() and out.stat().st_size > 0:
                continue
            try:
                rows = judge_cluster(cid, clusters[cid], families,
                                     pep_pool, cds_pool, workdir, cfg)
            except Exception as exc:            # one bad cluster must not
                logger.error("%s failed: %s", cid, exc)   # kill the campaign
                (verdict_dir / f"{cid}.FAILED").write_text(str(exc))
                continue
            out.write_text("".join("\t".join(r) + "\n" for r in rows))
        logger.info("judged: %d verdict files",
                    len(list(verdict_dir.glob("*.rows.tsv"))))

    if args.apply:
        parts = sorted(verdict_dir.glob("*.rows.tsv"))
        expect = args.expect if args.expect is not None else len(clusters)
        if len(parts) != expect:
            sys.exit(f"refusing to apply a partial campaign: {len(parts)} "
                     f"verdict files, expected {expect}. Re-run --judge; a "
                     f"partial merge is indistinguishable from a finished one "
                     f"in the output.")
        rows = []
        for p in parts:
            rows.extend(load_verdict_rows(p))
        combined = run_dir / "cluster_verdicts.tsv"
        with open(combined, "w") as f:
            f.write(VERDICT_HEADER + "\n")
            for r in rows:
                f.write("\t".join(r) + "\n")

        groups = merge_groups_from_rows(rows)
        merged, provenance = apply_merges(families, groups)

        before = sorted(g for m in families.values() for g in m)
        after = sorted(g for m in merged.values() for g in m)
        if after != before:
            sys.exit(f"gene set changed: {len(before)} -> {len(after)}; "
                     f"refusing to write")
        if len(after) != len(set(after)):
            sys.exit("duplicate genes after merge; refusing to write")

        cluster_of = {f: c for c, fams in clusters.items() for f in fams}
        out = Path(args.out) if args.out else run_dir / "summary_v3.tsv"
        write_summary_v3(merged, provenance, run_dir / "summary.tsv", out,
                         cluster_of)
        n_inter = sum(1 for r in rows if r[2] == "INTERLEAVED")
        n_mono = sum(1 for r in rows if r[2] == "MONOPHYLETIC")
        logger.info(
            "INTERLEAVED %d / MONOPHYLETIC %d (undecided) -> %d merge groups; "
            "%d -> %d families, %d genes preserved",
            n_inter, n_mono, len(groups), len(families), len(merged),
            len(after))
        logger.info("wrote %s and %s", out, combined)


if __name__ == "__main__":
    main()
