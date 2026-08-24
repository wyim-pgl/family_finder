#!/usr/bin/env python3
"""Tree-validate family boundaries and apply the merges (issue #47).

`steps/profile_assign.vote_edges` nominates: it says which families look like
pieces of one family. `steps/cluster_validate.fragment_verdict` decides, on
each cluster's own codon tree. Until this file existed, only the nominating
half had a caller - the deciding half ran from ad-hoc scripts on one cluster,
so the published family table could not be reproduced, re-run at another
threshold, or applied to a different dataset.

The verdict rule (see steps/cluster_validate for the reasoning):

  INTERLEAVED   the fragment does not form its own edge-defined clade. It
                merges only with other INTERLEAVED fragments that no edge can
                separate from it; a non-clade internal strip can have no such
                partner and therefore no merge group.
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
import hashlib
import inspect
import logging
import os
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
FINGERPRINT_PREFIX = "# judge_fingerprint\t"
LOGIC_PREFIX = "# judge_logic\t"
MEMBERSHIP_PREFIX = "# judge_membership\t"

# Verdicts that merge, and the full vocabulary. Kept as named sets because the
# apply path once hard-coded INTERLEAVED and silently ignored every outgroup
# verdict - a campaign could have judged thousands of clusters and merged
# nothing, with no error anywhere.
MERGING_STATUSES = {"INTERLEAVED", "SAME_FAMILY"}
VERDICT_STATUSES = MERGING_STATUSES | {"MONOPHYLETIC", "SEPARATE",
                                       "OUTGROUP_SPLIT", "ABSENT"}


# ---------------------------------------------------------------------------
# reading the run's own outputs
# ---------------------------------------------------------------------------

def load_families(summary_path: Path) -> Dict[str, List[str]]:
    """family_id -> member gene ids, from a run's summary.tsv."""
    families: Dict[str, List[str]] = {}
    owner: Dict[str, str] = {}
    with open(summary_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            fi = header.index("family_id")
            gi = header.index("gene_list")
        except ValueError as exc:
            raise ValueError(
                f"{summary_path}: expected family_id and gene_list columns"
            ) from exc
        ni = header.index("n_genes") if "n_genes" in header else None
        need = max(fi, gi, ni if ni is not None else 0)
        for lineno, line in enumerate(fh, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= need:
                raise ValueError(
                    f"{summary_path}:{lineno}: truncated family row"
                )
            family_id = parts[fi]
            genes = [g for g in parts[gi].split(",") if g]
            if not family_id or not genes:
                raise ValueError(
                    f"{summary_path}:{lineno}: empty family id or gene list"
                )
            if family_id in families:
                raise ValueError(
                    f"{summary_path}:{lineno}: duplicate family id {family_id!r}"
                )
            if len(genes) != len(set(genes)):
                raise ValueError(
                    f"{summary_path}:{lineno}: duplicate gene within {family_id}"
                )
            if ni is not None:
                try:
                    declared = int(parts[ni])
                except ValueError as exc:
                    raise ValueError(
                        f"{summary_path}:{lineno}: invalid n_genes {parts[ni]!r}"
                    ) from exc
                if declared != len(genes):
                    raise ValueError(
                        f"{summary_path}:{lineno}: {family_id} declares "
                        f"{declared} genes but lists {len(genes)}"
                    )
            for gene in genes:
                previous = owner.get(gene)
                if previous is not None:
                    raise ValueError(
                        f"{summary_path}:{lineno}: gene {gene!r} occurs in both "
                        f"{previous} and {family_id}"
                    )
                owner[gene] = family_id
            families[family_id] = genes
    return families


def load_clusters(clusters_path: Path) -> Dict[str, List[str]]:
    """cluster_id -> member family ids, from fragmentation_clusters.tsv."""
    clusters: Dict[str, List[str]] = {}
    family_owner: Dict[str, str] = {}
    with open(clusters_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ci = header.index("cluster_id")
            fi = header.index("families")
        except ValueError as exc:
            raise ValueError(
                f"{clusters_path}: expected cluster_id and families columns"
            ) from exc
        ni = header.index("n_families") if "n_families" in header else None
        need = max(ci, fi, ni if ni is not None else 0)
        for lineno, line in enumerate(fh, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= need:
                raise ValueError(
                    f"{clusters_path}:{lineno}: truncated cluster row"
                )
            cluster_id = parts[ci]
            members = [f for f in parts[fi].split(",") if f]
            if not cluster_id or not members:
                raise ValueError(
                    f"{clusters_path}:{lineno}: empty cluster id or family list"
                )
            if cluster_id in clusters:
                raise ValueError(
                    f"{clusters_path}:{lineno}: duplicate cluster id {cluster_id!r}"
                )
            if len(members) != len(set(members)):
                raise ValueError(
                    f"{clusters_path}:{lineno}: duplicate family in {cluster_id}"
                )
            if ni is not None:
                try:
                    declared = int(parts[ni])
                except ValueError as exc:
                    raise ValueError(
                        f"{clusters_path}:{lineno}: invalid n_families "
                        f"{parts[ni]!r}"
                    ) from exc
                if declared != len(members):
                    raise ValueError(
                        f"{clusters_path}:{lineno}: {cluster_id} declares "
                        f"{declared} families but lists {len(members)}"
                    )
            for family_id in members:
                previous = family_owner.get(family_id)
                if previous is not None:
                    raise ValueError(
                        f"{clusters_path}:{lineno}: family {family_id!r} occurs "
                        f"in both {previous} and {cluster_id}"
                    )
                family_owner[family_id] = cluster_id
            clusters[cluster_id] = members
    return clusters


def load_verdict_rows(path: Path) -> List[tuple]:
    """Verdict rows as tuples; comments and an optional header are skipped."""
    rows = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            if line.startswith("cluster_id\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 7:
                raise ValueError(
                    f"{path}:{lineno}: expected 7 verdict columns, got "
                    f"{len(parts)}"
                )
            rows.append(tuple(parts))
    return rows


def validate_cluster_membership(
    clusters: Dict[str, List[str]], families: Dict[str, List[str]],
) -> None:
    """Refuse clusters that name a family absent from the source summary."""
    unknown = sorted(
        (cluster_id, family_id)
        for cluster_id, members in clusters.items()
        for family_id in members
        if family_id not in families
    )
    if unknown:
        preview = "; ".join(f"{c}: {f}" for c, f in unknown[:10])
        raise ValueError(
            "fragmentation clusters name families absent from summary.tsv: "
            + preview
        )


def load_sequence_pool(directory: Path, label: str) -> Dict[str, str]:
    """Load FASTAs without silently overwriting duplicate sequence ids."""
    pool: Dict[str, str] = {}
    source: Dict[str, Path] = {}
    for fasta in sorted(Path(directory).glob("*.fa*")):
        seqs = read_fasta(str(fasta))
        overlap = sorted(set(pool) & set(seqs))
        if overlap:
            gene = overlap[0]
            raise ValueError(
                f"duplicate {label} id {gene!r} in {source[gene]} and {fasta}"
            )
        pool.update(seqs)
        source.update({gene: fasta for gene in seqs})
    return pool


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
        if status not in MERGING_STATUSES or not group:
            continue
        members = sorted(f for f in group.split("+") if f)
        if len(members) < 2:
            continue
        key = tuple(members)
        if key not in seen:
            seen.add(key)
            groups.append(members)
    return groups


def validate_verdict_coverage(
    rows: Sequence[tuple], clusters: Dict[str, List[str]],
    families: Dict[str, List[str]] = None,
) -> None:
    """Refuse incomplete, inconsistent, or cross-cluster verdict rows."""
    expected = {cid: set(fams) for cid, fams in clusters.items()}
    seen: Dict[str, set] = {cid: set() for cid in clusters}
    rows_by_fragment = {}

    for row in rows:
        if len(row) != 7:
            raise ValueError(
                f"malformed verdict row: expected 7 columns, got {len(row)}"
            )
        cluster_id, fam = row[0], row[1]
        if cluster_id not in expected:
            raise ValueError(
                f"verdict names unknown cluster {cluster_id!r}; verdicts and "
                f"fragmentation_clusters.tsv are from different runs"
            )
        if fam not in expected[cluster_id]:
            raise ValueError(
                f"verdict names {fam!r} in {cluster_id}, but that family is "
                f"not a member of the cluster"
            )
        if fam in seen[cluster_id]:
            raise ValueError(
                f"duplicate verdict row for {fam!r} in {cluster_id}"
            )
        status = row[2]
        if status not in VERDICT_STATUSES:
            raise ValueError(
                f"invalid verdict status {status!r} for {fam!r} in {cluster_id}"
            )
        try:
            n_members = int(row[3])
            n_missing = int(row[4])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"invalid member counts for {fam!r} in {cluster_id}"
            ) from exc
        if n_members < 1 or n_missing < 0:
            raise ValueError(
                f"invalid member counts for {fam!r} in {cluster_id}: "
                f"n_members={n_members}, n_missing={n_missing}"
            )
        if n_missing:
            raise ValueError(
                f"verdict for {fam!r} in {cluster_id} omitted {n_missing} "
                "expected tree member(s); refusing a partial-tree judgement"
            )
        if families is not None and n_members != len(families[fam]):
            raise ValueError(
                f"stale verdict for {fam!r} in {cluster_id}: records "
                f"{n_members} members, summary.tsv has {len(families[fam])}"
            )

        group = tuple(sorted(f for f in row[5].split("+") if f))
        if row[5] and (len(group) < 2 or fam not in group):
            raise ValueError(
                f"invalid merge group {row[5]!r} for {fam!r} in {cluster_id}"
            )
        outside = sorted(set(group) - expected[cluster_id])
        if outside:
            raise ValueError(
                f"merge group for {fam!r} in {cluster_id} crosses the cluster "
                f"boundary: {','.join(outside)}"
            )
        if status not in MERGING_STATUSES and group:
            raise ValueError(
                f"{status} fragment {fam!r} in {cluster_id} carries a merge "
                "group; only " + "/".join(sorted(MERGING_STATUSES)) +
                " may name one"
            )
        seen[cluster_id].add(fam)
        rows_by_fragment[(cluster_id, fam)] = (status, group)

    missing = []
    for cluster_id, fams in expected.items():
        absent = sorted(fams - seen[cluster_id])
        if absent:
            missing.append(f"{cluster_id}: {','.join(absent)}")
    if missing:
        raise ValueError(
            "incomplete verdict set: missing fragment rows for "
            + "; ".join(missing[:10])
        )

    # Every named merge group is a connected component emitted on every one of
    # its INTERLEAVED member rows. A truncated or hand-edited group must not
    # quietly merge fragments that did not carry the same decision.
    for (cluster_id, fam), (status, group) in rows_by_fragment.items():
        if not group:
            continue
        for member in group:
            other_status, other_group = rows_by_fragment[(cluster_id, member)]
            if (other_status != status or other_group != group):
                raise ValueError(
                    f"inconsistent merge group {row_group(group)!r} in "
                    f"{cluster_id}: {member!r} does not carry the same "
                    f"{status} decision"
                )


def row_group(group: Sequence[str]) -> str:
    """Stable display form for an already parsed verdict merge group."""
    return "+".join(group)


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


def load_vote_edges(path: Path) -> List[tuple]:
    """(from, to, votes, from_size, frac) rows from vote_edges.tsv."""
    edges = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                edges.append((p[0], p[1], int(p[2]), int(p[3]), float(p[4])))
    return edges


def pick_outgroup_families(cluster: Sequence[str], edges: Sequence[tuple],
                           n: int = 3) -> List[str]:
    """The nearest families OUTSIDE the cluster, by vote fraction.

    The tree can only answer "are these one family" when it has something
    from outside the family to compare against; without that, a subfamily is
    a clade and the verdict is permanently MONOPHYLETIC. vote_edges already
    ranks each family's neighbours, so the outgroup costs nothing extra: the
    cluster is the connected component, and the best-scoring families just
    outside it are its near-outside.

    Members of the cluster are never chosen. An outgroup drawn from inside
    the family makes the test circular.

    ⚠️ Pass the UNCUT edge dump, not the one the clusters were built from. A
    cluster is a connected component of those edges, so by construction none
    of them leaves it: measured on the shipped 15sp file, all four edges out
    of C0297's five PEPC fragments point back inside, and the file's minimum
    frac is exactly 0.600 - the component threshold. The outgroup lives in
    the edges BELOW that threshold, which only `vote_edges(..., min_frac=0)`
    keeps. Given a thresholded file this returns nothing, every cluster falls
    back to the topology rule, and the campaign quietly reverts to permanent
    MONOPHYLETIC - which is why the caller counts the empty ones and says so.
    """
    members = set(cluster)
    best: Dict[str, float] = {}
    for fam_from, fam_to, _votes, _size, frac in edges:
        if fam_from not in members or fam_to in members:
            continue
        if frac > best.get(fam_to, -1.0):
            best[fam_to] = frac
    ranked = sorted(best, key=lambda f: (-best[f], f))
    return ranked[:n]


def resolve_outgroup(cluster: Sequence[str], edges: Sequence[tuple],
                     n_from_edges: int = 0,
                     explicit: Sequence[str] = None) -> List[str]:
    """The outgroup for one cluster: an explicit panel, or the graph fallback.

    Prefer explicit. A below-cut graph neighbour is exactly where membership
    is most ambiguous, so it is the least defensible place to manufacture a
    known negative; an externally justified panel is the defensible one.

    PEPC has such a panel and it needs no curation beyond what is already
    here. Scored against the four Arabidopsis anchors using 15-species
    sequences alone, the five PTPC fragments sit 700-1000 bits nearer
    PPC1/PPC2/PPC3 while R1_OG0009826 sits 812 bits nearer PPC4 - one sign
    flip across 128 genes, not a gradient. PPC4 is the bacterial-type PEPC
    (BTPC), a class that diverged from plant-type before land plants, and
    steps' own anchors.tsv already labels ATH_AT1G68750.1 as bacterial. It is
    an independently justified negative for the PTPC family.
    """
    if explicit:
        overlap = sorted(set(explicit) & set(cluster))
        if overlap:
            raise ValueError(
                f"outgroup families {overlap} are cluster members - an "
                f"outgroup drawn from inside the family is circular")
        return list(explicit)
    if n_from_edges and edges:
        return pick_outgroup_families(cluster, edges, n_from_edges)
    return []


def judge_cluster(cluster_id: str, fam_ids: Sequence[str],
                  families: Dict[str, List[str]], pep_pool, cds_pool,
                  outdir: Path, cfg,
                  outgroup_families: Sequence[str] = ()) -> List[tuple]:
    """Align/tree/judge one cluster and return its verdict rows.

    `outgroup_families` name families OUTSIDE the cluster whose members go
    into the alignment and tree but get no verdict row of their own. They are
    what lets the tree say SAME_FAMILY or SEPARATE instead of the permanent
    MONOPHYLETIC silence - see steps/cluster_validate._outgroup_verdict.
    """
    workdir = outdir / cluster_id
    workdir.mkdir(parents=True, exist_ok=True)
    fragments = {f: families[f] for f in fam_ids}
    overlap = set(outgroup_families) & set(fam_ids)
    if overlap:
        raise ValueError(
            f"{cluster_id}: outgroup families {sorted(overlap)} are also "
            f"cluster members - the test would be circular")
    outgroup_genes = [g for f in outgroup_families for g in families[f]]
    if len(outgroup_genes) != len(set(outgroup_genes)):
        raise ValueError(f"{cluster_id}: an outgroup gene is listed twice")
    wanted = [g for m in fragments.values() for g in m]
    if len(wanted) != len(set(wanted)):
        raise ValueError(
            f"{cluster_id}: a gene belongs to more than one fragment"
        )

    missing_pep = sorted((set(wanted) | set(outgroup_genes)) - set(pep_pool))
    if missing_pep:
        raise ValueError(
            f"{cluster_id}: {len(missing_pep)} family member peptide(s) are "
            f"absent from the peptide pool: {','.join(missing_pep[:10])}"
        )

    in_tree = wanted + outgroup_genes
    pep = workdir / "cluster.pep.fa"
    write_fasta({g: pep_pool[g] for g in in_tree}, str(pep))
    cds = workdir / "cluster.cds.fa"
    have_cds = {g: cds_pool[g] for g in in_tree if g in cds_pool}
    use_cds = len(have_cds) == len(in_tree)
    if use_cds:
        write_fasta(have_cds, str(cds))

    tree = build_cluster_tree(pep, cds if use_cds else None, workdir, cfg)
    verdict = fragment_verdict(tree.read_text(), fragments,
                               outgroup=outgroup_genes)
    missing_tree = {
        fam: d["n_missing"]
        for fam, d in verdict["fragments"].items()
        if d["n_missing"]
    }
    if missing_tree or verdict["n_tree_leaves_unclaimed"]:
        detail = ", ".join(f"{fam}={n}" for fam, n in sorted(missing_tree.items()))
        raise RuntimeError(
            f"{cluster_id}: tree leaf set does not equal the requested gene "
            f"set (missing: {detail or 'none'}; unclaimed: "
            f"{verdict['n_tree_leaves_unclaimed']})"
        )

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


def _logic_digest() -> str:
    """Digest of the code whose result a resumable verdict file caches."""
    source = "\n".join(
        inspect.getsource(fn)
        for fn in (build_cluster_tree, judge_cluster, fragment_verdict)
    )
    return hashlib.sha256(source.encode("utf-8")).hexdigest()


def judge_fingerprint(cluster_id: str, fam_ids: Sequence[str],
                      families: Dict[str, List[str]], pep_pool, cds_pool,
                      cfg) -> str:
    """Fingerprint inputs, tool settings, and judgement logic for resume."""
    digest = hashlib.sha256()

    def add(value) -> None:
        data = str(value).encode("utf-8")
        digest.update(str(len(data)).encode("ascii") + b":" + data)

    add(_logic_digest())
    add(cluster_id)
    for attr in ("mafft_bin", "fasttree_bin", "pal2nal", "threads"):
        add(attr)
        add(getattr(cfg, attr))
    for family_id in sorted(fam_ids):
        add(family_id)
        for gene in sorted(families[family_id]):
            add(gene)
            add(pep_pool.get(gene, "<MISSING>"))
            add(cds_pool.get(gene, "<MISSING>"))
    return digest.hexdigest()


def cluster_membership_fingerprint(
    cluster_id: str, fam_ids: Sequence[str], families: Dict[str, List[str]],
) -> str:
    """Digest the exact family/gene membership a verdict will be applied to."""
    digest = hashlib.sha256()

    def add(value) -> None:
        data = str(value).encode("utf-8")
        digest.update(str(len(data)).encode("ascii") + b":" + data)

    add(cluster_id)
    for family_id in sorted(fam_ids):
        add(family_id)
        for gene in sorted(families[family_id]):
            add(gene)
    return digest.hexdigest()


def _read_stamp(path: Path, prefix: str) -> str:
    if not path.exists():
        return ""
    with open(path) as fh:
        for line in fh:
            if line.startswith(prefix):
                return line[len(prefix):].strip()
            if line and not line.startswith("#"):
                break
    return ""


def verdict_file_is_current(path: Path, fingerprint: str,
                            cluster_id: str, fam_ids: Sequence[str],
                            families: Dict[str, List[str]]) -> bool:
    """A cache hit needs both a matching stamp and a complete valid row set."""
    if not path.exists() or path.stat().st_size == 0:
        return False
    if _read_stamp(path, FINGERPRINT_PREFIX) != fingerprint:
        return False
    if _read_stamp(path, LOGIC_PREFIX) != _logic_digest():
        return False
    expected_membership = cluster_membership_fingerprint(
        cluster_id, fam_ids, families,
    )
    if _read_stamp(path, MEMBERSHIP_PREFIX) != expected_membership:
        return False
    try:
        rows = load_verdict_rows(path)
        validate_verdict_coverage(
            rows, {cluster_id: list(fam_ids)}, families=families,
        )
    except (KeyError, ValueError):
        return False
    return True


def write_verdict_rows(path: Path, rows: Sequence[tuple], fingerprint: str,
                       membership_fingerprint: str) -> None:
    """Atomically write a stamped verdict cache file."""
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with open(tmp, "w") as fh:
            fh.write(f"{LOGIC_PREFIX}{_logic_digest()}\n")
            fh.write(f"{FINGERPRINT_PREFIX}{fingerprint}\n")
            fh.write(f"{MEMBERSHIP_PREFIX}{membership_fingerprint}\n")
            fh.write("".join("\t".join(r) + "\n" for r in rows))
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def archive_stale_verdict(path: Path) -> None:
    """Move an invalid cache aside so --apply cannot consume it by accident."""
    if not path.exists():
        return
    suffix = _read_stamp(path, FINGERPRINT_PREFIX)[:12] or "unstamped"
    stale = path.with_name(f"{path.name}.stale-{suffix}")
    counter = 1
    while stale.exists():
        stale = path.with_name(f"{path.name}.stale-{suffix}-{counter}")
        counter += 1
    path.replace(stale)


def validate_verdict_file_stamps(
    path: Path, rows: Sequence[tuple], clusters: Dict[str, List[str]],
    families: Dict[str, List[str]], allow_unstamped: bool = False,
) -> None:
    """Verify that an apply input was judged for this logic and membership."""
    logic = _read_stamp(path, LOGIC_PREFIX)
    if not logic:
        if allow_unstamped:
            return
        raise ValueError(
            f"refusing unstamped legacy verdict {path}; re-run --judge or "
            "pass --allow-unstamped-verdicts explicitly"
        )
    if logic != _logic_digest():
        raise ValueError(f"refusing stale verdict logic in {path}; re-run --judge")

    cluster_ids = {row[0] for row in rows}
    if len(cluster_ids) != 1:
        raise ValueError(
            f"{path}: expected verdict rows for exactly one cluster, got "
            f"{sorted(cluster_ids)}"
        )
    cluster_id = next(iter(cluster_ids))
    if cluster_id not in clusters:
        raise ValueError(f"{path}: verdict names unknown cluster {cluster_id!r}")
    expected = cluster_membership_fingerprint(
        cluster_id, clusters[cluster_id], families,
    )
    recorded = _read_stamp(path, MEMBERSHIP_PREFIX)
    if recorded != expected:
        raise ValueError(
            f"refusing stale verdict membership in {path}; summary.tsv or "
            "fragmentation_clusters.tsv changed, so re-run --judge"
        )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def write_summary_v3(families: Dict[str, List[str]],
                     provenance: Dict[str, str],
                     original: Path, out: Path,
                     cluster_of: Dict[str, str],
                     species_delimiter: str = "_") -> None:
    """Write the merged family table, carrying merge provenance."""
    if not species_delimiter:
        raise ValueError("species delimiter must not be empty")
    meta = {}
    with open(original) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if "family_id" not in header or "round" not in header:
            raise ValueError(
                f"{original}: family_id and round metadata are required"
            )
        fi, ri = header.index("family_id"), header.index("round")
        for lineno, line in enumerate(fh, start=2):
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(fi, ri) or not p[ri]:
                raise ValueError(
                    f"{original}:{lineno}: missing family id or round metadata"
                )
            if p[fi] in meta:
                raise ValueError(
                    f"{original}:{lineno}: duplicate family id {p[fi]!r}"
                )
            meta[p[fi]] = p[ri]
    with open(out, "w") as f:
        f.write("family_id\tround\tn_genes\tn_species\tgene_list\t"
                "merged_from\tcluster_id\n")
        for fam in sorted(families):
            genes = families[fam]
            if fam not in meta:
                raise ValueError(
                    f"merged family {fam!r} has no round metadata in {original}"
                )
            missing_delimiter = [g for g in genes if species_delimiter not in g]
            if missing_delimiter:
                raise ValueError(
                    f"cannot count species for {fam}: gene id "
                    f"{missing_delimiter[0]!r} lacks delimiter "
                    f"{species_delimiter!r}"
                )
            species = {g.split(species_delimiter, 1)[0] for g in genes}
            f.write(f"{fam}\t{meta[fam]}\t{len(genes)}\t"
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
    ap.add_argument("--pep-dir", help="proteome dir (required with --judge)")
    ap.add_argument("--cds-dir", help="CDS dir")
    ap.add_argument("--species-delimiter", default="_",
                    help="gene-id species delimiter used by the source run")
    ap.add_argument("--workdir", help="default <run-dir>/validate")
    ap.add_argument("--out", help="apply: default <run-dir>/summary_v3.tsv")
    ap.add_argument("--mafft-bin", default="mafft")
    ap.add_argument("--fasttree-bin", default="FastTree")
    ap.add_argument("--pal2nal", default="pal2nal.pl")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument(
        "--outgroup-from-edges", metavar="N", type=int, default=0,
        help="judge: add the N nearest families OUTSIDE each cluster to its "
             "tree as an outgroup, taken from vote_edges.tsv. Without one, a "
             "subfamily is a clade and the verdict is permanently "
             "MONOPHYLETIC; with one the tree can say SAME_FAMILY or "
             "SEPARATE. Use several (3+): a single outgroup leaf is "
             "separated from everything by its own pendant edge and so "
             "discriminates nothing.")
    ap.add_argument("--vote-edges",
                    help="judge: default <run-dir>/vote_edges.tsv")
    ap.add_argument(
        "--outgroup-families", nargs="+", metavar="FAMILY",
        help="judge: use these families as the outgroup for every cluster, "
             "overriding --outgroup-from-edges. This is the defensible form: "
             "an externally justified negative, not a below-cut graph "
             "neighbour. For PEPC that is the BTPC family (R1_OG0009826).")
    ap.add_argument("--expect", type=int,
                    help="apply: required cluster count; refuses a partial set")
    ap.add_argument(
        "--allow-unstamped-verdicts", action="store_true",
        help="apply legacy verdict files that predate resume fingerprints; "
             "never permits a stamped verdict from different logic",
    )
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    run_dir = Path(args.run_dir)
    workdir = Path(args.workdir) if args.workdir else run_dir / "validate"
    workdir.mkdir(parents=True, exist_ok=True)
    verdict_dir = workdir / "verdicts"
    verdict_dir.mkdir(exist_ok=True)

    if not (args.judge or args.apply):
        ap.error("nothing to do: pass --judge and/or --apply")

    try:
        families = load_families(run_dir / "summary.tsv")
        clusters = load_clusters(run_dir / "fragmentation_clusters.tsv")
        validate_cluster_membership(clusters, families)
    except ValueError as exc:
        sys.exit(str(exc))
    logger.info("%d families, %d clusters", len(families), len(clusters))

    if args.judge:
        if not args.pep_dir:
            sys.exit("--judge requires --pep-dir; an empty peptide pool would "
                     "turn missing genes into partial-tree verdicts")
        args.mafft_bin = args.mafft_bin
        cfg = argparse.Namespace(
            mafft_bin=args.mafft_bin, fasttree_bin=args.fasttree_bin,
            pal2nal=args.pal2nal, threads=args.threads)
        pep_dir = Path(args.pep_dir)
        try:
            pep_pool = load_sequence_pool(pep_dir, "peptide")
            cds_pool = (
                load_sequence_pool(Path(args.cds_dir), "CDS")
                if args.cds_dir else {}
            )
        except ValueError as exc:
            sys.exit(str(exc))

        if args.cluster and args.cluster not in clusters:
            sys.exit(f"unknown cluster {args.cluster!r}")
        vote_edges = []
        if args.outgroup_from_edges:
            edges_path = (Path(args.vote_edges) if args.vote_edges
                          else run_dir / "vote_edges.tsv")
            if not edges_path.exists():
                sys.exit(f"--outgroup-from-edges needs {edges_path}; run the "
                         f"merge scan first or pass --vote-edges")
            vote_edges = load_vote_edges(edges_path)
            logger.info("outgroup: up to %d families per cluster from %d "
                        "vote edges", args.outgroup_from_edges,
                        len(vote_edges))

        todo = [args.cluster] if args.cluster else sorted(clusters)
        failed = []
        no_outgroup = []
        for cid in todo:
            out = verdict_dir / f"{cid}.rows.tsv"
            try:
                fingerprint = judge_fingerprint(
                    cid, clusters[cid], families, pep_pool, cds_pool, cfg,
                )
                membership_fingerprint = cluster_membership_fingerprint(
                    cid, clusters[cid], families,
                )
                if verdict_file_is_current(
                    out, fingerprint, cid, clusters[cid], families,
                ):
                    continue
                archive_stale_verdict(out)
                og = resolve_outgroup(clusters[cid], vote_edges,
                                      args.outgroup_from_edges,
                                      args.outgroup_families)
                if args.outgroup_from_edges and not args.outgroup_families \
                        and not og:
                    no_outgroup.append(cid)
                rows = judge_cluster(cid, clusters[cid], families,
                                     pep_pool, cds_pool, workdir, cfg,
                                     outgroup_families=og)
            except Exception as exc:            # one bad cluster must not
                logger.error("%s failed: %s", cid, exc)   # kill the campaign
                (verdict_dir / f"{cid}.FAILED").write_text(str(exc))
                failed.append(cid)
                continue
            write_verdict_rows(
                out, rows, fingerprint, membership_fingerprint,
            )
            (verdict_dir / f"{cid}.FAILED").unlink(missing_ok=True)
        if args.outgroup_from_edges and no_outgroup:
            share = 100.0 * len(no_outgroup) / max(len(todo), 1)
            logger.error(
                "%d of %d clusters (%.1f%%) got NO outgroup and fell back to "
                "the topology rule - they can only return MONOPHYLETIC. A "
                "cluster is a connected component of the vote edges, so a "
                "thresholded edge file cannot supply one: re-run the merge "
                "scan so vote_edges.tsv keeps the sub-threshold edges.",
                len(no_outgroup), len(todo), share)
        logger.info("judged: %d verdict files",
                    len(list(verdict_dir.glob("*.rows.tsv"))))
        if failed:
            sys.exit(
                f"{len(failed)} cluster(s) failed judgement: "
                + ",".join(failed[:10])
            )

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
            try:
                file_rows = load_verdict_rows(p)
                validate_verdict_file_stamps(
                    p, file_rows, clusters, families,
                    allow_unstamped=args.allow_unstamped_verdicts,
                )
            except ValueError as exc:
                sys.exit(str(exc))
            rows.extend(file_rows)
        try:
            validate_verdict_coverage(rows, clusters, families=families)
        except ValueError as exc:
            sys.exit(str(exc))
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
                         cluster_of, args.species_delimiter)
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
