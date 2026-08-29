"""Freeze the anchor-less subfamily catalog (consultation 2026-08-25).

Per multi-member Possvm OG (codon123 principal split):
  - leaf-universe assertion: codon123 vs codon12 trees vs BOTH Possvm
    universes vs TreeCluster must match, else the family's OGs are
    NOT_EVALUATED (technical absence must not masquerade as biological
    disagreement)
  - codon12 exact-membership recovery (the independent replication axis)
  - TC_NESTED veto: 100% of members inside ONE TreeCluster cluster
    (t=1.0, s=0; TC singletons count as their own clusters); spill and
    coarseness fields kept so nesting under a giant TC cluster is visible
  - own-clade dual supports (SH-aLRT, UFBoot) read from the ORIGINAL
    IQ-TREE treefile via unrooted splits; an OG that is no topological
    split grades UNRESOLVED, an OG that IS a split but carries no
    parseable support label is a technical gap -> NOT_EVALUATED
Grades (pre-registered, #25 2026-08-25 — MUST NOT change between versions):
  HIGH         sh>=80 and uf>=95 and codon12 recovered and TC_NESTED
               and HOG nested (when evaluated)
  PROVISIONAL  HOG not evaluated (capped), or recovered with uf 70-94
  UNRESOLVED   evaluated-and-failed recovery, TC conflict, HOG conflict
               (any conflict OBSERVED among covered members counts, even
               under partial coverage), low support, no topological split
Membership ids (family.OG) are stable; labels come later as a separate
versioned layer and must not move memberships.

2026-08-27 hardening (Codex review, 8 findings — grading thresholds unchanged):
  manifest-driven iteration; missing inputs become persistent NOT_EVALUATED
  rows; atomic tmp-file write to an explicit --out; codon12 universe joined
  to the universe check; topology splits recorded independently of support
  labels; duplicate leaves/assignments are hard failures; iterative
  postorder (no recursion, no quadratic leaf recount); deterministic TC
  parent tie-break.
"""
import argparse
import csv
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path.home() / "subfam" / "pkg"))
from newick import parse_newick

DUAL = re.compile(r"^(\d+(?:\.\d+)?)/(\d+(?:\.\d+)?)$")


def die(msg):
    raise SystemExit(f"FATAL: {msg}")


def postorder(root):
    stack, out = [root], []
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(n.children)
    return reversed(out)


def leafsets(root):
    """node -> frozenset(descendant leaves), iterative; duplicate leaves fatal."""
    ls = {}
    names = []
    for n in postorder(root):
        if not n.children:
            ls[id(n)] = frozenset([n.name])
            names.append(n.name)
        else:
            acc = set()
            for c in n.children:
                acc |= ls[id(c)]
            ls[id(n)] = frozenset(acc)
    if len(names) != len(set(names)):
        die("duplicate leaf names in tree")
    return ls


def split_supports(root, where):
    """(total, topo_splits, side->(sh,uf)). Topology recorded independently
    of the support label so an unlabelled split stays distinguishable from
    a genuine non-split."""
    ls = leafsets(root)
    total = ls[id(root)]
    topo, out = set(), {}
    for n in postorder(root):
        for child in n.children:
            side = ls[id(child)]
            if not 1 < len(side) < len(total):
                continue
            topo.add(side)
            topo.add(total - side)
            m = DUAL.match(child.name or "")
            if m:
                sup = (float(m.group(1)), float(m.group(2)))
                out[side] = sup
                out[total - side] = sup
    return total, topo, out


def possvm_ogs(path):
    ogs, seen = {}, {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g, og = row["gene"], row["orthogroup"]
            if g in seen and seen[g] != og:
                die(f"{path}: gene {g} in two OGs ({seen[g]}, {og})")
            seen[g] = og
            ogs.setdefault(og, set()).add(g)
    return ogs


def tc_assignment(path):
    tc = {}
    singleton = 0
    for line in open(path).read().splitlines()[1:]:
        gene, cluster = line.split("\t")
        if gene in tc:
            die(f"{path}: duplicate gene {gene}")
        if cluster == "-1":
            singleton += 1
            tc[gene] = f"S{singleton}"
        else:
            tc[gene] = cluster
    return tc


def parse_hog(path):
    """gene -> HOG id from an N*.tsv (fixed 3 leading columns)."""
    hog_of = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                die(f"{path}: row width mismatch")
            for cell in parts[3:]:
                for gene in cell.split(","):
                    gene = gene.strip()
                    if gene:
                        if gene in hog_of and hog_of[gene] != parts[0]:
                            die(f"{path}: gene {gene} in two HOGs")
                        hog_of[gene] = parts[0]
    return hog_of


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True,
                    help="one family id per line — the authoritative freeze set")
    ap.add_argument("--out", required=True)
    ap.add_argument("--base", default=str(Path.home() / "subfam" / "subfam_trees"))
    ap.add_argument("--hog-dir", default=str(Path.home() / "subfam" / "hog"))
    args = ap.parse_args()

    base, hog_dir = Path(args.base), Path(args.hog_dir)
    fams = [x.split("\t")[0].strip() for x in open(args.manifest)
            if x.strip() and not x.startswith("#")]
    if len(fams) != len(set(fams)):
        die("manifest contains duplicate family ids")

    tmp = args.out + ".tmp"
    ok = False
    cat = open(tmp, "w")
    cat.write("og_id\tfamily\tn_members\tsh_alrt\tufboot\tcodon12_recovered\t"
              "tc_nested\ttc_clusters_touched\ttc_spill\ttc_parent_size\t"
              "tc_parent_og_count\thog_status\tgrade\tgrade_reason\tmembers\n")

    try:
        fam_rows = []
        counts = {}
        for fam in sorted(fams):
            d = base / fam
            need = [d / "codon123.treefile", d / "codon12.treefile",
                    d / "possvm_codon123.ortholog_groups.csv",
                    d / "possvm_codon12.ortholog_groups.csv", d / "treecluster.txt"]
            missing = [p.name for p in need if not p.exists()]
            if missing:
                fam_rows.append((fam, "MISSING_OUTPUTS"))
                counts["NOT_EVALUATED"] = counts.get("NOT_EVALUATED", 0) + 1
                cat.write("\t".join(map(str, [
                    f"{fam}.NA", fam, "", "", "", "", "", "", "", "", "",
                    "not_evaluated", "NOT_EVALUATED",
                    "missing inputs: " + ",".join(missing), ""])) + "\n")
                continue

            t123 = parse_newick(need[0].read_text().strip())
            t12 = parse_newick(need[1].read_text().strip())
            u123, topo123, sup123 = split_supports(t123, need[0])
            u12 = leafsets(t12)[id(t12)]
            ogs = possvm_ogs(need[2])
            og_universe = frozenset(g for m in ogs.values() for g in m)
            ogs12 = possvm_ogs(need[3])
            og12_universe = frozenset(g for m in ogs12.values() for g in m)
            rec12_sets = set(map(frozenset, ogs12.values()))
            tc = tc_assignment(need[4])

            universes_ok = (u123 == u12 == og_universe == og12_universe
                            and set(tc) == set(u123))
            fam_rows.append((fam, "OK" if universes_ok else "UNIVERSE_MISMATCH"))

            tc_members = {}
            for g, c in tc.items():
                tc_members.setdefault(c, set()).add(g)
            tc_og_count = {}
            for og, members in ogs.items():
                for c in {tc[g] for g in members if g in tc}:
                    tc_og_count[c] = tc_og_count.get(c, 0) + 1

            hog_path = hog_dir / fam / "N_root.tsv"
            hog_of = parse_hog(hog_path) if hog_path.exists() else {}

            for og in sorted(ogs, key=lambda o: int(o[2:]) if o[2:].isdigit() else 0):
                members = ogs[og]
                if len(members) < 2:
                    continue
                fs = frozenset(members)
                touched = {tc[g] for g in members if g in tc}
                parent = max(sorted(touched),
                             key=lambda c: (len(tc_members[c] & members), c),
                             default=None)
                spill = len(members) - (len(tc_members[parent] & members)
                                        if parent else 0)
                nested = (len(touched) == 1
                          and all(g in tc for g in members))
                sup = sup123.get(fs)
                sh, uf = sup if sup else (None, None)

                # A conflict OBSERVED among covered members is disagreement even
                # under partial coverage; only a coverage gap with no observed
                # conflict is a technical absence.
                covered = {hog_of[g] for g in members if g in hog_of}
                if len(covered) > 1:
                    hog_status = "conflict"
                elif covered and all(g in hog_of for g in members):
                    hog_status = "nested"
                else:
                    hog_status = "not_evaluated"   # coverage gap, not disagreement
                if not universes_ok:
                    grade, reason = "NOT_EVALUATED", "leaf universes disagree"
                    recovered = None
                else:
                    # Veto order: every BIOLOGICAL failure (non-split, low support,
                    # codon12, TC, HOG conflict) before the TECHNICAL unlabelled-
                    # split branch — a technical gap must never mask a real failure.
                    recovered = fs in rec12_sets
                    if sup is None and fs not in topo123:
                        grade, reason = "UNRESOLVED", "OG is not a split in its own tree"
                    elif sup is not None and (sh < 80 or uf < 70):
                        grade, reason = "UNRESOLVED", f"support {sh:g}/{uf:g} below floor"
                    elif not recovered:
                        grade, reason = "UNRESOLVED", "codon12 evaluated and did not recover"
                    elif not nested:
                        grade, reason = "UNRESOLVED", "TreeCluster cuts through the OG"
                    elif hog_status == "conflict":
                        grade, reason = "UNRESOLVED", "HOG cuts through the OG"
                    elif sup is None:
                        grade, reason = ("NOT_EVALUATED",
                                         "split present but support label "
                                         "missing/unparseable")
                    elif hog_status == "not_evaluated":
                        grade, reason = ("PROVISIONAL",
                                         "HOG axis not evaluated (leaf coverage gap)")
                    elif uf >= 95:
                        grade, reason = ("HIGH", f"{sh:g}/{uf:g}, codon12 recovered, "
                                         "TC nested, HOG nested")
                    else:
                        grade, reason = "PROVISIONAL", f"recovered but UFBoot {uf:g} in 70-94"
                counts[grade] = counts.get(grade, 0) + 1
                cat.write("\t".join(map(str, [
                    f"{fam}.{og}", fam, len(members),
                    "" if sh is None else sh, "" if uf is None else uf,
                    "" if recovered is None else int(recovered),
                    int(nested), len(touched), spill,
                    len(tc_members[parent]) if parent else "",
                    tc_og_count.get(parent, "") if parent else "",
                    hog_status, grade, reason, ",".join(sorted(members))])) + "\n")

        cat.close()
        os.replace(tmp, args.out)
        ok = True
    finally:
        if not ok:
            cat.close()
            if os.path.exists(tmp):
                os.unlink(tmp)
    bad = [f for f, s in fam_rows if s != "OK"]
    print("families:", len(fam_rows), "universe_ok:", len(fam_rows) - len(bad),
          "flagged:", len(bad), bad[:10], "..." if len(bad) > 10 else "")
    print("grades:", dict(sorted(counts.items())))


if __name__ == "__main__":
    main()
