#!/usr/bin/env python3
"""PLAZA/AFDB naming pilot — cohort, orthology QC gates, structural QC, report.

The pilot decides two things BEFORE any bulk spend (design review 2026-08-25,
issues #25/#47): whether the full-panel DIAMOND->PLAZA naming layer is safe
enough to run, and whether the structural branch (AFDB downloads) deserves
disk. Every acceptance is gated and every rejection carries a reason —
abstention is an answer, silent fallback is not.

Subcommands:
  cohort        build the pilot cohort (23,744 family representatives +
                <=100/species comparable unplaced), deterministic
  qc-orthology  apply the SAFE_ORTHOLOGY gate contract to DIAMOND results,
                grid-calibrate the tunable gates, measure conflicts
  qc-afdb       coverage of named structural transfer from a foldseek m8
                against a checksum-pinned UniProt idmapping snapshot
  report        fold both QC outputs into the measurement table and evaluate
                the pre-committed decision rules (a)-(d)

Gate contract (starting values; pident/cov/margin are calibrated by the
pilot, RBH and unique-AGI mapping are non-negotiable):
  forward E <= 1e-10, pident >= 30, qcovhsp >= 70, scovhsp >= 70,
  relative bits margin >= 0.10 vs the best hit to a DIFFERENT AGI locus,
  species-scoped reciprocal best hit, unique AGI mapping.
"""

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from classify_unplaced import assess_cds, is_comparable, species_of
from extract_plaza_orthology import strip_transcript
from utils.seqio import read_fasta, write_fasta

PILOT_SALT = "plaza-afdb-pilot-v1"
PER_SPECIES_CAP = 100

DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                  "qcovhsp", "scovhsp"]


# ---------------------------------------------------------------- cohort --

def pick_representatives(families: Dict[str, Sequence[str]],
                         lengths: Dict[str, int]) -> Dict[str, str]:
    """family_id -> representative gene: longest member, ties break to the
    lexicographically smallest id. Refuses a family whose members carry no
    length (a rep that is not in the panel would silently shrink the pilot).
    """
    reps: Dict[str, str] = {}
    for family_id in sorted(families):
        members = families[family_id]
        if not members:
            raise ValueError(f"family {family_id!r} has no members")
        missing = [g for g in members if g not in lengths]
        if missing:
            raise ValueError(
                f"family {family_id!r}: member(s) missing from the panel "
                f"peptide pool: {sorted(missing)[:3]}"
            )
        reps[family_id] = min(members, key=lambda g: (-lengths[g], g))
    return reps


def _sample_key(gene_id: str) -> str:
    return hashlib.sha256(f"{PILOT_SALT}\t{gene_id}".encode()).hexdigest()


def sample_comparable_unplaced(
    eligible_by_species: Dict[str, Sequence[str]],
    per_species: int = PER_SPECIES_CAP,
) -> List[str]:
    """Deterministic per-species sample: salted-hash order, first N."""
    picked: List[str] = []
    for species in sorted(eligible_by_species):
        ordered = sorted(eligible_by_species[species], key=_sample_key)
        picked.extend(ordered[:per_species])
    return picked


# ------------------------------------------------------ orthology gates --

@dataclass(frozen=True)
class Gates:
    evalue: float = 1e-10
    pident: float = 30.0
    qcov: float = 70.0
    scov: float = 70.0
    margin: float = 0.10   # relative bits margin vs best DIFFERENT-AGI hit


def parse_diamond_tsv(lines: Iterable[str]) -> Dict[str, List[dict]]:
    """outfmt-6 rows (DIAMOND_FIELDS order) -> query -> hits, input order."""
    hits: Dict[str, List[dict]] = {}
    for lineno, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < len(DIAMOND_FIELDS):
            raise ValueError(f"line {lineno}: {len(parts)} columns, "
                             f"expected {len(DIAMOND_FIELDS)}")
        row = dict(zip(DIAMOND_FIELDS, parts))
        for key in ("pident", "evalue", "bitscore", "qcovhsp", "scovhsp"):
            row[key] = float(row[key])
        hits.setdefault(row["qseqid"], []).append(row)
    return hits


def safe_orthology(
    query: str,
    query_hits: Sequence[dict],
    reverse_best: Dict[Tuple[str, str], str],
    gates: Gates,
) -> dict:
    """One query's verdict under the gate contract.

    reverse_best maps (species_prefix, AGI) -> the panel gene that AGI's
    proteome-scoped reverse search ranks first. The margin compares the best
    hit against the best hit to a DIFFERENT AGI locus: isoforms of one locus
    must not count as competition, and a lone-locus hit passes by absence of
    a competitor, not by an arbitrary margin over nothing.
    """
    def fail(reason: str, **extra) -> dict:
        return {"query": query, "status": reason, "agi": None, **extra}

    if not query_hits:
        return fail("NO_HIT")
    best = max(query_hits, key=lambda h: (h["bitscore"], -h["evalue"]))
    agi = strip_transcript(best["sseqid"])
    if best["evalue"] > gates.evalue:
        return fail("EVALUE", evalue=best["evalue"])
    if best["pident"] < gates.pident:
        return fail("PIDENT", pident=best["pident"])
    if best["qcovhsp"] < gates.qcov:
        return fail("QCOV", qcov=best["qcovhsp"])
    if best["scovhsp"] < gates.scov:
        return fail("SCOV", scov=best["scovhsp"])
    rivals = [h["bitscore"] for h in query_hits
              if strip_transcript(h["sseqid"]) != agi]
    if rivals:
        second = max(rivals)
        rel = (best["bitscore"] - second) / best["bitscore"]
        if rel < gates.margin:
            return fail("MARGIN", margin=round(rel, 4), rival_bits=second)
    partner = reverse_best.get((species_of(query), agi))
    if partner is None:
        return fail("NO_REVERSE_HIT", agi_candidate=agi)
    if partner != query:
        return fail("NOT_RBH", agi_candidate=agi, reverse_partner=partner)
    return {"query": query, "status": "SAFE", "agi": agi,
            "pident": best["pident"], "bits": best["bitscore"],
            "qcov": best["qcovhsp"], "scov": best["scovhsp"]}


def wilson_upper(k: int, n: int, z: float = 1.96) -> float:
    """Upper bound of the Wilson score interval for k/n."""
    if n == 0:
        return 1.0
    p = k / n
    denom = 1 + z * z / n
    centre = p + z * z / (2 * n)
    spread = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (centre + spread) / denom


def conflict_rate(safe_calls: Dict[str, str],
                  direct: Dict[str, str]) -> dict:
    """Exact-description disagreement between direct and gated-orthology
    labels, over the genes carrying both (family_naming's rule)."""
    overlap = sorted(set(safe_calls) & set(direct))
    conflicts = [g for g in overlap if safe_calls[g] != direct[g]]
    n, k = len(overlap), len(conflicts)
    return {"n_overlap": n, "n_conflict": k,
            "rate": (k / n) if n else None,
            "wilson_upper": wilson_upper(k, n) if n else None,
            "conflicting_genes": conflicts}


PIDENT_GRID = (25.0, 30.0, 35.0, 40.0)
COV_GRID = (60.0, 70.0, 80.0)
MARGIN_GRID = (0.05, 0.10, 0.20)


def grid_calibrate(
    hits: Dict[str, List[dict]],
    reverse_best: Dict[Tuple[str, str], str],
    direct: Dict[str, str],
    descriptions: Dict[str, str],
    base: Gates = Gates(),
) -> List[dict]:
    """Per gate combination: safe count, conflict rate, Wilson upper.

    RBH and unique-AGI are not on the grid — they are safety requirements,
    not tunables (design review 2026-08-25).
    """
    rows = []
    for pident in PIDENT_GRID:
        for cov in COV_GRID:
            for margin in MARGIN_GRID:
                gates = replace(base, pident=pident, qcov=cov, scov=cov,
                               margin=margin)
                calls = {}
                for q, qh in hits.items():
                    verdict = safe_orthology(q, qh, reverse_best, gates)
                    if verdict["status"] == "SAFE":
                        desc = descriptions.get(verdict["agi"])
                        if desc:
                            calls[q] = desc
                stats = conflict_rate(calls, direct)
                rows.append({"pident": pident, "cov": cov, "margin": margin,
                             "n_safe": len(calls), **{
                                 k: stats[k] for k in
                                 ("n_overlap", "n_conflict", "rate",
                                  "wilson_upper")}})
    return rows


# ------------------------------------------------------ uniprot mapping --

def load_idmapping(lines: Iterable[str]) -> Dict[str, dict]:
    """ARATH idmapping.dat -> uniprot accession -> {'agi': set, 'names': set}.

    Only TAIR (AGI loci) and Gene_Name rows are consumed; everything else in
    the file is other databases' ids.
    """
    out: Dict[str, dict] = {}
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 3:
            continue
        acc, kind, value = parts
        rec = out.setdefault(acc, {"agi": set(), "names": set()})
        if kind == "TAIR":
            rec["agi"].add(strip_transcript(value))
        elif kind == "Gene_Name":
            rec["names"].add(value)
    return out


def resolve_uniprot(acc: str, mapping: Dict[str, dict]) -> Tuple[Optional[str], str]:
    """(AGI, status) under the one-to-many policy: unique AGI accepted,
    zero -> NO_AGI_MAP, several -> AMBIGUOUS_AGI."""
    rec = mapping.get(acc)
    if rec is None or not rec["agi"]:
        return None, "NO_AGI_MAP"
    if len(rec["agi"]) > 1:
        return None, "AMBIGUOUS_AGI"
    return next(iter(rec["agi"])), "OK"


def afdb_accession(target: str) -> Optional[str]:
    """AF-P10490-F1-model_v6 -> P10490 (None when the shape is foreign)."""
    parts = target.split("-")
    if len(parts) >= 3 and parts[0] == "AF":
        return parts[1]
    return None


FOLDSEEK_FIELDS = ["query", "target", "fident", "alnlen", "qcov", "tcov",
                   "evalue", "bits"]


def afdb_named_calls(lines: Iterable[str], mapping: Dict[str, dict],
                     qcov_min: float = 0.70, tcov_min: float = 0.70) -> dict:
    """Per query: best AFDB hit meeting the coverage floors that resolves to
    a uniquely-mapped, named accession. Counts every abstention reason."""
    best: Dict[str, dict] = {}
    for lineno, line in enumerate(lines, start=1):
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < len(FOLDSEEK_FIELDS):
            raise ValueError(f"line {lineno}: truncated foldseek row")
        row = dict(zip(FOLDSEEK_FIELDS, parts))
        for key in ("fident", "qcov", "tcov", "evalue", "bits"):
            row[key] = float(row[key])
        cur = best.get(row["query"])
        if cur is None or row["bits"] > cur["bits"]:
            best[row["query"]] = row

    calls: Dict[str, dict] = {}
    reasons: Dict[str, int] = {}

    def count(reason: str) -> None:
        reasons[reason] = reasons.get(reason, 0) + 1

    for query, row in best.items():
        if row["qcov"] < qcov_min or row["tcov"] < tcov_min:
            count("LOW_COVERAGE")
            continue
        acc = afdb_accession(row["target"])
        if acc is None:
            count("FOREIGN_TARGET")
            continue
        agi, status = resolve_uniprot(acc, mapping)
        rec = mapping.get(acc) or {"names": set()}
        names = {n for n in rec["names"] if n}
        if len({n.upper() for n in names}) > 1:
            count("CONFLICTING_UNIPROT_NAMES")
            continue
        if not names:
            count("UNNAMED")
            continue
        calls[query] = {"accession": acc, "agi": agi, "agi_status": status,
                        "name": sorted(names)[0], "bits": row["bits"],
                        "qcov": row["qcov"], "tcov": row["tcov"]}
    return {"calls": calls, "abstained": reasons, "n_queried": len(best)}


# -------------------------------------------------------- decision rules --

def decision_rules(m: dict) -> dict:
    """The pre-committed (a)-(d) rules over the measurement dict. Keys:
    safe_pct, named_rep_pct, conflict_rate, conflict_wilson, n_overlap,
    worst_species_conflict, swissprot_named_pct, arath_addressable_pct.
    Missing metrics leave a rule 'UNDECIDED' instead of silently passing."""
    def have(*keys):
        return all(m.get(k) is not None for k in keys)

    out: Dict[str, dict] = {}

    if have("safe_pct", "named_rep_pct", "conflict_rate", "conflict_wilson",
            "n_overlap", "worst_species_conflict"):
        ok = (m["safe_pct"] >= 50.0 and m["named_rep_pct"] >= 30.0
              and m["conflict_rate"] <= 0.02 and m["conflict_wilson"] <= 0.05
              and m["n_overlap"] >= 100
              and m["worst_species_conflict"] <= 0.10)
        out["a_full_panel_diamond"] = {"go": ok}
    else:
        out["a_full_panel_diamond"] = {"go": None, "reason": "UNDECIDED"}

    if have("swissprot_named_pct", "arath_addressable_pct"):
        out["b_arath_download"] = {
            "go": (10.0 <= m["swissprot_named_pct"] < 50.0
                   and m["arath_addressable_pct"] >= 15.0)}
        out["d_abandon_structural"] = {
            "go": (m["swissprot_named_pct"] < 10.0
                   and m["arath_addressable_pct"] < 15.0)
            or (m.get("conflict_rate") is not None
                and m["conflict_rate"] > 0.10)}
    else:
        out["b_arath_download"] = {"go": None, "reason": "UNDECIDED"}
        out["d_abandon_structural"] = {"go": None, "reason": "UNDECIDED"}

    out["c_broader_proteomes"] = {
        "go": None,
        "reason": "decided only after (b) runs: needs >=10pp gain from ARATH",
    }
    return out


# ------------------------------------------------------------------ CLI --

def _cmd_cohort(args) -> int:
    families: Dict[str, List[str]] = {}
    with open(args.summary) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi = header.index("gene_list")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            families[parts[0]] = [g for g in parts[gi].split(",") if g]

    pep = {}
    for fa in sorted(Path(args.pep_dir).glob("*.fa")):
        pep.update(read_fasta(str(fa)))
    lengths = {g: len(s.rstrip("*")) for g, s in pep.items()}

    reps = pick_representatives(families, lengths)

    unplaced = read_fasta(args.unplaced)
    cds = {}
    for fa in sorted(Path(args.cds_dir).glob("*.fa")):
        cds.update(read_fasta(str(fa)))
    rep_set = set(reps.values())
    eligible: Dict[str, List[str]] = {}
    for gene, seq in unplaced.items():
        if gene in rep_set:
            continue
        plen = len(seq.rstrip("*"))
        integrity = assess_cds(cds.get(gene))
        if is_comparable(plen, integrity.complete):
            eligible.setdefault(species_of(gene), []).append(gene)
    sampled = sample_comparable_unplaced(eligible, args.per_species)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for family_id in sorted(reps):
        g = reps[family_id]
        rows.append((family_id, g, "representative"))
    for g in sampled:
        rows.append(("-", g, "comparable_unplaced"))
    with open(outdir / "cohort.tsv", "w") as f:
        f.write("family_id\tgene_id\tspecies\tclass\taa_length\t"
                "n_family_members\tsha256\n")
        for family_id, g, klass in rows:
            seq = pep.get(g) or unplaced.get(g) or ""
            f.write(f"{family_id}\t{g}\t{species_of(g)}\t{klass}\t"
                    f"{len(seq.rstrip('*'))}\t"
                    f"{len(families.get(family_id, []))}\t"
                    f"{hashlib.sha256(seq.encode()).hexdigest()}\n")
    write_fasta({g: (pep.get(g) or unplaced[g]) for _, g, _ in rows},
                str(outdir / "cohort.fa"))
    per_species = {sp: len(v) for sp, v in sorted(eligible.items())}
    print(f"representatives: {len(reps)}  comparable-unplaced sampled: "
          f"{len(sampled)} (eligible per species: {per_species})")
    if len(reps) != len(families):
        print("FATAL: representative count != family count", file=sys.stderr)
        return 1
    return 0


def _cmd_qc_orthology(args) -> int:
    hits = parse_diamond_tsv(open(args.forward))
    reverse_best: Dict[Tuple[str, str], str] = {}
    for path in args.reverse:
        species = Path(path).stem.split(".")[0]
        for agi_query, rows in parse_diamond_tsv(open(path)).items():
            best = max(rows, key=lambda h: (h["bitscore"], -h["evalue"]))
            reverse_best[(species, strip_transcript(agi_query))] = best["qseqid"]

    gates = Gates(evalue=args.evalue, pident=args.pident,
                  qcov=args.qcov, scov=args.scov, margin=args.margin)
    cohort = [line.split("\t")[1] for line in
              open(args.cohort).read().splitlines()[1:] if line.strip()]
    verdicts = [safe_orthology(q, hits.get(q, []), reverse_best, gates)
                for q in cohort]
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "orthology_verdicts.tsv", "w") as f:
        f.write("query\tstatus\tagi\n")
        for v in verdicts:
            f.write(f"{v['query']}\t{v['status']}\t{v['agi'] or ''}\n")

    summary = {"n": len(verdicts)}
    for v in verdicts:
        summary[v["status"]] = summary.get(v["status"], 0) + 1
    with_hit = sum(n for s, n in summary.items()
                   if s not in ("n", "NO_HIT"))
    summary["safe_pct_of_hits"] = (100.0 * summary.get("SAFE", 0) / with_hit
                                   if with_hit else None)
    if args.direct and args.descriptions:
        direct = {p[0]: p[1] for p in
                  (l.split("\t") for l in open(args.direct))
                  if len(p) >= 2}
        descriptions = {strip_transcript(p[0]): p[1].rstrip("\n") for p in
                        (l.split("\t") for l in open(args.descriptions))
                        if len(p) >= 2}
        grid = grid_calibrate(hits, reverse_best, direct, descriptions)
        with open(outdir / "gate_calibration.tsv", "w") as f:
            cols = list(grid[0])
            f.write("\t".join(cols) + "\n")
            for row in grid:
                f.write("\t".join(str(row[c]) for c in cols) + "\n")
    (outdir / "orthology_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


def _cmd_qc_afdb(args) -> int:
    mapping = load_idmapping(open(args.idmapping))
    result = afdb_named_calls(open(args.m8), mapping,
                              qcov_min=args.qcov, tcov_min=args.tcov)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "afdb_calls.tsv", "w") as f:
        f.write("query\taccession\tagi\tagi_status\tname\tbits\tqcov\ttcov\n")
        for q in sorted(result["calls"]):
            c = result["calls"][q]
            f.write(f"{q}\t{c['accession']}\t{c['agi'] or ''}\t"
                    f"{c['agi_status']}\t{c['name']}\t{c['bits']}\t"
                    f"{c['qcov']}\t{c['tcov']}\n")
    summary = {"n_queried": result["n_queried"],
               "n_named": len(result["calls"]),
               "abstained": result["abstained"]}
    (outdir / "afdb_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


def _cmd_report(args) -> int:
    metrics = json.loads(Path(args.metrics).read_text())
    verdicts = decision_rules(metrics)
    print(json.dumps({"metrics": metrics, "decisions": verdicts}, indent=2))
    out = Path(args.outdir) / "pilot_decisions.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"metrics": metrics, "decisions": verdicts},
                              indent=2))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("cohort")
    c.add_argument("--summary", required=True)
    c.add_argument("--pep-dir", required=True)
    c.add_argument("--cds-dir", required=True)
    c.add_argument("--unplaced", required=True)
    c.add_argument("--per-species", type=int, default=PER_SPECIES_CAP)
    c.add_argument("--outdir", required=True)

    q = sub.add_parser("qc-orthology")
    q.add_argument("--cohort", required=True)
    q.add_argument("--forward", required=True)
    q.add_argument("--reverse", nargs="+", required=True,
                   help="per-species reverse DIAMOND TSVs (ATH as query)")
    q.add_argument("--direct", help="gene<TAB>description curated annotations")
    q.add_argument("--descriptions",
                   help="AGI<TAB>description table for calibration")
    q.add_argument("--evalue", type=float, default=Gates.evalue)
    q.add_argument("--pident", type=float, default=Gates.pident)
    q.add_argument("--qcov", type=float, default=Gates.qcov)
    q.add_argument("--scov", type=float, default=Gates.scov)
    q.add_argument("--margin", type=float, default=Gates.margin)
    q.add_argument("--outdir", required=True)

    a = sub.add_parser("qc-afdb")
    a.add_argument("--m8", required=True)
    a.add_argument("--idmapping", required=True,
                   help="decompressed ARATH_3702_idmapping.dat snapshot")
    a.add_argument("--qcov", type=float, default=0.70)
    a.add_argument("--tcov", type=float, default=0.70)
    a.add_argument("--outdir", required=True)

    r = sub.add_parser("report")
    r.add_argument("--metrics", required=True)
    r.add_argument("--outdir", required=True)

    args = ap.parse_args(argv)
    return {"cohort": _cmd_cohort, "qc-orthology": _cmd_qc_orthology,
            "qc-afdb": _cmd_qc_afdb, "report": _cmd_report}[args.cmd](args)


if __name__ == "__main__":
    raise SystemExit(main())
