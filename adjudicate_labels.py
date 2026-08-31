#!/usr/bin/env python3
"""Thin CLI for the multi-symbol label adjudication (design §11).

All judgment logic lives in steps/label_symbols.py; this script only does
file I/O, schema validation and the run manifest. Axis-specific quality
gates (coverage/margin thresholds of the ATH DIAMOND or AFDB foldseek
searches) are applied where those axes are computed — each call arrives
here already carrying its ``gate_passed`` verdict. Boolean columns accept
an explicit vocabulary only; anything else is an error, never a silent
pass (fail-open booleans were review finding C3).

Inputs
------
--families    TSV with header containing family_id and gene_list columns
              (arbiter summary format; gene_list comma-separated).
--calls       SOURCE=PATH, repeatable. SOURCE must be one of the known
              axes (mcry_curated | ath_diamond | afdb_swissprot). TSV with
              header: family_id, gene, symbol, gate_passed, gate_reason,
              calibrated (last three optional; defaults: True, "", the
              axis's registered calibration status).
--overrides   Versioned symbol_stem_overrides TSV (keys must be normalized).
--catalog     Optional TSV: og_id, family_id, grade, gene_list.
--tree-results Optional TSV: stem, og_id, transferable, grade, coverage
              (stem stored under its normalized key).
--min-label-coverage  Required, preregistered — no built-in default.

Outputs (into --out DIR)
------------------------
labels.tsv      v1 12 columns + verdict reason_code sources support coverage
                evidence_id suffix (design §9 + an explicit suffix column so
                the §12.4 non-empty-suffix invariant is checkable on the
                artifact itself; metadata columns left blank here — the
                metadata join is a separate enrichment step).
evidence.tsv    audit sidecar: one verdict row per family (NO_EVIDENCE
                included — absence of a row never encodes a decision),
                plus failures/abstentions/variants/withheld/uncovered/
                stale-override records.
verdicts.tsv    verdict -> family count, reconcilable against labels.tsv.
manifest.json   input checksums + policy. An existing, different manifest
                in --out refuses to run (no silent config drift on resume).
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from steps.label_symbols import (  # noqa: E402
    KNOWN_SOURCES, Policy, adjudicate_family, load_aliases, load_overrides,
    stale_overrides)

PARSER_VERSION = "label_symbols_v2"

V1_COLUMNS = ["level", "target_id", "mcry_gene", "symbol", "stem", "desc", "ec",
              "agi", "uniprot", "module", "og_grade", "transfer_tier"]
APPENDED = ["verdict", "reason_code", "sources", "support", "coverage",
            "evidence_id", "suffix"]

_TRUE = ("true", "1", "yes")
_FALSE = ("false", "0", "no")


def parse_bool(raw, default, context):
    v = raw.strip().lower()
    if not v:
        return default
    if v in _TRUE:
        return True
    if v in _FALSE:
        return False
    raise SystemExit("%s: boolean must be one of %s, got %r"
                     % (context, "/".join(_TRUE + _FALSE), raw))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path, required):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        missing = [c for c in required if c not in header]
        if missing:
            raise SystemExit("%s: missing required columns %s" % (path, missing))
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))
    return rows


def load_families(path):
    fams = {}
    for row in read_tsv(path, ["family_id", "gene_list"]):
        fams[row["family_id"]] = {g for g in row["gene_list"].split(",") if g}
    if not fams:
        raise SystemExit("%s: no families" % path)
    return fams


def load_calls(specs):
    from steps.label_symbols import SOURCE_CALIBRATED
    by_family = {}
    for spec in specs:
        if "=" not in spec:
            raise SystemExit("--calls expects SOURCE=PATH, got %r" % spec)
        source, path = spec.split("=", 1)
        if source not in KNOWN_SOURCES:
            raise SystemExit("unknown call source %r (known: %s)"
                             % (source, ", ".join(KNOWN_SOURCES)))
        for row in read_tsv(path, ["family_id", "gene", "symbol"]):
            ctx = "%s: family %s gene %s" % (path, row["family_id"], row["gene"])
            by_family.setdefault(row["family_id"], []).append({
                "source": source,
                "gene": row["gene"],
                "symbol": row["symbol"],
                "gate_passed": parse_bool(row.get("gate_passed", ""), True,
                                          ctx + " gate_passed"),
                "gate_reason": row.get("gate_reason", ""),
                "calibrated": parse_bool(row.get("calibrated", ""),
                                         SOURCE_CALIBRATED[source],
                                         ctx + " calibrated"),
            })
    return by_family


def load_catalog(path):
    cat = {}
    for row in read_tsv(path, ["og_id", "family_id", "grade", "gene_list"]):
        cat.setdefault(row["family_id"], {})[row["og_id"]] = {
            "grade": row["grade"],
            "members": {g for g in row["gene_list"].split(",") if g},
        }
    return cat


def load_tree_results(path):
    from steps.label_symbols import symbol_key
    out = {}
    for row in read_tsv(path, ["stem", "og_id", "transferable", "grade", "coverage"]):
        ctx = "%s: stem %s og %s" % (path, row["stem"], row["og_id"])
        out[(symbol_key(row["stem"]), row["og_id"])] = {
            "transferable": parse_bool(row["transferable"], False,
                                       ctx + " transferable"),
            "grade": row["grade"],
            "coverage": float(row["coverage"] or 0.0),
        }
    return out


def write_manifest(out_dir, args, inputs):
    manifest = {
        "parser_version": PARSER_VERSION,
        "min_label_coverage": args.min_label_coverage,
        "inputs": {str(p): sha256(p) for p in sorted(inputs)},
    }
    blob = json.dumps(manifest, indent=1, sort_keys=True)
    target = out_dir / "manifest.json"
    if target.exists() and target.read_text() != blob:
        raise SystemExit("refusing to overwrite %s: inputs or policy changed "
                         "(delete the output directory to rerun)" % target)
    target.write_text(blob)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--families", required=True, type=Path)
    ap.add_argument("--calls", action="append", required=True, metavar="SOURCE=PATH")
    ap.add_argument("--overrides", required=True, type=Path)
    ap.add_argument("--aliases", type=Path,
                    help="evidence-backed stem equivalence TSV "
                         "(data/symbol_alias_v1.tsv)")
    ap.add_argument("--catalog", type=Path)
    ap.add_argument("--tree-results", type=Path)
    ap.add_argument("--min-label-coverage", required=True, type=float)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    inputs = [args.families, args.overrides]
    if args.aliases:
        inputs.append(args.aliases)
    inputs += [Path(s.split("=", 1)[1]) for s in args.calls]
    if args.catalog:
        inputs.append(args.catalog)
    if args.tree_results:
        inputs.append(args.tree_results)
    write_manifest(args.out, args, inputs)

    families = load_families(args.families)
    calls = load_calls(args.calls)
    overrides = load_overrides(args.overrides.read_text().splitlines())
    aliases = (load_aliases(args.aliases.read_text().splitlines())
               if args.aliases else None)
    catalog = load_catalog(args.catalog) if args.catalog else {}
    trees = load_tree_results(args.tree_results) if args.tree_results else None
    policy = Policy(min_label_coverage=args.min_label_coverage)

    unmatched = sorted(set(calls) - set(families))
    if unmatched:
        raise SystemExit("calls reference %d unknown families, e.g. %s"
                         % (len(unmatched), unmatched[:5]))

    verdict_counts = {}
    all_used_override_keys = set()
    with open(args.out / "labels.tsv", "w") as lab, \
         open(args.out / "evidence.tsv", "w") as ev:
        lab.write("\t".join(V1_COLUMNS + APPENDED) + "\n")
        ev.write("family_id\tkind\tdetail\n")
        for family_id in sorted(families):
            fam_calls = calls.get(family_id, [])
            res = adjudicate_family(
                family_id, families[family_id], fam_calls,
                overrides=overrides, policy=policy, aliases=aliases,
                catalog=catalog.get(family_id) or None, tree_results=trees)
            verdict_counts[res["verdict"]] = verdict_counts.get(res["verdict"], 0) + 1
            all_used_override_keys.update(tuple(k) for k in res["used_override_keys"])
            for row in res["rows"]:
                out = {c: "" for c in V1_COLUMNS + APPENDED}
                out.update({k: row[k] for k in
                            ("level", "target_id", "symbol", "stem", "suffix",
                             "transfer_tier", "verdict", "reason_code", "sources",
                             "evidence_id")})
                out["support"] = "%.4f" % row["support"]
                out["coverage"] = "%.4f" % row["coverage"]
                lab.write("\t".join(str(out[c]) for c in V1_COLUMNS + APPENDED) + "\n")
            # Absence of a row never encodes a decision: every family gets a
            # verdict record here, NO_EVIDENCE included (review finding H1).
            ev.write("%s\tverdict\t%s\n" % (family_id, json.dumps(
                {"verdict": res["verdict"], "support": res["support"],
                 "coverage": res["coverage"],
                 "coverage_all_sources": res["coverage_all_sources"]},
                sort_keys=True)))
            for kind, items in (("audit", res["audit"]),
                                ("failure", res["failures"]),
                                ("abstention", res["abstentions"]),
                                ("variant", res["variants"]),
                                ("withheld", res["withheld"]),
                                ("catalog_uncovered", res["uncovered_genes"])):
                for item in items:
                    ev.write("%s\t%s\t%s\n" % (family_id, kind,
                                               json.dumps(item, sort_keys=True, default=str)))
        for key in stale_overrides(overrides, all_used_override_keys):
            ev.write("-\tstale_override\t%s\n" % json.dumps(list(key)))

    with open(args.out / "verdicts.tsv", "w") as vc:
        vc.write("verdict\tn_families\n")
        for v in sorted(verdict_counts):
            vc.write("%s\t%d\n" % (v, verdict_counts[v]))

    print("families: %d" % len(families))
    for v in sorted(verdict_counts):
        print("  %-36s %d" % (v, verdict_counts[v]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
