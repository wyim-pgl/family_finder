"""Multi-symbol family label adjudication (MULTI_SYMBOL_LABEL_DESIGN.md).

Pure functions over plain dicts/dataclasses — no ete4, no external tools.
The CLI (file I/O, schema validation, manifests) lives in
adjudicate_labels.py; tree judgments come in pre-computed through the
``tree_results`` argument so this module never touches a tree library.

Three hard rules carried over from the design:

* Different stems are never merged by weighted vote. Source priority
  (curation > ATH > AFDB) only picks the display spelling and metadata of
  the SAME stem.
* A gate failure is an abstention, not a counter-vote; a parse failure in
  an accepted call forces review because its stem is unknown.
* Every non-emitted label leaves a positive needs_review/audit record —
  absence of a row is never the encoding of a decision.

Recorded deviations from the design doc (each mirrored in the tests):

* §9 vs §12.4: the family audit row of a MULTI_STEM_PARTITIONED family has
  transfer_tier=needs_review, so its per-row verdict is
  NEEDS_REVIEW_MULTI_STEM while the family-level verdict stays
  MULTI_STEM_PARTITIONED.
* §9: auto-confirmed subfamily rows transfer a STEM through the tree gate,
  so they carry transfer_tier=stem_auto (the suffix tier would violate the
  §12.4 non-empty-suffix invariant).
* §10's ``failures_affect_choice`` is deliberately conservative here: any
  parse failure in a gate-passed call forces review, because the failed
  symbol's stem is unknown and could always conflict.
* §3.3 collision rule: two primaries sharing a symbol key may parse to
  different stems only when EVERY participant went through an explicit
  override; any grammar-parsed participant makes the divergence a
  NormalizationCollision.
"""

import re
import unicodedata
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

KNOWN_SOURCES = ("mcry_curated", "ath_diamond", "afdb_swissprot")
SOURCE_PRIORITY = KNOWN_SOURCES
# Calibration status is an explicit allow-list, never an inequality against
# one literal (a typo'd source must fail loudly, not fail open).
SOURCE_CALIBRATED = {"mcry_curated": True, "ath_diamond": True, "afdb_swissprot": False}

# Diagnostic weights only — never used to pick between different stems.
DEFAULT_SOURCE_WEIGHTS = {"mcry_curated": 1.0, "ath_diamond": 0.5, "afdb_swissprot": 0.25}

AUTO_VERDICTS = frozenset({
    "SINGLE_STEM",
    "SINGLE_STEM_MULTI_SUFFIX",
    "ORTHOGRAPHIC_VARIANT",
    "MULTI_STEM_PARTITIONED",
})

_GREEK = {"α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta", "ε": "epsilon"}
_DASHES = "‐‑‒–—―−"
_ROMAN_RE = re.compile(r"^[IVXivx]+$")
_ALLELE_RE = re.compile(r";(\d+[A-Za-z]?)$")
_NUMERIC_RE = re.compile(r"(?<=[A-Za-z])(\d+(?:\.\d+)?[A-Za-z]?)$")
_BOUNDARY_SERIES_RE = re.compile(r"[-_]([A-Za-z]\d+)$")
_BOUNDARY_ROMAN_RE = re.compile(r"[-_ ]([IVXivx]+)$")
_SUFFIX_TOKEN_RE = re.compile(r"^\d+[A-Za-z]?(?:[.;/]\d+)?$")


class NormalizationCollision(Exception):
    """Two raw symbols share a key without an explaining rule or alias."""


class OverrideConflict(Exception):
    """Duplicate override keys with different stem/suffix rows."""


@dataclass(frozen=True)
class Norm:
    raw: str
    display: str
    primary: str
    aliases: Tuple[str, ...]
    key: str
    rules: Tuple[str, ...]
    status: str  # OK | AMBIGUOUS


@dataclass(frozen=True)
class Parse:
    stem: str
    suffix: str
    status: str  # PARSED | OVERRIDE | AMBIGUOUS | UNPARSED
    rule_id: str


@dataclass(frozen=True)
class Override:
    stem: str
    suffix: str
    scope_source: str
    reason: str


@dataclass(frozen=True)
class Policy:
    """Explicit thresholds — no invented defaults for preregistered gates."""

    min_label_coverage: float
    source_weights: Dict[str, float] = field(default_factory=lambda: dict(DEFAULT_SOURCE_WEIGHTS))


# ------------------------------------------------------------------- normalization


def symbol_key(text: str) -> str:
    """Comparison key: NFKC, dash/space tidy, greek named, casefolded."""
    s = unicodedata.normalize("NFKC", text)
    for d in _DASHES:
        s = s.replace(d, "-")
    s = " ".join(s.split())
    for g, name in _GREEK.items():
        s = s.replace(g, name)
        s = s.replace(g.upper(), name)
    return s.casefold()


def normalize_symbol(raw: str) -> Norm:
    """Canonicalize one raw symbol; never overwrite the original."""
    rules: List[str] = []
    s = unicodedata.normalize("NFKC", raw)
    if s != raw:
        rules.append("nfkc")
    for d in _DASHES:
        if d in s:
            s = s.replace(d, "-")
            rules.append("dash")
    collapsed = " ".join(s.split())
    if collapsed != s:
        rules.append("whitespace")
    s = collapsed

    if s.count("(") != s.count(")") or "((" in s.replace(" ", ""):
        return Norm(raw, s, s, (), symbol_key(s), tuple(rules), "AMBIGUOUS")

    aliases: List[str] = []
    def _take_alias(match):
        inner = match.group(1).strip()
        if inner:
            aliases.append(inner)
        rules.append("paren_alias")
        return ""

    primary = re.sub(r"\(([^()]*)\)", _take_alias, s).strip()
    if "(" in primary or ")" in primary:
        return Norm(raw, s, primary, tuple(aliases), symbol_key(primary), tuple(rules), "AMBIGUOUS")

    tokens = primary.split(" ")
    if len(tokens) == 2 and _SUFFIX_TOKEN_RE.match(tokens[1]):
        primary = tokens[0] + tokens[1]
        rules.append("suffix_glue")
    elif len(tokens) == 2 and _ROMAN_RE.match(tokens[1]):
        # §3.2: space IS the boundary the roman rule needs — keep it for the
        # boundary grammar instead of gluing it away.
        rules.append("roman_space_boundary")
    elif len(tokens) > 1:
        return Norm(raw, primary, primary, tuple(aliases), symbol_key(primary),
                    tuple(rules), "AMBIGUOUS")

    return Norm(raw, primary, primary, tuple(aliases), symbol_key(primary), tuple(rules), "OK")


# ----------------------------------------------------------------------- overrides


def load_overrides(lines: Iterable[str]) -> Dict[Tuple[str, str], Override]:
    """Versioned override TSV -> {(symbol_key, scope): Override}.

    Conflicting duplicates raise; a key that is not already in normalized
    form raises too (a raw-cased key would otherwise silently never match).
    """
    out: Dict[Tuple[str, str], Override] = {}
    header: Optional[List[str]] = None
    for line in lines:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if header is None:
            header = parts
            for col in ("symbol_key", "stem", "suffix", "scope_source", "reason"):
                if col not in header:
                    raise ValueError("override table missing column: %s" % col)
            continue
        row = dict(zip(header, parts))
        key = row["symbol_key"].strip()
        if key != symbol_key(key):
            raise ValueError("override key %r is not in normalized form (want %r)"
                             % (key, symbol_key(key)))
        scope = row.get("scope_source", "*") or "*"
        ov = Override(row["stem"], row.get("suffix", ""), scope, row.get("reason", ""))
        if (key, scope) in out and out[(key, scope)] != ov:
            raise OverrideConflict("conflicting override rows for key %r scope %r"
                                   % (key, scope))
        out[(key, scope)] = ov
    return out


def stale_overrides(overrides: Dict[Tuple[str, str], Override],
                    used_keys: Set[Tuple[str, str]]) -> List[Tuple[str, str]]:
    return sorted(k for k in overrides if k not in used_keys)


def load_aliases(lines: Iterable[str]) -> Dict[str, str]:
    """Versioned alias TSV -> {stem_key: class representative}.

    Rows declare evidence-backed stem equivalences (columns stem_key_a,
    stem_key_b, ..., source). Pairs are merged transitively (union-find);
    the class representative is the lexicographically smallest member, so
    the mapping is deterministic regardless of row order. Keys must be in
    normalized (symbol_key) form — raw-cased keys raise, they would
    otherwise never match anything.
    """
    parent: Dict[str, str] = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    header: Optional[List[str]] = None
    for line in lines:
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if header is None:
            header = parts
            for col in ("stem_key_a", "stem_key_b", "source"):
                if col not in header:
                    raise ValueError("alias table missing column: %s" % col)
            continue
        row = dict(zip(header, parts))
        a, b = row["stem_key_a"].strip(), row["stem_key_b"].strip()
        for k in (a, b):
            if k != symbol_key(k):
                raise ValueError("alias key %r is not in normalized form (want %r)"
                                 % (k, symbol_key(k)))
            parent.setdefault(k, k)
        ra, rb = find(a), find(b)
        if ra != rb:
            lo, hi = sorted((ra, rb))
            parent[hi] = lo

    return {k: find(k) for k in parent}


def _stem_class(stem_key: str, aliases: Optional[Dict[str, str]]) -> str:
    if aliases is None:
        return stem_key
    return aliases.get(stem_key, stem_key)


def _override_for(overrides, key: str, source: Optional[str]):
    """Source-scoped override wins over the wildcard scope."""
    if source is not None and (key, source) in overrides:
        return overrides[(key, source)], (key, source)
    if (key, "*") in overrides:
        return overrides[(key, "*")], (key, "*")
    return None, None


# -------------------------------------------------------------------------- parser


def _grammar_candidates(primary: str) -> List[Tuple[str, str, str]]:
    """All (stem, suffix, rule) decompositions the grammars allow."""
    allele = ""
    base = primary
    m = _ALLELE_RE.search(base)
    if m:
        allele = ";" + m.group(1)
        base = base[: m.start()]

    candidates: List[Tuple[str, str, str]] = []
    m = _NUMERIC_RE.search(base)
    if m and base[: m.start()]:
        candidates.append((base[: m.start()], m.group(1) + allele, "numeric_series"))
    m = _BOUNDARY_SERIES_RE.search(base)
    if m and base[: m.start()]:
        candidates.append((base[: m.start()], m.group(1) + allele, "boundary_series"))
    m = _BOUNDARY_ROMAN_RE.search(base)
    if m and base[: m.start()]:
        candidates.append((base[: m.start()], m.group(1).upper() + allele, "boundary_roman"))
    if not candidates and allele and base.isalpha():
        candidates.append((base, allele, "allele_only"))
    return candidates


def parse_symbol(norm: Norm, overrides, source: Optional[str] = None) -> Parse:
    """(stem, suffix, status, rule) for one normalized symbol.

    §4.1.5: exactly one grammar parse -> PARSED; several -> AMBIGUOUS.
    Zero parses are PARSED only for a purely alphabetic symbol (the bare
    stem, e.g. ``Acl``); anything with digits or separators that no grammar
    claimed is UNPARSED — an abstention, never an asserted stem.
    """
    if norm.status != "OK":
        return Parse("", "", "AMBIGUOUS", "normalization")
    primary = norm.primary
    if not re.search(r"[A-Za-z]", primary):
        return Parse("", "", "UNPARSED", "no_letters")

    ov, _ = _override_for(overrides, symbol_key(primary), source)
    if ov is not None:
        return Parse(ov.stem, ov.suffix, "OVERRIDE", "override")

    candidates = _grammar_candidates(primary)
    distinct = sorted({(s, x) for s, x, _ in candidates})
    if len(distinct) > 1:
        return Parse("", "", "AMBIGUOUS", ";".join(sorted(r for _, _, r in candidates)))
    if len(distinct) == 1:
        stem, suffix = distinct[0]
        rule = candidates[0][2]
        return Parse(stem, suffix, "PARSED", rule)
    if primary.isalpha():
        return Parse(primary, "", "PARSED", "bare")
    return Parse("", "", "UNPARSED", "no_grammar_match")


# ---------------------------------------------------------------------- axis gates


def different_stem_margin(hits: Sequence[Tuple[float, str]]):
    """(margin, status) where the competitor must carry a DIFFERENT stem key.

    §7: margin is the RATIO (best - rival) / best. Ties on the best bit
    score across different stems abstain regardless of input order; no
    rival stem is NA, never infinity; a non-positive best abstains.
    """
    if not hits:
        return None, "NA_no_hits"
    ordered = sorted(hits, key=lambda h: (-h[0], h[1]))
    best_bits, best_stem = ordered[0]
    rival = next(((b, s) for b, s in ordered if s != best_stem), None)
    if rival is None:
        return None, "NA_single_stem_candidate"
    if rival[0] == best_bits:
        return None, "ABSTAIN_tied_best"
    if best_bits <= 0:
        return None, "ABSTAIN_nonpositive_best"
    return (best_bits - rival[0]) / best_bits, "OK"


# -------------------------------------------------------------------- adjudication


def _display_for(class_key: str, units: List[dict],
                 aliases: Optional[Dict[str, str]] = None) -> Tuple[str, str]:
    """(symbol_display, stem_display) chosen by source priority, then lexicographic.

    ``class_key`` is an alias equivalence class; every unit whose stem
    belongs to the class competes, and source priority picks the spelling.
    """
    pool = [u for u in units if _stem_class(u["stem_key"], aliases) == class_key]
    pool.sort(key=lambda u: (SOURCE_PRIORITY.index(u["source"]), u["symbol_display"]))
    return pool[0]["symbol_display"], pool[0]["stem_display"]


def _evidence_units(calls: List[dict], overrides):
    """Normalize, parse and deduplicate accepted calls into evidence units.

    One vote per (source, gene, stem); a second distinct full symbol from
    the same gene keeps no extra vote but stays visible as a variant so its
    suffix still reaches the tree gate (§8) and the sidecar records it.
    """
    accepted: Dict[Tuple[str, str, str], dict] = {}
    variants: List[dict] = []
    failures: List[dict] = []
    abstentions: List[dict] = []
    used_override_keys: Set[Tuple[str, str]] = set()
    key_parse: Dict[str, Tuple[str, str, bool]] = {}  # full key -> (stem, suffix, via_override)

    for c in sorted(calls, key=lambda c: (c["source"], c["gene"], c["symbol"])):
        if c["source"] not in KNOWN_SOURCES:
            raise ValueError("unknown call source %r (known: %s)"
                             % (c["source"], ", ".join(KNOWN_SOURCES)))
        norm = normalize_symbol(c["symbol"])
        if not c.get("gate_passed", True):
            abstentions.append({"call": c, "reason": c.get("gate_reason", "gate_failed")})
            continue
        parsed = parse_symbol(norm, overrides, source=c["source"])
        if parsed.status == "OVERRIDE":
            _, okey = _override_for(overrides, symbol_key(norm.primary), c["source"])
            used_override_keys.add(okey)
        if parsed.status in ("AMBIGUOUS", "UNPARSED"):
            failures.append({"call": c, "status": parsed.status, "rule": parsed.rule_id})
            continue

        via_override = parsed.status == "OVERRIDE"
        this_parse = (symbol_key(parsed.stem), parsed.suffix)
        seen = key_parse.get(norm.key)
        if seen is not None and (seen[0], seen[1]) != this_parse:
            # §3.3: same key, different parse (case-insensitive on the stem —
            # display spelling is not a collision). Divergence is allowed only
            # when every participant is explicitly overridden (source-scoped
            # aliases); any grammar-parsed participant makes it a collision.
            if not (via_override and seen[2]):
                raise NormalizationCollision(
                    "symbol key %r parses to both %r and %r"
                    % (norm.key, (seen[0], seen[1]), this_parse))
        key_parse.setdefault(norm.key, (this_parse[0], this_parse[1], via_override))

        unit = {
            "source": c["source"],
            "gene": c["gene"],
            "stem_key": symbol_key(parsed.stem),
            "stem_display": parsed.stem,
            "suffix": parsed.suffix,
            "symbol_display": norm.display,
            "symbol_full_key": norm.key,
            "calibrated": bool(c.get("calibrated", SOURCE_CALIBRATED[c["source"]])),
            "aliases": norm.aliases,
        }
        unit_key = (c["source"], c["gene"], unit["stem_key"])
        if unit_key in accepted:
            if accepted[unit_key]["symbol_full_key"] != norm.key:
                variants.append(unit)
        else:
            accepted[unit_key] = unit
    return list(accepted.values()), variants, failures, abstentions, used_override_keys


def _row(level, target_id, symbol, stem, suffix, tier, verdict, reason_codes,
         sources, support, coverage, evidence_id):
    return {
        "level": level, "target_id": target_id, "symbol": symbol, "stem": stem,
        "suffix": suffix, "transfer_tier": tier, "verdict": verdict,
        "reason_code": ";".join(sorted({r for r in reason_codes if r})),
        "sources": ";".join(sorted(set(sources))),
        "support": round(support, 4), "coverage": round(coverage, 4),
        "evidence_id": evidence_id,
    }


def _result(family_id, verdict, rows, *, audit, failures, abstentions, variants,
            withheld, uncovered, support, coverage, coverage_all, used_override_keys):
    return {
        "family_id": family_id, "verdict": verdict, "rows": rows,
        "audit": sorted(audit), "failures": failures, "abstentions": abstentions,
        "variants": variants, "withheld": withheld,
        "uncovered_genes": sorted(uncovered),
        "support": support, "coverage": coverage, "coverage_all_sources": coverage_all,
        "used_override_keys": sorted(used_override_keys),
    }


def _review_rows(family_id, verdict, reason_codes, units, support, coverage):
    return [_row("family", family_id, "", "", "", "needs_review", verdict,
                 reason_codes, [u["source"] for u in units] or ["none"],
                 support, coverage, "%s:review" % family_id)]


def _partition_by_catalog(units, catalog, aliases=None):
    """class_key -> set(og_id) for evidence-bearing genes; uncovered genes kept."""
    gene_og = {}
    for og_id, info in catalog.items():
        for g in info["members"]:
            gene_og[g] = og_id
    part: Dict[str, Set[str]] = {}
    uncovered: Set[str] = set()
    for u in units:
        og = gene_og.get(u["gene"])
        if og is None:
            uncovered.add(u["gene"])
        else:
            part.setdefault(_stem_class(u["stem_key"], aliases), set()).add(og)
    return part, uncovered


def adjudicate_family(family_id: str, members: Set[str], calls: List[dict], *,
                      overrides, policy: Policy,
                      catalog: Optional[Dict[str, dict]] = None,
                      tree_results: Optional[Dict[Tuple[str, str], dict]] = None,
                      aliases: Optional[Dict[str, str]] = None) -> dict:
    """State machine of design §6 for one family. Membership is never altered."""
    for c in calls:
        if c["gene"] not in members:
            raise ValueError("evidence gene %r is not a member of %s" % (c["gene"], family_id))

    units, variants, failures, abstentions, used_keys = _evidence_units(calls, overrides)
    calibrated_units = [u for u in units if u["calibrated"]]
    # The coverage that unlocks a label counts calibrated evidence only —
    # an uncalibrated axis may corroborate but never unlock (§5).
    coverage = (len({u["gene"] for u in calibrated_units}) / len(members)
                if members else 0.0)
    coverage_all = (len({u["gene"] for u in units}) / len(members)
                    if members else 0.0)
    weights = policy.source_weights
    audit: List[str] = []
    for u in units:
        u["class_key"] = _stem_class(u["stem_key"], aliases)
    for u in units:
        if not u["calibrated"]:
            others = {v["class_key"] for v in calibrated_units}
            if others and u["class_key"] not in others:
                audit.append("AFDB_CONFLICT:%s:%s" % (u["gene"], u["stem_display"]))
    for v in variants:
        audit.append("VARIANT:%s:%s:%s" % (v["source"], v["gene"], v["symbol_display"]))

    stems = sorted({u["class_key"] for u in calibrated_units})
    total_w = sum(weights.get(u["source"], 0.0) for u in calibrated_units)
    if stems and total_w > 0:
        lead_w = max(sum(weights.get(u["source"], 0.0)
                         for u in calibrated_units if u["class_key"] == s) for s in stems)
        support = lead_w / total_w
    else:
        support = 0.0

    kw = dict(audit=audit, failures=failures, abstentions=abstentions,
              variants=variants, withheld=[], uncovered=set(),
              support=support, coverage=coverage, coverage_all=coverage_all,
              used_override_keys=used_keys)

    if failures:
        verdict = "NEEDS_REVIEW_PARSE"
        rows = _review_rows(family_id, verdict,
                            ["PARSE_%s" % f["status"] for f in failures],
                            units, support, coverage)
        return _result(family_id, verdict, rows, **kw)
    if not units:
        return _result(family_id, "NO_EVIDENCE", [], **kw)
    if not calibrated_units:
        verdict = "NEEDS_REVIEW_UNCALIBRATED_SOURCE"
        rows = _review_rows(family_id, verdict, ["UNCALIBRATED_SOURCE"],
                            units, support, coverage)
        return _result(family_id, verdict, rows, **kw)

    if len(stems) == 1:
        return _single_stem_result(family_id, stems[0], calibrated_units, variants,
                                   units, support, coverage, policy, kw, aliases)
    return _multi_stem_result(family_id, members, stems, calibrated_units, units,
                              support, coverage, policy, catalog, tree_results,
                              kw, aliases)


def _single_stem_result(family_id, class_key, calibrated_units, variants, units,
                        support, coverage, policy, kw, aliases=None):
    if coverage < policy.min_label_coverage:
        verdict = "NEEDS_REVIEW_LOW_COVERAGE"
        rows = _review_rows(family_id, verdict, ["LOW_COVERAGE"],
                            units, support, coverage)
        return _result(family_id, verdict, rows, **kw)
    _, stem_display = _display_for(class_key, calibrated_units, aliases)
    suffix_pool = calibrated_units + [v for v in variants if v["calibrated"]]
    full_keys = {u["symbol_full_key"] for u in calibrated_units}
    displays = {u["symbol_display"] for u in calibrated_units}
    suffixes = sorted({u["suffix"] for u in suffix_pool if u["suffix"]})
    sources = [u["source"] for u in calibrated_units]

    if len(full_keys) == 1 and len(displays) > 1:
        verdict = "ORTHOGRAPHIC_VARIANT"
    elif suffixes and (len(full_keys) > 1 or len(suffixes) > 1):
        verdict = "SINGLE_STEM_MULTI_SUFFIX"
    else:
        verdict = "SINGLE_STEM"

    # The stem row asserts a STEM transfer, so it shows the stem, not a
    # member symbol whose suffix could be misread as transferred (§1: the
    # 5.6% suffix-conflict figure is why).
    rows = [_row("family", family_id, stem_display, stem_display, "", "stem_auto",
                 verdict, [], sources, support, coverage, "%s:stem" % family_id)]
    for i, suffix in enumerate(suffixes):
        sym = next(u["symbol_display"] for u in sorted(
            suffix_pool, key=lambda u: (SOURCE_PRIORITY.index(u["source"]),
                                        u["symbol_display"]))
            if u["suffix"] == suffix)
        rows.append(_row("family", family_id, sym, stem_display, suffix,
                         "suffix_needs_tree_gate", verdict, [], sources,
                         support, coverage, "%s:suffix:%d" % (family_id, i)))
    return _result(family_id, verdict, rows, **kw)


def _multi_stem_result(family_id, members, stems, calibrated_units, units,
                       support, coverage, policy, catalog, tree_results, kw,
                       aliases=None):
    def review(reasons, withheld=(), uncovered=()):
        kw2 = dict(kw)
        kw2["withheld"] = list(withheld)
        kw2["uncovered"] = set(uncovered)
        verdict = "NEEDS_REVIEW_MULTI_STEM"
        rows = _review_rows(family_id, verdict, reasons, units, support, coverage)
        return _result(family_id, verdict, rows, **kw2)

    if not catalog or tree_results is None:
        return review(["STEM_CONFLICT"])

    part, uncovered = _partition_by_catalog(calibrated_units, catalog, aliases)
    og_stems: Dict[str, Set[str]] = {}
    for stem_key, ogs in part.items():
        for og in ogs:
            og_stems.setdefault(og, set()).add(stem_key)
    if any(len(s) > 1 for s in og_stems.values()):
        return review(["STEM_CONFLICT"], uncovered=uncovered)

    passed_rows, withheld = [], []
    for class_key in stems:
        if class_key not in part:
            withheld.append({"stem": class_key, "og": "", "reason": "CATALOG_UNCOVERED"})
            continue
        _, stem_display = _display_for(class_key, calibrated_units, aliases)
        sources = [u["source"] for u in calibrated_units
                   if u["class_key"] == class_key]
        for og in sorted(part[class_key]):
            grade = catalog[og].get("grade")
            if grade != "HIGH":
                withheld.append({"stem": stem_display, "og": og, "reason": "CATALOG_NOT_HIGH"})
                continue
            if not catalog[og]["members"] <= members:
                # §6.1.1 / §12.4: a label must never leak past the family.
                withheld.append({"stem": stem_display, "og": og,
                                 "reason": "CATALOG_MEMBER_MISMATCH"})
                continue
            tres = tree_results.get((symbol_key(stem_display), og))
            if tres is None:
                tres = tree_results.get((stem_display, og))
            if tres is None:
                withheld.append({"stem": stem_display, "og": og, "reason": "TREE_MISSING"})
                continue
            if not tres.get("transferable"):
                withheld.append({"stem": stem_display, "og": og,
                                 "reason": "TREE_NOT_TRANSFERABLE"})
                continue
            if tres.get("grade") != "HIGH":
                withheld.append({"stem": stem_display, "og": og, "reason": "CATALOG_NOT_HIGH"})
                continue
            if tres.get("coverage", 0.0) < policy.min_label_coverage:
                withheld.append({"stem": stem_display, "og": og, "reason": "LOW_COVERAGE"})
                continue
            passed_rows.append(_row(
                "subfamily", og, stem_display, stem_display, "",
                "stem_auto", "MULTI_STEM_PARTITIONED", [], sources,
                support, tres.get("coverage", 0.0), "%s:%s" % (family_id, og)))

    if not passed_rows:
        return review(["STEM_CONFLICT"] + [w["reason"] for w in withheld],
                      withheld=withheld, uncovered=uncovered)

    family_reasons = sorted({w["reason"] for w in withheld})
    high_covered = set()
    for og, info in catalog.items():
        if info.get("grade") == "HIGH" and info["members"] <= members:
            high_covered |= info["members"]
    if members - high_covered or uncovered:
        family_reasons.append("PARTIAL_HIGH_COVERAGE")

    family_row = _row("family", family_id, "", "", "", "needs_review",
                      "NEEDS_REVIEW_MULTI_STEM",
                      family_reasons or ["MULTI_STEM_FAMILY_LEVEL_UNNAMED"],
                      [u["source"] for u in calibrated_units], support, coverage,
                      "%s:family_audit" % family_id)
    kw2 = dict(kw)
    kw2["withheld"] = withheld
    kw2["uncovered"] = uncovered
    return _result(family_id, "MULTI_STEM_PARTITIONED", [family_row] + passed_rows, **kw2)
