"""Sub- vs neo-functionalization verdict + narrative for a focal subfamily.

The report must not stop at diagnostics (SDP columns, monophyly, structure)
— it must EXPLAIN how the subfamily diverged. This module integrates the
selection-regime evidence with the functional-partition evidence and writes
that explanation.

Evidence axes (all optional; the verdict degrades honestly):
  RELAX      relaxation/intensification of selection on the focal clade
             (K < 1 with small p = relaxed constraint)
  MEME       episodic positive selection per site — restricted to the
             signal region when one is known (a subfunctionalized signal
             region should show NONE)
  aBSREL     per-branch positive selection, split into stem/internal vs
             terminal branches (stem-positive = neofunctionalization
             signature; terminal-only = lineage fine-tuning)
  expression the focal subfamily's share of family expression (dominance
             = partitioned ancestral role)
  signals    a subfamily-specific localization-signal region (DeepLoc
             attention windows) and Fitch retargeting events
  disorder   pLDDT inside that region vs the whole protein, CONTROLLED
             against the same alignment columns in non-focal members —
             protein termini are floppy in every structure, so the
             uncontrolled number means nothing

Classification logic (the PEPC SF3 case is the reference pattern):
  subfunctionalization  relaxed selection (K<1, p<alpha) AND no positive
                        selection on the stem AND an ancestral function
                        partitioned (expression dominance and/or a
                        subfamily-specific signal region)
  neofunctionalization  positive selection on the stem OR retargeting
                        events gained (a NEW state, not a partitioned one)
  conserved/inconclusive selection tests present but neither pattern
  insufficient evidence  no selection tests at all

JSON shapes follow HyPhy 2.5 output (verified against relax_sf3.json /
meme_sf3.json / absrel_sf3.json from the PEPC pilot).
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from utils.alignment import alignment_id, translate_columns

RELAX_ALPHA = 0.05
BRANCH_ALPHA = 0.05
DISORDER_ALPHA = 0.05
EXPRESSION_DOMINANCE = 0.5   # focal subfamily carries >half the family's expression


# ---------------------------------------------------------------------------
# HyPhy JSON parsers
# ---------------------------------------------------------------------------

def parse_relax(path: Path) -> Dict:
    """RELAX json -> {k, p, lrt, direction}."""
    tr = json.loads(Path(path).read_text())["test results"]
    k = float(tr["relaxation or intensification parameter"])
    return {
        "k": k,
        "p": float(tr["p-value"]),
        "lrt": float(tr["LRT"]),
        "direction": "relaxed" if k < 1 else "intensified",
    }


def parse_meme(path: Path, max_p: float = 0.05) -> List[Tuple[int, float]]:
    """MEME json -> [(site_1based, p)] for episodic sites with p <= max_p."""
    mle = json.loads(Path(path).read_text())["MLE"]
    p_idx = next(
        i for i, h in enumerate(mle["headers"]) if h[0].lower() == "p-value"
    )
    sites = []
    for i, row in enumerate(mle["content"]["0"], start=1):
        p = float(row[p_idx])
        if p <= max_p:
            sites.append((i, p))
    return sites


def parse_absrel(path: Path, max_p: float = 0.05) -> List[Tuple[str, float]]:
    """aBSREL json -> [(branch, corrected_p)] significant branches."""
    branches = json.loads(Path(path).read_text())["branch attributes"]["0"]
    out = []
    for name, attrs in branches.items():
        if not isinstance(attrs, dict):
            continue
        p = attrs.get("Corrected P-value")
        if p is not None and float(p) <= max_p:
            out.append((name, float(p)))
    return sorted(out, key=lambda x: x[1])


def apply_branch_names(branches: List[Tuple[str, float]],
                       name_map: Dict[str, str]) -> List[Tuple[str, float]]:
    """Map HyPhy's renamed leaf labels (g###) back to real gene ids."""
    return [(name_map.get(b, b), p) for b, p in branches]


# ---------------------------------------------------------------------------
# Coordinate systems (issue #42)
#
# MEME and BEB site numbers are columns of the alignment codeml/HyPhy was
# given — for the pipeline path that is the family codon alignment
# (steps/codeml.py, five species, untrimmed). The signal region comes from
# extract_signal_windows.py, whose columns belong to the clan protein
# alignment: a different taxon set, a different width, sometimes trimmed. The
# same residue therefore has two different column numbers.
#
# Counting sites inside a region is a plain interval test, so it returns a
# number whatever the two inputs mean, and a mismatch returns zero. classify()
# spends a zero as evidence AGAINST neofunctionalization, so the bug does not
# surface as an error — it surfaces as a specific scientific conclusion.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RegionCount:
    """Sites inside a region, and whether the two agreed on a coordinate system."""
    n_in_region: int
    n_untranslatable: int
    coordinates_verified: bool
    reason: str


def count_sites_in_region(sites: Iterable[int], lo: int, hi: int, *,
                          site_alignment: Optional[Dict[str, str]] = None,
                          region_alignment: Optional[Dict[str, str]] = None,
                          bridge: Optional[str] = None) -> RegionCount:
    """Count `sites` falling inside alignment columns [lo, hi].

    Give both alignments to have the answer verified. Identical stamps mean
    the numbers are directly comparable; different stamps mean they are not,
    and the sites are carried across through `bridge` — a sequence present in
    both with the same residues — rather than compared as if they were. A site
    the bridge cannot carry (a gap there) is counted in `n_untranslatable`,
    never as "outside the region".

    Called without alignments the count is still produced, because a number
    the caller can see is better than none, but `coordinates_verified` is
    False and `classify` will refuse to spend it.
    """
    sites = list(sites)
    if site_alignment is None or region_alignment is None:
        return RegionCount(
            n_in_region=sum(1 for s in sites if lo <= s <= hi),
            n_untranslatable=0,
            coordinates_verified=False,
            reason="no alignment was supplied for the site numbers, the region "
                   "bounds, or both, so the two coordinate systems were never "
                   "shown to be the same one",
        )

    if alignment_id(site_alignment) == alignment_id(region_alignment):
        return RegionCount(
            n_in_region=sum(1 for s in sites if lo <= s <= hi),
            n_untranslatable=0,
            coordinates_verified=True,
            reason="site numbers and region bounds come from the same alignment",
        )

    if bridge is None:
        raise ValueError(
            "The site numbers and the region bounds come from different "
            "alignments, so they cannot be compared directly. Supply a bridge "
            "sequence present in both (bridge=...) to translate, or the count "
            "is meaningless."
        )
    moved = translate_columns(sites, site_alignment, region_alignment,
                             via=bridge)
    return RegionCount(
        n_in_region=sum(1 for m in moved if m is not None and lo <= m <= hi),
        n_untranslatable=sum(1 for m in moved if m is None),
        coordinates_verified=True,
        reason=f"site numbers translated into the region alignment via {bridge!r}",
    )


# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------

def classify(evidence: Dict) -> Dict:
    """Integrate the axes into a verdict with explicit for/against lists."""
    relax = evidence.get("relax")
    stem = evidence.get("absrel_stem", [])
    terminal = evidence.get("absrel_terminal", [])
    expr = evidence.get("expression_share")
    signal = evidence.get("signal_partition", "")
    retarget = evidence.get("retargeting_events", 0)
    meme_region = evidence.get("meme_sites_in_region")
    # Absent means unverified. A coordinate system nobody checked is exactly
    # the case this guard exists for, so it must never default to True.
    coords_ok = bool(evidence.get("coordinates_verified"))

    if relax is None and not stem and not terminal and meme_region is None:
        return {"verdict": "insufficient evidence",
                "evidence_for": [], "evidence_against": [],
                "coordinates_verified": coords_ok, "cannot_judge": []}

    ev_for: List[str] = []
    ev_against: List[str] = []
    cannot_judge: List[str] = []

    def region_is_clean(count, method: str, label: str) -> bool:
        """Is "zero sites in the region" usable as evidence at all?

        Only when the site numbers and the region bounds were shown to belong
        to the same alignment. Otherwise the zero is recorded as something
        that could not be judged — never as a silent negative.
        """
        if count != 0 or not signal:
            return False
        if coords_ok:
            return True
        cannot_judge.append(
            f"{method} reports 0 sites inside the signal region, but the site "
            "numbers and the region bounds were never shown to share a "
            f"coordinate system — cannot judge ({label}). Cross them through "
            "count_sites_in_region with both alignments and a bridge sequence."
        )
        return False

    relaxed = bool(relax and relax["direction"] == "relaxed"
                   and relax["p"] <= RELAX_ALPHA)
    if relaxed:
        ev_for.append(f"relaxed selection (K = {relax['k']:g}, p = {relax['p']:.2g})")
    stem_positive = bool(stem)
    if stem_positive:
        ev_for.append(
            "stem positive selection: "
            + ", ".join(f"{b} (p = {p:.2g})" for b, p in stem)
        )
    partitioned = bool(
        (expr is not None and expr >= EXPRESSION_DOMINANCE) or signal
    )
    if expr is not None and expr >= EXPRESSION_DOMINANCE:
        ev_for.append(f"expression dominance ({expr:.0%} of family expression)")
    if signal:
        ev_for.append(f"signal-region partition: {signal}")
    if region_is_clean(meme_region, "MEME", "meme_sites_in_region"):
        ev_for.append("no episodic positive selection inside the signal region "
                      "(MEME: 0 sites)")
    if terminal:
        ev_against.append(
            "terminal-branch positive selection only: "
            + ", ".join(f"{b} (p = {p:.2g})" for b, p in terminal)
        )

    bs = evidence.get("branchsite")
    if bs:
        # A clade-wide foreground makes branch-site significance ambiguous:
        # it fires whether the stem or any descendant branch carries the
        # signal. aBSREL is what localizes it — so a significant test with
        # no stem branch implicated is context, not a neofunctionalization
        # call.
        if bs["p"] <= RELAX_ALPHA:
            (ev_against if not stem_positive else ev_for).append(
                f"branch-site test significant (LRT = {bs['lrt']:g}, "
                f"p = {bs['p']:.3g}), {bs.get('beb_sites_total', 0)} BEB sites "
                "P>=0.95 — clade-wide foreground, not localized to the stem"
                if not stem_positive else
                f"branch-site test significant (LRT = {bs['lrt']:g}, "
                f"p = {bs['p']:.3g})"
            )
        if region_is_clean(bs.get("beb_sites_in_region"), "The branch-site BEB scan",
                           "beb_sites_in_region"):
            ev_for.append("no BEB positively-selected site inside the signal "
                          "region either (branch-site: 0 sites)")

    dis = evidence.get("disorder")
    if dis and dis.get("p") is not None and dis["p"] <= DISORDER_ALPHA \
            and dis.get("delta", 0) < 0:
        ev_for.append(
            f"signal region is structurally disordered and specifically so "
            f"(mean pLDDT {dis['delta']:+.3f} vs the whole protein, "
            f"p = {dis['p']:.3g} against {dis['n_other']} non-focal members "
            f"at the same alignment columns)"
        )

    if stem_positive or retarget:
        verdict = "neofunctionalization"
        if retarget:
            ev_for.append(f"retargeting events on the family tree: {retarget}")
    elif relaxed and partitioned:
        verdict = "subfunctionalization"
    else:
        verdict = "conserved/inconclusive"

    return {"verdict": verdict, "evidence_for": ev_for,
            "evidence_against": ev_against,
            "coordinates_verified": coords_ok,
            "cannot_judge": cannot_judge}


# ---------------------------------------------------------------------------
# Narrative
# ---------------------------------------------------------------------------

def _relax_point(subfamily: str, relax: Optional[Dict]) -> Optional[str]:
    """Point 1: the selection regime on the focal clade."""
    if not relax:
        return None
    return (
        f"1. **Selection regime.** Selection on the {subfamily} clade is "
        f"RELAXED, not positive (RELAX K = {relax['k']:g}, "
        f"p = {relax['p']:.2g}, LRT = {relax.get('lrt', float('nan')):g}). "
        "Relaxation is what division of labor predicts: once the "
        "ancestral function is split, each copy tolerates change the "
        "full-role ancestor could not afford — without any burst of "
        "adaptive substitutions."
    )


def _meme_point(signal: str, meme_region: Optional[int],
                meme_total: Optional[int], coords_ok: bool) -> Optional[str]:
    """Point 2: MEME found nothing in the region — only sayable if the site
    numbers and the region bounds are known to share a coordinate system."""
    if not (meme_region == 0 and signal and coords_ok):
        return None
    return (
        f"2. **No positive selection where it would matter.** MEME finds "
        f"zero episodic positively-selected sites inside the "
        f"subfamily-specific region ({signal}); "
        + (f"{meme_total} episodic sites exist family-wide but none in "
           "the region. " if meme_total else "")
        + "The region diverged under relaxed constraint, not adaptation."
    )


def _expression_point(subfamily: str, expr: Optional[float]) -> Optional[str]:
    """Point 3: the partition itself, as a share of family expression."""
    if expr is None:
        return None
    return (
        f"3. **The partition itself.** {subfamily} carries {expr:.0%} of "
        "the family's expression — the ancestral expression domain was "
        "divided, and this clade inherited the dominant share."
    )


def _branchsite_point(bs: Optional[Dict], coords_ok: bool) -> Optional[str]:
    """Point 4: what the clade-wide branch-site test can and cannot localize.

    The LRT and the BEB total hold whatever coordinate system the region
    uses; only the "inside the region" clause depends on the two matching,
    so only that clause is withheld.
    """
    if not bs:
        return None
    para = (
        f"4. **What the branch-site test adds — and what it does not.** "
        f"With the whole clade as foreground, codeml branch-site Model A "
        f"rejects the null (LRT = {bs['lrt']:g}, p = {bs['p']:.3g}) and "
        f"BEB identifies {bs.get('beb_sites_total', 0)} sites at "
        f"P >= 0.95. Positive selection therefore exists somewhere in "
        "this clade — but a clade-wide foreground cannot say WHERE, and "
        "aBSREL's per-branch scan places it on a terminal branch, not "
        "the stem."
    )
    if coords_ok:
        return para + (
            " Decisively, "
            f"{bs.get('beb_sites_in_region', 0)} of those BEB sites fall "
            "inside the subfamily-specific signal region: the two "
            "site-level methods (MEME and BEB) independently agree that "
            "the region defining this subfamily is free of positive "
            "selection, while adaptive change sits elsewhere in the "
            "protein body."
        )
    return para + (
        " Where those sites sit relative to the signal region "
        "CANNOT BE JUDGED: the BEB site numbers and the region "
        "bounds were never shown to belong to the same alignment, "
        "and an unverified overlap of zero is indistinguishable "
        "from a coordinate mismatch."
    )


def _disorder_point(dis: Optional[Dict]) -> Optional[str]:
    """Point 5: the region is disordered, and — against the control — only here."""
    if not (dis and dis.get("p") is not None and dis["p"] <= DISORDER_ALPHA
            and dis.get("delta", 0) < 0):
        return None
    return (
        f"5. **The region is disordered, and only here.** Mean pLDDT "
        f"inside the region runs {dis['delta']:+.3f} below the rest of "
        f"the protein across {dis['n_focal']}/{dis['n_focal']} focal "
        "members. Protein termini are floppy in every structure, so "
        "that alone would prove nothing — the control is what counts: "
        f"the same alignment columns in {dis['n_other']} non-focal "
        "clan members show no such drop "
        f"(Mann-Whitney p = {dis['p']:.3g}). A disordered, "
        "subfamily-specific segment is what relaxed constraint "
        "produces; adaptive remodelling would leave an ordered one."
    )


def _terminal_only_point(stem: List, terminal: List) -> Optional[str]:
    """Closing paragraph: positive selection AFTER the split, not its cause."""
    if not (stem == [] and terminal):
        return None
    return (
        "The only positive selection detected is on terminal branches ("
        + ", ".join(f"{b}, p = {p:.2g}" for b, p in terminal)
        + ") — species-level fine-tuning AFTER the split, not the cause "
        "of it. The stem branch that defines the subfamily shows no "
        "positive selection, which rules out classic "
        "neofunctionalization for the split itself."
    )


def _subfunctionalization_body(subfamily: str, evidence: Dict,
                               coords_ok: bool) -> List[str]:
    """The numbered case for a partitioned ancestral role."""
    lines = [
        f"{subfamily} shows the subfunctionalization pattern: the ancestral "
        "role was PARTITIONED among daughter copies rather than a new "
        "function being gained. Three independent lines support this:",
        "",
    ]
    points = [
        _relax_point(subfamily, evidence.get("relax")),
        _meme_point(evidence.get("signal_partition", ""),
                    evidence.get("meme_sites_in_region"),
                    evidence.get("meme_sites_total"), coords_ok),
        _expression_point(subfamily, evidence.get("expression_share")),
        _branchsite_point(evidence.get("branchsite"), coords_ok),
        _disorder_point(evidence.get("disorder")),
        _terminal_only_point(evidence.get("absrel_stem", []),
                             evidence.get("absrel_terminal", [])),
    ]
    return lines + [p for p in points if p is not None]


def _neofunctionalization_body(subfamily: str, verdict: Dict) -> List[str]:
    """A NEW state was gained; the verdict's own evidence carries the case."""
    return [
        f"{subfamily} shows a neofunctionalization signature: a NEW state "
        "was gained rather than an ancestral role divided."
    ] + [f"- {item}" for item in verdict["evidence_for"]]


def _undecided_body(verdict: Dict) -> List[str]:
    """Neither pattern: say so, and show both sides rather than nothing."""
    return (
        ["The available evidence does not distinguish a functional shift "
         f"({verdict['verdict']})."]
        + [f"- for: {item}" for item in verdict["evidence_for"]]
        + [f"- against: {item}" for item in verdict["evidence_against"]]
    )


def _evidence_summary(verdict: Dict) -> List[str]:
    """The trailing ledger.

    Withheld evidence is listed, not omitted: an axis that could not be
    judged has to be visible, or its absence reads as an absent signal.
    """
    return (
        ["", "Evidence summary:"]
        + [f"- [supports] {item}" for item in verdict["evidence_for"]]
        + [f"- [context]  {item}" for item in verdict["evidence_against"]]
        + [f"- [cannot judge] {item}" for item in verdict.get("cannot_judge", [])]
    )


def narrative(family: str, subfamily: str, verdict: Dict, evidence: Dict) -> str:
    """A self-contained explanation of HOW the subfamily diverged."""
    call = verdict["verdict"]
    lines = [f"## {family} / {subfamily}: {call}", ""]

    if call == "subfunctionalization":
        lines += _subfunctionalization_body(
            subfamily, evidence, bool(verdict.get("coordinates_verified")))
    elif call == "neofunctionalization":
        lines += _neofunctionalization_body(subfamily, verdict)
    else:
        lines += _undecided_body(verdict)

    return "\n".join(lines + _evidence_summary(verdict))
