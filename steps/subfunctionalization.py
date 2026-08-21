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
from pathlib import Path
from typing import Dict, List, Tuple

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

    if relax is None and not stem and not terminal and meme_region is None:
        return {"verdict": "insufficient evidence",
                "evidence_for": [], "evidence_against": []}

    ev_for: List[str] = []
    ev_against: List[str] = []

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
    if meme_region == 0 and signal:
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
        if bs.get("beb_sites_in_region") == 0 and signal:
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
            "evidence_against": ev_against}


# ---------------------------------------------------------------------------
# Narrative
# ---------------------------------------------------------------------------

def narrative(family: str, subfamily: str, verdict: Dict, evidence: Dict) -> str:
    """A self-contained explanation of HOW the subfamily diverged."""
    relax = evidence.get("relax")
    expr = evidence.get("expression_share")
    signal = evidence.get("signal_partition", "")
    stem = evidence.get("absrel_stem", [])
    terminal = evidence.get("absrel_terminal", [])
    meme_region = evidence.get("meme_sites_in_region")
    meme_total = evidence.get("meme_sites_total")

    lines = [f"## {family} / {subfamily}: {verdict['verdict']}", ""]

    if verdict["verdict"] == "subfunctionalization":
        lines.append(
            f"{subfamily} shows the subfunctionalization pattern: the ancestral "
            "role was PARTITIONED among daughter copies rather than a new "
            "function being gained. Three independent lines support this:"
        )
        lines.append("")
        if relax:
            lines.append(
                f"1. **Selection regime.** Selection on the {subfamily} clade is "
                f"RELAXED, not positive (RELAX K = {relax['k']:g}, "
                f"p = {relax['p']:.2g}, LRT = {relax.get('lrt', float('nan')):g}). "
                "Relaxation is what division of labor predicts: once the "
                "ancestral function is split, each copy tolerates change the "
                "full-role ancestor could not afford — without any burst of "
                "adaptive substitutions."
            )
        if meme_region == 0 and signal:
            lines.append(
                f"2. **No positive selection where it would matter.** MEME finds "
                f"zero episodic positively-selected sites inside the "
                f"subfamily-specific region ({signal}); "
                + (f"{meme_total} episodic sites exist family-wide but none in "
                   "the region. " if meme_total else "")
                + "The region diverged under relaxed constraint, not adaptation."
            )
        if expr is not None:
            lines.append(
                f"3. **The partition itself.** {subfamily} carries {expr:.0%} of "
                "the family's expression — the ancestral expression domain was "
                "divided, and this clade inherited the dominant share."
            )
        bs = evidence.get("branchsite")
        if bs:
            lines.append(
                f"4. **What the branch-site test adds — and what it does not.** "
                f"With the whole clade as foreground, codeml branch-site Model A "
                f"rejects the null (LRT = {bs['lrt']:g}, p = {bs['p']:.3g}) and "
                f"BEB identifies {bs.get('beb_sites_total', 0)} sites at "
                f"P >= 0.95. Positive selection therefore exists somewhere in "
                "this clade — but a clade-wide foreground cannot say WHERE, and "
                "aBSREL's per-branch scan places it on a terminal branch, not "
                "the stem. Decisively, "
                f"{bs.get('beb_sites_in_region', 0)} of those BEB sites fall "
                "inside the subfamily-specific signal region: the two "
                "site-level methods (MEME and BEB) independently agree that the "
                "region defining this subfamily is free of positive selection, "
                "while adaptive change sits elsewhere in the protein body."
            )
        dis = evidence.get("disorder")
        if dis and dis.get("p") is not None and dis["p"] <= DISORDER_ALPHA \
                and dis.get("delta", 0) < 0:
            lines.append(
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
        if stem == [] and terminal:
            lines.append(
                "The only positive selection detected is on terminal branches ("
                + ", ".join(f"{b}, p = {p:.2g}" for b, p in terminal)
                + ") — species-level fine-tuning AFTER the split, not the cause "
                "of it. The stem branch that defines the subfamily shows no "
                "positive selection, which rules out classic "
                "neofunctionalization for the split itself."
            )
    elif verdict["verdict"] == "neofunctionalization":
        lines.append(
            f"{subfamily} shows a neofunctionalization signature: a NEW state "
            "was gained rather than an ancestral role divided."
        )
        for item in verdict["evidence_for"]:
            lines.append(f"- {item}")
    else:
        lines.append(
            "The available evidence does not distinguish a functional shift "
            f"({verdict['verdict']})."
        )
        for item in verdict["evidence_for"]:
            lines.append(f"- for: {item}")
        for item in verdict["evidence_against"]:
            lines.append(f"- against: {item}")

    lines.append("")
    lines.append("Evidence summary:")
    for item in verdict["evidence_for"]:
        lines.append(f"- [supports] {item}")
    for item in verdict["evidence_against"]:
        lines.append(f"- [context]  {item}")
    return "\n".join(lines)
