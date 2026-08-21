"""EC prediction sources: eggNOG-mapper + CLEAN parsers (issue #28).

ESM-ECForest was rejected on its known-answer gate (issue #20: EC 4.1.1.31
never predicted, anchor misclassification, uniform-noise probabilities).
Its replacements — both PASSED the same gate on the PEPC clan (2026-08-20):

- eggNOG-mapper (orthology transfer): anchors 7/7 -> 4.1.1.31, PPC symbols
  recovered; no probability, so ``confidence`` stays 0.0.
- CLEAN (contrastive embedding, Science 2023): anchors 7/7 -> 4.1.1.31 at
  conf 0.986-0.996; always emits an EC, so the CONFIDENCE carries the
  signal (remote clan members drop to ~0).

Both parsers return the ``steps.ecforest`` prediction-dict shape
(``{gene: {is_enzyme, ec, confidence}}``) so the downstream EC layer
(``annotate_families.py --ec``: family consensus + Fitch EC-switch events,
save_cache/load_cache round-trip) is reused unchanged.

Shared caveat, measured on SF2 array genes (21430/21450/21460): NEITHER
source sees residue-level catalytic loss (missing PF00311 / catalytic His)
— all three still get 4.1.1.31. EC calls must be paired with the domain /
catalytic-residue evidence layer; they are annotations, never filters.
"""

import csv
import logging
from pathlib import Path
from typing import Dict

logger = logging.getLogger("family_finder")

_EMAPPER_ID_COL = "#query"
_EMAPPER_EC_COL = "EC"


def parse_emapper(path: Path) -> Dict[str, dict]:
    """eggNOG-mapper ``.emapper.annotations`` -> ecforest prediction dicts.

    ``##`` preamble/trailer lines are skipped; the ``#query`` header row
    names the columns. emapper's EC field is ``-`` when absent and may hold
    a comma-separated list — the first (most specific transfer) is kept.
    emapper reports no probability, so confidence is 0.0.
    """
    preds: Dict[str, dict] = {}
    header = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##") or not line:
                continue
            if line.startswith(_EMAPPER_ID_COL):
                header = line.split("\t")
                continue
            if header is None:
                continue
            row = dict(zip(header, line.split("\t")))
            ec_raw = (row.get(_EMAPPER_EC_COL) or "").strip()
            ec = "" if ec_raw in ("", "-") else ec_raw.split(",")[0].strip()
            preds[row[_EMAPPER_ID_COL]] = {
                "is_enzyme": bool(ec),
                "ec": ec,
                "confidence": 0.0,
            }
    if header is None:
        raise ValueError(f"No #query header found in emapper file: {path}")
    return preds


def parse_clean(path: Path) -> Dict[str, dict]:
    """CLEAN maxsep CSV (``gene,EC:x.x.x.x/conf[;EC:.../conf]``) -> dicts.

    CLEAN always predicts an EC — the GMM confidence is the real signal
    (near-zero for remote members). Multiple candidates keep the highest-
    confidence one. ``is_enzyme`` is True whenever an EC parses; consumers
    gate on confidence (ecforest._gene_ec min_conf), not on the flag.
    """
    preds: Dict[str, dict] = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            gene, _, rest = line.partition(",")
            best_ec, best_conf = "", -1.0
            for cand in rest.split(";"):
                cand = cand.strip()
                if not cand.startswith("EC:"):
                    continue
                body = cand[3:]
                ec, _, conf_s = body.partition("/")
                try:
                    conf = float(conf_s)
                except ValueError:
                    conf = 0.0
                if conf > best_conf:
                    best_ec, best_conf = ec.strip(), conf
            if best_ec:
                preds[gene] = {
                    "is_enzyme": True,
                    "ec": best_ec,
                    "confidence": best_conf,
                }
    return preds


def merge_ec_predictions(
    emapper: Dict[str, dict], clean: Dict[str, dict]
) -> Dict[str, dict]:
    """Merge the two sources per gene; emapper's EC is authoritative.

    Orthology transfer outranks embedding similarity (the SF2 lesson cuts
    both ways, but emapper's call is at least evidence-traceable to a seed
    ortholog). CLEAN contributes its confidence — the one number emapper
    lacks — and an ``agree`` flag marks source agreement. Extra keys
    (source/agree) survive save_cache round-trips by being dropped, so the
    downstream ecforest cache contract is unchanged.
    """
    merged: Dict[str, dict] = {}
    for gene in sorted(set(emapper) | set(clean)):
        e, c = emapper.get(gene), clean.get(gene)
        if e and c:
            ec = e["ec"] or c["ec"]
            merged[gene] = {
                "is_enzyme": bool(ec),
                "ec": ec,
                "confidence": c["confidence"],
                "source": "emapper+clean",
                "agree": bool(e["ec"]) and e["ec"] == c["ec"],
            }
        elif e:
            merged[gene] = {**e, "source": "emapper", "agree": None}
        else:
            merged[gene] = {**c, "source": "clean", "agree": None}
    return merged
