"""Retargeting event detection on family gene trees (issue #16).

Maps DeepLoc localizations onto confirmed gene trees, reconstructs ancestral
states with Fitch parsimony, and reports branches where the localization
switches — the neofunctionalization candidates fed to the branch-site test.

Safeguard: genes in the exclude set (truncated models from pseudogene
detection) are treated as state-less. An N-terminally truncated model loses
its transit peptide and would otherwise produce a FAKE retargeting call.
"""

import logging
from dataclasses import dataclass
from typing import Dict, FrozenSet, List, Optional, Set

from steps.deeploc import top_location
from utils.newick import fitch, parse_newick, state_switches

logger = logging.getLogger("family_finder")


@dataclass(frozen=True)
class RetargetingEvent:
    family_id: str
    parent_state: str
    child_state: str
    clade_leaves: FrozenSet[str]

    @property
    def clade_id(self) -> str:
        return ",".join(sorted(self.clade_leaves))


def find_retargeting_events(
    family_id: str,
    newick_text: str,
    predictions: Dict[str, dict],
    exclude: Set[str] = frozenset(),
    min_prob: float = 0.5,
) -> List[RetargetingEvent]:
    """Retargeting events for one family gene tree.

    Args:
        family_id: family label carried into the events.
        newick_text: the family's confirmed_tree.nwk contents.
        predictions: DeepLoc predictions (steps.deeploc.parse_deeploc/load_cache).
        exclude: gene IDs to treat as state-less (truncation filter — mandatory
            safeguard; pass the pseudogene-flagged set).
        min_prob: minimum class probability for a leaf to carry a state.
    """
    root = parse_newick(newick_text)

    states: Dict[str, Optional[str]] = {}
    for leaf in root.leaf_names():
        if leaf in exclude:
            states[leaf] = None
        else:
            states[leaf] = top_location(predictions.get(leaf), min_prob)

    informative = {s for s in states.values() if s}
    if len(informative) < 2:
        return []  # uniform or uninformative — no switch possible

    fitch(root, states)
    return [
        RetargetingEvent(family_id, parent, child, leaves)
        for parent, child, leaves in state_switches(root)
    ]


def write_events_tsv(events: List[RetargetingEvent], path) -> None:
    with open(path, "w") as f:
        f.write("family_id\tparent_state\tchild_state\tn_clade_genes\tclade_genes\n")
        for e in events:
            f.write(
                f"{e.family_id}\t{e.parent_state}\t{e.child_state}\t"
                f"{len(e.clade_leaves)}\t{e.clade_id}\n"
            )
