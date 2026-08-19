"""Minimal Newick parsing, writing, and Fitch parsimony (issue #16).

Self-contained (no ete4) so retargeting detection and its tests run outside
the pipeline conda environment. Handles the FastTree output this pipeline
produces: unquoted leaf names, internal support labels, branch lengths.
"""

from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Set


@dataclass
class Node:
    name: str = ""
    dist: float = 0.0
    children: List["Node"] = field(default_factory=list)
    mark: str = ""                 # codeml foreground tag, e.g. "#1"
    state: Optional[str] = None    # Fitch-resolved state

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self) -> Iterator["Node"]:
        if self.is_leaf:
            yield self
        else:
            for child in self.children:
                yield from child.leaves()

    def leaf_names(self) -> Set[str]:
        return {leaf.name for leaf in self.leaves()}

    def postorder(self) -> Iterator["Node"]:
        for child in self.children:
            yield from child.postorder()
        yield self


_DELIMS = ",():;"


def parse_newick(text: str) -> Node:
    """Parse a Newick string into a Node tree."""
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = 0

    def parse_node() -> Node:
        nonlocal pos
        node = Node()
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                node.children.append(parse_node())
                if pos < len(s) and s[pos] == ",":
                    pos += 1
                    continue
                if pos < len(s) and s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"Malformed Newick near position {pos}")
        # name (leaf name, or internal support label)
        start = pos
        while pos < len(s) and s[pos] not in _DELIMS:
            pos += 1
        node.name = s[start:pos]
        # branch length
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in _DELIMS:
                pos += 1
            node.dist = float(s[start:pos])
        return node

    root = parse_node()
    if pos != len(s):
        raise ValueError(f"Trailing characters in Newick at position {pos}")
    return root


def write_newick(
    root: Node, with_dist: bool = True, keep_internal_names: bool = True
) -> str:
    """Serialize a Node tree back to Newick. Node.mark tags (e.g. '#1') are
    emitted after the name — the placement codeml expects for foreground
    branches. For codeml trees use with_dist=False, keep_internal_names=False
    (PAML chokes on support labels and does not need branch lengths)."""

    def rec(n: Node) -> str:
        if n.is_leaf:
            base = n.name
        else:
            inner = ",".join(rec(c) for c in n.children)
            label = n.name if keep_internal_names else ""
            base = f"({inner}){label}"
        if n.mark:
            base += f" {n.mark}"
        if with_dist:
            base += f":{n.dist:g}"
        return base

    return rec(root) + ";"


def fitch(root: Node, leaf_states: Dict[str, Optional[str]]) -> None:
    """Fitch parsimony ancestral-state reconstruction, in place.

    Leaves whose state is None (no DeepLoc call, excluded, or ambiguous) are
    treated as uninformative: they contribute nothing and inherit downward.
    After the call every node has .state set (possibly None if the whole
    tree is uninformative).
    """
    statesets: Dict[int, Set[str]] = {}
    # bottom-up
    for node in root.postorder():
        if node.is_leaf:
            st = leaf_states.get(node.name)
            statesets[id(node)] = {st} if st else set()
        else:
            child_sets = [statesets[id(c)] for c in node.children if statesets[id(c)]]
            if not child_sets:
                statesets[id(node)] = set()
            else:
                inter = set.intersection(*child_sets)
                statesets[id(node)] = inter if inter else set.union(*child_sets)

    # top-down resolution (deterministic: sorted first when parent state absent)
    def assign(node: Node, parent_state: Optional[str]) -> None:
        ss = statesets[id(node)]
        if not ss:
            node.state = parent_state if node.is_leaf is False else None
            # uninformative leaves keep None; internal gaps inherit the parent
            if node.is_leaf:
                node.state = None
        elif parent_state in ss:
            node.state = parent_state
        else:
            node.state = sorted(ss)[0]
        for child in node.children:
            assign(child, node.state)

    assign(root, None)


def state_switches(root: Node) -> List[tuple]:
    """Branches where the Fitch state changes: (parent_state, child_state,
    frozenset of the child clade's leaf names)."""
    events = []

    def rec(node: Node) -> None:
        for child in node.children:
            if node.state and child.state and node.state != child.state:
                events.append((node.state, child.state, frozenset(child.leaf_names())))
            rec(child)

    rec(root)
    return events


def mark_clade(root: Node, clade_leaves: Set[str], mark: str = "#1") -> bool:
    """Tag the node whose subtree contains exactly clade_leaves. Returns
    False if no such clade exists in this tree."""
    target = set(clade_leaves)
    for node in root.postorder():
        if node.leaf_names() == target:
            node.mark = mark
            return True
    return False
