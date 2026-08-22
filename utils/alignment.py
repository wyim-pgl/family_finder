"""Alignment column selection, occupancy reporting and coordinate mapping.

Issue #42. Trimming exists to help tree inference; residue-level analysis must
not inherit the trimmed matrix. When a matrix does have to be reduced, the
column filter must be group-aware.

A **global** occupancy threshold cannot preserve a region unique to a minority
subfamily: if that subfamily is k of N sequences, the region's global occupancy
is capped at k/N, so any threshold above k/N deletes it. This is arithmetic,
not a tuning problem, and it removes exactly the columns that carry
subfamily-specific signal. Measured on the PEPC clan: ppc-1E2 is 17 of 87
sequences (20%), so anything unique to it tops out at 0.20.

The group-aware rule keeps a column when **any** group covers it, which is
scale-free — a five-member subfamily's own region is 100% occupied within that
group. Only columns empty in every group are genuinely uninformative.
"""
from typing import Dict, List, Optional, Sequence, Tuple

GAP = "-"


def _columns(alignment: Dict[str, str]) -> int:
    lengths = {len(s) for s in alignment.values()}
    if len(lengths) > 1:
        raise ValueError(f"Alignment has unequal sequence lengths: {sorted(lengths)}")
    return lengths.pop() if lengths else 0


def column_occupancy(alignment: Dict[str, str]) -> List[float]:
    """Fraction of sequences with a residue (not a gap) in each column."""
    width = _columns(alignment)
    n = len(alignment)
    if not n:
        return []
    return [sum(1 for s in alignment.values() if s[i] != GAP) / n
            for i in range(width)]


def group_occupancy(alignment: Dict[str, str],
                    groups: Dict[str, Sequence[str]]) -> Dict[str, List[float]]:
    """Per-column occupancy computed within each group separately.

    Group members absent from the alignment are ignored; a group with no
    member present yields all zeros rather than an error, so it can never
    rescue a column.
    """
    width = _columns(alignment)
    out: Dict[str, List[float]] = {}
    for name, members in groups.items():
        present = [m for m in members if m in alignment]
        if not present:
            out[name] = [0.0] * width
            continue
        out[name] = [sum(1 for m in present if alignment[m][i] != GAP) / len(present)
                     for i in range(width)]
    return out


def select_columns(alignment: Dict[str, str], *, threshold: float,
                   groups: Optional[Dict[str, Sequence[str]]] = None) -> List[int]:
    """Columns to keep, as 1-based positions.

    Without `groups` this is a global occupancy filter — use it only for the
    matrix that feeds tree inference. With `groups` a column survives when any
    group covers it at `threshold`, which is what preserves minority-subfamily
    regions. All-gap columns are dropped in both modes, at any threshold.
    """
    width = _columns(alignment)
    if not alignment:
        return []
    overall = column_occupancy(alignment)
    per_group = None
    if groups:
        per_group = group_occupancy(alignment, groups)
        if not any(m in alignment for members in groups.values() for m in members):
            raise ValueError(
                f"The {len(groups)} supplied group(s) have no members in the "
                "alignment — selecting on them would silently return an empty "
                f"matrix. Groups: {sorted(groups)}"
            )
    kept = []
    for i in range(width):
        if overall[i] == 0.0:
            continue  # nothing to keep, whatever the threshold says
        if per_group is None:
            passes = overall[i] >= threshold
        else:
            passes = any(occ[i] >= threshold for occ in per_group.values())
        if passes:
            kept.append(i + 1)
    return kept


def apply_columns(alignment: Dict[str, str], columns: Sequence[int]) -> Dict[str, str]:
    """Subset every sequence to `columns` (1-based)."""
    width = _columns(alignment)
    bad = [c for c in columns if c < 1 or c > width]
    if bad:
        raise ValueError(f"Columns outside the alignment (width {width}): {bad}")
    idx = [c - 1 for c in columns]
    return {name: "".join(seq[i] for i in idx) for name, seq in alignment.items()}


def column_map(columns: Sequence[int],
               reference: str) -> List[Tuple[int, int, Optional[int]]]:
    """`(position in the kept matrix, original column, reference residue)`.

    Reference numbering is ungapped, so a column where the reference itself is
    gapped maps to None. Emitting this alongside a trimmed matrix is what makes
    a diagnostic position citable later without redoing the alignment by hand.
    """
    ref_pos: Dict[int, int] = {}
    n = 0
    for col, ch in enumerate(reference, start=1):
        if ch != GAP:
            n += 1
            ref_pos[col] = n
    return [(k, col, ref_pos.get(col)) for k, col in enumerate(columns, start=1)]
