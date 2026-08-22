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
from collections import Counter
from typing import Dict, List, Optional, Sequence, Tuple

GAP = "-"


def _most_common(residues: Sequence[str]) -> Tuple[str, int]:
    return Counter(residues).most_common(1)[0]


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


def invariant_columns(alignment: Dict[str, str],
                      groups: Dict[str, Sequence[str]],
                      *, threshold: float = 0.95,
                      min_cover: float = 0.7) -> Dict[int, str]:
    """Columns every group holds at the same residue: `{column: residue}`.

    This is the reference-free stand-in for "the functionally constrained
    core". Saying that a diagnostic residue avoids the catalytic site requires
    knowing which motif is catalytic; saying it avoids the columns on which
    every subfamily independently agrees requires only the alignment, and is
    therefore available for every family a run produces.

    Agreement is judged **within each group separately** and then compared
    across groups, not pooled. Pooling would let a large group carry a column
    on its own — the same size bias that makes a global occupancy threshold
    unusable (see `select_columns`). A group covering the column below
    `min_cover` disqualifies it, because absence is not agreement.

    `threshold` below 1.0 is deliberate: on real data a functional motif
    carries the occasional substitution, and strict identity would miss most
    of one. Measured on the PEPC clan, strict identity keeps 255 columns and
    0.95 keeps 422, while the catalytic motif GYSDSGKDAGR itself contains a
    column at 0.53.
    """
    width = _columns(alignment)
    if not alignment:
        return {}
    if not any(m in alignment for members in groups.values() for m in members):
        raise ValueError(
            f"The {len(groups)} supplied group(s) have no members in the "
            f"alignment. Groups: {sorted(groups)}"
        )
    present = {name: [m for m in members if m in alignment]
               for name, members in groups.items()}
    out: Dict[int, str] = {}
    for i in range(width):
        agreed: Optional[str] = None
        for members in present.values():
            if not members:
                agreed = None
                break
            column = [alignment[m][i] for m in members]
            residues = [c for c in column if c != GAP]
            if not residues or len(residues) / len(members) < min_cover:
                agreed = None
                break
            top, count = _most_common(residues)
            if count / len(residues) < threshold:
                agreed = None
                break
            if agreed is None:
                agreed = top
            elif agreed != top:
                agreed = None
                break
        if agreed is not None:
            out[i + 1] = agreed
    return out


def choose_reference(alignment: Dict[str, str],
                     candidates: Optional[Sequence[str]] = None) -> str:
    """Pick the sequence to number diagnostic positions against.

    A family with no characterised member still needs a stable coordinate
    system, so the reference is chosen from the family itself: the most
    complete sequence, i.e. the one with the most non-gap residues. Coverage
    is the criterion rather than typicality (a medoid) because the reference's
    only job is to give as many columns as possible a citable residue number —
    a column where the reference is gapped cannot be cited at all.

    Ties break on the identifier so the same alignment always yields the same
    reference, and therefore the same published positions.
    """
    pool = list(candidates) if candidates is not None else list(alignment)
    missing = [c for c in pool if c not in alignment]
    if missing:
        raise ValueError(f"Candidate sequences not in the alignment: {missing}")
    if not pool:
        raise ValueError("Cannot choose a reference from an empty alignment")
    return min(pool, key=lambda name: (-sum(1 for c in alignment[name] if c != GAP),
                                       name))
