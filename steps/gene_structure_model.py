"""The one gene model both `steps.gene_structure` and `steps.gff_join` use.

Its own module only so the two can import it without a cycle.
"""

from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class GeneModel:
    """One gene's coding blocks and the programme that predicted them."""
    gene_id: str
    transcript_id: str
    source: str                       # GFF column 2 — the covariate
    strand: str
    blocks: Tuple[Tuple[int, int], ...]
