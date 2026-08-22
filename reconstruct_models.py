#!/usr/bin/env python3
"""Rebuild full-length gene models from characterised anchors, and check them.

Issue #40. Six sequences entered the PEPC clan set carrying an `_M` suffix,
each replacing a truncated released annotation (250-834 aa) with a full-length
protein (801-1,042 aa). No script produced them, and neither this repository,
the lab wiki, the manuscript directory nor `ppc_resolve/README.md` mentions the
convention — so models underpinning a published claim could not be regenerated.

Checking them against their own genomes split them cleanly:

    0.9968  0.9937  0.9892  0.9817  |  0.8550  0.7898
      genuinely that species' gene   |   not that species at all

The two on the right were wrong. The *O. cochenillifera* one is 98.8% to the
*O. basilaris* locus and absent from all four O. cochenillifera assemblies,
including the pre-purge one that carries 266 MB more sequence; it was the
evidence for a cactus ppc-1E2 member that does not exist.

The tool was never the problem — four of six were right. The missing piece was
a check that the output is encoded in the genome it claims to come from. That
check is the point of this module.

Procedure:

    miniprot --gff --trans --outn=10 --outs=0.5 <species genome> <anchors.faa>

Anchors should be characterised full-length proteins (for PEPC: SwissProt
P10490 / P16097 via Mcry_Mcr8G11630 / Mcry_Mcr7G08600), never the models being
repaired — seeding with those makes the check circular.

Then parse, keep the best prediction per locus, and accept only what clears the
gate. Identity is miniprot's, i.e. of the protein against the genome it was
predicted from; a real gene sits at 0.98-0.997, allowing for splice-boundary
and frameshift losses.
"""
import argparse
import re
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

# A real gene predicted from its own genome lands at 0.98-0.997; the two bad
# models sat at 0.855 and 0.790. 0.95 separates them with room on both sides.
MIN_IDENTITY = 0.95
MIN_LENGTH = 300

_ATTR = re.compile(r"(\w+)=([^;\s]+)")


@dataclass(frozen=True)
class Prediction:
    contig: str
    start: int
    end: int
    strand: str
    identity: float
    rank: int
    anchor: str
    protein: str
    frameshifts: int = 0


def parse_miniprot(text: str) -> List[Prediction]:
    """Parse `miniprot --gff --trans` output.

    Each mRNA row is preceded by its `##STA` translated protein. An mRNA with
    no preceding `##STA` is skipped rather than paired with a stale one — the
    protein is the whole product here, and guessing it would be worse than
    dropping the row.
    """
    preds: List[Prediction] = []
    protein: Optional[str] = None
    for line in text.splitlines():
        if line.startswith("##STA"):
            parts = line.split("\t", 1)
            protein = parts[1].strip() if len(parts) > 1 else None
            continue
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "mRNA":
            continue
        if protein is None:
            continue
        attrs = dict(_ATTR.findall(fields[8]))
        target = attrs.get("Target", "")
        preds.append(Prediction(
            contig=fields[0], start=int(fields[3]), end=int(fields[4]),
            strand=fields[6], identity=float(attrs.get("Identity", 0.0)),
            rank=int(attrs.get("Rank", 0)), anchor=target.split()[0] if target else "",
            protein=protein, frameshifts=int(attrs.get("Frameshift", 0)),
        ))
        protein = None
    return preds


def accept(prediction: Prediction, min_identity: float = MIN_IDENTITY,
           min_length: int = MIN_LENGTH,
           max_frameshifts: int = 0) -> Tuple[bool, str]:
    """Is this model supported by the genome it was predicted from?

    Frameshifted predictions are rejected by default. `--trans` skips
    frameshifts when writing the protein, so the GFF CDS intervals no longer
    splice back into it: all four models reconstructed for issue #40 missed
    their own protein by 3 to 14 residues and could not enter the codon
    alignment. Raise `max_frameshifts` only when the protein alone is wanted.
    """
    if prediction.frameshifts > max_frameshifts:
        return False, (
            f"{prediction.frameshifts} frameshift(s) — the translated protein "
            "skips them, so its CDS cannot be recovered from the GFF and the "
            "model is unusable in a codon alignment"
        )
    if prediction.identity < min_identity:
        return False, (
            f"identity {prediction.identity:.4f} below {min_identity} — the "
            "protein is not encoded in this genome; it most likely belongs to "
            "another species or accession"
        )
    if "*" in prediction.protein.rstrip("*"):
        return False, "internal stop codon in the predicted protein"
    if len(prediction.protein.rstrip("*")) < min_length:
        return False, (
            f"length {len(prediction.protein)} below {min_length} — the point "
            "of the repair is to replace a truncated model, not emit one"
        )
    return True, "accepted"


def _overlaps(a: Prediction, b: Prediction) -> bool:
    return a.contig == b.contig and a.start <= b.end and b.start <= a.end


def best_per_locus(predictions: Iterable[Prediction]) -> List[Prediction]:
    """One prediction per genomic locus, keeping the highest identity.

    Several anchors hit the same gene, so the raw output holds near-duplicates
    that would otherwise be counted as separate loci.
    """
    kept: List[Prediction] = []
    for p in sorted(predictions, key=lambda x: -x.identity):
        if not any(_overlaps(p, k) for k in kept):
            kept.append(p)
    return sorted(kept, key=lambda x: (x.contig, x.start))


def main(argv: Optional[Sequence[str]] = None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gff", help="miniprot --gff --trans output")
    ap.add_argument("--species", required=True, help="prefix for emitted ids")
    ap.add_argument("--min-identity", type=float, default=MIN_IDENTITY)
    ap.add_argument("--min-length", type=int, default=MIN_LENGTH)
    ap.add_argument("-o", "--out", required=True, help="accepted proteins (FASTA)")
    ap.add_argument("--rejected", help="rejected predictions with reasons (TSV)")
    args = ap.parse_args(argv)

    preds = best_per_locus(parse_miniprot(open(args.gff).read()))
    kept, dropped = [], []
    for p in preds:
        ok, reason = accept(p, args.min_identity, args.min_length)
        (kept if ok else dropped).append((p, reason))

    with open(args.out, "w") as f:
        for p, _ in kept:
            f.write(f">{args.species}_{p.contig}_{p.start} "
                    f"identity={p.identity:.4f} anchor={p.anchor} "
                    f"locus={p.contig}:{p.start}-{p.end}{p.strand}\n{p.protein}\n")
    if args.rejected:
        with open(args.rejected, "w") as f:
            f.write("contig\tstart\tend\tidentity\tanchor\tlength\treason\n")
            for p, reason in dropped:
                f.write(f"{p.contig}\t{p.start}\t{p.end}\t{p.identity:.4f}\t"
                        f"{p.anchor}\t{len(p.protein)}\t{reason}\n")

    print(f"loci: {len(preds)}  accepted: {len(kept)}  rejected: {len(dropped)}")
    for p, reason in dropped:
        print(f"  REJECT {p.contig}:{p.start} id={p.identity:.4f} — {reason}")


if __name__ == "__main__":
    main()
