"""Expression share per family member — with its own coverage on the label.

The finding this axis was built for is real: in *Opuntia cochenillifera* the
single SF3 copy `OcoChr03G21370` carries 4,777 TPM, 74% of all PEPC
expression in the family, while the nine-copy SF2 array together makes 1,491
and the four SF1 copies 178. That is expression subfunctionalization measured
rather than asserted.

What the axis must not do is imply the other species were measured. Of the
seventeen species in the PEPC clan (fifteen Caryophyllales plus the ATH and
Aco anchors), exactly TWO have an expression matrix:

    Mcry   iceplant_10032022_Lomas.featurecounts.tpm.tsv   114 samples
    Ococ   rnaseq_29405/clustering/RSEM_TPM.average.tsv       7 timepoints

The other fifteen have none, and a member of one of them must come back as
`expression unavailable`, never as a missing row and never as zero: a silent
drop makes an unmeasured gene indistinguishable from a silent one. The
coverage fraction ("2 of 17" on the clan) therefore travels on the summary,
counted from the family actually given rather than assumed.

Shares are always computed WITHIN one species. TPM is normalised against its
own transcriptome, so a total pooled across species is not a quantity.

Which matrix, exactly (issue #38): the published 74% comes from the 29,365-gene
diel reanalysis (`rnaseq_29405/clustering/RSEM_TPM.average.tsv`, the matrix
behind supplementary Table S14), NOT from the older 33,746-gene
`results/rnaseq/RSEM_TPM.tsv`. The two disagree about whole genes —
`OcoChr03G21430` is 709.6 TPM in the first and identically zero in the second
— so run on the older matrix the same family gives SF3 3,822 TPM and 82.8%,
not 4,777 and 74.1%. The matrix path is an input here precisely so an answer
stays attached to the matrix that produced it.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from utils.gene_ids import match_ids

logger = logging.getLogger("family_finder")

UNAVAILABLE = "expression unavailable"
MEASURED = "measured"
ABSENT = "absent from the species matrix"


@dataclass(frozen=True)
class Matrix:
    """One species' expression matrix: sample names and per-gene values."""
    species: str
    path: str
    samples: Tuple[str, ...]
    values: Dict[str, Tuple[float, ...]]


def read_matrix(path: Path, species: str) -> Matrix:
    """Read a TSV whose first column is the gene id and the rest are samples.

    Both matrices in hand are of this shape — featurecounts TPM for Mcry and
    RSEM TPM for Ococ — so no format flag is needed, only the species the
    matrix belongs to, which is what keeps shares from crossing species.
    """
    samples: Tuple[str, ...] = ()
    values: Dict[str, Tuple[float, ...]] = {}
    with open(path) as handle:
        for index, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            if index == 0:
                samples = tuple(fields[1:])
                continue
            if len(fields) < 2:
                continue
            values[fields[0]] = tuple(_number(v) for v in fields[1:])
    return Matrix(species=species, path=str(path), samples=samples,
                  values=values)


def _number(text: str) -> float:
    try:
        return float(text)
    except ValueError:
        return float("nan")


def gene_means(matrix: Matrix, gene_ids: Sequence[str],
               max_unmatched: Optional[float] = None
               ) -> Tuple[Dict[str, float], dict]:
    """Mean expression across all samples, per requested gene.

    Ids go through `utils.gene_ids.match_ids` because the pep ids carry a
    transcript suffix the matrices do not (`Ococ_OcoChr03G21370.t1` vs
    `OcoChr03G21370`), and because a join that loses members quietly turns a
    measured family into a smaller one without saying so (#42).
    """
    if not gene_ids:
        return {}, {"level": "exact", "n_requested": 0, "n_unmatched": 0,
                    "unmatched": [], "unmatched_fraction": 0.0}

    stripped = {gene: _strip_species(gene) for gene in gene_ids}
    matched = match_ids(list(stripped.values()), list(matrix.values))

    means: Dict[str, float] = {}
    for gene, key in stripped.items():
        reference = matched.mapping.get(key)
        if reference is None:
            continue
        row = matrix.values[reference]
        if row:
            means[gene] = sum(row) / len(row)

    unmatched = [gene for gene in gene_ids if gene not in means]
    fraction = len(unmatched) / len(gene_ids)
    report = {
        "level": matched.level,
        "n_requested": len(gene_ids),
        "n_matched": len(means),
        "n_unmatched": len(unmatched),
        "unmatched": unmatched,
        "unmatched_fraction": fraction,
        "n_samples": len(matrix.samples),
        "matrix": matrix.path,
    }
    if max_unmatched is not None and fraction > max_unmatched:
        raise ValueError(
            f"{len(unmatched)} of {len(gene_ids)} {matrix.species} gene ids "
            f"stayed unmatched against {matrix.path} ({fraction:.0%} > "
            f"{max_unmatched:.0%}): {unmatched[:5]} — refusing rather than "
            "reporting a share over a silently reduced family"
        )
    return means, report


def _strip_species(gene_id: str) -> str:
    _species, _, rest = gene_id.partition("_")
    return rest.lstrip("_")


def species_of(gene_id: str) -> str:
    return gene_id.split("_", 1)[0]


def species_coverage(gene_ids: Sequence[str],
                     matrices: Dict[str, Matrix]) -> dict:
    """How much of the family expression could speak to — stated, not implied."""
    species = sorted({species_of(gene) for gene in gene_ids})
    with_matrix = [s for s in species if s in matrices]
    without = [s for s in species if s not in matrices]
    return {
        "n_species": len(species),
        "n_species_with_matrix": len(with_matrix),
        "species_with_matrix": with_matrix,
        "species_without_matrix": without,
        "coverage": f"{len(with_matrix)} of {len(species)}",
        "note": (
            "Expression is measured in "
            f"{len(with_matrix)} of {len(species)} species. Members of the "
            f"other {len(without)} are reported as {UNAVAILABLE!r} — that is "
            "an absent measurement, not an absent transcript, and the two "
            "must not be read as the same thing."
        ),
    }


def subfamily_shares(means: Dict[str, float],
                     groups: Dict[str, List[str]]) -> List[dict]:
    """Summed mean TPM per group, and each group's share of the total.

    The share denominator is the sum over the genes that actually have a
    value, and `n_unmeasured` says how many members did not — a group whose
    members were never measured must not look quiet.
    """
    total = sum(means.values())
    rows = []
    for name in sorted(groups):
        members = groups[name]
        measured = [means[m] for m in members if m in means]
        subtotal = sum(measured)
        rows.append({
            "group": name,
            "n_members": len(members),
            "n_measured": len(measured),
            "n_unmeasured": len(members) - len(measured),
            "total_tpm": subtotal,
            "share": (subtotal / total) if total else None,
            "max_tpm": max(measured) if measured else None,
        })
    return rows


def family_expression(gene_ids: Sequence[str], matrices: Dict[str, Matrix],
                      groups: Optional[Dict[str, List[str]]] = None,
                      max_unmatched: Optional[float] = None) -> dict:
    """Per-gene mean TPM and within-species share, plus coverage and groups.

    There is deliberately no family-wide total in the result. Summing TPM
    across species would be a number with no unit behind it, and once printed
    somebody would quote it.
    """
    by_species: Dict[str, List[str]] = {}
    for gene in gene_ids:
        by_species.setdefault(species_of(gene), []).append(gene)

    rows: List[dict] = []
    per_species: Dict[str, dict] = {}

    for species in sorted(by_species):
        members = by_species[species]
        matrix = matrices.get(species)
        if matrix is None:
            for gene in members:
                rows.append(_row(gene, species, None, None, UNAVAILABLE,
                                 groups, matrix=""))
            continue

        means, report = gene_means(matrix, members, max_unmatched=max_unmatched)
        total = sum(means.values())
        for gene in members:
            if gene in means:
                rows.append(_row(gene, species, means[gene],
                                 (means[gene] / total) if total else None,
                                 MEASURED, groups, matrix=matrix.path,
                                 n_samples=len(matrix.samples)))
            else:
                rows.append(_row(gene, species, None, None, ABSENT, groups,
                                 matrix=matrix.path,
                                 n_samples=len(matrix.samples)))

        entry = {
            "matrix": matrix.path,
            "n_samples": len(matrix.samples),
            "n_members": len(members),
            "n_measured": len(means),
            "family_total_tpm": total,
            "id_join": report,
        }
        if groups:
            entry["groups"] = subfamily_shares(
                means, {name: [m for m in gs if species_of(m) == species]
                        for name, gs in groups.items()})
            # Members of the group in species with no matrix are counted as
            # unmeasured here too, so a group is never made to look small
            # because thirteen of its species have no RNA-seq.
            for row in entry["groups"]:
                row["n_members"] = len(groups[row["group"]])
                row["n_unmeasured"] = row["n_members"] - row["n_measured"]
        per_species[species] = entry

    return {
        "rows": sorted(rows, key=lambda r: r["gene"]),
        "by_species": per_species,
        "coverage": species_coverage(gene_ids, matrices),
    }


def _row(gene: str, species: str, mean: Optional[float],
         share: Optional[float], status: str,
         groups: Optional[Dict[str, List[str]]], matrix: str,
         n_samples: Optional[int] = None) -> dict:
    subfamily = ""
    for name, members in (groups or {}).items():
        if gene in members:
            subfamily = name
            break
    return {
        "gene": gene,
        "species": species,
        "subfamily": subfamily,
        "mean_tpm": mean,
        "share": share,
        "n_samples": n_samples,
        "matrix": matrix,
        "status": status,
    }


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

COLUMNS = ("gene", "species", "subfamily", "mean_tpm", "share", "n_samples",
           "matrix", "status")


def _cell(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def write_table(rows: Sequence[dict], path: Path) -> None:
    """TSV for `annotation_matrix.py --expression`."""
    with open(path, "w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(_cell(row.get(c)) for c in COLUMNS) + "\n")


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--members", required=True,
                        help="File of family gene ids, one per line")
    parser.add_argument("--matrix", action="append", default=[],
                        metavar="SPECIES=PATH",
                        help="Expression matrix for one species; repeatable")
    parser.add_argument("--groups", default=None,
                        help="TSV of gene<TAB>subfamily")
    parser.add_argument("-o", "--out", required=True, help="Output TSV")
    parser.add_argument("--json", default=None, help="Also write the summary")
    parser.add_argument("--max-unmatched", type=float, default=None,
                        help="Refuse when more than this fraction of a "
                             "species' members fail to reach its matrix")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")

    members = [line.strip() for line in open(args.members) if line.strip()]

    matrices: Dict[str, Matrix] = {}
    for spec in args.matrix:
        species, _, path = spec.partition("=")
        if not path:
            parser.error(f"--matrix expects SPECIES=PATH, got {spec!r}")
        matrices[species] = read_matrix(Path(path), species)

    groups: Optional[Dict[str, List[str]]] = None
    if args.groups:
        groups = {}
        for line in open(args.groups):
            if not line.strip():
                continue
            gene, _, label = line.rstrip("\n").partition("\t")
            groups.setdefault(label, []).append(gene)

    report = family_expression(members, matrices, groups=groups,
                               max_unmatched=args.max_unmatched)
    write_table(report["rows"], Path(args.out))
    if args.json:
        Path(args.json).write_text(json.dumps(
            {"coverage": report["coverage"], "by_species": report["by_species"]},
            indent=2))

    coverage = report["coverage"]
    logger.info("expression coverage: %s species (%s); no matrix for %s",
                coverage["coverage"],
                ", ".join(coverage["species_with_matrix"]) or "none",
                ", ".join(coverage["species_without_matrix"]) or "none")
    for species, entry in sorted(report["by_species"].items()):
        logger.info("  %s: %d/%d members measured over %d samples, "
                    "family total %.1f TPM", species, entry["n_measured"],
                    entry["n_members"], entry["n_samples"],
                    entry["family_total_tpm"])
        for row in entry.get("groups", []):
            logger.info("    %-8s %8.1f TPM  share %s  (%d measured of %d)",
                        row["group"], row["total_tpm"],
                        "NA" if row["share"] is None else f"{row['share']:.1%}",
                        row["n_measured"], row["n_members"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
