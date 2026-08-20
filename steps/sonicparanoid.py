"""SonicParanoid2 tier-1 clustering adapter (issue #22).

The v2 plan (issue #21) runs SonicParanoid2 head-to-head against
OrthoFinder on the 5sp panel; whichever wins on our metrics becomes the
tier-1 default (config.clustering_method). This module is the ADAPTER
half only: a template-driven runner, a parser for SonicParanoid2's
ortholog_groups.tsv, and a converter to OrthoFinder-style Orthogroups.tsv
so every existing consumer (steps.orthofinder.parse_orthogroups,
score_recluster.py, check_pairs.py) reads the output unchanged. Nothing
here redefines homology — groups still pass through align/tree/prune.

SonicParanoid2 group table layout — PROVISIONAL, unverified against a
real run (the project wiki was unreachable; revisit after the #22
head-to-head produces an actual ortholog_groups.tsv): a header row with
the group-id column followed by per-species columns (headers may be plain
prefixes or embedded file names like "Mcry.faa"); one group per row;
cells are comma-separated gene ids with "*" for an empty cell. parse_groups
is deliberately defensive: known metadata columns are skipped by header
name, short rows are tolerated, and cell contents are taken VERBATIM —
nothing is stripped, so any inline per-gene confidence annotation survives
in the parsed ids. If a real run reveals a dedicated per-gene confidence
COLUMN, extend the parser to carry it rather than dropping the column.
"""

import logging
from pathlib import Path
from typing import Dict, Iterable, List, Set

from config import Config
from steps.esm import _run_template

logger = logging.getLogger("family_finder")

# SonicParanoid2 marks an empty per-species cell with "*".
EMPTY_CELL = "*"

# Non-species metadata columns some SonicParanoid versions insert between
# the group id and the per-species columns; skipped defensively so the
# generic rule (column 0 = group id, the rest = species) stays safe.
_METADATA_COLUMNS = {"group_id", "group_size", "sp_in_grp", "seed_ortholog_cnt"}


def run_sonicparanoid(pep_dir: Path, outdir: Path, config: Config) -> Path:
    """Run SonicParanoid2 on a directory of per-species proteome FASTAs.

    Driven by the config.sonicparanoid_cmd TEMPLATE (cluster installs
    vary), e.g. "sonicparanoid -i {pep_dir} -o {outdir} -t 16". Raises
    ValueError when the template is not configured. Returns outdir; the
    groups table is located afterwards with find_groups_file.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    _run_template(
        config.sonicparanoid_cmd,
        {"pep_dir": str(pep_dir), "outdir": str(outdir)},
        "SonicParanoid2",
    )
    return outdir


def find_groups_file(sp_outdir: Path) -> Path:
    """Locate ortholog_groups.tsv under a SonicParanoid output tree.

    SonicParanoid2 nests results as runs/<run_name>/ortholog_groups/
    ortholog_groups.tsv; the run name is timestamped, so search
    recursively and take the newest match (repeat runs append).
    """
    candidates = sorted(
        Path(sp_outdir).rglob("ortholog_groups.tsv"),
        key=lambda p: p.stat().st_mtime,
    )
    if not candidates:
        raise FileNotFoundError(
            f"No ortholog_groups.tsv found under {sp_outdir}"
        )
    return candidates[-1]


def parse_groups(groups_file: Path) -> Dict[str, List[str]]:
    """Parse SonicParanoid2's ortholog_groups.tsv into group_id -> gene ids.

    Generic parse: first column is the group id; every remaining column is
    a per-species comma-separated gene list ("*" = empty), except known
    metadata columns (_METADATA_COLUMNS), which are skipped by header
    name. Headers may be plain species prefixes or file names — cells are
    parsed identically either way and gene ids are kept verbatim.
    Groups with no genes are dropped.
    """
    groups: Dict[str, List[str]] = {}
    with open(groups_file) as f:
        header_line = f.readline().rstrip("\n")
        if not header_line.strip():
            raise ValueError(f"Empty header in {groups_file}")
        header = header_line.split("\t")
        gene_cols = [
            i for i, name in enumerate(header)
            if i > 0 and name.strip().lower() not in _METADATA_COLUMNS
        ]
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            group_id = parts[0].strip()
            genes: List[str] = []
            for i in gene_cols:
                if i >= len(parts):
                    continue
                cell = parts[i].strip()
                if not cell or cell == EMPTY_CELL:
                    continue
                genes.extend(
                    g.strip() for g in cell.split(",")
                    if g.strip() and g.strip() != EMPTY_CELL
                )
            if genes:
                groups[group_id] = genes
    logger.info(f"Parsed {len(groups)} SonicParanoid groups from {groups_file}")
    return groups


def _group_sort_key(group_id: str):
    """Numeric SonicParanoid group ids sort numerically, others lexically."""
    return (0, int(group_id), "") if group_id.isdigit() else (1, 0, group_id)


def write_orthogroups_tsv(
    groups: Dict[str, List[str]],
    species_order: List[str],
    out_tsv: Path,
    unassigned_out: Path,
    all_genes: Iterable[str],
    species_delimiter: str = "_",
) -> None:
    """Convert groups to OrthoFinder-style Orthogroups*.tsv files.

    out_tsv gets the OrthoFinder layout (header "Orthogroup" + one column
    per species in species_order; genes joined with ", ") so
    steps.orthofinder.parse_orthogroups, score_recluster.py and
    check_pairs.py consume it unchanged. unassigned_out gets the same
    header with one singleton row per unassigned gene (= all_genes minus
    every grouped gene), ids UNG0000000... A gene's species is its
    SpeciesPrefix_GeneID prefix (house convention); genes whose prefix is
    not in species_order are counted as grouped but cannot be written into
    a column — logged as warnings and omitted from out_tsv.
    """
    out_tsv = Path(out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    unassigned_out = Path(unassigned_out)
    unassigned_out.parent.mkdir(parents=True, exist_ok=True)

    header = "Orthogroup\t" + "\t".join(species_order) + "\n"
    grouped: Set[str] = set()

    with open(out_tsv, "w") as f:
        f.write(header)
        for group_id in sorted(groups, key=_group_sort_key):
            by_species: Dict[str, List[str]] = {sp: [] for sp in species_order}
            for gene in groups[group_id]:
                grouped.add(gene)
                species = gene.split(species_delimiter, 1)[0]
                if species in by_species:
                    by_species[species].append(gene)
                else:
                    logger.warning(
                        f"Group {group_id}: gene {gene} has species prefix "
                        f"{species!r} not in species_order — omitted from "
                        f"{out_tsv.name}"
                    )
            f.write(
                group_id + "\t"
                + "\t".join(", ".join(by_species[sp]) for sp in species_order)
                + "\n"
            )

    unassigned = sorted(set(all_genes) - grouped)
    with open(unassigned_out, "w") as f:
        f.write(header)
        for i, gene in enumerate(unassigned):
            species = gene.split(species_delimiter, 1)[0]
            row = [gene if species == sp else "" for sp in species_order]
            f.write(f"UNG{i:07d}\t" + "\t".join(row) + "\n")

    logger.info(
        f"Wrote {len(groups)} groups to {out_tsv} and "
        f"{len(unassigned)} unassigned genes to {unassigned_out}"
    )
