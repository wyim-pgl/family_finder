#!/usr/bin/env python3
"""Subfamily diagnostics for one family: WHY does it split?

Standalone (no ete4, no external tools): takes a family's protein
alignment plus a gene->subfamily table (from Possvm / TreeCluster / any
splitter) and writes three TSVs into --outdir:

  sdp_residues.tsv          diagnostic residues per subfamily (with
                            reference numbering when --ref-seq is given)
  taxonomy_attribution.tsv  per subfamily: lineage-specific vs
                            paralog-split. Default evidence (issue #27)
                            is the species tree (--species-tree, the
                            same data/species_tree.nwk the pipeline
                            requires): a subfamily whose species set is
                            a clade strictly inside the family span is
                            lineage-specific; non-monophyletic or
                            family-span sets are paralog splits. A
                            taxonomy TSV (--taxonomy) is an optional
                            label layer (rank names); without a tree it
                            alone drives the legacy rank-purity verdict.
  structure_coherence.tsv   (with --pairs) foldseek all-vs-all within- vs
                            between-subfamily coherence
  anchor_transferability.tsv (with --family-tree and --anchors) whether the
                            reference subfamily names designate anything in
                            the query species, and how well established each
                            assignment is: HIGH / PROVISIONAL / UNRESOLVED
                            from the weaker of SH-aLRT/UFBoot plus membership
                            recovered from the independent datasets given with
                            --cross-tree (with none, that condition reads NOT
                            EVALUATED and nothing can grade HIGH)

With --focal-subfamily plus the selection-evidence inputs (--relax-json,
--meme-json/--meme-region, --absrel-json, --expression-share,
--signal-partition, --disorder-json), also writes subfunctionalization.md — a narrative
explaining HOW the subfamily diverged (sub- vs neo-functionalization
verdict with the per-axis evidence) — and subfunctionalization.tsv.

Example (PEPC clan):
  python subfamily_report.py \
    --alignment clan_anchor.aln --groups possvm_groups.tsv \
    --species-tree data/species_tree.nwk --taxonomy taxonomy.tsv \
    --pairs foldseek_allvsall.tsv \
    --ref-seq ATH_AT1G53310.2 --outdir subfam_report
"""

import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from steps.subfamily import (  # noqa: E402
    GRADE_MIN_SH_ALRT,
    GRADE_PROVISIONAL_MIN_UFBOOT,
    anchor_transferability,
    load_taxonomy,
    parse_pair_identities,
    parse_pairwise_scores,
    resolve_groups,
    resolve_reference,
    sdp_scan,
    sequence_confounding,
    structure_coherence,
    DEFAULT_PAIR_COLUMNS,
    UNCONTROLLED_WARNING,
    taxonomic_composition,
)
from steps.subfunctionalization import (  # noqa: E402
    apply_branch_names,
    classify,
    count_sites_in_region,
    narrative,
    parse_absrel,
    parse_meme,
    parse_relax,
)
from steps.codeml import lrt_pvalue  # noqa: E402
from utils.newick import parse_newick  # noqa: E402
from utils.seqio import read_fasta  # noqa: E402


def read_groups(path: Path) -> dict:
    groups: dict = {}
    seen: dict = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        gene, subfamily = line.rstrip("\n").split("\t")[:2]
        # A Possvm/TreeCluster export often carries a header; parsing it as
        # data invents a gene named "gene_id" in a phantom subfamily.
        if gene in ("gene", "gene_id"):
            continue
        if gene in seen and seen[gene] != subfamily:
            raise ValueError(
                f"{path}: gene {gene!r} is assigned to both {seen[gene]!r} "
                f"and {subfamily!r} — a combined groups file carrying two "
                "labels for one gene must be resolved, not split across "
                "subfamilies"
            )
        if gene in seen:
            continue
        seen[gene] = subfamily
        groups.setdefault(subfamily, []).append(gene)
    return groups


def write_tsv(rows, path: Path):
    if not rows:
        path.write_text("")
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _region_bounds(region: str):
    """'LO-HI' -> (lo, hi)."""
    lo, hi = region.split("-")
    return int(lo), int(hi)


def read_branch_name_map(path: Path) -> dict:
    """TSV short<TAB>real — HyPhy renames leaves to g### for its own runs."""
    name_map = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            short, real = line.rstrip("\n").split("\t")[:2]
            name_map[short] = real
    return name_map


def _coordinate_frame(args) -> dict:
    """The alignments the site numbers and the region bounds each belong to.

    MEME/BEB sites are columns of whatever alignment the selection test ran
    on; `--meme-region` bounds are columns of the alignment the region was
    read off, which for the PEPC work was the clan protein alignment and for
    the pipeline path is not the same matrix at all (#42). Both must be named
    before an overlap between them means anything.
    """
    site_aln = read_fasta(args.sites_alignment) if args.sites_alignment else None
    region_path = (getattr(args, "region_alignment", None)
                   or getattr(args, "alignment", None))
    region_aln = read_fasta(region_path) if (site_aln and region_path) else None
    return {"site_alignment": site_aln, "region_alignment": region_aln,
            "bridge": getattr(args, "coord_bridge", None)}


def collect_selection_evidence(args) -> dict:
    """Gather whatever selection/partition evidence the CLI was given.

    Every axis is optional; classify() degrades honestly on what is absent.
    """
    evidence: dict = {
        "expression_share": args.expression_share,
        "signal_partition": args.signal_partition or "",
        "retargeting_events": args.retargeting_events,
    }
    frame = _coordinate_frame(args)
    # Verified until an axis says otherwise; with no region at all there is
    # nothing to verify and nothing that depends on it.
    verified = True

    if args.relax_json:
        evidence["relax"] = parse_relax(Path(args.relax_json))

    if args.meme_json:
        sites = parse_meme(Path(args.meme_json))
        evidence["meme_sites_total"] = len(sites)
        if args.meme_region:
            lo, hi = _region_bounds(args.meme_region)
            count = count_sites_in_region([s for s, _p in sites], lo, hi, **frame)
            evidence["meme_sites_in_region"] = count.n_in_region
            evidence["meme_sites_untranslatable"] = count.n_untranslatable
            verified = verified and count.coordinates_verified

    if args.absrel_json:
        branches = parse_absrel(Path(args.absrel_json))
        if args.branch_name_map:
            branches = apply_branch_names(
                branches, read_branch_name_map(Path(args.branch_name_map))
            )
        # HyPhy names internal nodes Node###; anything else is a leaf. The
        # stem/terminal split is what separates neofunctionalization from
        # post-split lineage fine-tuning.
        evidence["absrel_stem"] = [b for b in branches if b[0].startswith("Node")]
        evidence["absrel_terminal"] = [
            b for b in branches if not b[0].startswith("Node")
        ]

    if args.branchsite_mlc:
        from beb_cross import parse_beb  # noqa: E402

        sig = [(site, p) for site, _aa, p in
               parse_beb(Path(args.branchsite_mlc).read_text()) if p >= 0.95]
        bs: dict = {"beb_sites_total": len(sig)}
        if args.meme_region:
            lo, hi = _region_bounds(args.meme_region)
            count = count_sites_in_region([s for s, _p in sig], lo, hi, **frame)
            bs["beb_sites_in_region"] = count.n_in_region
            bs["beb_sites_untranslatable"] = count.n_untranslatable
            verified = verified and count.coordinates_verified
        if args.branchsite_lnl:
            alt, null = (float(x) for x in args.branchsite_lnl.split(","))
            bs["lrt"] = 2 * (alt - null)
            bs["p"] = lrt_pvalue(alt, null)
        evidence["branchsite"] = bs

    if args.disorder_json:
        d = json.loads(Path(args.disorder_json).read_text())
        evidence["disorder"] = {
            k: d[k] for k in ("delta", "n_focal", "n_other", "p", "all_below")
            if k in d
        }

    # classify() refuses to spend "0 sites in the region" as counter-evidence
    # for neofunctionalization unless this is True — an unverified zero is
    # indistinguishable from a coordinate mismatch.
    evidence["coordinates_verified"] = bool(verified and args.meme_region)
    return evidence


def build_structure_coherence(args, groups: dict):
    """Coherence rows plus the sequence-identity confounding report.

    Foldseek scores track sequence identity, and subfamilies are groups of
    close relatives by construction, so the raw within/between ratio can be
    true and uninformative at the same time (#39). When the pair table carries
    `fident` the comparison is redone on identity-adjusted residuals; when it
    does not, every row carries an explicit UNCONTROLLED warning rather than a
    bare ratio somebody might quote.

    Returns `(rows, confounding)`, where `confounding` is None if the table
    had no identities to measure the coupling with.
    """
    pairs = Path(args.pairs)
    columns = getattr(args, "pair_columns", None)
    metric = getattr(args, "pair_metric", None) or "bits"
    scores = parse_pairwise_scores(pairs, metric=metric, columns=columns)
    try:
        identities = parse_pair_identities(pairs, columns=columns)
    except ValueError:
        identities, confounding = None, None
    else:
        confounding = sequence_confounding(scores, identities)
    rows = structure_coherence(groups, scores, identities=identities,
                               metric=metric)
    return rows, confounding


def read_anchors(path: Path) -> dict:
    """TSV gene_id<TAB>subfamily_label — reference genes that name groups."""
    anchors = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        gene, label = line.rstrip("\n").split("\t")[:2]
        anchors[gene] = label
    return anchors


def _cross_tree_memberships(specs, anchors, query_prefixes) -> dict:
    """`NAME=tree.nwk` -> {NAME: {label: members}}.

    The same clade rule is run on each independent dataset, so "consistent"
    means the membership was RECOVERED there, not asserted by hand — and
    "recovered" requires the replicate's OWN clade to be transferable and
    support-qualified. A replicate that finds the same leaf set at 40/60,
    or behind a blocking rival, did not establish the membership in that
    dataset; counting it as agreement would let it promote the main tree's
    grade to HIGH on evidence the second dataset never provided. Filtered
    labels are simply absent from that dataset's dict, which the agreement
    check already reads as "not recovered" (disagreeing).
    """
    out: dict = {}
    for spec in specs or []:
        name, _, path = spec.partition("=")
        if not path:
            raise ValueError(f"--cross-tree wants NAME=path, got {spec!r}")
        rows = anchor_transferability(Path(path).read_text(), anchors,
                                      query_prefixes=query_prefixes)
        out[name] = {
            r["label"]: r["members"] for r in rows
            if r["label"] and r["transferable"]
            and r["sh_alrt"] is not None and r["ufboot"] is not None
            and r["sh_alrt"] >= GRADE_MIN_SH_ALRT
            and r["ufboot"] >= GRADE_PROVISIONAL_MIN_UFBOOT
        }
    return out


def build_transferability(args) -> list:
    """anchor_transferability rows, flattened for the TSV.

    Empty when no --family-tree/--anchors were given: the question was not
    asked, which is not the same as every name being fine.
    """
    if not (getattr(args, "family_tree", None)
            and getattr(args, "anchors", None)):
        return []
    anchors = read_anchors(Path(args.anchors))
    prefixes = (args.query_prefix.split(",") if args.query_prefix else None)
    rows = anchor_transferability(
        Path(args.family_tree).read_text(), anchors,
        query_prefixes=prefixes, min_support=args.anchor_min_support,
        cross_tree_members=_cross_tree_memberships(
            args.cross_tree, anchors, prefixes),
        tree_name=args.tree_name,
    )
    flat = []
    for r in rows:
        flat.append({k: (";".join(str(x) for x in v) if isinstance(v, list)
                         else "" if v is None else str(v))
                     for k, v in r.items()})
    return flat


def write_subfunctionalization(args, verdict: dict, evidence: dict,
                               outdir: Path) -> None:
    """subfunctionalization.md (narrative) + .tsv (machine-readable)."""
    text = narrative(args.family_name, args.focal_subfamily, verdict, evidence)
    (outdir / "subfunctionalization.md").write_text(text + "\n")
    with open(outdir / "subfunctionalization.tsv", "w") as f:
        f.write("subfamily\tverdict\tcoordinates_verified\tn_evidence_for\t"
                "n_evidence_against\tn_cannot_judge\tevidence_for\t"
                "evidence_against\tcannot_judge\n")
        f.write("\t".join([
            args.focal_subfamily, verdict["verdict"],
            str(verdict.get("coordinates_verified", False)),
            str(len(verdict["evidence_for"])),
            str(len(verdict["evidence_against"])),
            str(len(verdict.get("cannot_judge", []))),
            "; ".join(verdict["evidence_for"]),
            "; ".join(verdict["evidence_against"]),
            "; ".join(verdict.get("cannot_judge", [])),
        ]) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--alignment", required=True,
                    help="Family protein alignment (FASTA)")
    ap.add_argument("--groups", required=True,
                    help="TSV: gene_id<TAB>subfamily")
    ap.add_argument("--species-tree", default=None,
                    help="Species tree (Newick, leaves = species "
                         "prefixes) — default verdict evidence: clade "
                         "test on each subfamily's species set")
    ap.add_argument("--taxonomy", default=None,
                    help="TSV: species<TAB>genus<TAB>family<TAB>order — "
                         "label layer with --species-tree, sole verdict "
                         "evidence without it (legacy)")
    ap.add_argument("--pairs", default=None,
                    help="foldseek easy-search all-vs-all TSV. Include "
                         "fident in --format-output to get the "
                         "sequence-identity control; without it the coherence "
                         "table is emitted with an UNCONTROLLED warning")
    ap.add_argument("--pair-columns", default=None,
                    help="The foldseek --format-output string --pairs was "
                         f"written with (default: {DEFAULT_PAIR_COLUMNS})")
    ap.add_argument("--pair-metric", default="bits",
                    help="Pair-table column to score structural similarity on. "
                         "bits grows with alignment length AND identity; "
                         "alntmscore is length-normalised and less confounded "
                         "(default: bits)")
    ap.add_argument("--ref-seq", default=None,
                    help="Alignment sequence id for reference numbering. "
                         "Without it a characterised member is preferred, "
                         "then the family's own representative.")
    ap.add_argument("--characterised", default=None,
                    help="File of characterised gene ids, one per line "
                         "(SwissProt-backed anchors and the like). The most "
                         "complete one becomes the coordinate reference, so "
                         "positions can be compared with published ones.")
    # --- can the reference names be transferred, and how well established? ---
    ap.add_argument("--family-tree", default=None,
                    help="Family/clan gene tree (Newick, IQ-TREE "
                         "SH-aLRT/UFBoot node labels) — with --anchors, "
                         "writes anchor_transferability.tsv")
    ap.add_argument("--anchors", default=None,
                    help="TSV: reference_gene_id<TAB>subfamily_label")
    ap.add_argument("--anchor-min-support", type=float, default=None,
                    help="Support bar below which a clade is not called "
                         "transferable at all (the grade is reported "
                         "regardless, on its own thresholds)")
    ap.add_argument("--query-prefix", default=None,
                    help="Comma-separated species prefixes that count as "
                         "query genes (default: every non-anchor leaf)")
    ap.add_argument("--cross-tree", action="append", default=None,
                    metavar="NAME=TREE",
                    help="Independent dataset to test membership against "
                         "(amino acid / codon 1+2 / codon 1+2+3), repeatable. "
                         "Without any, cross-tree consistency is reported as "
                         "NOT EVALUATED and no assignment can grade HIGH")
    ap.add_argument("--tree-name", default="this tree",
                    help="Name of the dataset --family-tree came from, used "
                         "in the consistency report")
    ap.add_argument("--min-group", type=int, default=5)
    ap.add_argument("--max-unmatched", type=float, default=None,
                    help="Refuse when more than this fraction of --groups "
                         "members fail to match an alignment id. Without it "
                         "the loss is reported and the run continues.")
    ap.add_argument("--delimiter", default="_")
    ap.add_argument("--outdir", required=True)
    # --- subfunctionalization narrative (optional evidence axes) ---
    ap.add_argument("--focal-subfamily", default=None,
                    help="Subfamily to explain (writes subfunctionalization.md)")
    ap.add_argument("--family-name", default="family",
                    help="Family label used in the narrative")
    ap.add_argument("--relax-json", default=None, help="HyPhy RELAX json")
    ap.add_argument("--meme-json", default=None, help="HyPhy MEME json")
    ap.add_argument("--absrel-json", default=None, help="HyPhy aBSREL json")
    ap.add_argument("--branch-name-map", default=None,
                    help="TSV short<TAB>real to restore HyPhy-renamed leaves")
    ap.add_argument("--meme-region", default=None,
                    help="Alignment-column range 'LO-HI' of the subfamily-"
                         "specific signal region (MEME sites counted inside)")
    ap.add_argument("--sites-alignment", default=None,
                    help="Alignment whose columns the MEME/BEB site numbers "
                         "are (the matrix the selection test ran on, as "
                         "protein columns). Required before a count inside "
                         "--meme-region can be treated as evidence.")
    ap.add_argument("--region-alignment", default=None,
                    help="Alignment whose columns --meme-region refers to "
                         "(default: --alignment)")
    ap.add_argument("--coord-bridge", default=None,
                    help="Sequence id present in BOTH alignments, used to "
                         "translate site numbers when the two differ")
    ap.add_argument("--expression-share", type=float, default=None,
                    help="Focal subfamily's share of family expression (0-1)")
    ap.add_argument("--signal-partition", default=None,
                    help="One-line description of the signal-region partition")
    ap.add_argument("--branchsite-mlc", default=None,
                    help="codeml branch-site Model A output (BEB section)")
    ap.add_argument("--branchsite-lnl", default=None,
                    help="'ALT,NULL' lnL pair for the branch-site LRT")
    ap.add_argument("--disorder-json", default=None,
                    help="region_disorder.py output — pLDDT of the signal "
                         "region vs the rest, controlled against non-focal "
                         "members at the same alignment columns")
    ap.add_argument("--retargeting-events", type=int, default=0,
                    help="Fitch retargeting events gained by the focal clade")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    alignment = read_fasta(args.alignment)
    groups = read_groups(Path(args.groups))
    print(f"{len(alignment)} sequences, {len(groups)} subfamilies")

    # Resolve the group ids against the alignment ONCE, loudly. A '.t1' the
    # alignment does not carry used to drop those members from every layer
    # below without a word (#34), and a subfamily nobody looked at reads
    # exactly like a subfamily with no signal.
    groups, id_report = resolve_groups(alignment, groups,
                                       max_unmatched=args.max_unmatched)
    print(f"gene-id matching: {id_report['level']} level, "
          f"{id_report['n_unmatched']}/{id_report['n_requested']} members "
          f"unmatched ({id_report['unmatched_fraction']:.1%})")
    for name, info in sorted(id_report["groups"].items()):
        if info["unmatched"]:
            print(f"  {name}: {len(info['unmatched'])} unmatched, e.g. "
                  f"{info['unmatched'][:3]}")

    characterised = ([l.strip() for l in open(args.characterised) if l.strip()]
                     if args.characterised else None)
    reference = resolve_reference(alignment, ref_seq_id=args.ref_seq,
                                  characterised=characterised)
    print(f"coordinate reference: {reference.seq_id} "
          f"({reference.source} — {reference.reason})")
    if reference.unmatched_characterised:
        print(f"  {len(reference.unmatched_characterised)} characterised id(s) "
              f"not in the alignment, e.g. "
              f"{reference.unmatched_characterised[:3]}")

    sdp = sdp_scan(alignment, groups, min_group=args.min_group,
                   ref_seq_id=reference.seq_id)
    write_tsv(sdp, outdir / "sdp_residues.tsv")
    print(f"sdp_residues.tsv: {len(sdp)} diagnostic columns, positions in "
          f"{reference.seq_id}")

    taxonomy = load_taxonomy(Path(args.taxonomy)) if args.taxonomy else {}
    species_tree = (parse_newick(Path(args.species_tree).read_text())
                    if args.species_tree else None)
    tax_rows = taxonomic_composition(groups, taxonomy, args.delimiter,
                                     species_tree=species_tree)
    write_tsv(tax_rows, outdir / "taxonomy_attribution.tsv")
    source = "species tree" if species_tree else "taxonomy rank purity"
    print(f"taxonomy_attribution.tsv (verdict evidence: {source})")
    for r in tax_rows:
        label = f"  [{r['clade_label']}]" if r.get("clade_label") else ""
        print(f"  {r['subfamily']}: {r['verdict']}{label}")

    if args.pairs:
        coh, confounding = build_structure_coherence(args, groups)
        write_tsv(coh, outdir / "structure_coherence.tsv")
        n_coherent = sum(1 for r in coh if r["coherent"])
        print(f"structure_coherence.tsv: {len(coh)} subfamilies on "
              f"{args.pair_metric}, {n_coherent} coherent before any control")
        if confounding and confounding["r"] is not None:
            print(f"  sequence identity explains the score with r = "
                  f"{confounding['r']:.3f} over {confounding['n_pairs']} pairs")
            # An unsupported control must not read as a negative result, so
            # the three outcomes are counted separately (#39).
            for verdict in ("coherent", "not_coherent",
                            "no_interpretation_available"):
                n = sum(1 for r in coh if r["verdict"] == verdict)
                if n:
                    print(f"  identity-adjusted {verdict}: {n}")
            for r in coh:
                if r["verdict"] == "no_interpretation_available":
                    print(f"    {r['subfamily']}: {r['reason']}")
        else:
            print("  " + UNCONTROLLED_WARNING)

    transfer = build_transferability(args)
    if transfer:
        write_tsv(transfer, outdir / "anchor_transferability.tsv")
        print(f"anchor_transferability.tsv: {len(transfer)} rows")
        for r in transfer:
            if not r["label"]:
                continue
            print(f"  {r['label']}: {r['grade']} (support {r['support']}, "
                  f"{r['clade_size']} members) — {r['grade_reason']}")

    if args.focal_subfamily:
        evidence = collect_selection_evidence(args)
        verdict = classify(evidence)
        write_subfunctionalization(args, verdict, evidence, outdir)
        print(f"subfunctionalization.md: {args.focal_subfamily} -> "
              f"{verdict['verdict']}")
        for item in verdict.get("cannot_judge", []):
            print(f"  [cannot judge] {item}")


if __name__ == "__main__":
    main()
