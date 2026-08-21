#!/usr/bin/env python3
"""Foldseek AFDB-SwissProt hits -> per-gene function transfer TSV.

Takes `foldseek easy-search` output against the AlphaFold-SwissProt
database, keeps the best hit per query, and joins the target's UniProt
annotation (protein name, gene symbol, EC, GO) from the UniProt REST API.

Column order must match the `--format-output` used for the search:
  query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,alntmscore,
  qcov,tcov,taxname

Usage:
  python fs_transfer.py --hits hits.tsv -o structure_function_transfer.tsv

The output is what `annotation_matrix.py --foldseek-tsv` consumes.
Measured on the PEPC clan: 99 of 100 best hits carried EC 4.1.1.31, and
the two *M. crystallinum* queries recovered their own species' curated
entries (P10490/P16097) at sequence identity 1.000.
"""

import argparse
import re
import sys
import urllib.request
from typing import Dict, List

COLUMNS = ("query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,"
           "alntmscore,qcov,tcov,taxname").split(",")

# AFDB target ids look like AF-P10490-F1-model_v6, and isoform entries add a
# segment: AF-A0PJE2-2-F1-model_v6. Both must yield the bare accession.
ACCESSION_RE = re.compile(r"AF-([A-Z0-9]+)(?:-\d+)?-F1")

UNIPROT_FIELDS = "accession,protein_name,ec,go_id,gene_primary"
BATCH = 50


def best_hits(hits_path) -> Dict[str, dict]:
    """Highest-bitscore hit per query; '.pdb' stripped from query names."""
    best: Dict[str, dict] = {}
    with open(hits_path) as f:
        for line in f:
            if not line.strip():
                continue
            row = dict(zip(COLUMNS, line.rstrip("\n").split("\t")))
            q = row["query"].replace(".pdb", "")
            if q not in best or float(row["bits"]) > float(best[q]["bits"]):
                best[q] = row
    return best


def accession(target: str) -> str:
    m = ACCESSION_RE.match(target)
    if not m:
        raise ValueError(f"unrecognised AFDB target id: {target}")
    return m.group(1)


def fetch_uniprot(accessions: List[str]) -> Dict[str, dict]:
    """UniProt REST lookup, batched. Network failures degrade to empty."""
    anno: Dict[str, dict] = {}
    for i in range(0, len(accessions), BATCH):
        batch = accessions[i:i + BATCH]
        url = ("https://rest.uniprot.org/uniprotkb/search?query="
               + "+OR+".join(f"accession:{a}" for a in batch)
               + f"&fields={UNIPROT_FIELDS}&format=tsv&size=500")
        try:
            with urllib.request.urlopen(url, timeout=60) as h:
                rows = h.read().decode().splitlines()
        except Exception as e:                     # noqa: BLE001
            print(f"UniProt lookup failed for batch {i//BATCH}: {e}",
                  file=sys.stderr)
            continue
        if not rows:
            continue
        hdr = rows[0].split("\t")
        for row in rows[1:]:
            d = dict(zip(hdr, row.split("\t")))
            anno[d["Entry"]] = d
    return anno


def write_transfer(best: Dict[str, dict], anno: Dict[str, dict], out_path):
    with open(out_path, "w") as out:
        out.write("gene_id\tbest_hit\tuniprot\tprotein_name\tgene_symbol\tec\t"
                  "go_ids\tqtmscore\tfident\tevalue\tbits\ttarget_taxon\n")
        for q in sorted(best):
            r = best[q]
            a = anno.get(accession(r["target"]), {})
            out.write("\t".join([
                q, r["target"], accession(r["target"]),
                a.get("Protein names", ""), a.get("Gene Names (primary)", ""),
                a.get("EC number", ""), a.get("Gene Ontology IDs", "")[:200],
                r["qtmscore"], r["fident"], r["evalue"], r["bits"],
                r["taxname"],
            ]) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--hits", required=True, help="foldseek easy-search TSV")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--no-uniprot", action="store_true",
                    help="Skip the REST join (offline: ids and scores only)")
    args = ap.parse_args()

    best = best_hits(args.hits)
    accs = sorted({accession(r["target"]) for r in best.values()})
    print(f"{len(best)} queries with hits, {len(accs)} unique accessions",
          file=sys.stderr)
    anno = {} if args.no_uniprot else fetch_uniprot(accs)
    write_transfer(best, anno, args.out)
    print(f"-> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
