"""fs_transfer.py — foldseek AFDB hits -> function transfer TSV.

No network: UniProt join is skipped or stubbed.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fs_transfer import accession, best_hits, write_transfer

HITS = "\n".join([
    # query target fident alnlen evalue bits qtm ttm alntm qcov tcov taxname
    "GeneA.pdb\tAF-P10490-F1-model_v6\t1.000\t900\t1e-100\t2000\t0.95\t0.94\t0.94\t0.9\t0.9\tMesembryanthemum crystallinum",
    "GeneA.pdb\tAF-Q84VW9-F1-model_v6\t0.850\t900\t1e-90\t1500\t0.93\t0.92\t0.92\t0.9\t0.9\tArabidopsis thaliana",
    "GeneB.pdb\tAF-A0PJE2-2-F1-model_v6\t0.380\t300\t1e-20\t300\t0.80\t0.78\t0.78\t0.5\t0.6\tHomo sapiens",
]) + "\n"


def test_accession_plain_and_isoform():
    assert accession("AF-P10490-F1-model_v6") == "P10490"
    assert accession("AF-A0PJE2-2-F1-model_v6") == "A0PJE2"


def test_accession_rejects_unknown_shape():
    try:
        accession("not-an-afdb-id")
    except ValueError:
        return
    raise AssertionError("must raise on unrecognised target id")


def test_best_hits_keeps_highest_bitscore(tmp_path):
    p = tmp_path / "hits.tsv"
    p.write_text(HITS)
    best = best_hits(p)
    assert set(best) == {"GeneA", "GeneB"}          # .pdb stripped
    assert best["GeneA"]["target"] == "AF-P10490-F1-model_v6"   # 2000 > 1500


def test_write_transfer_emits_header_and_rows(tmp_path):
    p = tmp_path / "hits.tsv"
    p.write_text(HITS)
    out = tmp_path / "t.tsv"
    write_transfer(best_hits(p), {}, out)
    lines = out.read_text().splitlines()
    assert lines[0].split("\t")[:3] == ["gene_id", "best_hit", "uniprot"]
    assert len(lines) == 3
    assert lines[1].split("\t")[2] == "P10490"


def test_write_transfer_joins_annotation_when_present(tmp_path):
    p = tmp_path / "hits.tsv"
    p.write_text(HITS)
    out = tmp_path / "t.tsv"
    anno = {"P10490": {"Entry": "P10490",
                       "Protein names": "Phosphoenolpyruvate carboxylase 1",
                       "EC number": "4.1.1.31",
                       "Gene Names (primary)": "PPC1",
                       "Gene Ontology IDs": "GO:1"}}
    write_transfer(best_hits(p), anno, out)
    row = [l for l in out.read_text().splitlines() if l.startswith("GeneA")][0]
    f = row.split("\t")
    assert f[3] == "Phosphoenolpyruvate carboxylase 1"
    assert f[5] == "4.1.1.31"


def test_output_columns_match_annotation_matrix_reader(tmp_path):
    """annotation_matrix.load_axes reads these exact column names."""
    p = tmp_path / "hits.tsv"
    p.write_text(HITS)
    out = tmp_path / "t.tsv"
    write_transfer(best_hits(p), {}, out)
    from annotation_matrix import load_axes, build_matrix
    m = build_matrix(load_axes(foldseek_tsv=out))
    assert m["GeneA"]["fs_uniprot"] == "P10490"
    assert abs(m["GeneA"]["fs_qtm"] - 0.95) < 1e-9
