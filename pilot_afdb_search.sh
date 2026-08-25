#!/bin/bash
# PLAZA pilot — structural branch: cohort 3Di subDB vs AFDB-SwissProt (gpu).
#
# No new ProstT5 conversion: the cohort is a subset of the 459,398-protein
# panel whose 3Di DB already exists, so the sub-database is carved out of it
# by id (main, _ss, _h with one shared key file). Judgments downstream are
# bits/E-only — the query side has no CA coordinates.
set -euo pipefail

FOLDSEEK=${FOLDSEEK:-$HOME/scratch/bin/foldseek_src/build/src/foldseek}
PANEL_DB=${PANEL_DB:-$HOME/prostt5_15sp/merged/all_3di}
AFDB=${AFDB:-$HOME/foldseek_dbs/afdb_swissprot}
WORK=${WORK:-$HOME/plaza_pilot}
COHORT_TSV=${COHORT_TSV:-$WORK/cohort/cohort.tsv}
THREADS=${THREADS:-14}

mkdir -p "$WORK"/{subdb,search}

# gene ids -> internal keys via the panel lookup
awk 'NR>1 {print $2}' "$COHORT_TSV" | sort -u > "$WORK/subdb/cohort_ids.txt"
awk 'NR==FNR {want[$1]; next} ($2 in want) {print $1}' \
  "$WORK/subdb/cohort_ids.txt" "$PANEL_DB.lookup" > "$WORK/subdb/cohort_keys.txt"
n_ids=$(wc -l < "$WORK/subdb/cohort_ids.txt")
n_keys=$(wc -l < "$WORK/subdb/cohort_keys.txt")
echo "cohort ids: $n_ids -> panel keys: $n_keys"
if [ "$n_keys" -eq 0 ]; then echo "FATAL: no keys resolved" >&2; exit 1; fi

for suffix in "" _ss _h; do
  "$FOLDSEEK" createsubdb "$WORK/subdb/cohort_keys.txt" \
    "$PANEL_DB$suffix" "$WORK/subdb/cohort_3di$suffix" --subdb-mode 1
done

"$FOLDSEEK" search "$WORK/subdb/cohort_3di" "$AFDB" \
  "$WORK/search/res" "$WORK/search/tmp" --threads "$THREADS" -e 1e-5 --max-seqs 10
"$FOLDSEEK" convertalis "$WORK/subdb/cohort_3di" "$AFDB" \
  "$WORK/search/res" "$WORK/search/cohort_vs_afdb.m8" --threads "$THREADS" \
  --format-output query,target,fident,alnlen,qcov,tcov,evalue,bits
wc -l "$WORK/search/cohort_vs_afdb.m8"
echo "PILOT_AFDB_DONE"
