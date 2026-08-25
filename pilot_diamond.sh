#!/bin/bash
# PLAZA pilot — DIAMOND forward + per-species reverse searches (pronghorn).
#
# Reverse searches run against 15 SEPARATE per-species DBs, not one 459K
# concatenation: RBH across a 459K-vs-27K asymmetry is only meaningful
# species-scoped (design review 2026-08-25). No identity/coverage filters at
# search time — the QC step calibrates gates on the retained hits.
set -euo pipefail

DIAMOND=${DIAMOND:-/data/gpfs/assoc/pgl/bin/conda/conda_envs/orthofinder/bin/diamond}
ATH_FA=${ATH_FA:-/data/gpfs/assoc/pgl/data/sequence_data/plaza/protein/selected/proteome.selected_transcript.ath.fasta.gz}
PEP_DIR=${PEP_DIR:-$HOME/scratch/bin/family_finder/data_15sp/pep}
WORK=${WORK:-$HOME/scratch/plaza_pilot}
COHORT_FA=${COHORT_FA:-$WORK/cohort/cohort.fa}
THREADS=${THREADS:-8}

FLAGS=(--very-sensitive --evalue 1e-10 --max-target-seqs 25 --max-hsps 1
       --outfmt 6 qseqid sseqid pident length qlen slen qstart qend
       sstart send evalue bitscore qcovhsp scovhsp)

mkdir -p "$WORK"/{db,forward,reverse}

if [ ! -s "$WORK/db/ath.dmnd" ]; then
  zcat -f "$ATH_FA" | "$DIAMOND" makedb --db "$WORK/db/ath" --threads "$THREADS"
fi

echo "== forward: cohort vs ATH =="
"$DIAMOND" blastp --db "$WORK/db/ath" --query "$COHORT_FA" \
  --out "$WORK/forward/cohort_vs_ath.tsv" --threads "$THREADS" "${FLAGS[@]}"

echo "== reverse: ATH vs each species =="
for fa in "$PEP_DIR"/*.fa; do
  sp=$(basename "$fa" | cut -d. -f1)
  db="$WORK/db/$sp"
  [ -s "$db.dmnd" ] || "$DIAMOND" makedb --in "$fa" --db "$db" --threads "$THREADS"
  out="$WORK/reverse/$sp.tsv"
  if [ ! -s "$out.done" ]; then
    zcat -f "$ATH_FA" | "$DIAMOND" blastp --db "$db" --query - \
      --out "$out" --threads "$THREADS" "${FLAGS[@]}"
    date > "$out.done"
  fi
done
echo "PILOT_DIAMOND_DONE"
