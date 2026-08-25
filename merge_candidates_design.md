# `merge-candidate` design

This design lives at repo root because there is no `docs/` directory in this checkout (`rg --files steps tests docs` reported `docs: No such file or directory` during verification), so a root-level markdown file is the least surprising standalone-analysis document.

## Repo verification

The requested pre-design checks were verified against the repo, not assumed:

- `CLAUDE.md` confirms the architecture, per-round resume artifact, and assignment-tier framing: the main control flow is `family_finder.py -> pipeline.run()`, resume reads `round_NN/confirmed_families.tsv`, and the planned fourth tier is structural via ESM/Foldseek ([`CLAUDE.md:38-65`](CLAUDE.md), [`CLAUDE.md:125-126`](CLAUDE.md)).
- The prompt's summary of `resume.md` is **not** what the current file says. `resume.md` now says its older "confirmed cause is MCL/graph fragmentation" conclusion was overturned for the residual unplaced pool on 2026-08-22, warns not to start graph intervention from that assumption, and attributes most of the Mcry gap to annotation quality plus taxon sampling instead ([`resume.md:17-21`](resume.md), [`resume.md:32-33`](resume.md), [`resume.md:90-99`](resume.md)). `CLAUDE.md` still carries the older summary in its investigation-record section ([`CLAUDE.md:148-155`](CLAUDE.md)).
- Commit `2852ecad9f165bf813520a81b74b83b72850a4a9` is real and its message states the measured finding that structure cannot draw the PTPC/BTPC family boundary. The same conclusion is recorded in `results_15sp.md`: TM-score separation is tiny, ProstT5 3Di within-vs-cross separation is only 1.11, and structure is explicitly framed as unsuitable for drawing the family boundary by itself ([`git show 2852eca`](.git), [`results_15sp.md:495-544`](results_15sp.md)).
- Existing step conventions match a thin-wrapper/pure-logic split. `steps/esm.py` states that external tool calls are module-level mockable wrappers and decision logic is pure Python ([`steps/esm.py:8-21`](steps/esm.py), [`steps/esm.py:56-137`](steps/esm.py)); `steps/prostt5_chunks.py` and `steps/hmm_chunks.py` use the same pattern for foldseek/HMMER subprocess wrappers plus pure parsers/mergers ([`steps/prostt5_chunks.py:123-170`](steps/prostt5_chunks.py), [`steps/hmm_chunks.py:1-24`](steps/hmm_chunks.py), [`steps/hmm_chunks.py:156-240`](steps/hmm_chunks.py)).
- The required 3Di validator already exists and hard-fails malformed/incomplete DBs: `verify_3di_db()` checks for `_ss`, `_ss.index`, empty 3Di, byte-identical amino-acid/3Di payloads, and index-count mismatches ([`steps/prostt5_chunks.py:123-160`](steps/prostt5_chunks.py)).
- The round artifact this step must consume is exactly `family_id<TAB>comma-joined gene ids` with no header, written by `_write_round_families()` in `pipeline.py`; `_load_round_families()` enforces uniqueness and round ownership when it re-reads them ([`pipeline.py:309-318`](pipeline.py), [`pipeline.py:321-427`](pipeline.py)).
- `Config` field style is additive dataclass fields with short trailing comments or short measured comment blocks near the relevant section. Structural knobs already live near `foldseek_bin` / tier-3 settings ([`config.py:49-60`](config.py), [`config.py:198-211`](config.py)).
- Test conventions are exactly as requested: no real external tools, module-level monkeypatching of wrappers, pure functions fed synthetic fixtures, and local stubbing for module-level imports like `ete4` when needed ([`CLAUDE.md:139-146`](CLAUDE.md), [`tests/test_size_gate.py:20-27`](tests/test_size_gate.py), [`tests/test_esm_modules.py:1-5`](tests/test_esm_modules.py)).
- The "nominate only, never merge automatically" rule is already local repo policy for fragmentation scans: `tests/test_merge_scan.py` says candidates still need tree validation and nothing is merged automatically ([`tests/test_merge_scan.py:14-17`](tests/test_merge_scan.py), [`tests/test_merge_scan.py:61-67`](tests/test_merge_scan.py)). README also says EPA-ng placement is the later ambiguity arbiter and annotation layers are never membership filters ([`README.md:23-27`](README.md)).

## Scope

`merge-candidate` is a standalone structural nomination step. It does **not** change `pipeline.py` and does **not** declare any family merge. It answers only:

> "Which structurally co-clustered sets contain genes from multiple confirmed OrthoFinder families in this round, plus which unplaced genes sit with them?"

The later arbiter remains the planned EPA-ng + external-anchor phylogenetic judge already reflected in repo docs/config (`README.md:23-27`, `config.py:213-235`).

## CLI / step interface

Implemented as `python -m steps.merge_candidates`:

```bash
python -m steps.merge_candidates \
  --run-dir output_15sp_v2 \
  --db /path/to/all_proteins_prostt5_db \
  --outdir output_15sp_v2/round_01/merge_candidates \
  --config config_15sp_v2.json \
  --threads 32 \
  --verbose
```

Inputs:

- `--run-dir`: the pipeline output directory; the family table is the UNION of every `round_*/confirmed_families.tsv` (each file is that round's newly confirmed families only — review finding, 2026-08-25).
- `--db`: existing ProstT5 3Di database prefix, built elsewhere by `steps/prostt5_chunks.py`.
- `--outdir`: output directory for cluster TSVs and summaries.
- `--config`: optional JSON config; omitted means `Config()`.
- `--threads`: optional override for `config.n_workers`.

Outputs:

- `structural_merge_candidates.tsv`: candidate clusters only.
- `structural_merge_overflow.tsv`: oversized foldseek clusters intentionally excluded from nomination.
- `structural_merge_summary.tsv`: aggregate counts and the verified 3Di DB entry count.
- `foldseek_clusters.tsv`: raw `createtsv` rep/member pairs kept for provenance/debugging.

## Exact `foldseek` commands

The step runs exactly two foldseek commands after validating the DB:

```bash
foldseek cluster <db> <outdir>/foldseek_clusters <outdir>/foldseek_tmp \
  --threads <config.n_workers> \
  -e <config.merge_candidate_evalue> \
  -c <config.merge_candidate_min_coverage> \
  --min-seq-id <config.merge_candidate_min_seq_id>

foldseek createtsv <db> <db> <outdir>/foldseek_clusters <outdir>/foldseek_clusters.tsv
```

Chosen parameters and why:

- `--threads`: follows existing config practice; no new tool-path field is needed because `foldseek_bin` already exists.
- `-e merge_candidate_evalue` (default `1e-5`): keeps weak structural chaining out of the nomination set.
- `-c merge_candidate_min_coverage` (default `0.8`): requires broad structural agreement instead of short local overlap, because PTPC/BTPC already showed that global structural similarity over-merges.
- `--min-seq-id merge_candidate_min_seq_id` (default `0.3`): adds another conservative floor so a cluster cannot expand through extremely weak 3Di/state similarity alone.

Runaway-cluster guard:

- Primary guard is `merge_candidate_max_cluster_size` (default `150`).
- Any foldseek cluster with more members than the cap is written to `structural_merge_overflow.tsv` and **excluded** from `structural_merge_candidates.tsv`.
- Rationale: the measured PEPC PTPC family that eventually needed adjudication is 111 genes, so `150` leaves room for that scale while still treating larger greedy set-cover expansions as suspicious overflow instead of evidence.

This guard is mandatory because commit `2852eca` and `results_15sp.md` show that structural similarity alone over-merges PTPC/BTPC; the safe action here is to nominate conservatively and skip suspiciously large chains, not to redraw a family boundary.

## Mapping algorithm

Resolution from repo evidence:

- `confirmed_families.tsv` is the authoritative per-round mapping, but it is only two columns and has no header (`pipeline.py:309-318`).
- Therefore this step builds its own `gene_id -> family_id` index from that file.

Algorithm:

1. Call `verify_3di_db(db)` from `steps.prostt5_chunks.py`. If validation fails, abort before any foldseek run.
2. Run `foldseek cluster` on the existing 3Di DB.
3. Run `foldseek createtsv` to obtain rep/member cluster rows.
4. Parse only the first two TSV columns as `(representative_gene_id, member_gene_id)`.
5. Group rows by representative. Reject malformed output if the same member appears under two different representatives.
6. Cross-reference every member against `gene_to_family` built from the union of all `round_*/confirmed_families.tsv`.
7. For each structural cluster, compute:
   - `family_gene_counts`: how many members belong to each confirmed family.
   - `unplaced_gene_ids`: genes absent from `confirmed_families.tsv`.
   - `species_counts`: counts by the prefix before the first underscore.
8. Assign a status:
   - `candidate`: cluster spans at least two confirmed families and is within the size cap.
   - `overflow`: cluster exceeds `merge_candidate_max_cluster_size`; recorded separately, never nominated.
   - `single_family_with_unplaced`: one confirmed family plus unplaced genes; informative but not a family-merge candidate.
   - `single_family`: only one confirmed family.
   - `unplaced_only`: no confirmed-family members.

Unplaced/singleton handling:

- "Unplaced" means "present in the structural cluster but absent from every round's `confirmed_families.tsv`".
- Unplaced genes are written explicitly in candidate rows when they co-cluster with multi-family candidates.
- Clusters that are only one family plus unplaced genes, or unplaced-only, are counted in the summary but are **not** reported as merge candidates.

## Output schema

### `structural_merge_candidates.tsv`

One row per structural cluster with `status == candidate`:

- `cluster_id`: synthetic stable ID `FSCLUST_000001`, `FSCLUST_000002`, ...
- `representative_gene_id`
- `n_genes`
- `n_species`
- `n_families`
- `family_ids`: comma-joined confirmed family IDs
- `family_gene_counts`: comma-joined `family_id:n_members_in_cluster`
- `n_unplaced`
- `unplaced_gene_ids`: comma-joined gene IDs, or `-`
- `species_counts`: comma-joined `SpeciesPrefix:n_members`

### `structural_merge_overflow.tsv`

Same idea, but only for `status == overflow`. This is the audit trail for skipped runaway clusters.

### `structural_merge_summary.tsv`

Key/value TSV with:

- `verified_3di_entries`
- `candidate_clusters`
- `candidate_genes`
- `candidate_unplaced_genes`
- `overflow_clusters`
- `single_family_clusters`
- `single_family_with_unplaced_clusters`
- `unplaced_only_clusters`
- `total_clusters`
- `max_cluster_size_seen`

## New `Config` fields

Added next to the existing structural tier knobs in `config.py`:

- `merge_candidate_evalue: float = 1e-5`
  Structural cluster significance floor (`foldseek cluster -e`).
- `merge_candidate_min_coverage: float = 0.8`
  Coverage floor for structural clustering (`foldseek cluster -c`).
- `merge_candidate_min_seq_id: float = 0.3`
  Additional similarity floor to resist weak chaining.
- `merge_candidate_max_cluster_size: int = 150`
  Hard overflow guard: larger clusters are recorded but not nominated.

## Test plan

Unit tests are synthetic and do not invoke foldseek:

1. Parse `confirmed_families.tsv` exactly as `pipeline.py` writes it.
2. Parse synthetic foldseek cluster TSV fixtures using only the first two columns.
3. Reject malformed cluster fixtures where one member appears under two representatives.
4. Cross-reference logic:
   - cluster spanning two confirmed families becomes `candidate`;
   - family counts are correct;
   - unplaced genes are retained in `unplaced_gene_ids`.
5. Oversized-cluster guard:
   - `n_genes > merge_candidate_max_cluster_size` becomes `overflow`;
   - overflow rows are written to the separate TSV;
   - such clusters never appear in `structural_merge_candidates.tsv`.
6. Species-prefix handling:
   - counts come from the prefix before the first underscore;
   - candidate rows carry `species_counts`.
7. End-to-end module test:
   - monkeypatch `verify_3di_db`, `run_foldseek_cluster`, and `run_foldseek_createtsv`;
   - write synthetic cluster TSV output;
   - assert that candidate, overflow, and summary files are emitted with the expected rows.

## Ambiguities resolved from the repo

- I did **not** add pipeline wiring because the prompt explicitly forbids it and the repo already distinguishes standalone analyses from the main round loop (`CLAUDE.md:11-13`, `find_pseudogenes.py` CLI).
- I did **not** consume `summary.tsv`; the round artifact is the verified source for this step, and `pipeline.py` explicitly treats `summary.tsv` as end-of-run output while resume rebuilds from per-round `confirmed_families.tsv` (`CLAUDE.md:45-53`, `pipeline.py:321-427`).
- I treated structure as nomination-only because the repo's own measured PTPC/BTPC result says structure cannot draw the family boundary (`results_15sp.md:495-544`) and because the existing merge-scan tests already encode "candidate only, never merge automatically" (`tests/test_merge_scan.py:14-17`, `tests/test_merge_scan.py:61-67`).
