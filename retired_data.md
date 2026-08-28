# 격리된 데이터 대장 (RETIRED_DO_NOT_USE)

> 2026-08-24. 폐기 판정을 받은 산출물·입력을 **삭제하지 않고 격리**했다.
> 목적: 낡은 데이터 위에서 재분석이 다시 시작되는 것(측정 순환, #47)을 물리적으로 차단한다.
> **이 목록의 어떤 것도 새 분석의 입력으로 쓰지 말 것.** 필요한 건 전부 대체물이 있다.

## pronghorn: `~/scratch/bin/family_finder/RETIRED_DO_NOT_USE/`

| 격리물 | 왜 폐기됐나 | 대체물 |
|---|---|---|
| `data_12sp/` `data_14sp/` `data_17sp/` | **더미 종트리(전 가지 1.0)로 돌았음** — 종 인지 프루닝이 비활성이었다 (resume.md 폐기 결정) | `data_15sp/` + `species_tree_estimated.nwk` |
| `output_12sp_portulacineae/` (+`.pre_final`) | 더미 트리 런. **원고의 "Ppc2 = R1_OG0000168" 오기의 출처** — 실제 PEPC는 R1_OG0000165였다 (#43) | `output_15sp_v2/` |
| `output_15sp/` (2026-05 v1) | v1 파이프라인 런. 종트리 마커 추출용 1패스로 재사용된 뒤 소임 종료 (#41) | `output_15sp_v2/` |
| `wd_I1.05`–`wd_I1.3`, `wd15_I1.*`, `wd_noCgigH` | MCL inflation 스윕 하드링크 사본 — **inflation은 죽은 노브로 판정** (#10) | 불필요 |
| `pepc_pilot_seltest/{seltest,seltest_fast}` | **폐기 branch-site 정본** — LRT 15.3631이 나온 곳. `_M` 오염 서열 집합(102-taxon) 위 계산 (#40, #44). 사본은 `~/scratch/pepc_branchsite_20260823/retired_2026-08-20/`에도 있음 | `~/scratch/pepc_branchsite_20260823/` (교정 세트 재생성, 진행 중) |
| `codon95.codon.aln.NOGAP_599codons` | `pal2nal -nogap`이 **코돈 58% 삭제**(1,428→599). 선택압·BEB 좌표에 치명적 (#44, resume.md §5) | `pepc_branchsite_20260823/inputs/`의 1,428코돈 재정렬 |

| `fragmentation_clusters.tsv` `vote_edges.tsv` (2026-08-24 초판) | **vote-inflation 버그 산출물** — domtblout의 프로파일 단위 정렬을 유전자 단위 연속으로 오인해 유전자당 다중 투표. 간선 1.02M(상한 23,744), 거대 컴포넌트 13,347 families 전부 인공물 (resume.md §5) | `merge_analyze2` 재계산본 (`output_15sp_v2/`) |

## gpu: `~/RETIRED_DO_NOT_USE/`

| 격리물 | 왜 | 대체물 |
|---|---|---|
| `ppc_resolve/` (구 `~/pepc_pilot/ppc_resolve`) | `_M` 오염 102-taxon 세트의 트리들(T8/T10/T14) — #32 시대 산출물 | `pronghorn:~/scratch/mM_repair/clan_corrected.*` + `clan_corrected_tree.*` |

## 격리하지 않은 것 (혼동 방지)

| 유지물 | 이유 |
|---|---|
| `data/` (5종 MAKER Cgig 포함) | **의도적 유지** — 한 유전체의 두 주석이 §2.X.8 길이-정책 대조군의 유일한 장치 (#43) |
| `data_15sp/species_tree.nwk` (더미) | provenance. 실행은 `species_tree_estimated.nwk`가 담당 |
| `~/scratch/mM_repair/` 나머지 | **PEPC 정본** (`clan_corrected.faa` 108서열) |
| gpu `~/pepc_pilot/` 나머지 | 파일럿 정본 (signal_windows, allvall.tsv, DeepLoc 등) |
| 로컬 `/mnt/f/.../opuntia_analysis/family_finder/output_15sp/` | 다른 프로젝트 디렉터리라 이동하지 않음. **사용 금지** — v2로 대체됨 |

## ⚠️ 사건 기록 (2026-08-24): 격리가 15sp 입력을 끊었다

`data_15sp/{pep,cds}` 24개가 `data_12sp` 실파일로의 심볼릭 링크였다 — 12sp **패널**(더미 트리 런 구성)은
폐기 대상이 맞지만 그 **FASTA 파일들은 15sp가 공유하는 살아있는 입력**이었다. 격리 직후 깨진 링크 검사를
하지 않아 다음 분석이 실패할 때야 드러났다. 수리: 24개 전부 실파일 복사로 교체(자급화), pep 총 484,752
서열 재검증. **교훈: 격리/이동 직후 `find <활성 트리> -xtype l` 필수.**

## 규칙

- 격리물이 다시 필요하면(재현 검증 등) **읽기만** 하고, 결과에 인용할 때는 "폐기 산출물 재파싱"임을 명시한다
- 새 폐기가 생기면 이 파일과 현장 `README.md`에 **같이** 기록한다 — 한 곳에만 적힌 사실이 이 프로젝트 오류의 공통 패턴이었다

## 2026-08-28 추가 — 세대 구획(partition) 완결

**무효 이관 2건**: `output_14sp/`, `output_17sp/` (+`output_14sp_launch.log`) →
`RETIRED_DO_NOT_USE/`. 사유는 data_14sp/17sp와 동일 — **더미 종트리(전 가지 1.0) 런**.
data만 격리되고 output이 활성 트리에 남아 있던 잔무. 이관 직후 `find -L` 깨진 심링크 0 확인.

**상태 마커 체계 도입**: 활성 트리의 모든 분석 세대 디렉터리에 `_CANON_STATUS.txt`
(CANON / SUPERSEDED / MIXED + 대체물 포인터) 배치 — pronghorn 11곳
(output_15sp_v2·output_5sp_v2 = CANON, output_5sp·output_11sp* = SUPERSEDED,
arbiter_v3·subfam_trees·hog_pilot·vote_rescan·mM_repair·mcry_direct = CANON),
gpu 4곳 (subfam·prostt5_15sp·annot_panel = CANON, **pepc_pilot = MIXED** — canonical
파일과 _M 오염 시대 파일이 공존, 개별 파일은 resume.md 대조).
최상위 규칙 파일: `pronghorn:~/scratch/bin/family_finder/ANALYSIS_PARTITIONS.txt`.

**구획 3계층**: ① `RETIRED_DO_NOT_USE/` = 무효(입력 금지, 이 대장이 사유 기록)
② `SUPERSEDED` 마커 = 구세대(비교 기준선으로만, 새 분석 입력 금지)
③ `CANON` 마커 = 현행 정본. 새 분석은 CANON만 소비하고, 산출 디렉터리에 자기 마커를 만든다.
