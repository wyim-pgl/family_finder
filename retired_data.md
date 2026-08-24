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

## 규칙

- 격리물이 다시 필요하면(재현 검증 등) **읽기만** 하고, 결과에 인용할 때는 "폐기 산출물 재파싱"임을 명시한다
- 새 폐기가 생기면 이 파일과 현장 `README.md`에 **같이** 기록한다 — 한 곳에만 적힌 사실이 이 프로젝트 오류의 공통 패턴이었다
