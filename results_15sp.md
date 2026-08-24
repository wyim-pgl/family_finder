# 15종 v2 최종 family 결과 — **v2** (2026-08-24)

> **정본 테이블**: `pronghorn:~/scratch/bin/family_finder/output_15sp_v2/summary.tsv`
> (23,744 families / 459,398 배정 유전자)
> v1(`cb2c973`)의 잠정 항목 두 개 — rescue 층, family 경계 — 가 이 판에서 확정/주석화됐다.
> 남은 미결은 codeml 선택압 하나뿐이다(맨 아래).

## 런 정보

| | |
|---|---|
| SLURM | 6097190, **COMPLETED**, 1-12:43:29 |
| config | `config_15sp_v2.json` · 추정 종트리(135 loci, 전 내부노드 100/100) · Cgig=Helixer |
| 수렴 | R1 22,051 → R6–7 0 |
| 무결성 | 유실 0 · 중복 0 · rescue 청크 48/48 |

## 인용 가능 — 배치율

**comparable(≥100 aa + 완전 CDS)만 종간 비교에 쓴다.** headline: Mcry raw 10.8% vs baseline 14.9% **PASS**.

| 종 | comparable | raw | | 종 | comparable | raw |
|---|---|---|---|---|---|---|
| Ococ | **0.0%** | 0.0% | | Pole | 1.7% | 3.6% |
| Cjam | 0.3% | 1.6% | | Cgig | 1.7% | 5.9% |
| Ccac | 0.7% | 11.2% | | Sund | 1.9% | 3.0% |
| Obas | 0.8% | 2.8% | | Bvul | 2.7% | 5.3% |
| Ahyp | 0.9% | 4.0% | | Sof | 2.7% | 3.1% |
| Tfru | 1.0% | 12.8% | | **Mcry** | **2.8%** | 10.8% |
| Dcar | 1.6% | 2.5% | | Pami | 5.2% | 7.2% |
| Sole | 1.6% | 4.1% | | | | |

raw 상위 3종(Tfru·Ccac·Mcry)은 전부 주석 정책 산물. **"Mcry 문제"는 종결** — comparable에서 Bvul·Sof 동급.

## ✅ 확정 — rescue 층 (v1에서 잠정이었던 것)

게이트 보정 판정(#47): 프로파일 커버리지 게이트는 rescue 맥락에서 **판별력 0**
(comparable 거부 96% vs 깨진 모델 99%; 건강한 전장 멤버의 profile_cov 중앙값 0.22 — 프로파일
길이가 갭으로 부풀기 때문). `rescue_min_profile_coverage=0.0` 기본(off), 쿼리 게이트 유지(`91dc96c`).

```
최종 rescue: HIGH 8,159 / PROVISIONAL 1,139 / UNRESOLVED 1,083  ->  배정 9,298
             (옛 E-value 판정 10,381 대비: 게이트 탈락 ~1,083, family 변경 28)
배정 목록: ~/scratch/rescue_dom_15sp/final_assign.json
```

## ✅ 주석화 — family 경계의 파편화 (v1의 두 번째 잠정)

전수 병합 스캔(459,398 유전자 × 23,744 프로파일, SLURM 6098723 + 6099002) 산출물 3종이
`output_15sp_v2/`에 있다. **자동 병합은 하지 않았다** — 전부 트리 검증 전 후보다.

**1층 — 상호 쌍 (고신뢰 하한):** `merge_candidates_cov0.5.tsv` — **3,614쌍**
(min-reciprocal ≥0.6 양방향). 전부 검증·병합되면 23,744 → 20,130. 구조상 family당 최대 1쌍.

**2층 — 단방향 이웃 클러스터 (상한, 병합 권고 아님):** `fragmentation_clusters.tsv` —
frac ≥0.6 단방향 간선 14,842개의 연결 성분 **5,618개, family 16,061개(68%) / 유전자 320,721개** 연루.
최대 47 families/705 유전자. ⚠️ 단방향 간선은 "최근접 이웃 family에 표가 몰림"이라 **구성상
과병합**한다 — 진짜 파편화의 상한이지 병합 목록이 아니다. 전부 병합하면 13,301이 되지만
그 숫자를 인용하지 말 것.

**known-answer:** PEPC 6조각 중 **5개가 C0297 하나로 정확히**(119유전자, 외부 혼입 0) 묶인다.
쌍 규칙(1층)은 PEPC를 **하나도 못 잡았다** — 표가 사슬(0.94→0.94→1.00)을 이뤄 어느 쌍도
상호성을 못 채운다. 6번째 조각 `R2_OG0000359`(Mcry Ppc1 splinter)는 9멤버 중 4개만 유효
표를 내 어느 층에도 안 걸린다 — **검출 한계는 splinter의 작고 잡음 많은 쪽에 있다.**

**결론: clan 수준 분석은 `fragmentation_clusters.tsv`의 클러스터 단위로 서열을 모아 시작하고,
family 경계 자체는 트리 검증 없이 병합하지 않는다.** PEPC 정본은 여전히
`mM_repair/clan_corrected.faa`(108서열; ppc-1E1 62 [95/91] · ppc-1E2 15 [100/76], 둘 다 PROVISIONAL).

## 미결 — 하나

**codeml branch-site** (교정 세트 1,428코돈, SLURM 6097933-38): 2/8 lnL, 진행 중. 완료 시
methods.md에 반영한다. 감도 축(599코돈 -nogap)은 LRT 0.33/p 1.0이었으나 코돈 58%가 없는
정렬이라 정본이 아니다.

## 산출물 인덱스 (pronghorn)

```
output_15sp_v2/summary.tsv                     정본 family 테이블
output_15sp_v2/merge_candidates_cov{0.5,0.2}.tsv   1층: 상호 쌍
output_15sp_v2/fragmentation_clusters.tsv      2층: 이웃 클러스터 (C0297 = PEPC)
output_15sp_v2/vote_edges.tsv                  간선 원자료
~/scratch/rescue_dom_15sp/final_assign.json    확정 rescue 배정
RETIRED_DO_NOT_USE/                            폐기물 (대장: retired_data.md)
```
