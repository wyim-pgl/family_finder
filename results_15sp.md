# 15종 최종 family 결과 — **v3** (2026-08-24)

> **정본 테이블: `output_15sp_v2/summary_v3.tsv` — 20,133 families / 459,398 유전자.**
> v2의 "파편화 주석"이 이 판에서 **트리 검증된 병합**으로 확정됐다: 5,618개 클러스터 전수에
> 코돈 트리를 세워, 조각이 실제로 뒤섞인 경우만 병합했다. 남은 미결은 codeml 선택압 하나(맨 아래).
> v2 테이블(`summary.tsv`, 23,744)은 병합 전 상태로 보존된다.

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

## ✅ 확정 — family 경계 (트리 검증 완료)

5,618개 파편화 클러스터 **전수**에 대해 정렬(MAFFT)→코돈 사영(pal2nal)→트리(FastTree)를 세우고,
`steps/cluster_validate.py`의 규칙으로 판정했다 (**5,618/5,618, 실패 0**):

- **INTERLEAVED** (조각 9,228개) — 트리가 조각을 뒤섞는다 = 계통 축 파편화 = 한 family.
  **병합군 2,410개**(family 6,021개)로 묶여 v3에 반영: **23,744 → 20,133 families (−3,611)**
- **MONOPHYLETIC** (조각 6,833개) — 자기 clade를 이룬다. 같은 family의 서브패밀리인지 이웃
  family인지 **위상만으로는 못 정하므로 병합하지 않고** 그대로 둔다 (PPC4가 그 사례 —
  PEPC 클러스터 안에서 단계통이고 5sp에서 별도 family였다)

무결성: 유전자 459,398 전량 보존 · 중복 0 · family 수 = 23,744 − Σ(병합군−1) 단언 통과.
`merged_from`/`cluster_id` 컬럼이 병합 provenance를 담는다.

**known-answer**: PEPC 계통 조각 3개(0000440+0013449+0019668)가 84유전자 한 family로 병합,
1E2(0008467)·PPC4(0009826)는 미판정 유지 — 파일럿과 전수 결과 일치.

PEPC 서브패밀리 정본은 여전히 `mM_repair/clan_corrected`(1E1 62 [95/91] · 1E2 15 [100/76],
둘 다 PROVISIONAL). **서브패밀리는 v3에서도 family 테이블이 아니라 clan 트리에서 나온다.**

## 미결 — 하나

**codeml branch-site** (교정 세트 1,428코돈, SLURM 6097933-38): 2/8 lnL, 진행 중. 완료 시
methods.md에 반영한다. 감도 축(599코돈 -nogap)은 LRT 0.33/p 1.0이었으나 코돈 58%가 없는
정렬이라 정본이 아니다.

## 산출물 인덱스 (pronghorn)

```
output_15sp_v2/summary_v3.tsv                  정본 family 테이블 (v3, 병합 적용)
output_15sp_v2/summary.tsv                     v2 (병합 전, 보존)
output_15sp_v2/cluster_verdicts.tsv            클러스터 트리 판정 전수 (5,618)
output_15sp_v2/merge_candidates_cov{0.5,0.2}.tsv   1층: 상호 쌍
output_15sp_v2/fragmentation_clusters.tsv      2층: 이웃 클러스터 (C0297 = PEPC)
output_15sp_v2/vote_edges.tsv                  간선 원자료
~/scratch/rescue_dom_15sp/final_assign.json    확정 rescue 배정
RETIRED_DO_NOT_USE/                            폐기물 (대장: retired_data.md)
```
