# 15종 v2 최종 family 결과 — v1 (2026-08-23)

> **정본 테이블**: `pronghorn:~/scratch/bin/family_finder/output_15sp_v2/summary.tsv`
> (23,744 families / 459,398 배정 유전자, 12.7 MB)
> 이 문서는 그 테이블을 **어떻게 읽어야 하는지** — 무엇이 인용 가능하고 무엇이 잠정인지 — 를 정의한다.
> v2 갱신 조건은 맨 아래에 있다.

## 런 정보

| | |
|---|---|
| SLURM | 6097190, **COMPLETED**, 1-12:43:29 |
| config | `config_15sp_v2.json` (`profile_assign_per_round: true`, chunked rescue 500×48) |
| 종트리 | `data_15sp/species_tree_estimated.nwk` — 135 단일카피 loci / 250,329 nt / 전 내부노드 100/100 |
| Cgig | **Helixer** (27,583) — 5종 패널(MAKER 29,163)과 비교 불가 |
| 수렴 | R1 22,051 → R2 1,448 → R3 212 → R4 32 → R5 1 → R6–7 0 |

## 무결성 — 전부 통과

```
배정 459,398 + 미배치 25,354 = 484,752 = 입력 총계   (유실 0)
중복 배정 0  ·  rescue 청크 48/48 (# [ok] 전수 확인)
```

## 인용 가능한 배치율

**comparable = ≥100 aa AND 완전 CDS(시작·종결·3배수·내부종결 없음). 종간 비교는 이 열만 쓴다.**

| 종 | comparable | raw | 제외 (short / invalid CDS) |
|---|---|---|---|
| Ococ | **0.0%** | 0.0% | 806 / 271 |
| Cjam | 0.3% | 1.6% | 2,214 / 0 |
| Ccac | 0.7% | 11.2% | 1,104 / 2,845 |
| Obas | 0.8% | 2.8% | 2,383 / 36 |
| Ahyp | 0.9% | 4.0% | 2,634 / 2,547 |
| Tfru | 1.0% | 12.8% | 10,920 / 9 |
| Dcar | 1.6% | 2.5% | 2,358 / 880 |
| Sole | 1.6% | 4.1% | 3,408 / 9 |
| Cgig | 1.7% | 5.9% | 3,132 / 117 |
| Pole | 1.7% | 3.6% | 2,459 / 0 |
| Sund | 1.9% | 3.0% | 2,027 / 1,213 |
| Bvul | 2.7% | 5.3% | 2,012 / 325 |
| Sof | 2.7% | 3.1% | 2,566 / 285 |
| Mcry | **2.8%** | 10.8% | 3,374 / 2,237 |
| Pami | 5.2% | 7.2% | 3,621 / 385 |

- headline: **Mcry raw 10.8% vs v1 baseline 14.9% — PASS**
- raw 상위 3종(Tfru 12.8 / Ccac 11.2 / Mcry 10.8)은 **전부 주석 정책 산물** — comparable에서 1.0 / 0.7 / 2.8
- comparable 기준 Mcry는 Bvul·Sof와 동급, Pami보다 낮다 → **"Mcry 문제"는 종결**

## ⚠️ 잠정 — v2에서 바뀔 수 있는 것

**1. 대형 family의 경계.** CAM known-answer FAIL (`TOGETHER 0 / SPLIT 2`). PEPC clan이
**최소 6개 family, 128유전자로 파편화**: `R1_OG0000440`(63) · `R1_OG0008467`(18) ·
`R1_OG0009826`(17) · `R1_OG0013449`(15) · `R2_OG0000359`(9) · `R1_OG0019668`(6).
분할은 서브패밀리가 아니라 계통을 따른다(MCL 그래프 파편화). 전수 병합 스캔(6098723) 진행 중.
**clan 수준 분석은 family 경계를 믿지 말고 손으로 모을 것.**

**2. PEPC의 정본은 파이프라인 산출물이 아니다.** `pronghorn:~/scratch/mM_repair/clan_corrected.faa`
(108서열, `_M` 교정 완료)가 정본이다. 서브패밀리: **ppc-1E1** 62멤버(SH-aLRT/UFboot 95/91) ·
**ppc-1E2** 15멤버(100/76) — 등급 기준으로 **둘 다 PROVISIONAL**(UFboot < 95). 선인장 ppc-1E2
전장 사본은 *O. basilaris* 1개, Ococ·Cjam은 절단 좌위(synteny·재구성 근거).

**3. rescue 층의 10,381 유전자는 잠정 배정이다.** 이 배정은 옛 E-value 판정으로 나왔다.
bits 재판정으로 37건이 다른 family를 가리키고, tier-2 기본 커버리지 게이트를 적용하면
**325건만 통과(97% 거부)** 한다 — 게이트가 옳은지(잔여 풀의 69%가 short/broken 모델) 임계값이
안 맞는지(대형 발산 family의 HMM 길이 부풀림) 보정 전이다(#47). rescue 층 유전자 목록은
`output_15sp_v2/hmmer_rescue/rescue_summary.tsv`.

## v2 갱신 조건 (이 순서, 새 측정 스레드 없이)

1. 병합 스캔(6098723) 완료 → 후보 쌍 목록 → 상위 후보 트리 검증 → **family 경계 v2**
2. rescue 커버리지 임계값 보정(1건 확인: 거부가 short/broken 모델에 집중되는가) → **rescue 층 확정**
3. codeml 주 분석(1,428코돈) 완료 → SF3 선택압 재생성 값 → methods.md

작업 추적: 이슈 #47. 이 문서가 갱신되면 v2로 개판한다.
