# resume.md — Mcry(외군) 클러스터링 문제 조사 기록

> 📌 **정본**: §7 (재개 가이드, 2026-08-27) · §1.5 (Mcry 원인 판정) · §1.6 (#40 Ppc1/Ppc2) ·
> §5 (함정 목록 — 누적 유효). 그 밖의 절은 이력 보존용일 수 있다 — 각 절 머리의
> ❌/⚠️/✏️ 격리 마커를 확인할 것. 표기 규약은 위키 `guide/handoff-hygiene.md`.

작성일: 2026-08-19. 이 문서는 코드 조사 + 독립 검증(스트레스 테스트)을 거친 현재 상태 스냅샷이다.
모든 `file:line` 은 커밋 `f1226aa` 기준으로 원문 대조 확인했다.

## 1. 문제 정의

> ✏️ **PARTIAL (2026-08-22):** "원인은 클러스터링 상류" 추정은 §1.5에서 판정 완료 —
> 입력 주석 품질 약 2/3 + taxon sampling 약 1.7배, 파이프라인 결함 아님.
> 문제 정의 자체는 유효하나 아래 "핵심 결론" 문단은 인용 금지.

5종 런(Mcry, Cgig, CgigH, Ococ, Obas; 143,961 단백질)에서 외군 *Mesembryanthemum
crystallinum*(Mcry, Aizoaceae) 유전자가 **Round 1 OrthoFinder 클러스터링 단계에서부터**
진짜 cactus ortholog와 다른 orthogroup으로 배정된다. HMMER rescue 구제율이 낮은 것은
증상이지 원인이 아니다.

핵심 결론: **원인은 `steps/prune.py` 상류(클러스터링 · 크기 게이트 · 데이터 정합성)에 있을
가능성이 높고, 어느 것인지는 아직 측정되지 않았다.** 결정적 측정은 대부분 클러스터에 이미
있는 산출물로 가능하며 파이프라인 재실행이 필요 없다.

## 1.5 🔴 9차 라운드 (2026-08-22) — **이 문서의 기존 결론이 뒤집혔다**

**먼저 읽을 것.** 이 문서는 오랫동안 *"확정된 원인은 MCL/그래프 파편화"* 로 기록해 왔다.
v2 산출물에서 미배치 유전자를 전수 분류한 결과 **그 결론은 잔여 미배치 풀에 대해 성립하지
않는다.** 재진단을 시작하기 전에 아래를 읽고, MCL 그래프 개입에 착수하지 말 것.

### 측정 — `classify_unplaced.py` (신규, 테스트 14건), 4종 전부

| 종 | 미배치 | 그래프 절단 | lineage-specific | 진짜 고아 | 프루닝 |
|---|---|---|---|---|---|
| **Mcry** | 3,383 | **15.8%** ← 최저 | 10.2% | **73.9%** ← 최고 | 0.1% |
| Cgig | 738 | 33.5% | 11.7% | 52.8% | 2.0% |
| Ococ | 1,354 | 47.6% | 9.7% | 41.6% | 1.2% |
| Obas | 1,611 | 49.5% | 7.3% | 41.4% | 1.8% |

**Mcry는 그래프 파편화가 가장 적은 종이다.** 프루닝 기여도는 전 종에서 0.1–2.0%로,
"프루닝 ✗" 판정이 v2 산출물에서 재확인된다.

#### 🔴 2026-08-23 추가 — **위 표의 그림 자체가 상당 부분 모델 품질 산물이다**

CDS 무결성(시작·종결 코돈 / 3의 배수 / 내부 종결 없음)을 축으로 추가해 재측정했다
(`567e461`). **위 판정 비율은 소수점까지 그대로다** — 플래그는 가산적이다. 바뀌는 것은
**비교 가능한 모델만 볼 때의 그림**이다.

| 종 | 총계 | comparable(≥100 aa + 완전 CDS) | **comparable율** | 원시율 | 제외 (short / invalid CDS) |
|---|---|---|---|---|---|
| **Mcry** | 25,226 | 20,446 | **5.0%** | 13.4% | 3,374 / **2,237 (8.9%)** |
| Cgig | 29,163 | 29,163 | **2.5%** | 2.5% | 0 / **0** |
| CgigH | 27,583 | 24,360 | **4.4%** | 9.8% | 3,132 / 117 |
| Ococ | 33,745 | 31,529 | **2.6%** | 4.0% | 1,730 / 580 |
| Obas | 28,244 | 25,828 | **3.2%** | 5.7% | 2,383 / 36 |

Mcry:Cgig 격차가 **원시 5.4배 → 길이만 2.5배 → 길이+CDS 2.0배**로 줄어든다.

**그리고 comparable 부분집합 안에서는 판정 분포가 크게 움직인다:**

| | 전체 | comparable만 |
|---|---|---|
| Mcry 진짜 고아 | 73.9% | **53.0%** |
| Mcry 그래프 절단 | 15.8% | **35.0%** |
| Cgig 진짜 고아 | 52.8% | **52.8%** (미배치의 100%가 comparable) |
| Ococ / Obas 진짜 고아 | 41.6 / 41.4% | 32.4 / 22.7% |

**비교 가능한 모델만 놓으면 Mcry와 Cgig의 고아 비율은 53.0% vs 52.8%로 사실상 동률이다.**
Mcry 미배치 3,383개 중 comparable은 1,028개(30.4%)뿐이다. 즉 *"Mcry는 진짜 고아가
압도적"* 이라는 위 표의 인상은 **주로 비교 불가능한 모델이 만든 것**이다.

독립 확인: Mcry 미배치 중 `INVALID_CDS` **867개(25.6%)** 는 유전자 구조 축이 별도로 낸
"구조 파손 831개(24.6%)"와 거의 같은 크기다. 서로 다른 두 축이 같은 집단을 가리킨다.

⚠️ **이 표들을 인용할 때 원시율 단독으로 종간 비교하지 말 것.** `measure_v2.py`가 이제
raw 표 밑에 무조건 경고를 찍고 `comparable_unplaced_rate`를 병기한다.

### #8과 모순되지 않는다 — 서로 다른 집합을 잰 것이다

`forensics_r8.py`의 *"실패 사례의 90%가 강한 DIAMOND 히트 또는 splinter OG를 가짐"* 은
**HMMER rescue로 구제된 유전자**를 truth set으로 잰 값이다. 그 집합은 정의상 *HMM으로 family
소속이 검출되는* 유전자만 모은 것이라 히트가 있는 게 당연하다. 이번 측정 대상은 그
**여집합** — rescue조차 못 붙인 잔여물이다.

### 🔴 더 큰 문제 — 종간 미배치율 비교 자체가 교란돼 있었다

| 종 | 총계 | **최소 단백질 길이** | <100 aa | 원시 미배치율 | **≥100 aa 미배치율** |
|---|---|---|---|---|---|
| **Mcry** | 25,226 | **3 aa** | 3,374 | 13.4% | **6.2%** |
| Cgig | 29,163 | **151 aa** | **0** | 2.5% | 2.5% |
| Ococ | 33,745 | 46 aa | 1,730 | 4.0% | 2.9% |
| Obas | 28,244 | 33 aa | 2,383 | 5.7% | 3.2% |

**`Cgig.pep.fa`(MAKER)에는 100 aa 미만이 하나도 없다 — 최단 151 aa.** 즉 우리는 주석 정책이
서로 다른 프로테옴들의 미배치율을 비교해 왔다. `measure_v2.py`가 이제 길이 층화하고
`[!] CONFOUND` 경고를 낸다(`12d651f`).

### 확정된 결론

1. **Mcry 문제는 파이프라인 결함이 아니다.** 겉보기 격차 2.4–5.4배 중 **약 2/3가 입력 주석
   품질**(구조 파손 모델 + 짧은 모델 정책), 나머지 **약 1.7배는 Mcry가 패널의 유일한
   Aizoaceae라는 환원 불가능한 taxon sampling 성질**이다.
2. **15sp는 이것을 고치지 못한다.** 추정 종트리에서 Mcry는 여전히 Portulacineae+Cactaceae
   전체의 자매이고, 새로 들어오는 4종(Tfru/Ccac/Pami/Pole)은 전부 Mcry로부터 선인장과 같은
   거리다. **브리징 taxon이 하나도 늘지 않는다.**
3. **MCL inflation은 죽은 노브**(#10 스윕). **그래프 개입도 착수 근거가 없다.**

### 아직 갈리는 것 — 두 증거가 충돌한다

온전한 짧은 Mcry 모델 2,543개 중 몇 개가 진짜인가:
- **엑손 구조**: 단일엑손 58.6% vs CgigH 대조 28.3% → **초과 약 770개**
- **발현**(가뭄×일주기 TPM, 25,219/25,219 조인): **침묵 비율이 전 구간 9–13%로 평평.**
  770개가 spurious ORF라면 침묵 클래스에 몰려야 하는데 농축이 **없다**

발현은 *침묵 DNA 위의 ORF*만 배제하고 *실재 전사체 위의 오호출*은 배제하지 못한다.
**개별 식별은 미달성.** 구조 파손 831개(24.6%)만 가정 없이 방어 불가로 확정된다.

관련 이슈: #36(측정 전체), #43(주석 세대 불일치), #34(tier-3 — 표적이 짧아 기대 수익 낮음)

## 1.6 ✅ #40 종결 (2026-08-22) — Ppc1/Ppc2는 분리된다, 단 선인장 ppc-1E2는 1개뿐

**요구사항이었던 "선인장 유전자를 Ppc1/Ppc2에 배정 가능"이 충족됐다.** 결정적 요인은 앵커도
모델도 아니고 **clan 서열 집합 자체가 틀렸던 것**이었다.

### 교정 — `_M` 6개는 출처 불명의 대체 모델이었다

배포 유전자 모델(250–834 aa)을 전장 단백질(801–1,042 aa)로 갈아끼운 서열 6개가 recruitment
단계에 들어 있었고 **생성 절차가 기록돼 있지 않았다.** 자기 종 유전체에 대조한 결과 4개는
0.982–0.997로 진짜, **2개(Ococ 0.855 / Pami 0.790)는 그 종의 유전자가 아니었다.**

`reconstruct_models.py`(8373fc8)로 재구성: 비-`_M` ≥800 aa 앵커 94개 → 5개 선인장 유전체
miniprot 재예측 → **"자기 유전체 identity ≥0.95" 게이트**. 49좌위 중 44 채택 / 5 탈락(0.83–0.90).
`_M` 5개 교체 + 1개 제거 = **교정 세트 108**.

| 단계 | 서열 수 | 탈락 사유 |
|---|---|---|
| base | 109 | — |
| 교정 | 108 | `_M` 5 교체, 1 제거(게이트 0.9035) |
| 점유율 <50% | 99 | 9개 전부 575–683 aa 절단 모델 (8/9가 Ococ·Sund) |
| 코돈 정렬 | 95 | miniprot `--trans`가 frameshift를 건너뛰어 GFF CDS와 1:1 불일치 (교체본 4개) |

### 재추정 결과 — 분할 유지, 지지도 상승

```
ppc-1E1  n=62  (선인장 35)  SH-aLRT/UFboot 95/91   transferable=True
ppc-1E2  n=15  (선인장  1)  100/76                 transferable=True   overlap 없음
```

**교정 전 86.2/73 → 교정 후 100/76.** 틀린 서열을 빼니 노드가 더 견고해졌다.

### 잔기 축도 교정 세트에서 재계산 — 질적 결론 전부 유지

untrimmed 1,428컬럼 / plant-type 77서열 기준. **trim 행렬은 트리 전용**(#42 규칙).

| 지표 | 교정 전 | 교정 후 |
|---|---|---|
| SDP 컬럼 (1E1 / 1E2) | 31 / 58 | **31 / 56**, 상호 26 |
| 외군 복원 시 | 25 / 45 | **24 / 51** |
| S11 인산화 세린 | 17/17, 30/33 | **15/15, 29/32** |
| 잔기 9 (1E2 Met 고정) | 17/17 | **15/15** (1E1은 Leu 26/32) |
| 잔기 10 (1E2 Ala 고정) | 17/17 | **15/15** (1E1은 Ala 13 / Ser 9) |
| 촉매 모티프 내 SDP | 0, 최근접 611 | **동일** (1E1 Q 94% / 1E2 A 100%) |
| 크기·폭 매칭 null 500회 | 중앙값 0, 최대 1 | **동일** (p < 0.002) |
| 단계통 subclade null | 19개, 9–65 | **23개, 7–79** — 폭 논거 유지 (≥9종 최대 37 < 56) |

재현: `steps.subfamily.anchor_transferability` + `sdp_scan`을 `clan_corrected.aln` +
`clan_corrected_tree.treefile`에 적용하면 위 수치가 그대로 나온다 (순수 파이썬, 도구 불필요).

### 선인장 ppc-1E2 — 절단 정도와 프로파일 마진이 단조

| 종 | 상태 | 길이 | 마진 | 근거 |
|---|---|---|---|---|
| **Obas** | ✅ 전장, 트리 멤버 | 967 aa | **+251** | 배포 주석과 동일, occ 0.677 |
| Ococ | ⚠️ 절단 | 683 aa | **+156** | Chr10:36.21–36.23 Mb 재구성(0.985), occ 0.478 → 탈락 |
| Cjam | ⚠️ 절단 | 411 aa | **+61** | synteny(JBOBQE010000004.1 16.03 Mb), 점유율 탈락 |
| Sund | ❌ 부재 | — | — | 이웃이 chr04 anchor, 그 염색체에 PEPC 0개 |
| Cgig | ❓ 판정 불가 | — | — | 2,589 scaffold, anchor율 52% |

**967→+251, 683→+156, 411→+61.** 하나의 좌위가 계통마다 다른 잔존 정도로 남은 것으로 읽힌다.
트리 기반 멤버는 Obas 1개, 나머지는 **트리 밖 증거**(synteny·마진·좌위 재구성)로 서술한다.

산출물: `pronghorn:~/scratch/mM_repair/`. methods.md §2.X.6에 전부 반영(교정 절차 문단 신설).


## 2. 증상 지표

> ❌ **SUPERSEDED (2026-08-21):** 아래 표는 v1 런 수치다. 정본은 §6 7차 라운드의
> 5sp v2 (Mcry 13.4% 등, 전 종 개선), 15sp는 `results_15sp.md`(arbiter_v3 13,431).
> 표 자체도 rescue 전/후 혼합 가능성이 있다(아래 provenance caveat).

| 종 | 유전자 | R10 후 미배치 | 미배치율 | HMMER 구제 | 구제율 | 최종 배치율 |
|---|---|---|---|---|---|---|
| Cgig | 29,163 | 1,698 | 5.8% | 1,099 | 64.7% | 97.9% |
| Ococ | 33,745 | 2,639 | 7.8% | 1,692 | 64.1% | 97.2% |
| Obas | 28,244 | 2,364 | 8.4% | 1,262 | 53.4% | 96.1% |
| CgigH | 27,583 | 3,497 | 12.7% | 1,573 | 45.0% | 93.0% |
| **Mcry** | 25,226 | **3,765** | **14.9%** | 732 | **19.4%** | **88.0%** |

**출처 주의(provenance caveat):** 미배치(Orphan) 열은 `README.md:716-720`
(`species_comparison.tsv`)에서 왔고 합계가 정확히 13,963 = **rescue 이전** 값이다. 그런데 이
파일을 생성하는 코드 경로(`pipeline.py` → `write_species_comparison`)는 rescue **이후**에
실행되며 `steps/pseudogene.py:281` 은 해당 유전자를 `unplaced_after_rescue` 로 라벨링한다.
즉 이 표는 서로 다른 시점의 수치를 이어붙인 것일 수 있다. 다만 어느 해석으로 계산해도
(rescue 전 기준 Mcry 17.8% 미배치 / 16.3% 구제율) **Mcry가 압도적 최악이라는 정성적 결론은
동일**하다. 정확한 값은 클러스터의 원본 산출물에서 재계산할 것.

15종 런 참고치(`opuntia_analysis/family_finder/output_15sp/Orthogroups.GeneCount.tsv`):
Mcry 22,517 유전자 / 16,643 families. 단, 아래 §5의 "기각된 주장" 참조 — 15종 지표 다수는
교란 변수 때문에 증거로 쓸 수 없다.

## 3. 코드에서 확인된 사실 (전부 원문 대조 확인)

### 클러스터링 단계 (핵심 후보)

- `steps/orthofinder.py:30-35` — OrthoFinder에 `-f -t -o` 만 전달. `-I`(MCL inflation),
  `-S`(검색 프로그램), `-og` 모두 미지정. 유일한 우회로 `orthofinder_extra_args`
  (`config.py:38`, 차단 목록은 `orthofinder.py:36-42`)는 커밋된 모든 config에서 빈 문자열.
  **참고: OrthoFinder v3.1.3의 기본 inflation은 1.2다** (`src/orthofinder/__init__.py:16`,
  v2.x의 1.5에서 하향됨). DIAMOND 기본 모드는 `--more-sensitive`, `--max-target-seqs` 미지정
  (DIAMOND 기본 k=25) — v3.1.3 `run/config.json` 기준.
- `steps/orthofinder.py:24-28` — 매 라운드 출력 디렉터리를 무조건 `rmtree`. 같은 라운드를
  다른 `-I` 로 재실행하려면 DIAMOND를 다시 돌려야 한다(`-b` 재사용 불가). 라운드 간에는
  입력 서열이 달라지므로 `-b` 는 원래 라운드 간 재사용용이 아니라 **동일 라운드 스윕용**이다.
- Cgig(MAKER)와 CgigH(Helixer)는 같은 *Carnegiea gigantea* 게놈의 두 어노테이션인데 별개
  종으로 들어간다. MCL 관점에서 외군이 붙으려는 clique에 최대 가중치 이웃을 하나 더 얹는다.
- 5종 세트에는 Aizoaceae↔Cactaceae 간극을 메울 중간 분류군이 없다.

### 크기 게이트 / 데이터 정합성 (조용한 유전자 손실 경로)

- `pipeline.py:42`, `:53-55`, `:120-121` — `min_orthogroup_size`(4) 게이트 3곳. `:120-121` 은
  프루닝을 **통과한** 유전자까지 포함해 OG 전체를 outlier pool로 되돌린다.
- **`pipeline.py:53-57` 의 더 심한 누수(스트레스 테스트에서 발견):** `common_ids >= 4` 로 OG가
  **살아남는 경우에도**, CDS가 없는 유전자는 `prot_seqs` 에서 조용히 빠져 트리 leaf가 되지
  못하고, `confirmed` 에도 `outliers` 에도 없다. `:253-259` 의 미배정 보충은 "어떤 OG에도 없는"
  유전자만 회수하므로 이들은 **영구 소실**된다. Mcry의 pep/CDS ID 규칙이 조금이라도 다르면
  이것이 바로 보고된 증상의 기제다 — 크기 게이트에 걸리지 않아도 발생한다.
- `steps/align.py:72` — `_filter_internal_stops` 가 **정렬된** 단백질에 `rstrip("*")`. 정렬에서
  `"...PROT*----"` 로 끝나면 말단 stop이 internal로 오판된다. 제거된 유전자는 같은 경로로
  영구 소실. (5종 런 보고 규모는 Ococ 161 + Mcry 15 = 176개로 작음 — `methods.md:15`.
  버그는 실재하나 3,765개 Mcry 미배치의 설명은 아니다.)
- `pipeline.py:375-381` — `_write_final_output` 이 `remaining_pool`/`cds_pool` 을 받고 안 쓴다.
  최종 미배치 FASTA가 기록되지 않는다.
- `pipeline.py:169-191` + `utils/checkpoint.py:38-46` — resume은 `summary.tsv`(맨 끝에만 기록)
  에서 확정 family를 읽고, 체크포인트 status가 `completed` 가 아니면 라운드 1부터 전체 pool로
  재시작한다. 중간에 죽은 런에 대해 resume은 사실상 무의미하다.

### 프루닝 단계 (2차 — 기여도 미측정)

- `steps/prune.py:180` — S(i) = median_j(d_obs/d_exp) 를 **절대 임계값 5.0**과 비교.
- `prune.py:103-127` — TreeShrink 실패가 전부 `logger.debug`(`:109-110` 은 로그조차 없음)로
  삼켜짐. 5종 런에서 TreeShrink는 실행되지 않았고(`methods.md:70`), 로그상 "미실행"과
  "0개 발견"이 구분되지 않는다.
- `prune.py:148` — `tree.prune(surviving)` 에 `preserve_branch_length=True` 누락.
- `prune.py:44` — 종 수가 아니라 **leaf 수**를 `min_species_for_pruning` 과 비교.
- `prune.py:166-168` — 종 트리에 없는 종 쌍을 조용히 건너뜀. 전부 불일치하면 `:176` 에서
  점수 0.0 → 절대 프루닝 안 됨. 종 트리 leaf 이름과 데이터 접두어의 일치 검증이 어디에도 없다.
- `Ratio outlier:` 로그는 `prune.py:182-185`, DEBUG 레벨.

### HMMER rescue 단계

- 수렴 후 1회, 미배치 유전자에만 실행(`pipeline.py:305`). **이미 잘못된 family에 들어간
  유전자는 영원히 재검토되지 않는다.**
- full-sequence E-value 단독 판정(`hmmer_rescue.py:101`), coverage/margin/domain E-value 없음.
- **동점 편향(실증 확인됨):** `:106` 의 strict `<` 는 동점 시 tblout 첫 줄 승리. HMM DB 순서는
  `sorted(hmm_dir.glob("*.hmm"))`(`:42`) = 파일명 사전순인데 `'0' < '_'` 이므로 순서는
  `R10_*, R11_*, …, R1_*, R2_*` — **동점 시 `R10_` 계열이 이긴다(R1이 아님)**. HMMER 3.3.2
  실증: E-value 둘 다 `0` 출력 시 bit score 3795.7 짜리가 3135.2 짜리에게 진다. 동점은
  underflow(E<~1e-308 → `0`)뿐 아니라 tblout의 유효숫자 2자리 반올림으로도 흔히 발생한다.
- 배정 후 tree 재검증 없음(`:229-242`, 소속 확정이 재정렬보다 먼저), 프로파일은 프루닝 이후
  정렬에서 생성(`:250-268`) — 프루닝 오류를 스스로 강화.

### 종 트리 — 위상은 정상, 가지 길이가 문제

`((Mcry:0.5,(Cgig:0.01,CgigH:0.01):0.49):0.3,(Ococ:0.2,Obas:0.2):0.3)` (`methods.md:7`)

unroot 시 비자명 split은 `{Cgig,CgigH}` 와 `{Ococ,Obas}` 뿐 — "Mcry 외군" 트리와 위상 동일.
문제는 가지 길이: d(Mcry,Cgig)=1.00 < d(Cgig,Ococ)=1.30. 외군이 내군 한쪽에 더 가깝다고
주장하는 셈이며, 값 자체가 손으로 쓴 반올림 수다. 참고로 `family_finder.py:44` 와
`README.md:375`, `:436` 은 ASTRAL로 트리를 만들라고 안내하는데 `CLAUDE.md:73` 은 정확히
그것을 금지한다(coalescent 단위). **실제 런에 쓰인 트리 파일을 클러스터에서 확인할 것.**

### 문서 오류

`README.md:664` 가 Mcry를 *Mammillaria cristata*(선인장)로 표기 — `:324`, `:329` 도 동일 오류.
실제는 *Mesembryanthemum crystallinum*(Aizoaceae). 근거: `methods.md:7`,
`/mnt/f/scratch/develop/opuntia_analysis/scripts/cam_evolution_analysis.py:30`.

## 4. 확정 / 유력 가설 / 기각·정정된 주장

**확정 (코드 원문 또는 실증으로 검증):** §3 전체. 추가로 rescue 표의 산술
(1692+1573+1262+1099+732=6,358=45.5%), README:771 의 "18/14,712" 가 v1 개발 단계 기록이라는
독해, HMMER 동점 편향(실증), OrthoFinder v3.1.3의 `-I`/`-S`/`-og`/`-b` 지원
(`run/process_args.py:672/775/861/406`).

**유력 가설 (측정 필요) — 2026-08-22 기준 전부 해소됨, §1.5 참조:**
> ❌ **SUPERSEDED (2026-08-22):** 아래는 8차 라운드까지의 상태다. 9차 라운드에서
> 미배치 풀 전수 분류로 결론이 뒤집혔다 — 정본은 §1.5. 그래프 포렌식 3분할·
> BUSCO OJR·inflation 스윕은 **다시 돌리지 말 것.**

- pep/CDS ID 불일치로 인한 Mcry 유전자의 조용한 소실 (§3 크기 게이트 절 — 최우선 확인)
- Round-1 클러스터링에서의 실제 Mcry 오배정 (BUSCO OJR로 측정)
- MCL inflation / OrthoFinder 유전자별 포함 임계값 / DIAMOND 민감도 중 어느 것이 병목인지
  (그래프 포렌식 3분할로 판별)
- CgigH 중복 게놈의 MCL 왜곡, taxon sampling 효과 (15종 Round-1과 비교)

**기각·정정된 주장 (이슈/문서에 옮기지 말 것):**
- ~~"MCL 기본 inflation 1.5"~~ → **v3.1.3 기본은 1.2.** 스윕은 {1.05, 1.1, 1.15, 1.2(기준), 1.3}.
- ~~"동점 시 R1 계열 승리"~~ → **R10 계열 승리** (위 실증 참조).
- ~~"라운드 5-10에서 10번 더 돌고 772개 감소"~~ → 라운드 6-10, **5번** 더 돌아
  14,735→13,963 (772개, 5.2%) 감소하며 family 93/17,596개(0.5%) 추가.
- "Mcry의 family당 평균 1.35가 15종 중 최저" — 프로테옴 크기 + cactus WGD 교란. Bvul 1.39,
  Ahyp 1.42와 사실상 동일. 지표로 쓰지 말 것.
- "Mcry 포함·Cactaceae 없는 OG 935개" — 비-cactus에게 정상 (Sole 1,661, Dcar 1,725).
- "`-S diamond_ultra_sens` 가 해결책" — 기본이 이미 `--more-sensitive` 이고 Aizoaceae↔
  Cactaceae ortholog는 고정체성 구간으로 **예상**되므로 효과가 없을 것으로 **예측**된다.
  이것은 예측이지 측정이 아니다 — 그래프 포렌식이 판정한다.

## 5. 알려진 함정 (다음 작업자 주의)

- **`-og` 는 `orthofinder -h` 에 안 나온다** (v3.1.3 `helpinfo.py:274-278` 주석 처리) —
  `-h` 파싱으로 지원 여부를 검증하면 거짓 음성. `process_args.py:861` 에는 살아 있다.
- **`-og` 는 `Comparative_Genomics_Statistics/` 생성을 건너뛴다** (`main.py:203-206` 조기
  반환) — `-og` 를 채택하면 종별 통계 진단 파일이 미래 런에서 사라진다. `Orthogroups.tsv`
  계열은 유지됨(`tools/mcl.py:294-296`).
- `Ratio outlier:` 집계는 DEBUG 레벨이므로 `--verbose` 런에서만 가능.
- `opuntia_analysis` 는 이 리포에 없다 — `/mnt/f/scratch/develop/opuntia_analysis/`.
  `run_busco_species.sh` 는 `viridiplantae_odb12` 사용(`:29`).
- `config_11sp.json`, `data/`, `output_5sp/` 도 리포에 없고 클러스터에만 존재.

### 🔴 조용한 실패 — 에러 없이 그럴듯한 틀린 결과를 내는 것들 (8차 라운드에서 전부 실측)
| 증상 | 진짜 원인 | 판정법 |
| 🔴 **격리가 살아있는 입력을 부순다** | `data_15sp/{pep,cds}` 24개 파일이 **`data_12sp`로의 심볼릭 링크**였고, data_12sp를 RETIRED로 옮기자 전부 끊겼다. `measure_v2` 등은 격리 **전에** 돌아 통과했고, 끊긴 것은 다음 분석(파일럿)이 `FileNotFoundError`를 낼 때야 드러났다 — 대장에 경고해둔 "반대 방향 오류"를 그대로 저질렀다 | 격리/이동 **직후** `find <활성 트리> -xtype l` 로 깨진 링크 0을 확인할 것. 수리: 실파일 복사로 자급화(링크를 RETIRED로 되가리키지 말 것) |
| **domtblout을 유전자 단위 연속으로 가정** | hmmsearch domtblout은 **프로파일 단위 정렬**이라 같은 유전자가 파일 곳곳에 흩어진다. 유전자 변경 시 조기 flush하면 **유전자 하나가 프로파일마다 투표** — 간선 1,020,292개(상한 23,744), 13,347-family 거대 컴포넌트라는 그럴듯한 쓰레기가 나왔다 | 전 스트림 누적 후 마지막에 1회 투표. **불변식을 단언할 것**: family당 ≥0.6 표 대상은 최대 1개이므로 간선 ≤ family 수. 위반 시 산출 대신 즉사 (`analyze2.py`) |
| 🔴 **`find ~/scratch ...` 가 조용히 아무것도 안 찾는다** | **`~/scratch` 는 심볼릭 링크**(`→ /data/gpfs/assoc/pgl`)이고 `find` 는 **시작 경로가 링크면 내려가지 않는다.** 출력 없이 exit 0. 이것 때문에 "정본 branch-site 산출물이 두 호스트 어디에도 없다"는 **틀린 결론을 이슈에 게시**했다 — 실제로는 `~/scratch/bin/family_finder/pepc_pilot/seltest_fast/` 에 전부 있었다 | `find -L ~/scratch` 또는 `find ~/scratch/` (후행 슬래시), 또는 **실제 디렉터리에서 시작**. ⚠️ 링크를 지나 시작한 검색(`~/scratch/data/...`)은 정상 동작하므로 **다른 검색이 되는 걸 보고 믿으면 안 된다** |
| **코돈 정렬의 58%가 없다** | `mM_repair/p2n.sbatch` 가 pal2nal을 **`-nogap`** 으로 돌려 gap 포함 코돈 컬럼을 전부 버렸다. `codon95.codon.aln` = **599 코돈**(gap 문자 0), 원래 1,428. `cleandata=1` 함정과 결과가 같고 원인만 다르다 | 코돈 정렬은 **컬럼 수를 단백질 정렬과 대조**할 것. 선택압 검정·BEB 좌표에 치명적이다 |
| **`ldd` 로 CUDA 유무를 판정** | CUDA가 정적 링크라 **CUDA 빌드에서도 `ldd \| grep cuda` 는 0건**이다. 이 근거로 "CUDA 없음"을 판정했다가 상류에 오보를 올릴 뻔했다 | `strings <bin> \| grep -c cudaMalloc` (CUDA 39 / CPU 0) 또는 실행 |
| **주석 프로그램이 서브패밀리로 오독됨** | PEPC 인트론 이탈률이 Helixer 0.52 / EVM 0.39 / AUGUSTUS·UNR 0.00 — 서브패밀리가 아니라 **주석 세대**를 따라간다. 통제 없이 쓰면 주석 차이가 생물학으로 읽힌다 | `steps/gene_structure.py` 가 프로그램 미상이면 인트론 수치를 **아예 출력하지 않는다** |
| **발현 매트릭스를 파일명으로 지목** | SF3 74%는 `results/rnaseq/RSEM_TPM.tsv` 가 아니라 **29,405 diel 재분석**(`revision/rnaseq_29405/clustering/RSEM_TPM.average.tsv`) 값이다. 전자로는 82.8%가 나온다. `OcoChr03G21430` 이 한쪽 709.6 TPM / 다른 쪽 0.0 | 매트릭스 경로는 **CLI 입력**으로 두고 결과에 어느 매트릭스인지 남길 것 |
| **커버리지를 한 방향으로만 셈** | `region_disorder` 85/102를 "손실 17"로 읽어왔으나, 정렬 102서열 중 **7개(ATH·Aco 앵커)는 접힌 적이 없다**. 공유 가능 최대치가 애초에 95 | 양쪽 방향을 다 셀 것: `111 = 16 + 10 + 85`, `102 = 95 + 7` |
| **바이트 동일 증명이 통과했는데 경로가 증명 안 됨** | 기준선 픽스처가 분기를 안 태운다. 실측 2건: (1) rescue 기준선이 `_concat_hmms`를 덮었지만 프로파일이 ~50바이트라 8192 read 루프가 **1회만** 돌아 멀티청크 경로 미증명 (2) `narrative()` 픽스처 `_sf3_evidence()`가 `coordinates_verified=True`뿐이라 `CANNOT BE JUDGED` 문단과 `[cannot judge]` 요약줄이 **한 번도 실행 안 됨** | 체크섬 전에 **문장 커버리지**를 볼 것. `narrative()`를 만지면 픽스처에 `coordinates_verified=False`를 반드시 넣을 것 (`test_narrative_says_it_cannot_judge_instead_of_reporting_a_clean_region`, `test_unverified_coordinates_produce_an_explicit_cannot_judge_entry`가 그 분기를 탄다) |
| **통제 플래그가 통제 성립을 뜻하지 않음** | `sequence_controlled`가 "identities 인자를 받았는가"였다. reason이 `UNCONTROLLED`인 행이 True로 집계된다 | 플래그는 verdict와 같은 기준으로 정의할 것 (b73bde9). 회귀 테스트 `test_the_flag_and_the_verdict_never_disagree` |
| **fident 선형 통제가 가짜 음성을 만든다** | `alntmscore`는 1.0에서 포화하는데 `fident`는 계속 오른다 → 고identity 쌍(=within 쌍 대부분)이 구조와 무관하게 회귀선 아래. 실측 r2 0.43–0.50, 십분위 잔차 비단조(최악 13σ) | 잔차 비교 전에 **적합 적정성**과 **identity 중첩 지지도**를 검사. 미달이면 `no_interpretation_available` |
| **`.hmm` glob이 hmmpress 인덱스를 삼킴** | `*.hmm*`로 쓰면 `.h3i`/`.h3m`까지 연결된다. 그리고 정렬이 사전순이라 **`R10_`이 `R1_`보다 앞선다**(`'0' < '_'`) | E-value 타이브레이크(CLAUDE.md)와 같은 함정이 연결 순서에도 있다. 테스트로 고정할 것 |
|---|---|---|
| foldseek `createdb`가 **0.00초**에 끝남 | 3Di 없이 **아미노산 DB**를 뱉는다(로그엔 경로가 정상 출력됨). ⚠️ **2026-08-22 정정: 원인이 `.gguf` 파일 경로가 아니다.** 같은 바이너리·가중치·머신에서 `.gguf`를 줘도 72.75초 동안 정상 3Di를 만든다. 원 관측은 **불완전한 `.gguf` 다운로드**를 잡았을 가능성이 크다(소급 확인 불가) | 원인과 무관하게 **상태**로 판정: `steps/prostt5_chunks.py`의 `verify_3di_db()` — `_ss` 없음/빈 파일/인덱스 없음/**AA DB와 바이트 동일**/엔트리 수 불일치를 전부 하드 실패 |
| **foldseek GPU 여부를 로그로 판정하지 말 것 — 시계로 판정하라 (2026-08-24 확정)** | ~~`--gpu`는 파싱되지만 ProstT5 경로에 닿지 않는다~~ **2026-08-25 정정: 재빌드 바이너리에서는 플래그가 경로를 결정한다** — `--gpu 1` 없이는 CPU ProstT5ForkRunner(fork 워커가 각자 ~3GB 모델 로드 → 62GB×16 전멸 → 부모가 SysV 큐 msgsnd에 영원히 블록, CPU 0초·에러 0). `Use GPU 0`은 **다른 mmseqs 노브**이고, 실제 구동은 ggml이 등록하는 `CUDA0` 백엔드다. 이 로그를 믿고 **다섯 번 오진**했다 | **실행 시간으로 판정.** 112서열, 동일 커밋(`18478813`)·동일 심볼(39)·둘 다 `Use GPU 0`: 기존 `~/bin/foldseek-cuda` **497초** vs `-DENABLE_CUDA=1` 재빌드 **12초** = **41배**. 도구 표의 37.6배는 맞았고 우리가 못 쓰고 있었을 뿐. 전수 459K: 23일 → **약 14시간** |
| **재빌드 함정 2건** | ① `rust`가 위키 conda 목록에 없는데 mmseqs corrosion이 요구 → 없으면 cmake `Configuring incomplete`. ② 이 호스트에 CUDA 런타임 3벌(11.5 `/usr/lib`, 12.6 `/usr/local`, 12.6 conda)이고 맨 `nvcc`는 **11.5**가 잡힌다 → `CUDA_CUDART`를 12.x로 명시 안 하면 `undefined reference to cudaGetDeviceProperties_v2`(CUDA 12 심볼)로 링크 실패 | 레시피는 `steps/prostt5_chunks.py` 상단 |
| ~~`--gpu 1`을 줬는데 GPU 0%~~ ~~기각 2026-08-22~~ **← 이 기각이 틀렸다 (위 행 참조)** | ~~`No GPU devices found`로 exit 1, 빌드에 CUDA 없음~~ | 재측정에서 exit 0이고 CUDA 심볼 39건. `ldd` 대신 `strings`를 쓰라는 §5의 다른 규칙은 유효하나, 그 규칙으로 내린 결론이 반대로 갔다 |
| codeml 잡이 SLURM에 **FAILED** | `error: end of tree file` — 결과를 다 쓴 **뒤** exit 1. **3회 재현** | SLURM 상태 말고 `lnL` 줄 + BEB 블록 존재로 판정 |
| BEB 확률이 파싱값과 다름 | results.txt에 **NEB 블록이 BEB보다 먼저** 나오고 값이 다르다 (site 568: NEB 0.931 vs BEB 0.970) | `Bayes Empirical Bayes` 헤더 **이후만** 파싱 |
| 구조 관련 대조가 이상하게 적게 매칭됨 | 구조 파일명은 유전자 ID의 `.`을 `_`로 바꾼다. 정규화 없이 대조하면 **111개 중 29개만** 매칭되고 에러가 없다 | 양쪽 ID를 canonical form으로 |
| 서브패밀리 이름 support가 **1.00**인데 미덥지 않음 | support = 이긴 가중치 / **투표한** 가중치. 20멤버 중 2개만 주석이면 1.00이 나온다 | **`coverage`** 컬럼을 볼 것 (0703796) |
| 청크 hmmsearch 결과가 단일 실행과 `diff`에서 다름 | HMMER가 파일 내 최장 target 이름에 맞춰 **컬럼 padding**을 맞춘다. 값은 동일 | 공백 정규화 후 비교 (`awk '{$1=$1};1'`) |
| tblout 병합이 조용히 적은 결과를 냄 | 죽은 태스크가 **겉보기 멀쩡한 잘린 tblout**을 남긴다 | `merge_tblouts`가 `# [ok]` 검사로 차단 — 우회하지 말 것 |
| 원격 축 결과를 로컬 merge에 넣으면 실패 | 축은 원격에서 돌아 경로가 원격이다 | `annotate_stack.py`가 fetch 후 로컬 경로로 명령을 출력 (8c5cdc4) |
| `~/dir` 가 아니라 `~`라는 **디렉터리**가 생김 | `shlex.quote('~/x')` → `'~/x'` (리터럴) | `qpath()` 사용 |
| 종트리가 검증을 통과했는데 프루닝이 무력 | 가지 길이가 **전부 `1.0`인 위상 전용 트리**. 양수라 `Non-positive` 검사를 통과하고, 최대 쌍거리가 **정확히 10.0**인데 게이트가 `> 10.0`이라 딱 빗나간다. `data_12sp`·`14sp`·`15sp`·`17sp` **네 패널이 이 상태로 돌았다** | `validate_species_tree`가 이제 **모든 가지 길이가 동일**하면 실패시킨다 (`ce05e53`). 임계값 조정은 근본 수정이 아니다 — 전부 `0.1`이면 최대 0.4로 어떤 범위 검사도 통과한다 |
| 종간 미배치율 비교가 그럴듯한데 틀림 | 프로테옴마다 **주석 길이 하한이 다르다**. Cgig(MAKER) 최단 151 aa / Mcry 최단 3 aa | `measure_v2.py --pep-dir`로 길이 층화. 하한이 다르면 `[!] CONFOUND` 경고 (`12d651f`) |
| SDP 스캔이 "신호 없음"을 냄 | 트리용으로 **trim한 행렬을 잔기 분석에 재사용**했다. `-gappyout`은 gappy 컬럼을 자르는데 **N말단 조절 확장부는 어느 패밀리에서나 gappy하다** — 가장 볼 가치가 있는 구간이 먼저 잘린다. PEPC에서 CAM 인산화 모티프 SIDAQL(잔기 11–16)이 통째로 사라졌다 | **SDP는 항상 untrimmed 정렬에서 스캔.** 줄여야 하면 전역이 아니라 **그룹별 최대 점유율**(`select_columns(groups=...)`, `e31372d`) — 전역 임계값은 소수 서브패밀리 고유 구간을 산술적으로 지운다(천장 = k/N) |
| 그룹에 SDP가 없다고 보고됨 | `min_cover`(기본 0.7)가 그 그룹의 컬럼을 **통째로 조용히 스킵**했다. PEPC ppc-1E1은 N말단을 70개 중 33개(47%)만 덮어 **1468컬럼 중 539개(37%)가 판정조차 안 됐다** | `coverage_suppressed()`를 `sdp_scan` 옆에서 같이 호출해 보고할 것 (`e31372d`). ⚠️ 불변 코어 지표에는 아직 같은 보고가 없다 |
| TE 중첩 비율이 종간에 깨끗하게 갈림 | EDTA 라이브러리가 ***Opuntia*에서 만들어졌다.** 선인장 대조군(CgigH 6.1%)이 Aizoaceae(Mcry 0.9%)보다 라이브러리에 가깝다 — **TE 함량이 아니라 계통 거리를 재고 있다** | 종별 repeat 라이브러리 없이 종간 TE 비교를 하지 말 것 |
| 유전체에 없는 서열이 clan에 들어와 있음 | `_M` 접미사 6개는 절단된 배포 모델을 대체한 전장 PEPC인데, **자기 종 유전체 확인이 없었다.** 4개는 진짜(0.98–0.997), **2개는 그 종의 유전자가 아니다**(Ococ 0.855, Pami 0.790). Ococ 것은 4개 어셈블리(purge 이전 포함) 어디에도 없고 Obas 유전체에 0.9884 | `reconstruct_models.py` (`8373fc8`) — miniprot 재예측 + **"자기 종 어셈블리 최고 히트" 게이트**. identity 임계값이 아니라 *어느 유전체에 붙는가*로 판정. 탈락분은 사유와 함께 TSV 출력 |
| miniprot CDS가 단백질과 안 맞음 | `--trans`는 **frameshift를 건너뛰고** 단백질을 출력하므로 GFF의 CDS 구간과 1:1 대응하지 않는다. 4/4 불일치(±3~14 aa) | 코돈 정렬에 쓰려면 frameshift 없는 예측(`fs:i:0`)만 쓰거나 코돈 인지 재예측. 단백질 정렬에는 그대로 사용 가능 |
| 같은 이름의 두 서열이 **다른 단백질** | #40 clan 트리는 `Ococ_OcoChr10G09070_M`(**949 aa**, SIDAQL 있음)을, synteny 채점과 배포 주석은 `Ococ_OcoChr10G09070.t1`(**660 aa**, N말단 없음)을 썼다. `_M`은 recruitment 단계의 재예측 모델이고 6개 있다(Ococ/Pami/Sund×2/Tfru×2) | 좌위 이름으로 조인하기 전에 **서열 길이·N말단을 대조**할 것. `_M`은 isoform이 아니므로 `utils/gene_ids.py`가 정당하게 매칭 실패시킨다 |
| ID 불일치로 커버리지가 조용히 깎임 | 모듈마다 정규화 규칙이 달랐다. 구조 파일명은 `.`→`_`, transcript 접미사는 한쪽만 붙는다. 과거 111개 중 29개 매칭, `region_disorder`에서 16/111 손실 | `utils/gene_ids.match_ids()` — exact→canonical→isoform 순으로 필요한 만큼만 느슨해지고 **어느 레벨을 썼는지·무엇이 안 맞았는지·비율**을 반환. `max_unmatched` 상한 초과 시 실패 |
| 잔기 좌표를 다른 정렬끼리 교차 | `beb_cross`가 codeml 사이트(패밀리 코돈 정렬)와 signal window(clan 단백질 정렬)를 **구간 겹침만으로** 교차했다. 불일치 시 0 overlap이 나오고 `classify()`가 그것을 neofunctionalization **반증**으로 승격시킨다 | `alignment_id()`로 stamp를 달고, 다르면 `translate_columns()`로 다리 서열 경유 번역. `cross_windows`가 이제 요구한다 (`9a986f9`). **PEPC 실측 대조 결과 세 축 전부 1,371컬럼/102 taxa로 일치했으므로 기존 SF3 결론은 안전** |
| 미배치 유전자가 100% pseudogene 후보로 나옴 | 파이프라인이 **모든 미배치 유전자를 `is_orphan`으로 자동 플래그**한다. 순환 지표다 | `pseudogene_candidates.tsv`의 `is_orphan`·테이블 존재 여부를 미배치 질문의 증거로 쓰지 말 것. CDS를 직접 읽을 것 |

## 6. 다음 단계 — 열린 이슈 (2026-08-22 기준)

> ⚠️ **CAUTION (2026-08-27):** 이슈 상태·다음 단계는 §7이 정본이다. 아래 이슈 표는
> 2026-08-22 시점 (이후 #34 닫힘, #44 압축, #43 결정 완료 등 — §7 참조). 이 절의
> "N차 라운드" 블록들은 라운드별 이력 기록으로 유효하나, 각 라운드 안의 결론은 후속
> 라운드·§5에서 정정됐을 수 있다 — 수치 인용 전 §5·§7 대조 필수.

`wyim-pgl/family_finder` (PUBLIC). **#1~#31 중 열린 것 없음** — 아래가 전부.

| 이슈 | 우선순위 | 내용 | 막힘 |
|---|---|---|---|
| [#47](https://github.com/wyim-pgl/family_finder/issues/47) | **P0** | 15종 v2→v3 완결 경로 — 순환 중단 규칙 + 진행 로그 | 캠페인·codeml 실행 중 |
| [#35](https://github.com/wyim-pgl/family_finder/issues/35) | P0 | 15sp v2 — **완주**, 수용 판정 완료(배치율 PASS / CAM쌍 FAIL→검증 캠페인으로 이관) | — |
| [#44](https://github.com/wyim-pgl/family_finder/issues/44) | **P1** | codeml 미설정 + **정본 branch-site 산출물 부재** | — (코드 완료 `b82fee3`, 재생성 미실행) |
| [#42](https://github.com/wyim-pgl/family_finder/issues/42) | P1 | SDP/잔기 해석 일반화 | 임계값 보정만 #35 대기 |
| [#43](https://github.com/wyim-pgl/family_finder/issues/43) | P1 | 원고 불일치 3건 | — |
| [#33](https://github.com/wyim-pgl/family_finder/issues/33) | P1 | ATH MYB 앵커로 서브클레이드 분리 | 앵커 확보는 지금 가능, 본작업은 #35 |
| [#34](https://github.com/wyim-pgl/family_finder/issues/34) | P1 | tier-3 구조 스크린 | **함정 2건 기각됨** — 아래 §5 |
| [#36](https://github.com/wyim-pgl/family_finder/issues/36) | P1 | Mcry 미배치 | #35 재측정 |
| [#38](https://github.com/wyim-pgl/family_finder/issues/38) | P2 | 주석 축 — **두 축 구현 완료** (`e1a794f`) | cis-element만 남음 |
| ~~#39~~ ~~#40~~ ~~#41~~ ~~#32~~ ~~#37~~ | — | **전부 닫힘** (2026-08-21~22) | — |
| [#25](https://github.com/wyim-pgl/family_finder/issues/25) [#26](https://github.com/wyim-pgl/family_finder/issues/26) | — | SPEC / 트래커 (살아있는 문서) | — |

**3차 라운드 — v2 파이프라인 + 스윕 판정 (2026-08-19~20):**
- **v2 tiered 아키텍처 확정·구현 완료** (#21): Tier1 클러스터링(OrthoFinder/SonicParanoid2 분기) →
  Tier2 프로파일(#13) → Tier3 Foldseek(#17) → EPA-ng 심판(#23, 병합·모호 판정) → relative 프루닝(#14).
  전 계층 코드 병합·push, 테스트 165건 통과. 병렬 앙상블 기각(합의 규칙이 winner-take-all 재생산).
- **인플레이션 스윕 최종 판정: 죽은 노브.** 5sp {1.05–1.3} 전 구간 ALL 4.4–5.1%로 평탄,
  재구축 I=1.2 대조군이 원본과 완전 동일(그래프 구축 결정적). 15sp {1.05–1.2} 전 구간에서
  **Mcry Ppc1 미배정 유지** — edge가 MCL 이전 포함 임계값에서 잘리므로 원리적으로 닿지 않음.
- **노브 서열 (실측 co-clustering)**: CgigH 제외 9.0% (3–13배) ≫ 인플레이션/재구축 (±0.5pt) ≈
  taxon sampling (+0.6pt). → CgigH 제외 채택 권고 (#12 닫힘).
- **flagship 수용 기준**: Ppc1 쌍 (Mcry_Mcr8G11630↔Ococ_OcoChr03G21370, 74.7% id) —
  5sp TOGETHER / 15sp SPLIT(Mcry 미배정, hit 22개 bitscore ≤1701에도). Ppc2 쌍은 양 패널 SPLIT.
  cam_pairs.tsv + check_pairs.py 로 어느 Results든 즉시 판정. PEPC clan은 ≥5 family로 파편화 —
  단일 점수 재배정 불가(nbits 마진 0.015), 병합+트리 심판 필수 (#13 사례 연구).
- **종 트리 추정 완료** (#14 닫힘): 500 loci × 901,281 sites IQ-TREE GTR+G —
  외군 등거리 복원(0.278–0.284), Ococ–Obas 22배 과대추정 교정(0.40→0.018).
- **GPU 파일럿 진행 중** (pgl-gpu RTX 4090): PEPC clan 95개 ESMFold 접기 + ProstT5 3Di 벤치 +
  SonicParanoid2 5sp/15sp(--graph-only, d2v ZeroDivisionError 우회). 문헌 dossier 요점은 #21 코멘트.
- 신규 모듈: steps/{epa,sonicparanoid,esm,plant_glm,ecforest}.py + annotate_families.py

**8차 라운드 — Ppc1/Ppc2 분리 성공 + 앵커 전이 가능성 진단 (2026-08-21~22):**
- 🎯 **Ppc1/Ppc2 분리 해결 (#40).** 결정적 요인은 **뉴클레오타이드(코돈) 데이터** — 아미노산
  트리(PMSF, LG+G, 둘 다 AA-109)는 **여전히 실패**한다. 재구축 `N1_codon123`:
  - **ppc-1E1 (Ppc1형): 70잎, SH-aLRT 96.2 / UFBoot 94** — 식물형 선인장 PEPC **45개**
  - **ppc-1E2 (Ppc2형): 17잎, 86.2 / 73** — 선인장 **2개** (`Obas__…_000494`,
    `Ococ_OcoChr10G09070`, 둘 다 Opuntioideae, 서로 자매 100/100)
  - 45+2 = **47 = 식물형 선인장 PEPC 전부**, 두 집합 **disjoint**, 합집합 = 식물형 87잎 @100/100
  - 멤버 구성이 codon123·codon12·AA102 **세 트리에서 동일**
  - 베이스라인 `clan_tree.treefile`에서는 **Ppc1이 Ppc2를 제외하는 clade를 아예 갖지 못했다**
    (자기 자신 1잎, 최소 정보 clade 24잎 @ 지지도 27).
  - ⚠️ **ppc-1E2는 채택이지 증명이 아님** — UFBoot 73으로 기준(≥95) 미달. `min_support=80`
    게이트에서 False로 떨어진다. ppc-1E1(94)은 견고.
  - ⚠️ Ppc2 선인장 증거 얇음: 47개 중 2개, 둘 다 Opuntia. 유일한 Cactoideae 후보
    `Cjam__JBOBQE010000004.1_002232`는 점유율 <50% 필터에 제거됨. **다음 단계는 정렬에
    단편을 더 넣는 게 아니라 해당 좌위의 synteny를 세 유전체에서 직접 보는 것.**
  - 산출물: gpu `~/pepc_pilot/ppc_resolve/rebuild/`
- **아이스플랜트에 Ppc3 없음 (사용자 확정).** Mcry PEPC는 `Mcr8G11630`(Ppc1, P10490) /
  `Mcr7G08600`(Ppc2, P16097) / `Mcr6G09730`(세균형 PPC4) **3개뿐**. Mcry에 PPC3 라벨이
  붙으면 **오전이(false transfer) 신호**로 판정할 것.
- **신규 진단 도구 `anchor_transferability()`** (`steps/subfamily.py`, 228c984→**3fc5305**):
  레퍼런스 종의 서브패밀리 **이름이 질의 종에 전이 가능한지**를 사전 판정. 라이벌 라벨을
  포함하지 않는 **최대** clade + `min_support` 게이트(IQ-TREE `aLRT/UFBoot`는 **약한 쪽**으로
  판정) + 중복 감지(clade 바깥이 무엇을 더하는지).
  실측: **ATH PPC1/PPC3은 지지도 100의 자매(배추과 자체 중복) → 전이 불가**, PPC2도 미해상,
  **PPC4(세균형)만 전이 가능**. anchor-free clade 6개 = lineage-specific expansion.
  ⚠️ 초판(228c984)은 **최근접 이웃 clade**를 돌려주는 버그가 있어 선인장 0개로 오판했다 —
  서브에이전트가 잡음. 3fc5305에서 수정.
- **ATH 심볼로 명명하면 어떻게 되는지 실측**: 서로 다른 **6개 서브패밀리가 전부 "PPC1"**,
  **PPC2는 아무도 못 받음**. SF3(=OG2)은 유전자 **1개** 히트로 "PPC1"이 되는데 support는
  **1.00**으로 표시된다.
  → **명명 코드 결함 수정 (0703796)**: `support = 이긴 가중치 / **투표한** 가중치`라 20멤버 중
  2개만 주석을 가진 clade가 support 1.00으로 나왔다. **`coverage = 주석보유/전체멤버`** 필드
  추가 + 서브패밀리 이름이 멤버 절반 미만에 근거하면 경고. v2 런에서 정확히 5개에 발동.
- **Musa MYB 논문 방식(#33)은 우리 파이프라인에 이미 있다** — 앵커 합류·ML 트리·앵커 clade
  소속 분류·lineage-specific clade 보고·모티프(SDP)까지. 없는 것은 intron/exon 구조,
  cis-element, 체계화된 발현(#38). **단 논문은 "레퍼런스 서브패밀리가 전이 가능하다"를
  가정한다** — MYB는 고대 clade라 타당하지만 PEPC는 아니었다. `anchor_transferability()`가
  이 차이를 사전 판정한다.
- **#37 청크 rescue 실경로 검증 완료** — 클러스터에서 실제 Config로 `_run_hmmsearch_chunked`
  호출, **SLURM 배열 6097076_0..3** 제출·완료, 단일 실행과 **배정 31/31 일치**, 잘린 청크·
  누락 청크 **양쪽 차단 확인**. `sbatch --wait` 분기가 실제로 동작한다.
- ⚠️ **15sp 발사 선결 3건 (#41)**: (1) `data_15sp/species_tree.nwk`의 **가지 길이가 전부 1.0**
  (더미) — 종 인지 프루닝이 무력화된다. (2) `config_15sp.json`이 v1 구성 —
  `profile_assign_per_round`·`hmmer_chunk_size`가 빠져 기본값(false/0)으로 돈다.
  (3) 입력 실측 **484,752 vs 143,961 = 3.37배** → hmmsearch 약 11.3배 ≈ **65시간**.
  **결정 대기: 종트리를 (a) 2패스 재추정 / (b) 대체 트리 중 무엇으로 할지.** 권고는 (a).
- ⚠️ **ProstT5 조용한 실패 2건 (#34)**: `--prostt5-model`에 **`.gguf` 파일 경로**를 주면 3Di를
  만들지 않고 **0.00초에 아미노산 DB**를 뱉는다(에러 없음, 로그엔 경로가 정상 출력됨).
  **디렉터리**를 줘야 `_ss` DB가 생긴다. 그리고 **`--gpu 1`은 createdb에서 무시된다**
  (`Use GPU 0`). 실측 4.3 s/서열(16스레드) → 미배치 9,786개 ≈ 12시간, CPU 전용, 체크포인트
  없음 → 분할 필요.
- **`annotate_stack.py` 실행 검증에서 버그 수정 (8c5cdc4)**: 축은 원격에서 도는데 출력된 merge
  명령이 **원격 경로를 로컬 명령에 박아넣고** 있었다. fetch 단계 추가. `qpath()`(물결표 인용)
  수정도 실전 확인 — `~`라는 이름의 디렉터리가 생기지 않음.
- 문서: README 전면(micromamba 통일, 주석 스택 설치, **GPU/CPU 요구사항 실측**: RTX 4090 24GB /
  드라이버 560.35.03 / CUDA 12.6, DeepLoc 5.7GB VRAM, ESMFold가 진짜 제약, 디스크 ~75GB,
  조용한 CPU 폴백 3사례), methods.md §2.X.6·§2.X.7 신설.
- **테스트 352 그린.** 커밋: 496cfb4 4371899 14e5ed6 22920df 6ece01a 850b9ce 02b0181 cb09673
  e66d2bb 0afb62e 135b621 dbe6beb 5a38495 8c5cdc4 0703796 228c984 3fc5305 + 문서 다수.

**7차 라운드 — 주석 스택 완성·SF3 판정 확정 (2026-08-21):**
- **SF3 = subfunctionalization, 6축 확정.** RELAX 이완(K=0.636, p=1.1e-16) · 발현 74% 독식 ·
  신호영역 분할(col 165-209) · MEME 영역 내 0 · **BEB 영역 내 0** · **구조적 무질서**.
  branch-site **LRT=15.3631, p=8.87e-05** — ω 0.5/1.5/3.0 세 시작값 **전부 lnL −49044.147942
  (편차 0)** 로 수렴해 국소최적 유보 해제. 정본 BEB는 fast-track 46-taxa(컬럼 무손상 46×4113
  vs 102×4113, 좌표 독립 검증: site 102→D/GAC 등). full 102-taxa는 23h/100 iteration 미수렴 —
  정본 아님.
- **신규 6번째 축 — region_disorder.py**: SF3 Δ pLDDT **−0.304 (9/9)**, 비-SF3 76개 동일 컬럼
  대조 −0.071(중앙값 **+0.009**), **Mann-Whitney p=1.95e-05**. 대조군이 곧 발견 — N말단은
  어느 단백질이나 낮으므로 focal 수치 단독은 무의미.
- **구조 응집도 최초 실측**: 6개 중 5개 응집(비 중앙값 1.24). ⚠️ foldseek bits는 서열 유사도와
  교란 — 통제 전 인용 금지. OG4는 비응집(0.914).
- **ProstT5는 subfunctionalization에 직접 기여 못 함** — DB에 좌표가 없어 pLDDT·TM 불가.
  역할은 tier-3 소속 배정(값싼 유전체 스크린)에 한정.
- **이슈 #27·#28·#29·#30·#31 닫힘.** EC 축은 ECForest 기각 → eggNOG-mapper + CLEAN(앵커 7/7
  통과)로 교체. 구조 기능 전이 99/100이 EC 4.1.1.31, Mcry flagship 2건이 자기 종
  SwissProt(P10490/P16097) fident 1.000 회수.
- **공통 한계(스펙 명문화)**: 서열·orthology·구조 어느 축도 **잔기 수준 촉매 결손을 못 봄** —
  SF2 array를 셋 다 고신뢰 PEPC로 호출. EC·구조는 도메인 증거와 병용하는 주석 레이어.
- **5sp v2 완주** (6091850, 1-05:45:06 COMPLETED, 14,949 family / 8 라운드).
  **수용 지표 PASS — Mcry 미배치 13.4% vs v1 14.9%.** 5종 전부 개선:
  Mcry 14.9→**13.4%**(−1.5pt) / Cgig 5.8→2.5 / CgigH 12.7→9.8 / Ococ 7.8→4.0 / Obas 8.4→5.7.
  종마다 `placed + unplaced = total` **정확히 일치**(유전자 유실 0).
  ⚠️ **Mcry 문제는 미해결**: 개선폭이 가장 작고(−1.5 vs 나머지 −2.7~−3.8), 절대값도
  선인장 3종의 2~5배로 남음. MCL/그래프 파편화 진단(§ 위)은 여전히 유효한 과제.
  ✅ **CgigH 제외 근거 강화**: 클러스터링에서 빼고 프로파일/구제로만 배정했는데도
  자신의 미배정률이 12.7→9.8%로 개선 — 다른 종 이득이 CgigH 손해를 대가로 한 게 아님.
  HMMER rescue **5h43m**에 29,943 유전자 → 13,567 family.
- **#31 청크 hmmsearch**: 실데이터 등가성 통과(배정 4,379건 전건 일치; raw는 컬럼 padding만
  상이). 15sp는 프로파일×서열이라 2.7–4.1일 → **3일 제한 초과 위험, 옵션 필수**.
- ⚠️ **codeml FAILED는 허위 신호일 수 있음** — `error: end of tree file`로 결과를 다 쓴 뒤
  exit 1. 3회 재현. `lnL`+BEB 블록 존재로 판정할 것. tblout/results 파싱 시 **NEB 블록이
  BEB보다 먼저** 나오고 값이 다름(568: NEB 0.931 vs BEB 0.970).

**4차 라운드 — GPU 파일럿 완결, v2 구성 확정 (2026-08-20):**
- **구조 판정 (flagship)**: ESMFold(111/111, 두 flagship pLDDT 89.3) + foldseek —
  **Ppc1 쌍 qtmscore 0.90/0.93, E~1e-100, 한 클러스터** ✅. Ppc2 쌍도 동일 클러스터.
  97/111이 단일 구조 클러스터 = 구조는 clan 전체를 한 family로 봄(과병합 주의 → 소속 채널로만).
- **ProstT5 (접기 생략)**: 같은 쌍을 E=0.0 / bits ~6,020으로 재현 → 게놈 스크린은 ProstT5,
  접기는 최종 확인용. 주의: ProstT5 DB엔 CA 좌표 없음 → TM류 출력 불가(bits/E만).
- **SonicParanoid2 기각** (#22 닫힘): --graph-only(2.0.9 d2v ZeroDivisionError 우회) 실측 —
  배정률 5sp Mcry 65.1%/Cgig 68.2%(OF 86.8/95.9 대비 붕괴), 15sp Mcry 74.7%.
  Ppc1 양 패널 SPLIT(15sp에선 배정은 됨, Ppc2는 15sp TOGETHER — 부분 점수). BBH 문헌 예측 적중.
- **15sp 인플레이션 스윕 완결**: I∈{1.05,1.1,1.15,1.2} 전부 Mcry Ppc1 미배정, 지표 평탄 —
  포함 임계값이 MCL 이전에 자르므로 원리적으로 무효(5sp 대조군은 원본과 완전 동일 = 결정적 구축).
- **ATH 앵커 검증** (PLAZA selected ath): PPC1/2/3/4 4개 앵커 vs clan ProstT5 —
  두 flagship 모두 plant-type 확정(세 앵커에 E=0 동급 hit, bacterial-type과 ~2,500 bits 분리).
  plant-type 내부 구분은 마진 2–25 bits로 불가 → 하위 계보는 트리(EPA-ng) 소관.
  ATH 앵커는 clan 참조 트리의 라벨링 leaf로 재사용 가능.
- **v2 최종 구성 (전부 실측 근거)**: Tier1 = OrthoFinder + clustering_species_exclude(CgigH,
  co-cl 9.0%) · I=1.2 유지 / Tier2 = per-round 프로파일(#13) / Tier3 = ProstT5-foldseek
  (접기는 확인용) / 심판 = EPA-ng + ATH 앵커(#23) / 프루닝 = relative + 추정 종 트리(#14)
- flagship 스코어보드: 서열 그래프 ❌(모든 인플레이션·양 패널) · taxon ❌ · CgigH 제외 ❌(clan) ·
  BBH ❌ · 프로파일 단독 ❌(모호) → **구조 ✅ + 외부 앵커 ✅ + 트리 심판(설계)**

**다음 실행 항목:**
1. `clustering_species_exclude=[CgigH]` 소형 기능 구현 + 5sp v2 재실행(프로파일 per-round on)
2. ProstT5 게놈 전체 스크린: 최종 미배치 pool + family 대표 → 3Di 검색 → tier3_assign
3. EPA-ng LWR 보정: PEWO prune-and-replace를 PEPC clan에서(ATH 앵커 포함 IQ-TREE 참조)
4. clan 병합 실행: PEPC 5개 fragment → 구조 소속 + EPA 트리 심판으로 재구성 (Ppc1 쌍 TOGETHER 확인)
5. Mcry 미배치율 14.9% 대비 v2 최종 측정

**5차 라운드 — 서브패밀리 방법 확정 + 기능 축 + 이슈 재편 (2026-08-20 후반):**
- **서브패밀리 top-3 비교 완료**: Possvm(-r 앵커, min_support 95, bacterial 루팅; numpy<1.24 핀 +
  -skipprint 필요) 30 OG / TreeCluster(s=0,t=1.0) 9클러스터 / 기존 greedy 3 SF —
  **완전한 계층적 중첩, 경계 충돌 0** (SF3 ≡ Possvm OG2 정확히 10/10). 방법론 채택:
  HOG 정의 + Possvm 주분할 + TreeCluster 교차검증 (dossier 전문은 세션 로그, 요약은 #21 코멘트)
- **subfunctionalization 4축**: 발현 ✅ 양성(SF3 단독 Ococ 사본 = 전체 PEPC 발현 74%, 4,777 TPM) /
  위치 ✅ 음성(전원 세포질; NES/NLS 상관 → #24) / EC 🔄(ECForest 실행 중) / 선택압 ⏳(미착수)
- 기능 잔기 519/780/890 (ATH 좌표): 3 SF 모두 불변 E/T/R — 촉매 분화 없음 (purslane 좌표
  재매핑 필요, #26)
- **이슈 재편**: #25 = 전체 파이프라인 SPEC(정본) / #26 = 실시간 트래커 / #21 닫힘(증거 로그).
  오픈: #24, #25, #26 (등록 26건 중 23건 해결; #18은 2026-08-20 구조적 차단으로 닫힘)
- DeepLoc 2.1 설치 완료(라이선스 tarball, transgenic env, setuptools<81 + sentencepiece 함정),
  wiki guide/installs.md에 전체 설치 기록 push됨

**6차 라운드 — 선택압 검정 가동 + ECForest 음성 판정 (2026-08-20 오후, /clear 후):**
- **ESM-ECForest 최종 판정: 도구 신뢰 불가** ❌ — known-answer 실패 (102서열 중 EC 4.1.1.31
  예측 0건, ATH PPC4가 non-enzyme, EC 확률 0.03–0.10 uniform noise). SF2 array 판정(21430
  non-enzyme / 21450·21460 low-conf enzyme)은 인용 불가. **EC 축은 도메인/촉매잔기 증거
  (PF00311 결손, catalytic His 부재)로 대체** — #26 코멘트에 전체 데이터, CSV는 gpu
  `~/pepc_pilot/clan_anchor_predictions.csv`. 함정: fasta의 terminal `*`가 ESM 알파벳에
  없어 KeyError — `clan_anchor_esmclean.fa` 사본으로 우회.
- **선택압 검정 (대기열 1번) 가동**: clan 코돈 정렬 완성 (pronghorn
  `$FF/pepc_pilot/clan_codon.aln`, 102×4,113nt, 제외 0건) → SF3(=Possvm OG2 10-leaf,
  단일 clade 검증) `#1` 마킹 → branch-site alt/null SLURM 6091847/6091848 실행 중
  (`$FF/pepc_pilot/seltest/`, 짧은 ID g001–g102 매핑 `id_map.tsv`).
  함정 2건: **cleandata=1이 clan 정렬에서 전 사이트 제거**(1차 제출 즉사 → 0으로 수정),
  **Biopython phylip 변환의 10자 ID 절단**(짧은 ID로 우회).
- **RELAX 완료 — SF3 선택 이완 확정: K=0.636, LRT=68.3, p=1.1e-16** (gpu, HyPhy
  2.5.100 `hyphy_new` env, test=SF3 clade 19 branch, 산출 `~/pepc_pilot/relax_sf3.json`).
  subfunctionalization 정합, 강화/양성선택 주도 신기능화 배제 방향.
  hyphy 채널 순서 함정(bioconda 우선→2.5.29 구빌드)은 wiki installs.md 기록됨.
- **기능 잔기 재매핑 완료**: ATH 좌표 유효 — E519/T780/R890 = **AT1G53310.2** 기준(정렬
  컬럼 861/1152/1289), purslane 자체 번호로는 성립 안 함 → ATH 번호로 인용.
  **invariance 예외 2건**: `Obas__JBFLFP010000003.1_000523` E→D (**SF3 구성원**, 941aa
  완장 — RELAX 이완과 정합), `Ococ_OcoChr03G21480.t1` T→D+R결손 (OG5, SF2 array 아님).
  21430은 T/R 컬럼 gap(절단) — catalytic incompleteness 부합. 5차 라운드의 "3 SF 모두
  불변" 문장은 이 예외 명시로 정정됨. 스크립트 gpu `verify_residue_numbering.py`.
- **MEME·aBSREL 완료 (gpu, hyphy 2.5.100)**: MEME — 1,371코돈 중 FDR 통과 0
  (명목 p≤0.1 26개; attention 영역 col165–209 내 0개; 촉매 E/R 인접 codon 857/858·
  1285/1290은 명목 수준 기록만). aBSREL — 19 test 가지 중 보정 유의 1개 =
  **Tfru 말단**(p=6.0e-13, 5.7% 사이트 ω≈40; 조건적 CAM 종의 종 특이 후행 에피소드),
  SF3 stem 비유의. 산출 `~/pepc_pilot/{meme,absrel}_sf3.{json,log}`.
  함정: hyphy_new env는 scratch 쪽 micromamba-root에 있음 — `-n` 말고 바이너리 직접 경로.
- **fast-track codeml BEB** (pronghorn seltest_fast/, 6091859 alt/6091860 null):
  46서열 서브셋(SF3 10 + 앵커 7 + OG별 최장 대표 2) — **정렬 컬럼 무손상**이라 BEB
  사이트 번호가 signal_windows 컬럼과 직접 호환. full(102-taxa)은 robustness 확인용으로
  병행 유지(과학적 필수는 아님 — HyPhy 삼중으로 가설 판정 완료).
- **15sp v2 예상 소요** (5sp 실측 외삽, 입력 ~48만 = 4.1배): 8코어 ~2일
  (R1: DIAMOND ~11h + per-OG 12–15h + 프로파일 5–8h), **32코어 ~10–14시간**.
  반복 실험은 `-b` BLAST 재사용(#10)으로 R1 DIAMOND 비용 회피. 15sp 패널엔 CgigH 원래 없음.
- README 갱신(eeec64c): JSON 예시에 v2 키 + "Key v2 parameters" 표(실측 근거 포함).
- **PEPC clan 병합 실행 완료** (round-4 실행 항목 4번): 15sp fragment 정확히 5개 →
  단일 PEPC_clan(서브패밀리 라벨 유지). **Ppc1 쌍 TOGETHER ✔ + Ppc2 쌍 TOGETHER ✔**,
  bacterial 침입 0/95 (plant vs bacterial anchor 분리 1,978 bits). flagship 분열 실체:
  15sp에서 Mcry_Mcr8G11630=R2_OG0000399 vs Ococ_21370=R1_OG0000440.
  레코드 gpu `~/pepc_pilot/clan_merge.tsv` (95행) + merge_clan.py.
  꼬리: truncation 16개는 레코드 밖 — v2 조각 정책(epa_min_query_len)과 함께.
- **EPA 엘라스틱 게이트 구현**(454069a): `epa_lwr_aggregate="family"`(질량 합 판정,
  opt-in) + `epa_min_query_len=150`(조각 사전 거부, 기본 on). 테스트 184.
- 멀티세션 협업: cluster-diag(코돈 정렬·hyphy·RELAX·DeepLoc -p·signal_windows) +
  ecforest-runner(ECForest 완주·잔기 재매핑). ec_classifier 다운로드 중복 프로세스 정리.
- **EPA-ng LWR 보정 완료** (PEWO PAC, PEPC clan 30 prunings/189 queries):
  **min_lwr 0.8→0.2, margin 0.3→0.0(죽은 노브)** — 구 0.8은 recall 9% 손실·precision
  이득 0. 오배치 대부분 LWR=1.0 인접 가지(nd≤2), 파국적 오배치는 전원 80–129aa 조각
  → 길이 게이트가 실질 방어선. 코드 반영 0fc9f89+ee97a2d, 데이터 gpu
  `~/pepc_pilot/lwr_calibration.tsv`, 런 `~/pewo_pepc/`. PEWO 함정: snakemake<8 핀,
  --recursive 클론+ant 빌드, 잔여 트리 <4 가드 java 패치.
- **clustering_species_exclude 구현·병합**(0fdcfde) + **5sp v2 재실행 제출**(pronghorn
  6091850, output_5sp_v2/): CgigH 27,583개 제외 확인(풀 116,378), 추정 종 트리 +
  per-round 프로파일 + relative 프루닝. 완주 후 Mcry 미배치율 vs 14.9% 대조.
- **DeepLoc -p 완료**(#24 1단계): per-residue attribution 102서열,
  `~/pepc_pilot/deeploc_p_out/alpha_*.csv` (+ results_20260820-114221.csv).
- **signal_windows 완료 — SF3 특이 N말단 영역 발견**: `signal_windows.tsv`
  (86윈도우/74서열). SF별 공통 컬럼 스캔에서 **SF3만 통과(col 174–177+198–202)**;
  SF3 10멤버 중 9개가 N말단 13–35aa 윈도우(col 165–209 수렴), 7멤버 NLS 예측
  (flagship 포함, aa22–30=col174–202). 단 **K/R-rich 고전 NLS 조성 아님** — 비고전
  NLS 클래스 대조 필요. NES 40/40 소수성 통과. RELAX 이완·발현 74%와 같은 clade.
  남은 것: BEB 교차(branch-site 대기). 상세 #24 코멘트.

---

## 7. /clear 후 세션 재개 가이드 (2026-08-27 기준 — 카탈로그 완주 직전 · 전수 확장 준비 세션)

이 섹션이 최신이다. 아래 7.05(2026-08-26)·7.1(2026-08-24)은 이전 세션 기록으로 보존한다.
⚠️ **모니터·백그라운드 waiter는 /clear와 함께 전부 죽는다** — 재개 즉시 아래 "돌고 있는 것"을 손으로 확인할 것.

### 끝난 것 (전부 이슈/커밋에 기록됨)

- **① Mcry direct 주석 테이블 확보 + PLAZA 게이트 보정 완결** (#47 코멘트):
  정본은 `pronghorn:~/scratch/data/opuntia_coche/revision/cam_results/<종>/iceplant_*.tsv`
  (16 CAM 모듈, 507행/467 유전자, 컬럼 McrID·심볼·설명·EC·AGI·UniProt·IPR·위치; 패널 매칭
  460/467, 결측 7=수작업 모델). 보정 결과: **margin이 conflict 지배**(pident/cov 평평),
  운영셀 25/60/0.1에서 SAFE 178, conflict 10/178=5.62%, Wilson 10.03% → 사전 등록 기준
  (≤2%/5%)에 전 셀 미달 = **(a) 전수 DIAMOND 형식 NO-GO 확정**. 충돌 10건 중 7건은 같은
  family 파랄로그 번호 교체, 교차-family는 3건(≈1.7%) → **어간(family) 전이는 안전, 접미
  (subfamily)는 트리 심판 소관**이 정량 확정. 산출물 `pronghorn:~/scratch/mcry_direct/`.
  버그: `plaza_pilot.py --direct`가 개행을 안 벗겨 conflict 100% 아티팩트 → `97ffc95` 수정.
- **② sda1 복구 종결** (wiki installs.md `2492205`): 포맷 금지(사용자 지시), NTFS 유지로
  `/data/sda1` ntfs3 rw 마운트 + fstab(UUID 6AF1D3192ABDD7FE, nofail). 콜드 294G 이관
  (Deschampsia 174/opuntia 60/eggnog_data 48/temp 12) — rsync 검증→원본 삭제→심링크,
  **루트 99%→82%(여유 322G)**. 심링크는 각 사본의 `SYMLINKS.tsv`로 보존(19/2,439).
  구 NTFS 잔존물은 `/data/sda1/LEGACY_DO_NOT_USE/`로 격리, 신·구 혼합 금지 규칙을
  README+위키에 명시. 추가 회수 후보: jlomas 288G(타 사용자, 조율 필요).
- **subfam_tree.py 결함 2건 수정·실증** (#25 코멘트): pal2nal이 translate 사전 검증(≤max(2,1%))
  보다 엄격 → **stderr ERROR 블록에서 문제 ID 파싱해 격리·재정렬 재시도(≤6회)**,
  `MISSING.txt`에 `pal2nal_rejected` 기록. 실패 시 LOCK 잔류 → `atexit` 해제.
  `SUBFAM_THREADS` env 오버라이드 추가. OG87 백필로 실증(316잎, `Sole_..._575` 정확히
  격리), OG42도 완료. 스크립트는 `pronghorn:~/scratch/subfam_trees/subfam_tree.py`(repo 밖).
- **대형 473 트리 — 사실상 완주**: 3-어레이 협동(전방 6102594 %60 / 역순 6104456 %40 /
  **중앙 s1 6104681 %16** — cpu-48/49 idle 활용, 57→84 슬롯). /clear 시점 **471/473**,
  잔여 `R1_OG0000383`(UFBoot ~29분)·`R1_OG0000403`(~42분). gpu 증분 **967/1,014 축 완료**
  (전송→convert_supports→run_subfam_axes 루프, 실패 0).
- **이슈 라운드** (Codex 트리아지 + 내 gh 분업 — Codex 샌드박스 gh 토큰 무효라 로컬 판정
  파일 산출 후 게시만 내가): **#34 닫음**(tier-3 판정 완결). #44 → **Sund 복구 1건으로 압축**
  (fast46 완주·HyPhy 채택 판정, methods PAML 4.10.10/HyPhy 2.5.100 행 `1cce79d`).
  #47 사유 해소(results_15sp.md 최상단 정본 블록 `164d2cb` — arbiter_v3 13,431이 정본,
  20,133/20,726 인용 금지). #38 사유 해소(methods §2.X.8 인트론=주석 프로그램 교란 명문화
  `afa33a1`, cis-element만 남음). **#43 실측 3건**: 15sp Cgig=Helixer(27,583) 확인 ·
  🔴 **"Ppc2" 라벨 판정 (b)** — 12sp `R1_OG0000168`(Ococ 17, Chr3 tandem, `Ccac_g11025`)은
  교정 clan 트리의 **ppc-1E1(Ppc1형)**이고, 진짜 ppc-1E2(`OcoChr10G09070`·`Ccac_g9957`)는
  `R1_OG0014645`에 있다 → **원고의 "Cactaceae Ppc2 tandem 확장"은 Ppc1형 확장** (개명 권고) ·
  더미 종트리 서술은 이미 정정돼 있었음(잔여: (Tfru,Ccac) 의존 하류 목록화+추정트리 채택 결정).

### 돌고 있는 것 (/clear 시점)

- 트리 마지막 2개(383·403) — 러닝 태스크 3개가 LOCK 잡고 마무리 중. 확인:
  `ssh pronghorn 'cd ~/scratch/subfam_trees; d=0; for f in $(cut -f1 slate2.txt); do [ -e "$f/DONE" ] && d=$((d+1)); done; echo $d/473'`
- 실패 감시 기준선: FAILED 고유줄 **4** (OG42/OG87 구건, 전부 해소됨). 5 이상이면 신규.

### ✅ 재개 계약 1–5 완료 (2026-08-27 오후) — 현재 상태

- **473/473 완주** (383·403은 UFBoot 미수렴 연장 끝에 완료, FAILED 신규 0) → 증분 전송 →
  gpu **1,014/1,014 축 완료** (신규 47 전부 AXES_OK).
- **freeze_catalog.py Codex 2라운드 검토·수정** (#25 코멘트 2건): 8건 수정(매니페스트
  필수화·원자적 쓰기·HOG 관측 충돌 강등 등) + 수정이 만든 거부권 순서 회귀 1건 재수정.
  등급 임계값 불변, v2 대비 멤버십 0·HIGH 0·등급 이동 정확히 80(전부 규칙 교정).
  경화판은 `gpu:~/subfam/freeze_catalog.py` (구판 `freeze_catalog_v2impl.py.bak`).
- **v3 동결**: `gpu:~/subfam/subfamily_catalog_v3.tsv` — 1,014/1,014 universe 통과,
  14,580 OG. **HIGH 2,767(19.0%) / PROVISIONAL 1,294 / UNRESOLVED 10,515 / NOT_EVALUATED 4.**
  관전점 답: 대형 473 HIGH율 19.0% = 소형 18.9% **동률**. 행 없는 family 2
  (`R1_OG0000202`·`R1_OG0003665`, Possvm 전부 싱글턴). 보고는 #25.
- **slate3 발사됨**: 6111378(fwd %60)·6111379(rev %40)·6111380(s1 %16)·6111381(big6)·
  6111382(hog_full3 %40) — 12,411+6 family, 제출 직후 PENDING. 완주 시 v4 동결(버전 분리).

### ✅ v4 동결 + iceplant 라벨 레이어 v1 완결 (2026-08-29)

- **big6 6/6 완주** (R1_OG0003670은 codon123에만 19.5h — 최난물, 실패 0) → 증분·축 →
  **v4 동결**: `gpu:~/subfam/subfamily_catalog_v4.tsv` — 13,431/13,431, 36,150 OG,
  플래그 32 = 트리 적용불가 목록과 정확 일치. **v3 대비 공통 1,014에서 등급·멤버십
  diff 0** (결정론 실증). 코호트: v3 밴드 HIGH 19.0% / slate3 12.5%(HOG 적용불가가
  PROVISIONAL 14.0%를 만듦) / big6 2.7%(3/110 — >500유전자는 codon12 재현 거의 불가,
  30–500 밴드가 sweet spot이라는 사후 검증). 행 없는 family 1,101 = Possvm 전부-싱글턴.
  버전 체계: v2(구현 구판 기준선)·v3(사전 등록 밴드 정본)·v4(전 커버리지). 보고 #25.
- **iceplant 매핑 (사용자 지시 08-28) 완결**: `gpu:~/subfam/labels_mcry_direct_v1.tsv`
  — 큐레이션 457/467 배치, family 366 라벨(stem_auto, 모듈 16/16), OG 294
  suffix_needs_tree_gate(자동 전이 금지). **PEPC 검증: Ppc1+Ppc2 → R1_OG0000440(병합
  PTPC 111), Ppc4 → R1_OG0009826(BTPC 분리)** — 채택 구조 합치. 다심볼 family 17개는
  검토 목록. 조인 함정: summary_v3.tsv는 7컬럼 — 멤버는 `gene_list`(5번째), p[-1]은
  cluster_id.
- 남은 선택지: ATH 앵커·AFDB protein_name 축 추가(#33 연동), 다심볼 17 검토.

### (구) v4 진행 상태 (2026-08-28 아침 — big6만 남음)

> ❌ **SUPERSEDED (2026-08-29):** 위 블록으로 대체 — v4 동결·매핑까지 전부 완료.

- **slate3 완전 해결 (12,411 수지 확정)**: DONE **12,379** + **적용불가 32**
  (`slate3_inapplicable.txt`, #25 기록). 32건은 fwd·rev 각 2회 시도 후 동일 에러 —
  28건 "bootstrap with less than 4 sequences"(중복 제거 후 고유 서열 <4, 예:
  R1_OG0003873은 Sof 동일 파랄로그 다수·정렬 6nt) + 4건 "at least 3 sequences".
  UFBoot이 원리적으로 미정의 → 재실행 대상 아님, freeze가 기술 결손 흡수.
- **최종 증분 완료**: 16개 추가 전송·축 → gpu `axes.done` **13,393** = 1,014 + 12,379
  전체. v4에 남은 입력은 big6 6개뿐.
- **big6 진행 — 4/6 DONE (08-28 저녁)**: R2_OG0000000·R1_OG0000002·0000003·0000004
  완료. 잔여 2: **R1_OG0003670**(codon123에 19.5h 소요 후 codon12 진입, 최대 난물)
  + **R1_OG0000005**(진행 중). s1 원본(6111381)·S2 클론(6120178) 양 어레이 LOCK 협동
  정상, 실패 0. 밤사이 6/6 유력.
- ⚠️ 감시 함정 추가: slate3 실패 감시 grep이 `t3*.out`(없는 경로)을 봐서 **실패 0
  오보** — 실제는 `logs/t3f_*.log`의 `FAILED` 마커. 마커 문자열·로그 경로는 sbatch
  원문에서 확인할 것.
- **HOG 축 완결**: 산출 12,642 + **적용불가 790**(1종 501/2종 286/pep부재 2 — 종간
  orthology 미정의라 OrthoFinder가 exit 0·테이블 없이 종료, §7.05 침묵사망 계열.
  재실행 대상 아님, freeze가 기술 결손 처리) = 13,431 전체 커버, 갭 0 (#25 기록).
- **gpu 증분 완료**: 트리 12,363 전송 + convert_supports + **8-way 병렬 축**
  (`run_axes_list.sh`, /tmp/axes_chunk_*) 12,363/12,363 실패 0 — 직렬 10–20h를 ~2h로.
  HOG N_root 전송 완료. ⚠️ 함정: 감시에서 `pgrep -fc run_axes_list.sh`가 **감시 자신의
  ssh 명령줄을 매칭**해 영원히 2를 반환 (위키 pgrep -f 자기오살 규칙 그대로).
- **v4 스모크 통과**: `manifest_13431.txt`(= 1,014+12,411+6, arbiter_v3 정본 수와 일치)
  → 13,377 처리 / 플래그 54 = 미전송분과 일치. 인터림 HIGH 5,440 / PROV 4,286 /
  UNRES 26,056.
- **big6 재배선**: s1 원본 6111381 예상 시작 08-30(Priority 대기) → **S2 클론 6120178
  즉시 가동**(`subfam_big6_s2.sbatch`, 태스크 0·1 cpu-76/78 RUNNING, %2). S2 상한도
  14일이라 walltime 문제 없음. s1 원본은 LOCK/DONE 협동 백업으로 유지.
- **남은 것**: big6 완주 → 최종 증분+축 → `freeze_catalog.py --manifest
  manifest_13431.txt --out ~/subfam/subfamily_catalog_v4.tsv` 정식 동결 → 코호트 분포
  보고(#25).
- **v4 후속 계약 (사용자 지시 2026-08-28)**: v4 동결 후 **iceplant(Mcry) gene ID 매핑**
  — Mcry direct 큐레이션 테이블(467 유전자, `opuntia_coche/revision/cam_results/`)을
  v4 카탈로그의 Mcry 멤버십으로 조인해 family/subfamily 라벨 레이어(버전드, 멤버십
  불변)를 만든다. 게이트 보정 결과대로 **어간(family) 전이는 자동**(오류 1.7%),
  **접미(subfamily)는 `anchor_transferability()` 트리 게이트 경유**(conflict 5.6%).
- **세대 구획 완결 (2026-08-28)**: output_14sp/17sp → RETIRED 이관, 전 분석 디렉터리에
  `_CANON_STATUS.txt`(CANON/SUPERSEDED/MIXED) 15곳 배치, 규칙은
  `ANALYSIS_PARTITIONS.txt` + retired_data.md 2026-08-28 절.
- #43 원고 결정 기록(§ 사용자 입력 대기 참조), 격리 마커 규약은 위키
  `guide/handoff-hygiene.md`로 병합·이 문서에 적용.
- **이슈 정리 (08-28)**: **#47 닫음**(완결 경로 3건 소진 — arbiter_v3 정본·게이트 보정
  NO-GO·codeml→HyPhy/#44 이관) · **#43 닫음**(저장소 측 판정 3건 완결, 원고 반영은
  기록된 계약대로 추후 — Cgig 이원 운용이 공식 결정으로 명문화됨). 잔여 오픈:
  #44(Sund 복구) · #38(cis-element) · #33(MYB 본작업) · #26 · #25.

### (구) 재개 계약 — 473/473 확인 후 순서 (사용자 승인 완료: "전체 진행" + "완주 후 나머지 전체")

> ❌ **SUPERSEDED (2026-08-27):** 아래 1–5는 위 블록대로 전부 실행 완료. 6(v4 동결)만 남음.

1. 최종 증분: `ssh pronghorn 'cd ~/scratch/subfam_trees && for f in $(cut -f1 slate2.txt); do [ -e "$f/DONE" ] && printf "%s/codon123.treefile\n%s/codon12.treefile\n%s/DONE\n" "$f" "$f" "$f"; done | tar cf - -T -' | ssh gpu 'tar xf - -C ~/subfam/subfam_trees && cd ~/subfam && python3 convert_supports.py && bash run_subfam_axes.sh'`
   → gpu 1,014/1,014 확인
2. **Codex 검토**(사용자 지시): `freeze_catalog.py` 구현이 #25 사전 등록 등급 규칙과 일치하는지 + 입력 무결성 — 문제 있으면 먼저 수정
3. `freeze_catalog.py` 출력명 `subfamily_catalog_v3.tsv`로 **1,014 전체 동결** (등급 규칙 변경 금지)
4. 분포 보고(v2 대비 **대형 family HIGH율 변화**가 관전점) + #25 기록
5. **slate3 발사** (전부 작성·검증 완료, 제출만 하면 됨):
   `sbatch ~/scratch/subfam_trees/subfam_array3_fwd.sbatch` (S2 %60) + `..._rev.sbatch` (S2 %40)
   + `..._s1.sbatch` (s1 %16, 13일) + `subfam_big6.sbatch` (s1 16코어/96g, >500 6개)
   + `hog_full3.sbatch` (S2 %40, HOG 축 12,417). slate3.txt 12,411 + slate3_big.txt 6,
   크기 내림차순. 이후 같은 증분 루프(모니터 재구성 필요 — 잡명 subfam_trees3*, 로그 t3f/t3r/t3s/t3big).
   48h 한도 시 S2 어레이 재제출이 곧 재개(DONE/LOCK 스킵).
6. slate3 완주 시 **전 커버리지 확장판(v4) 동결** — v3(1,014)와 버전 분리 유지(v2↔v3 비교 보존)

### 이번 세션의 신규 함정

- **pronghorn `~/scratch/bin/ff_arbiter`는 GitHub 클론이 아니다** — origin이 로컬
  `~/scratch/bin/family_finder`(구판 0fdcfde)라 `git pull`이 "Already up-to-date"를 내며
  **아무것도 안 가져온다.** 수정 반영은 직접 패치/scp.
- **family ID는 런 스코프다**: `R1_OG0000168`을 11sp 표에서 찾으면 **다른 family**(Chr7/8/11)가
  나온다. 원고가 인용하는 12sp ID는 반드시 `RETIRED_DO_NOT_USE/output_12sp_portulacineae/`에서.
- **plaza_pilot `--direct` 값 개행** → conflict 100% 아티팩트 (`97ffc95`+회귀 테스트).
- **sda1은 NTFS**: 심볼릭 링크 못 실음(SYMLINKS.tsv 매니페스트로 보존), conda env는 아카이브
  전용(실행 불가). ntfs3 마운트 옵션은 fstab 참조.
- **iqtree2는 GPU 옵션이 없다**(전 버전, ML 부트스트랩 계열 공통) — 가속 축은 파티션 폭
  (s1 idle 노드). BEAGLE류는 지지도 통계가 달라 등급 규칙과 비교 불가.
- Codex 샌드박스 gh 토큰 무효 → **로컬 트리아지 파일 + 내 gh 게시** 분업이 작동 패턴.

### 사용자 입력 대기

- **앵커 라벨 레이어 착수 판단(구 ③)** — 선행 조건(게이트 보정) 이번에 충족됨. 어간 전이
  1.7% 오류/접미 5.6%라는 보정 결과가 설계 근거.
- ✅ ~~#43 원고 측 결정 2건~~ → **결정됨 (2026-08-27, #43 코멘트 게시)**: 원고 정정(Chr3
  tandem Ppc1형 개명 + 추정 종트리 채택)은 **분석 전부 완료 후 착수** — 473 트리 완주 →
  카탈로그 v3 동결 → slate3 전수 확장(v4)까지. 그 전에는 이슈에 사실만 축적, 원고 파일
  불가침. 착수 시 OG id는 최종 summary 기준으로 재확인(`0000168` vs `0000165` 이력 참조).
- 5sp Cgig MAKER 유지의 공식 결정 기록(#43 ①, 현행 이원 운용이 사실상 결정)

## 7.05 (구) /clear 후 세션 재개 가이드 (2026-08-26 기준 — v3 테이블·subfamily 카탈로그 세션)

> ❌ **SUPERSEDED (2026-08-27):** 2026-08-26 스냅샷, 기록 보존용. 이 절의 재개 계약·
> 진행 수치·대기 항목을 실행하지 말 것 — 정본은 §7.

### 끝난 것 (전부 이슈에 기록됨 — #47, #25, #34, #38)

- **v3 family 테이블**: uncut `vote_edges.tsv`(50,961 엣지, 21,260/23,744 family 커버)
  → `nominate_clusters`(54760a9; 클러스터 2,774·pairwise 12,394·deferred 221, cap 500)
  → 심판 배치(SLURM 6101783, 2,774/2,774 실패 0) → `validate_families --apply` =
  **23,744 → 13,431 families**, 유전자 전량 보존. **PEPC 5조각+양 flagship이 111유전자
  단일 family, PPC4(BTPC)는 분리 유지** — 트리 심판이 서열·구조가 못 긋던 경계를 그었다.
  산출: pronghorn `~/scratch/arbiter_v3/{summary_v3.tsv,cluster_verdicts.tsv}`.
- **구조 채널 판정**: 15종 전수 ProstT5 3Di(46청크 병합 DB, gpu `~/prostt5_15sp/merged/all_3di`)
  → foldseek cluster 153초 → 후보 10,070(`~/prostt5_15sp/merge_out/structural_merge_*.tsv`).
  사전 등록 rubric 판정(#47): **OrthoFinder 대체 기각 · Tier-0 기각 · 지명/주석 전용 채택**
  (PTPC/BTPC co-cluster + 최대 1,165 초과).
- **PLAZA/AFDB 파일럿 판정**(#47): (b) ARATH tar 다운로드 **NO-GO**(addressable 12.2%<15%),
  comparable 고아의 구조 명명 전이 **0/1,387**('고아는 iteration으로' 확정), SAFE orthology
  22.7%(QCOV 지배 — 게이트 보정 전 전수 DIAMOND 무의미). RBH 버그 수정 47800ce
  (역방향 맵이 qseqid 저장 → SAFE=0 전멸했었다).
- **subfamily 카탈로그 v2** (gpu `~/subfam/subfamily_catalog_v2.tsv`): 소형 541 family(≤200
  유전자)의 IQ-TREE 쌍(codon123+codon12, `-alrt 1000 -B 1000`) → Possvm×2 + TreeCluster
  → HOG(파일럿 20 → 전수 1,014) 5중 검증 = **HIGH 690 / PROVISIONAL 390 / UNRESOLVED 2,575**
  (전 건 사유 기록, 무앵커 — 라벨은 후속 버전드 레이어). 등급 규칙·판정 임계는 #25에 사전
  등록돼 있고 결과보다 먼저 커밋됐다.
- **주석 캐시 3축 전수** (gpu `~/annot_panel/<axis>/<종>/`): SignalP(GPU 변환, 전수 ~30분)
  · emapper(~1시간) · DeepLoc Accurate(~7시간). DeepLoc 캐시는 neofunctionalization 모듈
  입력 겸용.
- 코드 커밋 2c950c9..47800ce (전부 push): nominate_clusters(+테스트), structural merge
  candidates, plaza_pilot(게이트 계약·abstention·사전 결정규칙), prostt5 GPU 노브
  (`config.prostt5_gpu`), subfamily/naming 리뷰 수정 2라운드, manifest 신필드 호환.
  테스트 1,075+ green.

### 돌고 있는 것 (/clear 시점)

- **대형 473 family(201–500 유전자) 트리 배치**: 어레이 6102594(%60) + 역순 보조 6104456(%40),
  family당 `LOCK/`(mkdir 원자 락)으로 협동. 진행 ~572/1,014, 실패 1
  (`R1_OG0000042` — pal2nal 불일치, 디렉터리 정리됨, 패치본으로 백필 대상).
  48h 한도 도달 시 **같은 sbatch 재제출이 곧 재개**(DONE 스킵).
- HOG 축은 **1,014 전수 완료**(pronghorn `~/scratch/hog_pilot/*/N_root.tsv`, 실패 0).

### 재개 계약 — 트리 완주(또는 한도) 신호 후 자동 진행 순서

1. 백필: `sbatch ~/scratch/subfam_trees/subfam_array2.sbatch` 재제출(DONE/LOCK 스킵이 잔여만 잡는다)
2. 신규 treefile → gpu `~/subfam/subfam_trees/` 전송(tar 파이프) + `convert_supports.py`
   (이중 지지도 → min 단일값 `.possvm.nwk`)
3. `~/subfam/run_subfam_axes.sh` (axes.done 스킵 — 신규 family만 돈다)
4. 신규 `N_root.tsv` → gpu `~/subfam/hog/<fam>/` 전송
5. `freeze_catalog.py` — 출력 파일명을 `subfamily_catalog_v3.tsv`로 바꿔 1,014 전체 동결
   (등급 규칙 변경 금지 — 사전 등록분 그대로)
6. 분포 보고(v2 대비 대형 family의 HIGH율 변화가 관전점) + #25 기록

### 이번 세션의 신규 함정 (코드/위키에 반영됨)

- **pal2nal은 all-or-nothing**: 불일치 유전자 1개가 family 전체를 죽인다 →
  `subfam_tree.py`가 사전 translate 검증으로 `cds_mismatch`만 격리 (그리고 `-nogap`은
  다서열 정렬에서 전 컬럼을 지운다 — 빈 서열 dict는 truthy)
- **IQ-TREE 이중 지지도(`96.9/100`)를 ete3가 파싱 불가** → min(aLRT,UFBoot) 변환본 사용
- **단일 family OrthoFinder는 종 트리 추론이 exit 0으로 침묵 사망** → `-s`로 루팅 종
  트리 공급(무근이면 거부; Mcry 외군 루팅 + family별 가지치기), 이때 HOG 레벨은 **N1부터**
- **`--gpu 1` 없이는 CUDA foldseek도 CPU ProstT5ForkRunner**로 가서 msgsnd에 영원히
  블록될 수 있다(§5의 시계 판정 원칙 유지) — installs.md 2026-08-25
- **`vote_rescan/members.fa`(459,398)는 placed 유전자만**이다. 패널 총계는 484,752
  (+unplaced 25,354). '전체 패널'로 읽지 말 것(#47에 종결 기록)
- SignalP 6도 GPU가 된다(env torch가 +cpu였을 뿐) — 변환 레시피 installs.md

### 사용자 입력 대기 (자동 불가)

- **큐레이션 direct 주석 테이블**(Mcry 기능 주석) 위치 — conflict 축 정식 계산 + PLAZA
  게이트 보정 + 라벨 레이어의 선행 조건
- gpu 박스 `~/scratch` 7.3TB(sda1) 복구 마운트 — sudo 필요, 안전 절차 installs.md
- 앵커 라벨 레이어 착수 판단 (멤버십은 라벨과 무관하게 이미 동결)

## 7.1 (구) /clear 후 세션 재개 가이드 (2026-08-24 밤 기준)

> ❌ **SUPERSEDED (2026-08-26):** 2026-08-24 스냅샷, 기록 보존용 — 정본은 §7. 이 절의
> "돌고 있는 것"(vote_rescan 등)은 전부 종료됐고, family 수 20,133/20,726 논의는
> arbiter_v3 **13,431**로 종결됐다(`results_15sp.md` 정본 블록). 아래 "정본" 줄은 당시 기준이다.

**정본**: 이 파일 · `results_15sp.md` · `AGENTS.md` · `retired_data.md` · #47 · #26.
⚠️ **모니터·백그라운드 waiter는 /clear와 함께 전부 죽는다** — 손으로 확인할 것.

### 지금 돌고 있는 것 — 한 갈래

```bash
# 전수 재스캔 6099988 (S2, 24×8코어) — 459,398 유전자 × 23,744 프로파일, cov 0.2
ssh pronghorn 'ls ~/scratch/vote_rescan/dom/*.domtblout | wc -l'   # /clear 시점 23/24
# 24개가 차면 build_edges.py 가 자동 실행된다 (run_when_done.sh, setsid로 분리)
ssh pronghorn 'cat ~/scratch/vote_rescan/edges_run.log'
```

`build_edges.py`가 내는 것: uncut 간선 + **컷별 클러스터 프로파일**(0.0/0.2/0.4/0.6의
간선 수·클러스터 수·최대 성분·커버리지). **여기서 임계값을 골라야 한다** — 커버리지 단독은
목적함수가 아니다(cut 0이면 거대 성분 하나로 100%가 된다). 최대 성분이 판정 비용을
지배한다(119잎 클러스터 40–53초).

**확인할 것**: `R2_OG0000359`(Mcry flagship, frac 2/9=0.22)가 **cut 0.2에서 들어오는가.**
그게 재스캔이 메우려던 커버리지 구멍의 직접 증거다.

### 이 세션에서 뒤집힌 것 — 먼저 읽을 것

1. **🔴 무근 트리 결함 (#51)** — `fragment_verdict`가 rooted 질문을 했다. 캠페인 트리
   5,618개 재판정 시 **3,536개(62.9%) 판정이 바뀌고 전부 한 방향**(INTERLEAVED →
   MONOPHYLETIC). 옛 규칙은 **과병합만** 할 수 있었다. **두 번 고쳤다**: `5231b2f`가 split
   양쪽을 보게 했고(status는 불변이 됨), `ca2507f`가 "가장 작은 포함 split"이 여전히 루트
   의존임을 잡아 pairwise 규칙으로 바꿨다. family **20,133 → 20,726**.
2. **🔴 커버리지 구멍** — 트리를 본 family는 **16,061/23,744 = 67.6%**뿐. 나머지 7,683개는
   "별개로 판정된" 게 아니라 **심판대에 오르지 못했다**(frac 0.600 컷).
3. **`merge_min_profile_cov` 0.5 → 0.2** (`5c50a3e`) — 발행된 `vote_edges.tsv`가 실은 0.2로
   만들어졌다(스윕 191/200 vs 0.5의 84/200). 게이트는 표를 **버리는 게 아니라 돌린다**.
4. **20,133도 20,726도 확정이 아니다** — 두 결함이 **반대 방향**이라 스칼라로 상쇄하면
   안 된다. `results_15sp.md`의 절차를 따를 것(원본 23,744에서 새로, 델타 덧붙이기 금지).

### PEPC — 목표는 달성됐고 형태가 바뀌었다

- **파편화 확정**: clan이 **6조각 128유전자**. 15종 **전부가 3–5조각에 분산** → 파편화 축이
  계통이 아니라 **서브패밀리**다. 서브패밀리는 정의상 단계통이라 **트리 규칙이 원리적으로
  병합 못 한다**(6조각 중 2개만).
- **외군으로 해결**: PPC4(BTPC)를 외군으로 주자 5개 PTPC 조각이 `SAME_FAMILY` →
  **111유전자 한 family, 40초**. **두 CAM flagship이 처음으로 같은 family**
  (`Ococ_OcoChr03G21370` + `Mcry_Mcr8G11630`).
- **외군 정당화는 15sp만으로 된다**: 128유전자를 ATH 앵커 4개에 점수 매기면 **부호가 딱 한 번**
  뒤집힌다(5조각 PPC1/2/3 쪽 700–1,000 bits, `R1_OG0009826` PPC4 쪽 812).
- **채택된 구조 (사용자 결정)**: PTPC/BTPC를 **가르지 않고 한 family로 묶고 내부에서 나눈다.**
  128잎 트리에서 **간선 하나로 깨끗이 분리**됨을 확인(둘 다 단일 clade). 그러면 경계 판정
  축을 찾는 문제가 사라지고 **커버리지 문제만** 남는다.

```
family     PEPC (6조각 128)
  ├ class  PTPC (111)   ← 트리 간선 하나
  │   └ subfamily  1E1 / 1E2 / 2 / 3   ← Possvm + TreeCluster + HOG
  └ class  BTPC (17)
```

### 기각된 축 — 재제안 금지 (전부 실측)

| 축 | 기각 근거 |
|---|---|
| **Pfam 도메인** | PPC1/2/3/4가 **아키텍처 동일**(PEPcase + PEPcase_2). family 경계는 흔히 도메인 층위 **아래**에 있다 |
| **구조 — AFDB TM** | PTPC−BTPC 중앙값 **+0.033**(최대 0.079). 절단 모델에 오염(411 aa ppc-1E2가 BTPC에 TM 0.943) |
| **구조 — ProstT5 3Di** | within-PTPC 4,290 vs cross 3,882, **비 1.11**. 독립 두 표현이 같은 결론 |
| **Ks/Ka 거리로 경계** | 단일카피는 종분화에서 갈라지고 family는 더 오래된 파랄로그를 담는다 → **하한만** 준다. Ka/Ks는 전 쌍 0.311로 평탄(거리가 아님) |
| **아래-컷 간선에서 외군** | 클러스터가 **연결성분**이라 나가는 간선이 없다. C0297의 4간선 전부 내부 |
| **HOG를 family 계층으로** | OG 28,150 → N0 HOG 28,161, **여러 OG를 담은 HOG 0건**. 세분만 하고 라운드별이다 |
| **`prune.py` 통계량 승격** | family별 중앙값 정규화 + **동종 파랄로그 스킵** → 척도가 family 내부에서만 의미 |

### 남은 축 하나 — "general하게 이용"

**참조 DB 라벨을 자동으로 읽는 것.** 사람이 "PPC4는 BTPC"임을 안 게 아니라
**AFDB 헤더에 적혀 있었다**(`Phosphoenolpyruvate carboxylase 4`). 파일을 읽는 건
큐레이션이 아니라 자동화다. 그리고 **명명 규약이 계층을 담는다** — 어간이 family,
접미 숫자가 서브패밀리.

- 1·2단계(AFDB 검색 → UniProt `protein_name` 전이)는 **`fs_transfer.py`에 이미 있다**
  (PEPC clan 100개 중 99개 EC 4.1.1.31)
- 3·4단계(어간=family, 접미=subfamily)만 새로 필요
- **상한은 참조 커버리지**이고 잴 수 있다 — `~/scratch/vote_rescan/family_reps.py`가
  family당 대표(가장 긴 멤버)를 뽑아 둔다. AFDB 검색해 `protein_name` 커버리지를 보면
  이 방법의 가부가 결정된다

### ⚠️ foldseek — 로그가 아니라 시계로 판정할 것

~~`--gpu`는 파싱되지만 ProstT5 경로에 닿지 않는다~~ — **2026-08-25 정정: 재빌드 바이너리에서는
`--gpu 1`이 경로를 결정한다.** 플래그 없이는 CUDA 빌드도 CPU ProstT5ForkRunner로 들어가고, 그
러너는 스레드 수만큼 fork해 **워커마다 ~3GB 모델을 로드**한다(62GB 박스 × 16 = 워커 전멸,
부모는 가득 찬 SysV 큐에 msgsnd로 영원히 블록 — CPU 0초, 에러 없음, `Use GPU 0`은 그대로 찍힘).
진단: `ps` TIME 0초 + wchan `do_msgsnd` + GPU 0%; 정리: `ipcrm -q`. `--gpu 1` 실측
**0.034 s/서열**(10K 청크 340초). `Use GPU 0`은 **다른 mmseqs 노브**이고 실제 구동은 ggml의
`CUDA0` 백엔드다. **이 로그를 믿고 다섯 번 오진했다** — 시계 판정 원칙은 그대로다.

```
112서열, 동일 커밋 18478813, 동일 심볼(39), 둘 다 "Use GPU 0"
  기존 ~/bin/foldseek-cuda        497초
  -DENABLE_CUDA=1 재빌드           12초    ← 41배
```

재빌드본: `~/scratch/bin/foldseek_src/build/src/foldseek`.
전수 459K: **23일 → 실측 약 4.4시간**(config.prostt5_gpu, 5ab7f4c). 재빌드 함정 2건(`rust` 누락, CUDA 11.5/12.6 혼재)은
`steps/prostt5_chunks.py` 상단.

### codeml vs HyPhy — 충돌, HyPhy 채택

fast46(40 taxa, 이름과 달리 46 아님) 완주: LRT **27.11**, p 9.6e-08, BEB 5사이트.
**실행은 유효**(`#1`이 depth-10 내부 노드, 잎 10개가 fgA와 일치). 그런데 MEME은 FDR 통과 0,
aBSREL은 stem 비유의, RELAX는 **이완**(p 5.9e-13). branch-site Model A의 ω₂≥1이 이완을
흡수한 것으로 읽고 **HyPhy를 채택**. **BEB 사이트는 인용하지 않는다.**
gapped 6체인(6097933-38)은 정책 위반이라 scancel.

### 닫은 이슈 5건 / 남은 것

닫음: **#48**(resume 유전자 소실 3경로) · **#35**(15sp v2 프로덕션 런) · **#49**(실패가
completed) · **#51**(무근 판정 등 3건) · **#50**(캠페인 3건).
남음: #47 #44 #43 #38 #34 #33 #26 #25.

### 테스트 881 → 1,019, 커밋 23건

새 모듈: `validate_families.py`(캠페인 CLI) · `steps/hog.py`(서브패밀리 축) ·
`steps/domain.py`(Pfam 음성) · `measure_ks.py` · `data/ks_baseline_15sp.tsv`(105 종쌍).

### 이 세션의 방법 교훈

**재보지 않고 단언해서 다섯 번 틀렸다** — `Dcar`를 약어로 종 추정(실제 *Dianthus*,
종트리가 이미 반증하고 있었다), Ks 포화 강도 두 번, ESMFold chunk_size 미설정 가정
(이미 설정돼 있었다), foldseek GPU를 로그로 판정. **AGENTS.md의 "조용한 실패" 규칙이
도구 사용에도 그대로 적용된다.**

## 8. 재현 정보

- 5종 런: `python family_finder.py --protein-dir data/pep --cds-dir data/cds
  --species-tree data/species_tree.nwk --outdir output_5sp --config config_5sp.json
  --threads 8 --verbose` (conda env `orthofinder`)
- 클러스터 산출물: `/data/gpfs/assoc/pgl/bin/family_finder/` 또는
  `~/scratch/bin/family_finder/output_5sp/`
- 로컬 참고 산출물: `/mnt/f/scratch/develop/opuntia_analysis/family_finder/output_15sp/`
