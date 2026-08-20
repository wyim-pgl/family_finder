# resume.md — Mcry(외군) 클러스터링 문제 조사 기록

작성일: 2026-08-19. 이 문서는 코드 조사 + 독립 검증(스트레스 테스트)을 거친 현재 상태 스냅샷이다.
모든 `file:line` 은 커밋 `f1226aa` 기준으로 원문 대조 확인했다.

## 1. 문제 정의

5종 런(Mcry, Cgig, CgigH, Ococ, Obas; 143,961 단백질)에서 외군 *Mesembryanthemum
crystallinum*(Mcry, Aizoaceae) 유전자가 **Round 1 OrthoFinder 클러스터링 단계에서부터**
진짜 cactus ortholog와 다른 orthogroup으로 배정된다. HMMER rescue 구제율이 낮은 것은
증상이지 원인이 아니다.

핵심 결론: **원인은 `steps/prune.py` 상류(클러스터링 · 크기 게이트 · 데이터 정합성)에 있을
가능성이 높고, 어느 것인지는 아직 측정되지 않았다.** 결정적 측정은 대부분 클러스터에 이미
있는 산출물로 가능하며 파이프라인 재실행이 필요 없다.

## 2. 증상 지표

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

**유력 가설 (측정 필요):**
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

## 6. 다음 단계

GitHub 이슈로 등록됨 (`wyim-pgl/family_finder`, 2026-08-19):

**P0 진단 (이것부터 — 코드 수정 없음, 대부분 기존 산출물로 즉시 가능):**
- [#5](https://github.com/wyim-pgl/family_finder/issues/5) pep/CDS ID 정합성 + within-OG dropout 감사 (최유력 원인, 몇 분 소요, 여기서 종료될 수 있음)
- [#6](https://github.com/wyim-pgl/family_finder/issues/6) Round-1 OrthoFinder 종별 통계 (디스크에 이미 있음)
- [#7](https://github.com/wyim-pgl/family_finder/issues/7) BUSCO Outgroup Join Rate 벤치마크 + 실패 5분류
- [#8](https://github.com/wyim-pgl/family_finder/issues/8) DIAMOND/MCL 그래프 포렌식 — 실패 원인 3분할 (수정 방향 결정)
- [#9](https://github.com/wyim-pgl/family_finder/issues/9) 실제 프루닝 발생률 측정

**P1 클러스터링 (측정 결과에 따라 채택):**
- [#10](https://github.com/wyim-pgl/family_finder/issues/10) OrthoFinder -og/-I/-S 타입 필드 + -b 재사용 (인플레이션 스윕 {1.05,1.1,1.15,1.2,1.3})
- [#11](https://github.com/wyim-pgl/family_finder/issues/11) 크기 게이트 분리, OG 전체 해체 제거, max_recycles
- [#12](https://github.com/wyim-pgl/family_finder/issues/12) CgigH 제외 실험 + 15종 taxon sampling 비교 (정보량 최대 실험)
- [#13](https://github.com/wyim-pgl/family_finder/issues/13) HMM 프로파일 배정 1급 승격 + 재배정 + family 병합

**P2:**
- [#14](https://github.com/wyim-pgl/family_finder/issues/14) 프루닝 보정 기준 + 종 트리 데이터 추정 + ASTRAL 문서 모순
- [#15](https://github.com/wyim-pgl/family_finder/issues/15) 조용한 데이터 손실·resume·TreeShrink 로그·종명 오기

**진행 상황 (2026-08-19):**
- #7 닫힘 (BUSCO 벤치마크 — 유지자 결정으로 제외)
- **#8 실행 완료** — HMMER 구제 유전자 6,358개를 truth set으로 Round-1 운명 분류.
  판정: **B+C ≈ 90%가 MCL/그래프 파편화** (C: 강한 DIAMOND hit이 있는데도 미배정 43-59%,
  B: 파편 OG 39-47%, 전원 4개 미만 OG). D(민감도) 0.8-6.8%뿐 → `-S` 기각.
  결과: 이슈 #8 코멘트 + 클러스터 `forensics_r8/per_gene.tsv`
- **#10 구현 완료** — `build_orthofinder_cmd()` 분리, `-I`/`-S`/`-og`/`-b` 타입 필드,
  `-b` 모드 rmtree 금지. 테스트 7건 통과 (tests/test_orthofinder_cmd.py)
- **#11 핵심 구현 완료** — `min_family_size=2` 출력 게이트 분리 (`pipeline.py`),
  테스트 4건 통과 (tests/test_size_gate.py). max_recycles는 미구현
- [#16](https://github.com/wyim-pgl/family_finder/issues/16) 신규 — DeepLoc retargeting 가지 +
  branch-site dN/dS로 신기능화 탐지 (codeml.py에 branch-site Model A 추가 필요; #10·#13 선행)

**진단·구현 2차 라운드 (2026-08-19, 서브에이전트 4개):**
- **#5 기각** — pep/CDS ID 완전 일치(5종 모두 pep_only=0), round_01 누수 0. 데이터 정합성은 원인이 아님
- **#6 확정** — Round-1 원시 클러스터링에서 Mcry 최하위: 86.8% 배정(cactus 90.2–95.9%),
  outlier pool이 R05에 14.9–15.0%로 정체(cactus는 계속 빠짐). 원인은 클러스터링 단계로 확정
- **#9 판정** — 1,918개 트리 재채점: 임계값 5.0에서 **Mcry 0개 프루닝**(면책).
  단 **Cgig 10.8%·CgigH 4.9%가 임계 초과** — d_exp(Cgig,CgigH)=0.02 병리가 실측됨.
  #14가 "방법론 방어"에서 "Cgig/CgigH 실질 정확도 수정"으로 격상.
  실제 사용 트리 확인: data/species_tree.nwk = 손으로 쓴 값 맞음
- **#19 신규** — ProcessPoolExecutor 워커 로그가 메인 로그로 전파되지 않아 per-OG 진단이
  전부 유실됨 (--verbose여도 pipeline.log에 워커 줄 0개). QueueHandler 패턴 필요
- **#12 진행** — CgigH 제외 재클러스터링 SLURM job **6088858** 제출됨 (완료 후 채점)
- **#3 차단** — 5종 모두 isoform-level 단백질 없음(Helixer는 설계상 단일 isoform). Sof만 다중
- **#13 구현·병합** — steps/profile_assign.py (+1,179줄, 테스트 24건). per-round 훅은
  `profile_assign_per_round: true` opt-in
- **#14 코드 파트 구현·병합** — relative 프루닝 기준(terminal branch 제외 + family 정규화 +
  이중 조건), validate_species_tree, ASTRAL 문서 수정. 종 트리 추정(클러스터)만 남음
- **#15·#2·#4 구현·병합 후 닫힘** — 조용한 손실 경로 폐쇄, ID agreement 리포트,
  #4 원인 = 깨진 resume가 R1_* 재생성(수정됨)
- **#17·#18 신규** — ESM-2/ESMFold/Foldseek 3차 구조 배정 / 식물 gLM(AlphaGenome은
  human/mouse 전용이라 기각 기록)
- 테스트 73건 통과, master push됨 (ff4002a)

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
- **서브패밀리 top-3 비교 완료**: Possvm(‑r 앵커, min_support 95, bacterial 루팅; numpy<1.24 핀 +
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

## /clear 후 세션 재개 가이드 (2026-08-20 기준)

**정본 문서**: 이 파일(조사 기록) + 이슈 **#25**(파이프라인 스펙) + **#26**(작업 트래커).
서브에이전트·모니터는 세션 종료와 함께 사라짐 — 아래 "확인할 것"부터 시작.

### 즉시 확인할 것
1. ~~ESM-ECForest~~ **완료 (6차 라운드): known-answer 실패로 도구 기각** — #26/#20 기록됨.
2. **branch-site codeml 결과**: fast-track(46-taxa, **BEB 정본**) pronghorn 6091859 alt/
   6091860 null `$FF/pepc_pilot/seltest_fast/` + full(102-taxa, robustness) 6091847/48
   `seltest/`. 완료 시 `parse_lnl`+`lrt_pvalue`(df=1) → alt의 results.txt BEB 섹션
   사이트를 signal_windows 컬럼과 교차(#24 마지막 단계) → #24 닫기.
   full이 내일까지 미수렴이면 kill해도 무방(HyPhy 삼중으로 판정 완료).
3. **5sp v2 런** (6091850, output_5sp_v2/): 완주 시 Mcry 미배치율 vs 14.9% 대조 → #26.
   R1 실측: 클러스터링 39분 + per-OG 1h54 + 프로파일 배정(hmmsearch 장시간 정상).
4. ~~RELAX·MEME·aBSREL~~ **완료 — SF3 순수 subfunctionalization 확정** (#26 기록됨).
5. `gh issue view 26 --repo wyim-pgl/family_finder` — 대기열 순서대로 진행.
   (#18은 구조적 차단으로 닫힘; 오픈 #24/#25/#26)

### 다음 작업 레시피 (#26 대기열 1번: 선택압 검정)
```
# clan 코돈 정렬: 멤버 CDS는 pronghorn data_15sp/cds (+ 5sp data/cds),
#   앵커 CDS는 /data/gpfs/assoc/pgl/data/sequence_data/plaza/cds/selected (ath, aco)
# pal2nal을 gpu:~/pepc_pilot/clan_anchor.aln 에 맞춰 실행 → clan_codon.aln
# SF3 stem(10-leaf flux clade, 지지도 100)에 #1 마킹: steps/codeml.py write_marked_tree
# branch-site alt/null: steps/codeml.py generate_branch_site_ctl + lrt_pvalue
# RELAX: hyphy 설치 필요 (micromamba -n phylo -c bioconda hyphy)
```

### 자산 위치
- **pgl-gpu** (`ssh gpu`): `~/pepc_pilot/` — clan_anchor.{fa,aln}, clan_tree.treefile(+rooted),
  pdb/ (111 ESMFold 구조), plddt.tsv, aln_esm.tsv, aln_p5.tsv, db_p5/db_ath/db_aco,
  possvm_mcl/pepc.ortholog_groups.csv, tc_*.txt, deeploc_out/results_*.csv, anchors.tsv,
  excluded.txt(truncation 16개). envs: phylo(mafft/iqtree3/epa-ng/raxml-ng),
  subfam(possvm deps/treecluster), sonic, transgenic(torch+DeepLoc2.1). foldseek ~/bin/foldseek,
  ProstT5 ~/pepc_pilot/prostt5_weights. SonicParanoid 결과 ~/sonic{5,15}sp/
- **pronghorn** (`ssh pronghorn`, FF=~/scratch/bin/family_finder): output_5sp/, output_15sp/,
  forensics_r8/per_gene.tsv, score_recluster.py, check_pairs.py, cam_pairs.tsv, ath_ppc.fa,
  aco_ppc.fa, species_tree_est/species_tree_busco.treefile(추정 종 트리 — 채택 대기),
  sonic{5,15}sp_results/, wd_I*/wd15_I*/wd_noCgigH(스윕 산출물)
- **로컬 repo**: master 64f2323+ (테스트 165), 신규 모듈 steps/{profile_assign,epa,esm,
  plant_glm,ecforest,deeploc,retargeting,sonicparanoid}.py, utils/newick.py,
  find_neofunctionalization.py, annotate_families.py
- **wiki**: pgl-wiki guide/installs.md 2026-08-20 항목 (설치·함정 전체)

### 핵심 결론 요약 (재론 금지)
Round-1 문제 원인 = MCL 파편화(포함 임계값이 MCL 이전에 edge 절단) — 인플레이션·taxon
sampling·BBH·pep/CDS·프루닝 전부 실측 기각. 해법 = CgigH 제외(9.0%) + 프로파일(#13) +
구조 tier(ProstT5 재현, Ppc1 쌍 TM 0.90+) + EPA 심판. 두 Ppc1은 한 family, 다른 서브패밀리
(SF1 vs SF3). SF3 = 발현 주도 subfunctionalization(74%). AlphaGenome 기각(human/mouse 전용).
ESM-ECForest 기각(known-answer 실패) — EC 축은 도메인/촉매잔기 증거로.
**선택압 축 (HyPhy 삼중, 2026-08-20)**: RELAX 이완 K=0.636 p=1.1e-16 + MEME 사이트 FDR
통과 0 + aBSREL stem 비유의 → **SF3 = 순수 발현 주도 subfunctionalization**.
유일한 양성선택은 Tfru 말단 가지(p=6e-13, 종 특이 후행) — 별도 관찰 가치.
codeml BEB(사이트 목록)만 남음: full(6091847/48) + 46-taxa fast(6091859/60, 컬럼 보존).

## 7. 재현 정보

- 5종 런: `python family_finder.py --protein-dir data/pep --cds-dir data/cds
  --species-tree data/species_tree.nwk --outdir output_5sp --config config_5sp.json
  --threads 8 --verbose` (conda env `orthofinder`)
- 클러스터 산출물: `/data/gpfs/assoc/pgl/bin/family_finder/` 또는
  `~/scratch/bin/family_finder/output_5sp/`
- 로컬 참고 산출물: `/mnt/f/scratch/develop/opuntia_analysis/family_finder/output_15sp/`
