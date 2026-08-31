# 다심볼 family 라벨 자동 판정 설계

상태: 설계안  
범위: 13,431개 정본 family의 **멤버십은 바꾸지 않고**, 그 위의 버전드 라벨
레이어만 판정한다.

## 1. 결정 요약

1. 심볼은 `raw_symbol`, 표시용 `canonical_symbol`, 비교용 `symbol_key`를
   구분한다. 원문을 덮어쓰지 않는다.
2. stem 파서는 규칙의 적용 근거와 상태를 함께 반환한다. 확실하지 않은 문자열은
   `AMBIGUOUS` 또는 `UNPARSED`이며 자동 판정에 쓰지 않는다.
3. 정규화 후 stem이 하나이면 family stem은 `stem_auto`가 될 수 있다. 숫자·allele·
   서브패밀리 접미가 붙은 전체 symbol은 계속 `suffix_needs_tree_gate`이다.
4. 서로 다른 stem은 가중 다수결로 하나를 고르지 않는다. HIGH 서브패밀리 분할과
   라벨 clade가 완전히 정합할 때만 해당 HIGH subfamily 행들을 자동 확정한다. family
   전체를 대표하는 단일 stem은 만들지 않는다.
5. 큐레이션 > ATH > AFDB 우선순위는 같은 stem의 표시 철자와 부가 필드 선택에만
   쓴다. 서로 다른 stem의 충돌을 덮는 우선순위가 아니다.
6. 미파싱, 동점, 낮은 coverage, 미보정 AFDB, PROVISIONAL/UNRESOLVED 분할, 축간
   다-stem 충돌은 모두 `needs_review`라는 양성 상태로 출력한다. 행 부재로 표현하지
   않는다.

이 정책은 stem 전이 오류 1.7%와 suffix 충돌 5.6%의 보정 결과를 그대로 사용한다.
전자는 동일 stem의 family-level 전이를 허용하지만, 후자는 suffix를 트리 없이
전이하지 못하게 한다. 두 수치는 다-stem 충돌에서 승자를 고르는 가중치로 재해석하지
않는다.

## 2. 용어와 증거 단위

- `symbol`: 출처가 제공한 전체 유전자 심볼. 예: `Pth1;3`, `Snrk2.3`.
- `stem`: 동일 family 이름으로 전이 가능한 부분. 예: `Pth`, `Snrk`.
- `suffix`: stem 뒤의 paralog/subfamily/allele 표기 전체. 예: `1;3`, `2.3`.
- `evidence unit`: `(source, gene, normalized stem)` 하나. 같은 유전자의 동일 축에서
  여러 isoform/hit이 나와도 한 표만 센다.
- `accepted call`: 해당 축 자체의 품질·coverage·margin 게이트를 통과한 evidence
  unit. 게이트 실패는 반대표가 아니라 abstention이다.
- `coverage`: family 또는 subfamily 전체 멤버 중 accepted call을 가진 **고유 유전자**
  비율. `support`는 투표자 안의 일치도이며 둘을 합치지 않는다.

## 3. 심볼 정규화

정규화는 두 단계다. 표시 문자열을 정돈하는 `canonicalize_symbol()`과 비교 key를
만드는 `symbol_key()`를 분리한다. 반환값은 항상 원문과 적용 규칙 목록을 포함한다.

### 3.1 결정적 규칙

| 순서 | 입력 | 출력/처리 | 규칙 |
|---:|---|---|---|
| 1 | Unicode | NFKC | 전각 문자 등 호환문자만 통일한다. |
| 2 | 앞뒤/연속 공백 | trim, 내부 공백 1개 | `Cry  2`를 결정적으로 만든다. |
| 3 | dash 변형 | ASCII `-` | U+2010~U+2015, minus를 통일한다. |
| 4 | stem과 접미 사이 공백 | 제거 | `Cry 2` → 표시형 `Cry2`. 단, 두 독립 토큰이면 거부한다. |
| 5 | 괄호 | primary와 alias로 분리 | `Cipk10 (Sip1)` → primary `Cipk10`, alias `Sip1`. alias는 `alias_of=Cipk10`인 audit evidence이며 독립 stem 후보나 표가 아니다. 중첩/닫히지 않은 괄호는 `AMBIGUOUS`. |
| 6 | 세미콜론 | 보존 | `Fba2;1`, `Pth1;4`에서 `;1`, `;4`를 버리지 않는다. 이는 allele/카피 정보다. |
| 7 | 그리스 문자 | key에서 이름으로 치환 | `α/β/γ/δ` → `alpha/beta/gamma/delta`; 표시형은 원 출처 철자를 보존할 수 있다. |
| 8 | 대소문자 | key만 Unicode casefold | `CRY2`와 `Cry2`는 같은 key. 표시 철자는 출처 우선순위로 고른다. |

구두점은 전역 삭제하지 않는다. `-`, `.`, `;`, `_`는 suffix 구조의 일부일 수 있다.
괄호 alias도 삭제하지 않고 provenance가 있는 audit 후보로 남기되, 명시적 alias 사전이
primary와의 동치 관계를 확인하기 전에는 독립 stem 투표로 승격하지 않는다. `At` 같은 종/족보
접두는 일반 규칙으로 제거하지 않는다. 출처별 allow-list에 등록된 경우에만
`prefix_alias` 규칙으로 바꾸며, 원문과 규칙 ID를 기록한다.

### 3.2 로마 숫자

로마 숫자는 `-II`, `_IV`, 공백 뒤 `III`처럼 **경계가 있는 말단 토큰**일 때만 suffix로
읽고 key에서는 대문자 표준형으로 만든다. 붙어 있는 말단 `I/V/X`는 보통 영문 심볼과
구별할 수 없으므로 자동 파싱하지 않는다. 필요한 족보는 override 표에 넣는다. 이 규칙은
`PHYC`의 `C`를 로마/일반 접미로 오인하는 것도 막는다.

### 3.3 충돌 안전성

서로 다른 raw symbol 두 개가 같은 `symbol_key`가 되었는데 primary/alias 관계나 허용된
표기 변형으로 설명되지 않으면 즉시 `NormalizationCollision`을 낸다. last-wins나
사전순 선택은 금지한다. 정규화 함수는 idempotent여야 한다.

## 4. stem 추출

`parse_symbol()`은 `(stem, suffix, parse_status, rule_id, aliases)`를 반환한다.
`parse_status`는 `PARSED`, `OVERRIDE`, `AMBIGUOUS`, `UNPARSED` 중 하나다.

### 4.1 적용 순서

1. **버전드 override 사전**: `symbol_key -> stem, suffix, reason`.
2. **allele/copy 문법**: 말단 `;<digits>[letters]?`를 앞의 suffix에 결합한다.
3. **숫자 계열 문법**: 영문 base 뒤의 `digits`, `digits.digits`, `digits-letter`,
   `digits+letter`를 suffix로 분리한다.
4. **경계 있는 계열 문법**: `-a3`, `_a3`, 경계 있는 그리스/로마 토큰을 suffix로
   분리한다. 단, `Vha-a3`처럼 계열 문자가 family stem의 일부인 것으로 override된
   족보는 override가 먼저 이긴다.
5. 하나의 유일한 parse만 남으면 `PARSED`; 0개면 `UNPARSED`, 2개 이상이면
   `AMBIGUOUS`.

일반 정규식은 “선행 비숫자 전부가 stem”으로 끝내지 않는다. 후보를 만들 뿐이며,
허용 문법과 override가 유일성을 보장해야 한다.

### 4.2 known-answer

| 입력 | stem | suffix | 근거 |
|---|---|---|---|
| `Tpt1`, `Tpt2` | `Tpt` | `1`, `2` | 숫자 계열 |
| `Bam1`, `Bam3`, `Bam5`, `Bam7` | `Bam` | 숫자 | 숫자 계열 |
| `Pth1;1` … `Pth1;4` | `Pth` | `1;1` … `1;4` | 숫자 + allele |
| `CKB1;1`, `CKB1;2` | `CKB` | `1;1`, `1;2` | 숫자 + allele |
| `Fba2;1` | `Fba` | `2;1` | 숫자 + allele |
| `Ppck1a`, `Ppck1b` | `Ppck` | `1a`, `1b` | 숫자+문자 계열 |
| `Snrk2.3` | `Snrk` | `2.3` | 점 숫자 계열 |
| `Vha-a1`, `Vha-a2`, `Vha-a3` | `Vha-a` | `1`, `2`, `3` | override: `a`는 subunit stem 일부 |
| `PHYC`, `PHYD` | `PHY` | `C`, `D` | override: 숫자 없는 족보 |
| `Cry 2`, `CRY2` | `Cry` | `2` | 공백/대소문자 정규화 후 숫자 계열 |
| `Cipk10 (Sip1)` | `Cipk` | `10` | primary parse; `Sip1`은 alias evidence |

`Acl`, `ZFN1`, `JMJ14`의 말단 문자를 일반적으로 떼지 않는다. `Acl`은 숫자도 경계도
없으므로 stem `Acl` 전체다. `PHYC/PHYD`처럼 숫자 없는 접미는 반드시 override가 있어야
한다. `AtPPa`의 `At` 제거도 source-scoped override 없이는 금지한다.

override 파일은 코드 상수가 아니라 버전드 TSV로 둔다. 최소 열은
`symbol_key, stem, suffix, scope_source, reason, added_in`이다. 중복 key의 상충 행과
사용되지 않은(stale) override를 검증 보고서에 낸다.

## 5. 축별 call 자격

| 축 | accepted call 조건 | 자동 판정에서 역할 |
|---|---|---|
| Mcry 큐레이션 | 정본 467 유전자 집합과 ID가 유일하게 join되고 symbol parse 성공 | primary. 동일 stem의 표시·desc·EC·AGI·UniProt 선택 우선권 |
| ATH | DIAMOND best locus, query/subject coverage 게이트, species-scoped RBH/현행 안전 게이트, **다른 stem**인 최선 경쟁자 대비 bit margin 통과 | 보정된 stem 전이 call. 같은 stem paralog/isoform은 margin 경쟁자가 아님 |
| AFDB-SwissProt | 최고 bit score, 구조 coverage/TM 품질, 다른 stem 경쟁자 대비 margin, 유일 UniProt symbol join, known-answer 보정 통과 | 보정 완료 전에는 corroboration 또는 conflict flag만; 단독 자동 stem 금지 |

E-value는 검색 후보를 자르는 데 쓸 수 있지만 순위·margin·동점 해소에 쓰지 않는다.
동일 bit score는 입력 순서와 무관하게 abstain한다. AFDB 축의 구체 임계값은 진행 중인
보정 결과로 preregister해야 하며 이 문서가 임의 숫자를 만들지 않는다.

축별 coverage 분모는 family의 정본 멤버 수다. subfamily call에서는 v4 catalog의 해당 OG
멤버 수를 쓴다. `support=1.00, coverage=0.10`은 높은 신뢰가 아니라 낮은 관측 범위다.

## 6. family 판정 상태기계

먼저 각 accepted call을 정규화하고, alias-equivalence/override까지 적용한 stem 집합 `S`를
만든다. parse 실패 call은 `S`에 넣지 않되 실패 수와 원문을 반드시 남긴다.

| 조건 | family verdict | 출력 행동 |
|---|---|---|
| accepted call 없음 | `NO_EVIDENCE` | label을 만들지 않고 audit 행은 남긴다. |
| parse 실패/충돌이 하나라도 판정에 영향을 줌 | `NEEDS_REVIEW_PARSE` | `transfer_tier=needs_review`; 후보와 오류 기록 |
| `|S| = 1`, 축간 stem 충돌 없음, 축별 게이트 통과 | `SINGLE_STEM` | family stem 행 `stem_auto`; 전체 symbol은 suffix gate로 보냄 |
| `|S| > 1`, 아래 HIGH partition 조건 전부 통과 | `MULTI_STEM_PARTITIONED` | 통과한 HIGH subfamily에만 라벨 행; family 단일 stem은 만들지 않음 |
| `|S| > 1`, partition 불완전/충돌/낮은 coverage | `NEEDS_REVIEW_MULTI_STEM` | 승자 없음; 모든 후보와 출처를 출력 |
| 같은 stem, 서로 다른 suffix | `SINGLE_STEM_MULTI_SUFFIX` | stem은 `stem_auto`; 각 suffix는 별도 tree gate |
| stem은 같지만 표시 철자만 다름 | `ORTHOGRAPHIC_VARIANT` | source 우선순위로 표시형 선택, 모든 raw variant 보존 |

### 6.1 HIGH partition 자동화 조건

다-stem family는 다음 조건을 **모두** 만족해야 subfamily 행만 자동 확정한다.

1. 사용한 catalog가 v4 manifest와 정확히 일치하고 family 멤버가 정본 family의 부분집합을
   넘어가지 않는다.
2. 각 자동 라벨 대상 OG가 `grade=HIGH`이다. PROVISIONAL/UNRESOLVED/NA는 증거 없음이지
   반대 증거가 아니다.
3. 서로 다른 stem의 evidence-bearing gene들이 서로 겹치지 않는 HIGH OG 집합에 놓인다.
   한 HIGH OG 안에 두 stem이 있으면 즉시 review다.
4. 각 symbol/suffix에 대해 `anchor_transferability()`가 rival label을 포함하지 않는 최대
   clade를 찾고 `transferable=True`를 반환한다. SH-aLRT/UFBoot는 약한 쪽이 기준이며,
   자동 suffix에는 `grade=HIGH`와 설정된 최소 coverage도 요구한다.
5. `members`와 v4 OG 멤버의 join coverage를 검사한다. unmatched gene이 허용 상한을
   넘으면 예외이며, 0으로 축소된 clade를 음성 결과로 읽지 않는다.
6. family 전체가 HIGH OG로 덮이지 않으면 결과는 `PARTIAL`이다. HIGH OG에는 안전한
   subfamily 행을 낼 수 있지만 나머지 멤버나 family 전체에 그 이름을 전이하지 않는다.

즉 HIGH 분할은 다-stem을 “표로 이긴 stem 하나”로 축약하는 도구가 아니라, 서로 다른
이름이 family 내부의 어느 검증된 조각에 적용되는지를 제한하는 도구다.

## 7. 축간 충돌 규칙

출처 우선순위는 다음처럼 제한적으로 적용한다.

| 큐레이션 | ATH | AFDB | 판정 |
|---|---|---|---|
| A | A | A/없음 | A `stem_auto`; 표시·metadata는 큐레이션 우선 |
| A | A | B | AFDB가 보정 전이면 A 유지 + `AFDB_CONFLICT` audit; 보정 후 accepted B면 review/partition 검사 |
| A | B | 임의 | 단순 큐레이션 승리 금지; HIGH partition 검사, 실패하면 review |
| 없음 | A | A | ATH A; AFDB는 독립 corroboration으로 기록 |
| 없음 | A | B | HIGH partition 검사, 실패하면 review |
| 없음 | 없음 | A | AFDB 보정 완료 및 게이트 통과 후에만 자동; 그 전에는 `NEEDS_REVIEW_UNCALIBRATED_SOURCE` |

두 약한 축의 합이 큐레이션을 “2표 대 1표”로 뒤집지 않으며, 큐레이션 한 표가 독립적인
다-stem 신호를 덮지도 않는다. 가중치는 동일 stem 내부의 대표 desc/metadata를 정하거나
진단용 support를 계산하는 데만 쓴다. 최종 verdict는 위 상태기계의 논리 게이트로 정한다.

ATH stem margin은 같은 stem paralog를 rival로 세지 않는 현재 규약을 유지한다. best와
runner-up의 stem이 다른 경우에만 `(best_bits - rival_bits) / best_bits`를 계산한다.
경쟁 stem hit가 없으면 margin은 `NA(single_stem_candidate)`이지 무한대가 아니다.

## 8. suffix tree gate

family stem과 전체 symbol을 별개의 산출물로 취급한다.

- `Bam` family stem: 보정 오차 1.7% 규칙 아래 `stem_auto` 가능.
- `Bam1`, `Bam3`, `Bam5`, `Bam7`: 각각 `suffix_needs_tree_gate`.
- `PHYC`, `PHYD`: override로 동일 `PHY` stem임을 알아도 `C/D`는 suffix이므로 동일하게
  tree gate를 거친다.
- tree 부재, 약한 support, 낮은 coverage, reference-lineage-specific duplication,
  rival label 공유 clade, query gene 0은 전부 `needs_review`/abstention이다.
- `MONOPHYLETIC`이나 verdict 행 부재를 “다른 family” 또는 transferable로 바꾸지 않는다.

## 9. 출력 스키마

기존 v1 12열은 순서와 의미를 유지한다.

```text
level target_id mcry_gene symbol stem desc ec agi uniprot module og_grade transfer_tier
```

`transfer_tier`는 판정 상태가 아니라 적용한 전이 방법이다. 따라서 tree gate를 통과한
suffix도 `suffix_needs_tree_gate`이고, 실제 통과 여부는 `verdict`가 말한다. enum에는
`needs_review` 하나를 추가한다. 그리고 최소한 다음 6열을 뒤에
append한다.

```text
verdict reason_code sources support coverage evidence_id
```

- `verdict`: 위 상태기계의 machine-readable 값.
- `reason_code`: `STEM_CONFLICT`, `PARSE_AMBIGUOUS`, `LOW_COVERAGE`,
  `CATALOG_NOT_HIGH`, `TREE_NOT_TRANSFERABLE`, `UNCALIBRATED_SOURCE` 등의 정렬된
  `;` 목록.
- `sources`: `mcry_curated;ath_diamond;afdb_swissprot` 중 실제 accepted source.
- `support`: 투표자 안에서 선택 stem의 가중 일치율. 최종 게이트 대신 쓰지 않는다.
- `coverage`: target의 전체 고유 멤버 중 accepted evidence gene 비율.
- `evidence_id`: 상세 audit sidecar의 행/레코드 ID.

`symbol`은 선택된 표시 문자열이고 raw variants, alias, bits, rival bits, margin, gate
결과, parser rule, catalog version, tree 결과는 `label_evidence_vN.tsv` sidecar에 둔다.
다-stem review 행은 `symbol`/`stem`을 억지로 합치지 않고 비우며 후보는 sidecar로 연결한다.
빈칸이 모호해지지 않도록 이때 verdict는 반드시 `NEEDS_REVIEW_*`여야 한다.

권장 `level` 값은 기존 `family|subfamily`를 유지한다. `MULTI_STEM_PARTITIONED`에서
family audit 행은 `level=family`, `transfer_tier=needs_review`로 두고, 자동 확정된 각 HIGH
OG만 `level=subfamily` 행을 가진다. membership table은 어느 경우에도 수정하지 않는다.

## 10. 의사코드

```python
def adjudicate_family(family, raw_calls, catalog, trees, policy):
    assert set(raw_calls.genes) <= family.members

    calls = []
    failures = []
    for raw in stable_sort(raw_calls):
        norm = normalize_symbol(raw.symbol, policy.aliases)
        parsed = parse_symbol(norm.primary, policy.overrides)
        audit(raw, norm, parsed)
        if not raw.axis_gate_passed:
            failures.append(abstain(raw.axis_reason))
        elif parsed.status not in {"PARSED", "OVERRIDE"}:
            failures.append(review(parsed.status))
        else:
            calls.append(one_evidence_unit(raw, parsed))

    calls = deduplicate_by_source_gene_stem(calls)
    coverage = unique_gene_coverage(calls, family.members)
    stems = equivalence_collapsed_stems(calls, policy.aliases)

    if failures_affect_choice(failures, stems):
        return needs_review("PARSE_OR_GATE_AMBIGUITY", calls, failures)
    if not stems:
        return no_evidence(failures)
    if len(stems) == 1:
        stem = only(stems)
        family_row = emit_stem_auto(stem, choose_display_by_source(calls),
                                    support(calls), coverage)
        suffix_rows = []
        for symbol in distinct_full_symbols(calls):
            suffix_rows += tree_gate(symbol, catalog, trees, policy)
        assert_row_contracts(family_row, suffix_rows)
        return family_row, suffix_rows

    partition = map_calls_to_catalog_ogs(calls, catalog)
    if not high_partition_is_disjoint(partition):
        return needs_review("STEM_CONFLICT", calls, partition)

    passed, withheld = [], []
    for stem, ogs in stable_items(partition):
        for og in ogs:
            result = anchor_transferability_for_og(stem, og, trees, policy)
            if (result.transferable and result.grade == "HIGH"
                    and result.coverage >= policy.min_label_coverage):
                passed.append(emit_subfamily_label(stem, og, result))
            else:
                withheld.append(result)

    family_audit = partitioned_family_audit(
        partial=bool(withheld) or not catalog_high_covers_family(partition, family)
    )
    assert_row_contracts(family_audit, passed, withheld)
    return family_audit, passed, withheld
```

모든 iteration과 tie-break는 입력 순서가 아니라 명시적 stable key를 쓴다. 단, score 동점
자체는 사전순 승리가 아니라 abstention이다.

## 11. 구현 위치

라벨 판정은 파이프라인의 clustering step이 아니므로 `pipeline.py`에 넣지 않는다.

- `steps/label_symbols.py`: Unicode/alias 정규화, stem parser, evidence-unit 축약,
  family 상태기계. 전부 plain dataclass/dict 입출력의 순수 함수이며 ete4 불필요.
- `data/symbol_stem_overrides_v1.tsv`: 버전드 예외와 source-scoped prefix alias.
- `adjudicate_labels.py`: v1 labels, ATH, AFDB, v4 catalog, tree 판정물을 읽고 쓰는 얇은 CLI.
  파일 I/O, 스키마 검증, manifest/checksum만 담당한다.
- `steps/subfamily.py`: `anchor_transferability()`를 그대로 재사용한다. 심볼 파서나 축 병합
  로직을 이 큰 모듈에 추가하지 않는다.
- `utils/gene_ids.py`: gene ID만 담당한다. 심볼 정규화는 의미와 충돌 규칙이 달라 여기에
  섞지 않는다. join은 기존 `match_ids`/collision guard를 사용하고 coverage 상한을 건다.

CLI는 입력 파일 checksum, override version, policy thresholds, catalog manifest, parser
version을 출력 manifest에 기록한다. 동일 출력 디렉터리 resume 시 하나라도 달라지면 거부한다.

## 12. 회귀 테스트와 불변식

### 12.1 정규화/파서 단위 테스트

- NFKC, dash, 공백, casefold가 idempotent다.
- `CRY2 == Cry 2` key, 표시 문자열은 source 우선순위로 결정된다.
- `Cipk10 (Sip1)`은 primary/alias 두 audit 레코드이나 primary evidence unit만 투표한다.
- `Pth1;1`, `CKB1;2`, `Fba2;1`, `Snrk2.3`, `Ppck1a/b`, `Vha-a1/a3`,
  `PHYC/PHYD`가 위 known-answer 표와 일치한다.
- `Acl`의 `l`, 일반 `ABC`의 `C`, 임의 `At...` 접두를 떼지 않는다.
- 붙은 로마 문자와 닫히지 않은 괄호는 자동 parse하지 않는다.
- 서로 다른 raw symbol의 normalization collision은 예외다.
- override 중복 충돌은 예외이고 stale override는 보고된다.

### 12.2 축 게이트 테스트

- ATH margin의 competitor는 다른 AGI가 아니라 **다른 normalized stem**이다. 동일 stem
  paralog/isoform은 제외된다.
- best bit 동점은 파일 순서와 무관하게 abstain한다.
- E-value만 바꾸어도 winner/verdict가 변하지 않는다.
- gate 실패 call은 반대표가 아니며 support 분모에 들어가지 않는다.
- 같은 gene/source의 중복 hit은 한 evidence unit이다.
- support 1.00/coverage 0.10 fixture가 자동 coverage gate를 통과하지 않는다.

### 12.3 상태기계 테스트

- 동일 stem + 여러 suffix → family `stem_auto`, 모든 suffix는 tree gate 전까지 확정 아님.
- 큐레이션 A + ATH B → 가중치와 무관하게 review 또는 HIGH partition뿐이다.
- 큐레이션 A + ATH A + 미보정 AFDB B → A 행과 AFDB conflict audit이 함께 남는다.
- 다-stem이 서로 다른 HIGH OG에 완전히 놓이고 각 tree gate가 HIGH이면 subfamily 행만
  자동 확정된다.
- 한 HIGH OG에 두 stem, catalog PROVISIONAL, tree coverage 부족, rival label 포함,
  query 0 각각은 review/abstention이다.
- HIGH로 덮이지 않은 family 멤버가 있으면 `PARTIAL`이며 그 멤버에게 label이 번지지 않는다.
- 입력 행 순서를 무작위로 섞어도 byte-identical output이다.

### 12.4 출력 불변식

- 정본 membership checksum과 family 수 13,431은 전후 동일하다.
- 모든 evidence gene은 정확히 한 정본 family에 join되거나 명시적 unmatched 행이다.
- `placed + unmatched == input evidence units`; 중복 제거 수는 별도 보존한다.
- `transfer_tier=stem_auto` iff verdict가 허용된 단일-stem verdict다.
- `transfer_tier=suffix_needs_tree_gate`인 행은 suffix가 비어 있지 않다.
- `transfer_tier=needs_review`인 행은 비어 있지 않은 `reason_code`와 `evidence_id`를 가진다.
- `transferable=False` 또는 `grade!=HIGH`인 tree 결과는 자동 suffix 행을 만들 수 없다.
- `support`와 `coverage`는 각각 `[0,1]`; `coverage` 분모와 분자가 evidence sidecar와
  재계산되어 일치한다.
- family-level 단일 stem 행에는 서로 다른 accepted stem이 숨어 있지 않다.
- 자동 subfamily label의 멤버는 해당 v4 HIGH OG 멤버의 부분집합이며 family 밖으로 새지 않는다.
- candidate 수 합계가 verdict별 집계와 일치한다. 누락 verdict는 허용하지 않는다.
- flag와 verdict의 양방향 불변식을 한 표로 parameterize한다. 예:

```python
def test_the_flags_and_verdicts_never_disagree():
    for row in all_cases():
        assert row.needs_review is row.verdict.startswith("NEEDS_REVIEW")
        assert row.auto_accepted is (
            row.transfer_tier in {"stem_auto", "suffix_needs_tree_gate"}
            and row.verdict in AUTO_VERDICTS
        )
        assert (row.reason_code != "") is row.needs_review
```

### 12.5 실측 known-answer 회귀 세트

62개 Mcry 다심볼 family 전부를 fixture manifest로 고정하되 원자료 전체를 테스트에 복제하지
않고 `family_id, raw symbols, expected stems, expected coarse verdict`를 둔다. 최소 필수 사례는
Tpt, Bam, Vha-a, Pth, CKB, DiT, PHY, CRY, Cipk/Sip, Fba, Ppck, Nhx/Acl,
Cipk/Snrk/Cbl, ZFN/JMJ, Fbp/Sbp, Almt 충돌군이다. override나 parser 변경 시 이 fixture의
diff가 0인지 CI에서 확인한다.

## 13. 자동화 경계

이 설계에서 자동화의 목표는 모든 family에 한 이름을 강제로 채우는 것이 아니다. 안전하게
확정 가능한 stem과 HIGH subfamily만 확정하고, 나머지를 **왜 판정하지 않았는지 재계산 가능한
형태로 남기는 것**이다. 따라서 `needs_review` 비율은 실패율이 아니라 현재 증거의 적용 범위다.
