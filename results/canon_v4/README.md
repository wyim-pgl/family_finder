# CANON — v4 산출물 로컬 사본 (2026-08-29)

원격(pronghorn / gpu)에만 있던 **현행 정본 산출물**을 로컬로 가져온 것이다. 목적은
"문서로만 뒷받침되던 수치"를 파일로 검증 가능하게 만드는 것 — Codex 감사(2026-08-29)가
원격 접근 불가로 확인하지 못했던 항목들이 여기서 전부 재현된다.

> 📌 이 디렉터리는 **읽기 전용 사본**이다. 계산·갱신은 원격 정본에서 하고 다시 내려받는다.
> 상태 규약은 `ANALYSIS_PARTITIONS.txt`와 위키 `guide/handoff-hygiene.md`.

## 파일과 출처

| 파일 | 출처 | md5 |
|---|---|---|
| `subfamily_catalog_v4.tsv` | `gpu:~/subfam/` | `78f5b7d29029df772b86d969479b3221` |
| `subfamily_catalog_v3.tsv` | `gpu:~/subfam/` | `0f601b579b438f636569099920ff467c` |
| `subfamily_catalog_v2.tsv` | `gpu:~/subfam/` | (구현 구판 기준선) |
| `labels_mcry_direct_v1.tsv` | `gpu:~/subfam/` | `5e0c471ae741cd1f5c2520fc4e528bac` |
| `freeze_catalog.py` | `gpu:~/subfam/` | 경화판(Codex 2라운드 반영) |
| `manifest_1014.txt` / `manifest_13431.txt` | `gpu:~/subfam/` | v3 / v4 동결 매니페스트 |
| `summary_v3.tsv` | `pronghorn:~/scratch/arbiter_v3/` | family 멤버십 **정본** |
| `cluster_verdicts.tsv`, `subfamily_targets.tsv` | `pronghorn:~/scratch/arbiter_v3/` | 심판 판정 / v3 대상 밴드 |
| `slate3_inapplicable.txt` | `pronghorn:~/scratch/subfam_trees/` | 트리 축 적용불가 32 |
| `curated_full.tsv`, `direct.tsv` | `pronghorn:~/scratch/mcry_direct/` | Mcry 큐레이션 입력 |
| `ANALYSIS_PARTITIONS.txt` | `pronghorn:~/scratch/bin/family_finder/` | 세대 구획 규칙 |

체크섬은 전송 후 원격 값과 대조해 일치를 확인했다.

## 로컬 재검증 결과 (전부 이 사본에서 계산)

| 주장 | 검증 |
|---|---|
| arbiter_v3 = 13,431 families / 459,398 유전자 | ✅ 정확히 일치, 중복 유전자 0 |
| v4 = 13,431 매니페스트 전수, 36,150 OG 행 | ✅ 매니페스트 13,431 고유, 행 36,150 |
| v4 매니페스트 == arbiter_v3 family 집합 | ✅ 완전 일치 |
| 등급 분포 HIGH 5,445 / PROV 4,293 / UNRES 26,370 / NE 42 | ✅ 일치 |
| 트리 적용불가 32 == v4의 `.NA` 행 | ✅ **집합 동일** |
| v3↔v4 공통 14,580 OG에서 등급·멤버십 diff | ✅ **0 / 0** (결정론 재현) |
| 라벨 레이어 457 유전자 / 366 family / 294 OG 전부 게이트 대기 | ✅ 일치 |

## 재현 명령

```bash
# 등급 분포와 v3↔v4 동일성
python3 - <<'EOF'
import csv
from collections import Counter
rows = list(csv.DictReader(open("subfamily_catalog_v4.tsv"), delimiter="\t"))
print(dict(Counter(r["grade"] for r in rows)))
v3 = {r["og_id"]: (r["grade"], r["members"])
      for r in csv.DictReader(open("subfamily_catalog_v3.tsv"), delimiter="\t")}
v4 = {r["og_id"]: (r["grade"], r["members"]) for r in rows}
c = set(v3) & set(v4)
print(len(c), sum(1 for k in c if v3[k] != v4[k]))
EOF
```

`.tsv`는 크기(38 MB) 때문에 `.gitignore`로 추적하지 않는다 — 이 README와 마커만 커밋된다.
