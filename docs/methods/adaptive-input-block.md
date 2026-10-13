# 입력 파일에서 재분할 켜기 — `adaptive` 블록 (2026-09-20, 이름 변경 2026-09-23·2026-10-13)

입력 JSON에 `adaptive` 객체를 두면 `bin/main.out -f 입력.json`이 한 번 계산하는 대신 **적응형 재분할**을
돈다: 계산하고, 오차를 재고, 묶음(bundle)을 고쳐 다시 계산하기를 되풀이한다. 블록이 없으면 지금까지와
똑같이 한 번 계산한다.

```json
{
  "method": "cCCE", "order": 2,
  "cluster_mode": "bundle", "bundle": {"max_dim": 8, "merge_cutoff": 1e300},
  "...": "나머지는 보통의 CCEX 입력",

  "adaptive": {
    "policy": "merge",
    "output": "./adaptive_out",
    "atol": 1e-8,
    "rtol": 2e-6,
    "max_rounds": 8
  }
}
```

```
mpirun -n 4 bin/main.out -f 입력.json
```

- 결과는 `output` 폴더에 쌓인다: 단계별 폴더(`0001-plan`, `0002-probe`, …)와 전체 요약 `summary.json`.
- 랭크 수는 `mpirun -n`을 따른다. 계산은 같은 프로세스 안에서 돌므로 bath는 한 번만 읽는다.
- 재분할은 가변 클러스터링 위에서 돈다. 입력의 `cluster_mode`가 `bundle`(또는 `pcce`)이어야 한다(그 밖이면 시작 전에 거부).
- 이름: `cluster_mode=bundle`만 있으면 가변형(묶음을 한 번 만들어 계산), `adaptive` 블록이 더 있으면 적응형(오차로 되먹이는 반복)이다 (`clustering-terminology.md`). 2026-09-23까지 이 블록의 이름은 `Refine`이었다. 옛 이름은 읽지 않고 오류로 알린다. 2026-10-13까지는 `Adaptive`로 적었고(키는 대소문자를 가리지 않으므로 지금도 읽힌다), `cluster_mode`의 값은 `adaptive`, 묶음 키는 `bundle_max_dim` 같은 평평한 이름이었다(`docs/user/InputFile.md`의 옛 이름 표).
- 입력의 `#`·`!` 주석은 보통 입력과 같은 규칙으로 무시된다 (2026-09-23. 그 전에는 재분할 드라이버가 파일을 따로 읽어 주석이 있으면 실패했다). `seed`는 반드시 적어야 한다: 라운드마다 같은 계를 다시 계산하므로 bath 상태가 보존돼야 한다.

## 두 층의 구분

| | 하는 일 | 어디서 제어 |
|---|---|---|
| 가변 클러스터링 | 스핀을 묶음으로 묶어 **한 번** 계산 | 입력의 `cluster_mode`, `bundle` 블록 |
| 적응형 재분할 | 위 계산을 되풀이하며 묶음을 고침 | 입력의 `adaptive` 블록 |

## 키

이름은 명령줄 옵션(`--adaptive` 뒤에 쓰던 것)에서 `-`를 `_`로 바꾼 것이다. `true/false` 키는 `true`일 때만
켜진다. 모르는 키, 형식이 틀린 값, 빠진 필수 키는 시작할 때 `Error: adaptive.<키> …`로 알려 준다.

| 키 | 형식 | 기본값 | |
|---|---|---|---|
| `output` | 경로 | **필수** |
| `atol` | 실수 | **필수** |
| `rtol` | 실수 | **필수** |
| `max_rounds` | 정수 | `8` |
| `max_patches` | 정수 | `64` |
| `max_probes` | 정수 | `4096` |
| `max_candidates` | 정수 | — |
| `neighbors` | 정수 | `8` |
| `cache_mb` | 정수 | `256` |
| `timeout` | 실수 | `600` |
| `raw_cache` | 경로 | — |
| `probe_batch_clusters` | 정수 | `0` |
| `merge_priority` | 문자열 | `gain-per-cost` | 값: `gain-per-cost`, `worst-error`
| `policy` | 문자열 | `merge` | 값: `merge`뿐. `error-gated`, `budgeted`, `auto-partition`, `cost-preserving`은 2026-10-10에, `auto-order`는 2026-09-25에 제거했고, 쓰면 그 사실을 알리는 오류다
| `max_total_work` | 실수 | `8000000000000.0` |
| `prefer_physical` | true/false | — |
| `hub_candidates` | true/false | — |
| `regroup` | true/false | — | 방식 D(아래 절). 끄면 merge는 그대로다
| `regroup_max_regions` | 정수 | `64` | `regroup`과 함께만. 한 라운드에 검사하는 범인 영역 수 상한
| `regroup_max_candidates` | 정수 | `256` | `regroup`과 함께만. 한 영역의 후보 수 상한
| `regroup_factors` | 문자열 | `probes` | 값: `probes`, `factors`. `regroup`과 함께만. 범인 인자를 얻는 방법(아래 절)
| `regroup_whole` | 문자열 | `probes` | 값: `probes`, `except`. `regroup`과 함께만. 조건 (다)의 전체 곡선을 얻는 방법(아래 절)

2026-10-10에 제거한 정책들만 쓰던 키 20개(`plan_only`, `validation_fraction`, `max_order`, `max_total_seconds`,
`curve_change_budget` 등)는 이제 모르는 키로 오류다.

`input`, `binary`, `cwd`, `mpiexec`, `ranks`는 쓸 수 없다 — 입력은 이 파일이고, 바이너리는 이 프로세스이며,
랭크는 `mpirun`이 정한다. 블록에 쓰면 오류다.

## `regroup` — |L| > 1인 곳을 다시 묶기 (방식 D, 2026-10-11)

`"regroup": true`(명령줄 `--regroup`)이면 merge의 고르는 단계 대신 (internal notes)의
방식 D가 돈다(`src/clustering/algo_adaptive/policy_regroup.cpp`). 끄면(기본) merge의 결과는 바이트까지 그대로이고,
요약의 `settings`에도 세 키가 나타나지 않는다.

- 하는 일: 지금 분할로 계산해 모든 상태·시간점이 유한하고 |L| ≤ 1.01이면 끝. 아니면 최악점에서 연결된 bundle
  쌍마다 g = log|L_ij/(L_i L_j)|를 probe로 재고, g > 0인 쌍(범인) 둘레 세 bundle 영역 안에서 한 번의 이동(합치기·
  이분할·site 옮기기·site 바꾸기)을 찾는다. 후보는 영역의 국소 오차(상태 평균 곡선과 영역 직접해의 차, 시간 평균)가
  1e-4 넘게 줄고, 가장 가까운 바깥 bundle 셋 각각과의 국소 오차가 1e-4 넘게 커지지 않고, 전체 전개의 비물리 상태
  수가 늘지 않고 최악 상태의 초과량이 줄어야 받는다. 받은 분할은 다시 실행해 확인하고, 확인되지 않으면 되돌린 뒤
  막는다. 이미 가 본 분할로는 돌아가지 않는다. 받아들여 확인된 변경은 8번까지.
- 연결은 옛 bundle에서 옮기지 않고 site 거리(`rdip`, 안쪽 반지름 포함)로 매번 다시 정하며, 실행마다 CCEX가 계획한 연결과
  같은지 확인한다. 그래서 `response.physics.json`(site 좌표)을 첫 실행에서 받는다.
- 범위(시작 전에 이유와 함께 거부): `method` cCCE, `order` 2, `bundle.connections` distance(생략 시 distance)와
  `rdip`. gCCE, 다른 차수, `auto`·`checked`, 그리고 `prefer_physical`·`hub_candidates`와 함께 쓰면 거부한다.
- 비용: `max_total_work`가 실행 전체의 누적 한도다(D^3 × 상태 × 시간점; 첫 실행·범인 probe·영역 직접해·이웃 직접해·
  후보·바깥 항·재실행을 따로 센다). 첫 실행은 결과가 있어야 하므로 한도와 무관하게 돌고 그 비용도 장부에 적는다 —
  한도는 그 뒤의 D 탐색에 걸린다(첫 곡선이 physical이면 그대로 끝나고, 아니면서 첫 실행만으로 한도를 넘었으면 고칠 후보를 하나도 못 받아 `budget_limited`). 필요한 probe를 중복 없이 모아 실행 전에 예약하고, 받을 경우의 재실행 몫을
  남겨 둔다(한 번에 모은 후보들 가운데 가장 큰 재실행 몫). 감당 못 하는 후보는 검사를 건너뛰지 않고 거부한다(`budget_limited`). `max_probes`·`cache_mb`도 그대로
  적용된다. 실행 전체의 `bundle.max_clusters`·`bundle.max_work`를 넘는 후보도 CCEX에 묻기 전에 거부하므로,
  in-process에서 후보 하나 때문에 작업 전체가 끝나지 않는다.
- 끝(요약의 `status`): `physical`(종료 코드 0), `no_positive_culprit`, `no_candidate`, `budget_limited`,
  `change_limit`(2), `failed`(1). 결과는 마지막으로 확인된 분할이다. `physical`은 정확하다는 뜻이 아니다 — merge의
  `validated_full_bath_model`·`passed_convergence_checks`와 다른 이름을 쓰는 이유다.
- 요약의 `regroup` 블록: 알고리즘 판(`regroup-D-v1`), 상태 목록, ε·한도·이웃 수·변경 한도, 라운드마다 s*·t*·U·초과량,
  범인 쌍과 g·비유한 원인, 영역·이웃·후보 수·거부 사유별 개수·통과 수, 받은 분할과 새 연결, E_loc와 이웃별 오차
  (전·후), U·초과량의 예측과 재실행 실측, 바뀐 지수 항 수, 비용 장부(갈래별 work, 계산한 probe 수, 최대 차원).
- `regroup_factors`(명령줄 `--regroup-factors`, (internal notes)): 기본 `probes`는 매 라운드
  모든 bundle과 연결 쌍의 곡선을 모든 상태로 probe해 최악점의 값을 꺼낸다. `factors`는 CCEX의 `factors` 단계 하나로
  최악 상태 하나·시간점 하나의 값만 받는다(연결 하나 이상인 bundle과 연결 쌍; 값은 probe와 같은 함수로 계산해 비트까지
  같고, g는 드라이버가 같은 식으로 계산). 비용은 대략 실행 한 번의 1/상태 수이고 장부의 "culprit"에 한 상태 몫으로
  잡힌다. 받은 값은 probe 캐시에 넣지 않으므로 `max_probes`·`cache_mb`와 무관하다.
- `regroup_whole`(명령줄 `--regroup-whole`, (internal notes)): 기본 `probes`는 후보마다
  전체 전개의 모든 항을 probe해 다시 조립한다. `except`는 (가)를 통과한 후보가 있는 영역마다 CCEX의 `except` 단계를
  한 번 부른다 — 지금 분할을 그대로 실행하되 후보들이 지수를 바꾸는 항 X를 뺀 곱 P_K를 따로 내고, X의 곡선은 probe로
  받는다. 후보의 전체 곡선은 P_K × (후보의 새 항)이고, 새 항만 probe한다. except 단계의 전체 곡선이 직전 실행과,
  P_K × X가 그 곡선과 상대 1e-12 안에서 같아야 하며 아니면 `failed`다. 그 영역의 후보·바깥 항 곡선은 영역이 끝나면
  probe 캐시에서 버린다(장부의 work는 그대로; 요약 장부의 `curves_freed`). except 단계의 work는 실행 한 번 몫으로
  "outside"에 잡힌다.
- `except`는 long double이 x87 80비트인 빌드(x86의 wsl-gnu·idm4 등)에서만 쓸 수 있다. CCEX가 long double로 쌓은 곱을 같은 정밀도로 다시 조립해야 하기 때문이다(double로는 반올림이 1e-10까지 커졌다). long double이 64비트인 빌드(Apple Silicon)에서는 빌드는 되고 `except`만 시작 전에 거부한다.
- 둘 다 기본값이면 요약의 `settings`에 나타나지 않고 결과도 그대로다. `factors`·`except`는 `cluster_contributions`나
  `rdip_ladder`와 함께 쓰면 시작 전에 거부한다(CCEX가 그 단계를 거부하므로 in-process에서 작업 전체가 끝나지 않게).
  둘을 바꾼 결과 대조: (internal notes).
- merge의 `max_rounds`, `max_patches`, `max_candidates`, `neighbors`, `merge_priority`, `atol`, `rtol`은 D가 쓰지
  않는다(`atol`·`rtol`은 여전히 적어야 한다). ε = 1e-4, 이웃 3, 변경 8, |L| 한도 1.01은 저장 자료 시험의 재현용
  초기값이고 옵션이 아니다.

## 명령줄 방식과의 관계

같은 요청을 명령줄로도 할 수 있다(`bin/main.out --adaptive --input 입력.json --output … --atol …`). 둘은 같은
인자 목록이 되어 같은 코드로 들어가며, `test/cluster_verification/check_adaptive/inprocess_compare.py`가 merge 사례 다섯(merge-grow, gCCE, order 1, prefer_physical,
hub_candidates), regroup 사례 다섯(한 번 바꾸고 physical, rdip로 성긴 연결, 확장 한도에 걸린 후보 셋, 첫 사례의 `factors`판과
`except`판)과 `outfile` 없는 입력에서
두 방식의 산출물이 같음을 확인한다(2026-10-10 전에는 정책별 사례를 더해 여덟 사례, 2026-10-11 전에는 둘). 둘을 함께 쓰면 오류다.

재분할이 단계마다 쓰는 `input.json`에는 `adaptive` 블록이 복사되지 않는다 — 복사되면 단계 계산이 제 나름의
재분할을 시작하려 든다.
