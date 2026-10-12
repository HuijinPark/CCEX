# Adaptive 분할 모드

> The `archive/…` names below are dated records that were taken out of the tree on 2026-09-24. They are in git history up to commit 83fe451d5: `git show 83fe451d5:docs/archive/<name>`.

용어는 [clustering 용어 정의](clustering-terminology.md)를 따른다: site를 묶은 것이 **bundle**,
연결된 bundle `order`개의 합집합이 **cluster**, 전체 배정이 **partition**이다.
2026-09-18 이전 문서의 "그룹"이 bundle이다.

현재 입력 설명은 2026-09-18 구현 기준이다. bundle order 3 이상과 제한된 자동 차수 증가 (archive/adaptive-high-order-20260912.md)도 지원한다.
아래의 초기 예제·검증 기록은 각 검사 당시 조건을 유지하며, 옛 옵션 이름을 그대로 둔다.
일반 입력은 `rdip` 유무로 연결 방식을 선택한다. 생략하면 제한된 coherence 검사 후 약한 연결을 생략한다.
차원 상한은 **`bundle_max_dim` 하나**이며 필수다. cluster 차원은 이 상한과
`order`에서 따라 나오는 값이고 별도 옵션으로 받지 않는다.

**계산 결과로 bundle을 다시 묶으려면** 입력에 `"Adaptive": {"policy": "merge", …}` 블록을 넣는다
([`Adaptive` 블록](adaptive-input-block.md); `policy`는 생략해도 merge; 2026-09-20 전에는 `scripts/adaptive_refine.py`).
정책은 merge 하나다: bundle 쌍을 시험 계산해 국소 오차를 가장 많이 줄이는 합치기를 한 번에 하나씩 채택한다
(선택: `prefer_physical`, `hub_candidates`). 정확도로 분할·차수를 자동 선택하던 auto-partition과
error-gated, budgeted, cost-preserving은 2026-10-10에 뺐다((internal notes) B5).
일반 `cluster_mode=adaptive`만으로 이 탐색이 시작되지는 않는다.

## 후속 정확도 검증의 시간 범위

2026-09-13 사용자 결정: 목적은 T2뿐 아니라 시간에 따른 다른 거동도 확인하는 것이다.
따라서 `|L|=1/e`에 도달해도 검사를 종료하지 않는다. 입력의 전체 시간 격자
`t_j = j * deltat`, `j=0,...,nstep-1`를 후보 평가, 잔여 오차 검사,
차수·선별 임계값·`rbath`·`rdip` 민감도 비교에 공통으로 사용한다.
현재 검증 코드도 모든 입력 시간점의 복소 coherence 차이를 비교한다.
크기뿐 아니라 위상 변화, 늦은 시간의 진동과 revival도 평가 범위에 남긴다.
시간점 사이의 거동까지 확인하려면 시간 간격에 대한 별도 민감도 검사가 필요하다.

전체 곡선에 사용할 최종 허용오차는 `rbath`·`rdip` 민감도 결과를 본 뒤 사용자가 정한다.
기존 검증의 `atol=1e-7`, `rtol=1e-5`는 당시 검사값으로 유지하며,
후속 계산의 물리적 정확도 요구로 자동 확정하지 않는다. 기존 h-BN 검사는
50점, `deltat=0.003 ms`, 즉 0–0.147 ms 전체를 사용했다.

2026-09-15 사용자 결정(잠정): 후속 계산의 판정 허용오차는 절대 `atol=1e-3`, `rtol=0`이다.
기록만 했으며, 실행할 때 `Adaptive` 블록에 `"atol": 1e-3, "rtol": 0`을 명시한다. 그 필수 키와
rdip 생략 모드의 내부 상수(`include/adaptive_defaults.h`, `1e-7`/`1e-5`)는 바꾸지 않았다. 절단 민감도(`rbath`·`rdip`)가
이 값보다 크면 분할보다 절단 수렴을 먼저 본다. 최종값은 나중에 다시 정한다.

후속 구현에서 클러스터 중요도를 기록할 때도 동일한 전체 시간창을 사용한다.
`max |log(L_tilde)|`는 선별을 돕는 진단값이며 전체 coherence 오차의 상한이 아니다.
실제 영점에서 정의되지 않는 로그와 비유한 값은 별도로 기록하며 작은 상수나 0으로 치환하지 않는다.
2026-09-14 추가한 native `"cluster_contributions": true`는 모든 전파 클러스터의
진부분 그룹 보정을 계산해 전체 시간창의 `max |log L_tilde|`를 기록한다.
필요한 진단용 하위 항과 메모리·시간 비용이 추가될 수 있다.
검증된 전역 coherence 잔여 오차 추정과 이를 이용한 정확도 종료는 남아 있다.

추가 기능: 계산 결과를 이용한 분할 개선 도구 (archive/adaptive-feedback-20260912.md)를 통해
기존 초기 분할을 시험 계산으로 검사하고 병합할 수 있다(merge 정책). 일반 adaptive 실행 자체가
자동으로 이 반복 탐색을 수행하는 것은 아니며 입력에 `Adaptive` 블록을 넣어야 시작된다.

2026-09-12 첫 구현. `cluster_mode=adaptive`는 실제 결합과 스핀 양자수로 서로 다른
크기의 그룹을 만든 뒤 그룹 단위 CCE를 계산한다. 기존 고정 `sK` pCCE와 독립된 경로다.
현재 **group order 1·2·3·4·5 이상**, `method=cce/gcce`, `quantity=coherence`, ensemble/single,
기존 CPU/CUDA/auto 실행을 지원한다. gCCE의 propagator/evolution 제약은 기존과 같다.
고차 실행은 행렬 차원·전개 개수·탐색량 상한을 만족해야 한다.

## 입력: rdip를 쓰거나 생략하기

일반 사용에서는 **`rdip`의 유무로 연결 방식을 선택**한다. bath·qubit·자기장·펄스·시간
설정은 기존대로 사용하고, 가변 분할에는 `cluster_mode`, `order`, `bundle_max_dim`을 넣는다.

| 입력 | 연결 선택 | 추가 coherence 검사 |
|---|---|---|
| `rdip` 있음 | 지정한 거리 안의 연결을 사용 | 일반 실행에서는 없음 |
| `rdip` 생략 | 약한 연결을 생략할 후보를 만들고 부분계에서 확인 | 제한된 직접 계산과 비교, 실패하면 연결 복구 |

거리 연결 예:

```json
{
  "cluster_mode": "adaptive",
  "order": 2,
  "bundle_max_dim": 16,
  "rdip": 30
}
```

약한 연결의 자동 선택 예:

```json
{
  "cluster_mode": "adaptive",
  "order": 2,
  "bundle_max_dim": 16
}
```

두 예의 실행 가능한 전체 입력은 [거리 방식](../../test/cluster_verification/check_bundle/simple_ensemble.json),
[rdip 생략 방식](../../test/cluster_verification/check_bundle/checked_ensemble.json)이다. 다음처럼 직접 실행한다.

```bash
bash do_compile.sh wsl-gnu        # 또는 cmake -B build && cmake --build build; 둘 다 bin/main.out을 만든다
mpirun -n 2 bin/main.out -f test/cluster_verification/check_bundle/simple_ensemble.json -v
mpirun -n 2 bin/main.out -f test/cluster_verification/check_bundle/checked_ensemble.json -v
```

입력의 경로가 저장소 최상위 기준이므로 거기서 실행한다. 로그의 `[adaptive] algorithm=`이 거리 방식은
`coupling-bound-v1`, rdip 생략 방식은 `coherence-checked-links-v1`이다. `checked_ensemble.json`은
6 site가 bundle 하나에 다 들어가 `[adaptive-check] status=single_bundle_no_links`로 끝난다. 연결 생략 검사가
실제로 도는 모습은 `bundle_max_dim`을 2로 줄이면 보인다(이 bath에서는 `status=base_unconverged_links_restored`,
2026-09-25 확인).

`order`는 한 cluster에 들어가는 최대 **bundle 수**, `bundle_max_dim`은 **bundle 하나**의
최대 차원이다. 스핀 수를 고정하거나 차원 상한을 `order`로 균등 분할하지 않는다.
cluster 차원은 bundle 차원의 곱(gCCE는 중심스핀 차원도 포함)으로 따라 나오며, 감당할 수
없으면 실행 중에 거부된다. 이것은 차원 제한을 만족하는 분할이며 정확도 판정은 아니다.

`rbath`는 전체 bath를 선택한다. `rdip`는 양의 유한한 거리이며 `null`이나 `0` 대신 **키를 생략**한다.
거리 연결은 서로 다른 site의 실제 구성 스핀 중 `rdipcut < r < rdip`를 만족하는 쌍이 있으면
유지한다. defect 부가 스핀도 포함하므로 pCCE의 그룹 중심 거리와 의미가 다르다.
`rdipcut`만 지정하는 입력은 새 생략 방식에서 거절한다.

### rdip 생략 시 무엇을 확인하나

1. 결합과 차원 상한으로 초기 그룹을 만든다. 이 단계에는 coherence 계산이 없다.
2. 그룹 사이의 결합 합과 마지막 계산 시간을 사용해 약한 연결 후보를 고른다.
3. 후보 연결 주변을 부분계로 잡는다. 같은 스핀들의 직접 계산, 생략 전 전개, 생략 후 전개를
   비교한다. 계산 차원이 허용하면 강한 이웃 그룹도 포함하고, bath 전체가 들어가면 전체를 검사한다.
4. 요청한 **모든 시간점과 모든 single 상태**에서 검사를 통과할 때만 후보를 함께 채택한다.
   생략으로 오차가 커지거나 원래 차수의 부분계 오차가 크면 후보 연결을 모두 복구한다.
5. 최종 그룹·연결 하나로 본 계산을 실행한다. 검사에서 생성한 bath/defect 상태를 재사용하므로
   검사 때문에 난수 시퀀스나 평균에 사용되는 상태가 바뀌지 않는다.

첫 구현은 보수적으로 **최대 8개 연결 후보, 최대 8개 부분계, 상태당 64개 고유 raw 클러스터**만
검사한다. 상태 저장과 raw 곡선 버퍼 예산은 64 MiB이며, raw 곡선은 상태마다 해제한다. 이 예산은 프로세스나 GPU 전체 메모리 상한은 아니다.
검사하지 않은 연결은 유지한다. 차원 상한 때문에 이웃을 포함할 수 없으면 해당 제한을 보고한다.
부분계 허용오차는 `1e-7 + 1e-5 * |L_direct|`이며, 생략 전후 차이도 같은 형태의 기준으로 검사한다.
이는 전체 bath의 오차 보장이나 최적 분할 증명이 아니다. 초기에는 전체 연결로 분할과 전개 가능 여부를 확인하므로,
큰 bath의 초기 전개 한도 초과를 자동으로 해결하지 않는다. 이 버전은 큰 bath를 한 번에 강하게
희소화하는 기능도 아니며, 추가 검사가 소규모 계산의 전체 시간을 늘릴 수 있다.

결과와 함께 `<outfile>.adaptive.json`을 저장한다. 여기에 실제 그룹·최종 연결·후보·검사 부분계,
오차·시간·예산·통과 상태가 들어간다. 주요 상태는 다음과 같다.

| 상태 | 의미 |
|---|---|
| `local_checks_passed` | 선택한 부분계 검사 통과, 후보 생략 적용; 전체 bath 수렴은 미확인 |
| `pruning_failed_links_restored` | 생략 때문에 오차가 커져 전체 후보 복구 |
| `base_unconverged_links_restored` | 생략 전 차수/분할부터 부분계 기준에 미달, 후보 복구 |
| `*_budget_unverified_links_kept` | 검사 예산 부족, 생략하지 않고 원래 연결로 실행 |
| `single_group_no_links` | 한 그룹으로 직접 계산하므로 그룹 간 연결 선택 불필요 |
| `order_one_no_link_expansion_unchecked` | order 1의 여러 그룹은 그룹 간 전개가 없어 연결 검사 미실시 |

보고서의 `accuracy_certified`는 항상 false다. `order`는 여전히 필수이며 자동 차수 결정은 없다
(auto-order는 2026-09-25, 차수를 올리던 auto-partition은 2026-10-10에 뺐다). 이 기능은 기존 NaN/Inf 안정 곱셈이나 PL6 전체 bath의 미수렴 문제를 해결했다고
판정하지 않는다. 구현과 검증 기록 (archive/adaptive-checked-connections-20260913.md).

<details>
<summary>기존 입력 호환 및 고급 설정</summary>

일반적인 두 입력에서는 아래 옵션을 지정하지 않아도 된다. 기존 명시 입력은 유지한다.

| 옵션 | 의미 및 기본값 |
|---|---|
| `bundle_max_dim` | **필수** 양의 정수. bundle 하나의 차원 상한 |
| `bundle_merge_cutoff` | 병합 후보의 최소 결합 강도, 기본 0 rad/ms |
| `expansion_max_clusters` | 전개 목록 한도, 기본 200,000 |
| `expansion_max_work` | 전개/차원 탐색 작업량 한도, 기본 20,000,000 |
| `bundle_connections=distance` | 기존 거리 방식, `rdip` 필수 |
| `bundle_connections=auto` | 기존 수동 제외 예산 방식, `bundle_connection_budget` 필수 |

명시적인 기존 `auto`는 **새 coherence 검사 방식과 다르다**. 위 두 일반 입력은
`bundle_connections`를 생략한다. 내부 로그의 `checked`는 rdip 생략을 파서가 해석한 값이며
사용자가 넣을 새 옵션이 아니다. explicit legacy 방식에 대해서는 아래 별도 설명을 참고한다.

`adaptive_work_file`의 plan/probe/fixed-groups 프로토콜과 `Adaptive` 블록의 재분할 정책은
별도 검증 정책을 사용한다. 여기에는 `rdip` 또는 기존 연결 설정을 명시해야 한다.
새 rdip 생략 자동 검사는 CCEX를 직접 실행하면 된다.

`sK`, `max_trial`, `max_iter`, `kmeans_pp`, `iter_detail`, `pcce_partition`, `pcce_connections`, `pcce_cache_file`은
adaptive에서 거절하고 `pcce_cache_mode=off`, `nk=0`만 지원한다.

</details>

## 그룹 생성 알고리즘

기존 `distance`의 알고리즘 이름은 `coupling-bound-v1`이다. 아래 기준은 분할 휴리스틱이며 coherence 오차
추정치나 전역 최적 분할 알고리즘이 아니다.

1. 선택된 모든 bath site를 하나씩 시작한다. 주 스핀과 해당 defect 부가 스핀은
   불가분 site로 취급하고 실제 차원을 곱한다. 나머지 site를 버리지 않는다.
2. 거리 조건에 해당하는 구성원 쌍이 있으면 그 두 site를 연결한다. 서로 다른
   defect site 사이의 상호작용은 기존 Hamiltonian과 같은 point-dipole tensor를 쓴다.
3. 연결된 두 site의 모든 교차 구성원 쌍에 대해 `3 I_i I_j ||T_ij||_F`를 합한다.
   이는 결합 Hamiltonian 연산자 norm의 보수적인 상한이며 단위는 rad/ms다.
   동일한 거리에서도 gyro와 I에 따라 값이 달라진다. Frobenius norm은 일관되게
   좌표계를 회전했을 때 보존되지만, 이 기준 자체는 외부장 방향에 따른 공명 억제를
   구분하지 않는다.
4. 기준 이상의 강한 연결부터 병합한다. 동률은 원본 site ID로 결정한다.
   그룹 차원 및 병합 후 요청 order까지의 연결 조합 차원을 검사한 뒤에만 병합한다.
5. 병합 기준 아래의 연결도 그룹 간 CCE 연결로 유지한다. 차원 상한을 넘는다는
   이유로 생성된 CCE 항을 나중에 버리지 않는다. 최소 site/pair도 상한을 넘으면
   명확한 오류로 종료한다.

root가 그룹을 결정하고 그룹 목록/연결을 broadcast한다. 초기 분할은 bath-state RNG를
소비하지 않으며 MPI 수나 backend에 따라 바뀌지 않는다. 새 연결 검사의 판정은
수치 허용오차를 사용하므로 경계값에서는 backend 간 재확인이 필요하다. 한 시간 곡선 전체에 같은
분할을 사용한다. 분할 과정의 메모리는 site 수와 후보 연결 수에 비례하며,
모든 site 쌍의 거리 검사는 현재 O(N²)다.

group order 2에서는 그룹 쌍의 지수가 +1, 개별 그룹의 지수가 `1 - 연결된 그룹 수`다.
group order 1에서는 모든 개별 그룹의 지수가 +1이다. 이 지수를 실제 site 목록에
전달한 뒤 같은 길이 버킷에 함께 저장한다. spin 단위 tilde 보정을 다시 하지 않는다.
기존 실행기는 `L0 * product((L_C/L0)^power)`를 계산하므로 0차 항의 활성값은 1을
유지한다. 특히 gCCE에서 이 정규화를 독립 부분계 계산과 비교했다.

## 기존 명시 auto 방식: 수동 제외 예산

`test/cluster_verification/check_bundle/auto_ensemble.json`은 기존 명시 `bundle_connections=auto`의 예제다.
새 rdip 생략 자동 검사에는 위의 `checked_ensemble.json`을 사용한다. 기존 bath·qubit·시간·pulse
설정을 유지하면서 다음 옵션을 지정하고 `rdip`를 생략한다.

```json
{
  "cluster_mode": "adaptive",
  "method": "cce",
  "order": 2,
  "bundle_max_dim": 64,
  "bundle_merge_cutoff": 1.0,
  "bundle_connections": "auto",
  "bundle_connection_budget": 0.01
}
```

자동 모드의 이름은 `coupling-budget-v1`이다. 각 그룹 쌍 A, B의 결합 점수는 모든
교차 구성원 쌍의 `3 I_i I_j ||T_ij||_F` 합인 `w_AB`이며 단위는 rad/ms다.
계산의 마지막 시간 `tmax = deltat * (nstep - 1)`을 곱해 다음 조건을 적용한다.

```text
omitted_score = tmax * sum(w_AB for omitted group pairs)
omitted_score <= bundle_connection_budget
```

1. 모든 선택 site 쌍의 결합을 평가하고 강한 결합부터 차원 상한 안에서 병합한다.
   병합으로 이웃 그룹 쌍이 너무 커져 제외해야 하는 결합도 예산에 포함한다.
   그 때문에 새 예산 초과가 생기는 병합은 거절한다.
2. 완성된 그룹 사이의 결합을 합산한다. `rdip`/`rdipcut` 거리 조건 또는 행렬 차원
   상한 때문에 계산할 수 없는 그룹 쌍의 결합부터 제외 예산에 넣는다.
   그룹 쌍에 거리 조건을 만족하는 구성원 쌍이 하나라도 있으면 거리 후보로 인정한다.
   이 후보의 점수에는 거리 밖의 다른 교차 구성원 결합도 포함한다.
3. 남은 그룹 쌍은 약한 순서로 제외하되 전체 제외 점수의 합이 예산을 넘지 않게 한다.
   많은 약한 연결이 합쳐져 중요해지는 경우도 이 합산에 반영한다.
4. 연결 선택이 끝난 후 포함·배제 지수를 생성한다. 계산된 CCE 항을 사후 삭제하지 않는다.
   거리나 차원 제한만으로 이미 예산을 넘으면 이유를 출력하고 종료한다.

예산 0은 양의 결합을 제외하지 않는 엄격한 비교 기준이다. 차원이나 거리 제한으로
이를 만족할 수 없으면 오류다. `tmax=0`인 계산에서는 제외 점수도 0이다.
group order 1은 모든 잔여 그룹 간 연결을 제외하므로, 하나의 그룹으로 합쳐지지
않으면 작은 예산에서 거절될 수 있다. `order`, 전체 계산 차원 상한과 제외 예산은
이 기존 방식에서 사용자가 지정한다. 그룹 상한과 병합 cutoff는 생략할 수 있다. 어떤 bath에도 적합한 예산 기본값을 넣지는 않았다.

예제의 6스핀, pulse 1, 0–0.04 ms 조건은 3/2/1스핀 그룹을 만든다. 예산 0이면
그룹 간 연결 3개를 유지하고, 예산 0.01이면 가장 강한 연결 1개를 유지한다.
나머지 두 연결의 합산 제외 점수는 약 0.005956이다. 이 숫자들은 이 예제의 동작을
보이기 위한 설정이다. **예산 0.01은 coherence 오차 1%를 보장한다는 뜻이 아니다.**
group order 절단 오차와 곱·나눗셈의 수치 불안정성은 별개다.

같은 분할에서 유지할 연결 개수를 줄이는 선택이며 전역 최적 분할이나 최소 실행시간을
보장하지 않는다. 같은 bath와 마지막 시간이면 pulse 수·ensemble/single·backend에
따라 그래프가 바뀌지 않는다. pulse와 공명 특성은 기존 동역학 계산에서 반영한다.
실제 사용에서는 예산을 줄이면서 결과의 수렴을 확인해야 한다.

자동 모드는 `rdip`가 있더라도 제외 결합을 감사하기 위해 전체 site 쌍을 평가한다.
따라서 CPU 전처리의 결합 평가와 초기 그래프 메모리는 O(N²)이며, 작은 `rdip`를
주어도 이 부분의 비용은 줄지 않는다. 병합 과정에는 추가 그래프 갱신 비용이 있다.
큰 bath는 충분한 메모리의 서버에서 실행한다.

## 출력과 제한

`[adaptive]`는 알고리즘/방법, 그룹 order, site/그룹/pair 수, 병합 수, 상한으로
차단된 연결 수와 차원 제한을 출력한다. `[adaptive-group]`는 모든 그룹의 차원과
0부터 시작하는 원본 bath site IDs를 기록한다. `[adaptive-expansion]`은 실제로
전파하는 비영 지수 항 수와 최대 site 수를 출력한다. `[profile]`의 차원은 defect
및 gCCE qubit까지 반영한 실제 행렬 차원이다.

기존 명시 `bundle_connections=auto`에서는 `[adaptive-connections]`에 유지/제외 그룹 쌍 수, 제외 점수, 예산,
거리·차원으로 제외한 점수와 예산 때문에 차단한 병합 수를 기록한다. `-v` 실행은
`[adaptive-link]`에 선택된 그룹 ID 쌍과 결합 크기도 기록한다. 생략한 `rdip`는
설정 출력에서 `unlimited`로 표시한다. `backend=auto`는 CPU/GPU 선택 옵션이며,
`bundle_connections=auto`와 독립적이다.

초기 결합 기반 분할 자체는 아직 B·detuning·hyperfine·quadrupole에 따른 공명 특성을
분할 중요도에 반영하지 않는다. 해당 항들의 기존 동역학 계산은 그대로 수행하며,
새 오차 기반 드라이버는 그 계산 결과로 제한된 병합·재배치 후보를 평가한다.
새 rdip 생략 방식은 제한된 연결 생략을 coherence로 검사한다.
반복 분할 개선과 연결 선택을 함께 최적화하는 전역 탐색은 구현하지 않았다. 고차 전개와
제한된 자동 차수 증가는 별도 설명 (archive/adaptive-high-order-20260912.md)을 참조한다.
차원 제한만으로 동시에 실행하는 모든 MPI
프로세스의 전체 메모리 사용량을 보장하지는 않는다.

## 검증

`bundle_connections=auto`의 추가 검사, 실제 bath 크기·연결 수·실행시간 및
서버 작업 464/465 결과는 자동 연결 검증 보고서 (archive/adaptive-auto-connections-20260912.md)에 있다.
아래 기록은 첫 거리 기반 adaptive 구현의 검증이다.

- 독립 complex 포함·배제, 모든 작은 set partition, 혼합 차원, 입력 edge 순서,
  상한/overflow 등 13,257개 검사 및 동일 검사의 ASan/UBSan 통과.
- 최종 CPU/CUDA 통합 실행 117건 통과: singleton → 일반 CCE2, 하나/두 그룹 →
  전체 6스핀 Hamiltonian, 비균일 3/2/1 그룹의 독립 부분계 곱, MPI 1/2/4,
  ensemble/single, pulse 1/4, I=1/2·1·5/2 혼합, P1 부가 핵스핀,
  기존 pCCE 결과 및 잘못된 입력 처리.
- 비균일 그룹의 독립 부분계 곱과 차이는 저장된 출력 정밀도에서 최대
  `2.03e-10`이었다. 부분계들은 이전 실행파일에서 각각 전체 Hamiltonian으로
  계산하고 별도 Python 코드에서 group CCE 곱을 구성했다.

6스핀 예제의 0–0.04 ms, 9시간점에서 전체 bath Hamiltonian 대비 최대 절대 오차:

| 계산 | pulse | 일반 spin CCE2 | adaptive group order 2 |
|---|---:|---:|---:|
| cCCE ensemble | 1 | 0.033825 | 1.7e-9 |
| cCCE ensemble | 4 | 0.090813 | 1.0e-10 |
| gCCE ensemble | 1 | 0.033746 | 8.01e-9 |
| gCCE ensemble | 4 | 0.090716 | 2.42e-8 |

이는 강한 3스핀 그룹이 있는 특정 예제다. 일반 cCCE2의 최대 bath 차원은 4,
adaptive의 최대 bath 차원은 32이므로 같은 계산 비용의 비교가 아니다.
큰 bath에서의 오차 보장 또는 보편적인 속도 향상으로 해석하지 않는다.

실제 bath에서도 CPU/CUDA/auto 총 15회 실행이 통과했다. 공통으로 group order 2,
group cap=4, 전체 cluster cap=48, cutoff=1 rad/ms, ensemble, pulse 4,
33시간점이다. NV는 B=300 G, rbath/rdip=30/8 Å, 0–1.28 ms이고,
PL6는 B=57 G, rbath/rdip=40/10 Å, 0–2.3366 ms다.
CPU/auto MPI 4개, CUDA MPI 2개에서 그룹과 클러스터 목록이 동일했고 최대 출력
차이는 약 1e-10이었다. 이 실행들은 기능 비교이며 통제된 속도 벤치마크가 아니다.

| Bath/엔진 | 유지 site | 그룹 수 | 1-site / 2-site 그룹 | 전파 항 수 | 최대 행렬 차원 |
|---|---:|---:|---:|---:|---:|
| NV–13C / cCCE | 222 | 178 | 134 / 44 | 388 | 16 |
| NV–13C / gCCE | 222 | 178 | 134 / 44 | 388 | 48 |
| PL6 bath 1 / cCCE | 788 | 608 | 428 / 180 | 3,268 | 16 |
| PL6 bath 2 / cCCE | 765 | 584 | 403 / 181 | 3,093 | 16 |
| PL6 bath 3 / cCCE | 768 | 601 | 434 / 167 | 3,179 | 16 |

gCCE의 별도 L0 계산은 전파 항 수에서 제외했다. 신규 분할기와 MPI 전달/배열 변환
코드를 계측한 ASan/UBSan 실행 4건도 통과했다. MPI 실행에서는 라이브러리의 leak
검출을 끄고 주소/정수 오류를 검사했으며, 독립 분할기 검사는 leak 검출을 켰다.

구조화 검증 결과는 `test/cluster_verification/check_bundle/validation-20260912.json`,
통합 결과는 `integration-20260912.json`, 실제 bath 결과는 `baths-20260912.json`이다.
원본 입력/출력/소스/로그 809개 파일은 `test/cluster_verification/check_bundle/archive/adaptive-local-20260912.tar.gz`에
보존했다. SHA-256은 `06bfaf8032fa5410c9d00114984b9724858eea85e36bd1711d64b468b0d815b3`이다.

IDM/GPU 공유 저장소: `/home/dlwodnr95/tmp/ccex-adaptive-20260912`.
소스/의존 파일/입력/검증 도구 4,270개 파일을 전송 후 해시 확인했다.
GPU 서버 검증은 재실행 작업 **462**에서 **3분 10초, exit 0:0으로 통과**했다.
13,257개 분할 검사와 ASan/UBSan, 117회 통합 실행, 실제 bath 15회 비교를 완료했다.
실제 bath는 CPU/auto MPI 32개와 CUDA MPI 4개를 비교했다.
처음 작업 461은 중첩한 기준 소스의 legacy `../include` 경로가 상위 새 헤더를 읽어
링크 실패했다. 원본 기준 소스를 독립 위치에 그대로 복사해 빌드했고 해시 동일성을
확인했다. 물리 계산 코드와 분할 알고리즘은 바꾸지 않았다.
원래 소스 압축 파일과 `adaptive-repair.tar.gz`를 함께 보존했다.
전체 서버 결과와 기존 PL6의 미해결 NaN/Inf는 서버 결과 보고서 (archive/gpu-server-results-20260912.md)에 구분했다.

재현 도구 (2026-09-12 당시의 것이다. 아래 Python 검사는 그때의 옵션 이름과 바이너리 이름을 쓰므로 지금 빌드에서는
돌지 않는다. `g++` 줄의 단위 검사 셋은 지금도 돈다(2026-10-02: 13,256, 550, 9,363건 통과; 연결 검사의 줄에
`bundle_expansion.cpp`를 더했다). 그 밖의 지금 검사는 `test/cluster_verification/check_adaptive/verify.sh`이다):

```bash
g++ -std=c++17 -O2 -Iinclude test/cluster_verification/check_pcce/test_adaptive_partition.cpp src/clustering/bundle/bundle_partition.cpp src/clustering/bundle/bundle_expansion.cpp -o /tmp/adaptive-partition-test
/tmp/adaptive-partition-test
python3 test/cluster_verification/check_bundle/check_adaptive.py --baseline BEFORE_CPU --output NEW_CHECK_DIR
python3 test/cluster_verification/check_bundle/check_baths.py --fixtures FIXTURES --output NEW_BATH_CHECK_DIR
g++ -std=c++17 -O2 -Iinclude test/cluster_verification/check_pcce/test_adaptive_connections.cpp src/clustering/bundle/bundle_partition.cpp src/clustering/bundle/bundle_expansion.cpp -o /tmp/adaptive-connections-test
/tmp/adaptive-connections-test
g++ -std=c++17 -O2 -Iinclude test/cluster_verification/check_pcce/test_adaptive_high_order.cpp src/clustering/bundle/bundle_partition.cpp src/clustering/bundle/bundle_expansion.cpp -o /tmp/adaptive-high-order-test
/tmp/adaptive-high-order-test
python3 test/cluster_verification/check_bundle/check_auto_connections.py --baseline BEFORE_CPU --output NEW_AUTO_CHECK_DIR
```
