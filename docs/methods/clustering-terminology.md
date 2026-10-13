# Clustering 용어 정의

2026-09-18 확정. CCEX의 세 clustering 방법(표준 CCE, pCCE, 가변형/적응형)을 한 틀로
기술하기 위한 용어다. 이 문서의 이름을 코드·문서·논문에서 그대로 쓴다.

## 1. 네 단위

| 용어 | 뜻 | 차원 |
|---|---|---|
| **site** | bath spin 하나, 또는 나눌 수 없는 결함 하나 | `2S+1` (실측 2, 3, 4, 7) |
| **bundle** | 항상 통째로 cluster에 들어가는 site 묶음. 안쪽은 정확히 계산된다 | site 차원의 곱 |
| **cluster** | 연결된 bundle `order`개의 합집합. **실제로 행렬을 만들어 전파하는 단위** | 그 안 site 차원의 곱 (× 중심스핀) |
| **partition** | 모든 site를 각각 정확히 하나의 bundle에 배정한 것 **전체** | — |

`cluster`의 뜻은 원래 CCEX와 같다(`src/physical/bath.cpp:182` `BathArray_dim`이 그 차원을
site 차원의 곱으로 계산한다). 가변형이 새로 끼워넣은 것은 **bundle** 한 층뿐이다.

```
표준 CCE : site ────────────────> cluster        (order = site 개수)
pCCE     : site ──(위치 K개씩)──> center ──────> cluster  (order = center 개수)
가변형   : site ──(결합 세기)───> bundle ──────> cluster  (order = bundle 개수)
```

pCCE는 같은 자리의 층을 **center**라고 부른다(`src/clustering/cluster_pcce.cpp:24`,
`pCCE(N,K) indicates a clustering method of the N-clusters having K number of bath spins`).
혼동을 막기 위해 가변형의 층은 **bundle**로 부르고 `group`이라는 말은 쓰지 않는다.
`src/clustering/cluster_dsjitb.cpp`·`include/clustering/cluster_dsjitb.h`의 `GroupChecker`는 disjoint-set 알고리즘의 것으로 이 층과 무관하다.

## 2. 두 방법

**가변형 clustering (variable clustering)** — bundle을 만드는 방법. 한 bundle에 여러
spin을 넣는 것을 고려한다. 결합이 강한 site들을 한 bundle로 묶어 그 상호작용을 정확히
계산하고, 약한 결합만 bundle 경계에 남겨 `order`차까지의 전개로 근사한다.
bundle이 전부 1 site면 표준 CCE와 같다.

**적응형 clustering (adaptive clustering)** — 가변형이 만든 partition이 최선이 아닐 수
있으므로, **오차를 재서 bundle을 재구성하는 과정 전체**. 후보 partition을 만들고, 실제로
계산해 본 오차로 채택 여부를 판정하고, 필요하면 다시 만든다. 정책은 merge
하나다(국소 오차를 가장 많이 줄이는 bundle 합치기; 선택 기능 `prefer_physical`, `hub_candidates`).
error-gated, budgeted, auto-partition, cost-preserving은 2026-10-10에 뺐다: 적응형의 목표를
bundle에서 |L| > 1이 나오는 곳을 다시 묶는 것으로 정했고, 그 목표에는 정책 하나면 된다
((internal notes) B5). 그 전에 auto-order(분할을 고정하고 order만 올려 가며
budgeted를 반복)는 2026-09-25에 뺐다: 차수끼리 일치해도 정확도를 인증하지 않아, 입력의
`order`를 바꿔 따로 돌리는 것과 다른 판정을 하지 않았다.

즉 가변형은 **한 번의 구성**, 적응형은 **오차로 되먹이는 반복**이다.

## 3. order

`order`는 **cluster 하나에 들어가는 bundle의 개수**다. 표준 CCE의 "spin 개수"와 이름은
같지만 세는 대상이 다르다. CCEX 자신의 로그도 `bundle_order`로 구분해 찍는다
(`src/clustering/bundle/bundle_cluster.cpp:143`; 2026-09-18 전 기록에서는 `group_order`).

bundle 하나도 cluster다: 전개는 크기 1인 집합부터 만든다(`expandAdaptivePlan`, `src/clustering/bundle/bundle_expansion.cpp`).

## 4. 차원 상한

`bundle.max_dim` **하나만** 받는다. cluster 차원은 여기서 따라 나오는 값이며
(`bundle.max_dim^order`가 상한), 별도 옵션으로 받지 않는다. 2026-09-18 이전의
`adaptive_cluster_max_dim`은 제거했다 — cluster 상한을 받으면 `order`를 올릴 때마다
bundle이 조여들어, "어느 spin을 정확히 계산할 것인가"와 "상관을 어디까지 볼 것인가"라는
서로 다른 두 결정이 한 숫자에 묶였다.

## 5. 쓰지 않는 말

| 쓰지 않음 | 대신 |
|---|---|
| group | bundle |
| group order | bundle order |
| adaptive_group_max_dim, bundle_max_dim | bundle.max_dim |
| adaptive_cluster_max_dim | (제거) |

`test/cluster_verification/check_bundle/` 아래의 과거 실행 기록과 그것을 읽는 분석 스크립트는 옛 이름을 그대로
둔다. 그 기록은 작성 당시의 산출물이고, 새 코드로 다시 돌지 않는다.
