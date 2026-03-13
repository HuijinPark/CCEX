# CCEX 새로운 방법들 — 프레젠테이션 정리

---

## 1. 전체 방법 비교표

| 방법 | 클러스터링 | 해밀토니안 | 소산(Dissipation) | 구현 파일 |
|------|-----------|-----------|-------------------|----------|
| **CCE** | Hash (r_dipole 내 페어) | 단순 곱 공식 (secular) | ✗ | `simulator_cce.cpp` |
| **gCCE** | Hash (r_dipole 내 페어) | 풀 해밀토니안 | ✗ | `simulator_gcce.cpp` |
| **pCCE** | k-means (sK개/그룹) | 단순 곱 공식 (secular) | ✗ | `cluster_pcce.cpp` + cce 커널 |
| **pgCCE** | k-means (sK개/그룹) | 풀 해밀토니안 | ✗ | `cluster_pcce.cpp` + gcce 커널 |
| **ME-CCE** | Hash (r_dipole 내 페어) | 단순 곱 공식 (secular) | ✓ Lindblad | `simulator_mecce.cpp` |
| **ME-gCCE** | Hash (r_dipole 내 페어) | 풀 해밀토니안 | ✓ Lindblad | `simulator_megcce.cpp` |
| **pmeCCE** | k-means (sK개/그룹) | 단순 곱 공식 (secular) | ✓ Lindblad | `cluster_pcce.cpp` + mecce 커널 |
| **pmegCCE** | k-means (sK개/그룹) | 풀 해밀토니안 | ✓ Lindblad | `cluster_pcce.cpp` + megcce 커널 |

> 새로 추가된 방법: **ME-CCE, ME-gCCE, pCCE, pgCCE, pmeCCE, pmegCCE**

---

## 2. CCE와 gCCE의 차이

### CCE (secular approximation)
- 큐빗 상태(α, β)에 projection된 해밀토니안을 사용
- H_α = H_bath + ⟨α|H_qb|α⟩,  H_β = H_bath + ⟨β|H_qb|β⟩
- 코히런스: L_k(t) = Tr[ ρ_bath · e^{iH_β t} · e^{-iH_α t} ]
- 계산이 빠르고 메모리 절약

### gCCE (full Hamiltonian)
- 큐빗 + 배스 전체 공간에서 해밀토니안을 구성
- H_total = H_q ⊗ I + I ⊗ H_b + H_qb (전체 상호작용 포함)
- 코히런스: L_k(t) = ⟨ψ_α|Tr_bath[e^{-iH_tot t} ρ_tot e^{iH_tot t}]|ψ_β⟩ / L_0(t)
- 더 정확하지만 계산 비용이 큼
- **0번째 차수 클러스터 (큐빗만)를 계산해서 normalization 필요**

---

## 3. Lindblad Master Equation (새로 추가된 핵심)

### 수식

$$\frac{d\rho}{dt} = \mathcal{L}[\rho] = -i[H, \rho] + \sum_k \left( C_k \rho C_k^\dagger - \frac{1}{2}C_k^\dagger C_k \rho - \frac{1}{2}\rho C_k^\dagger C_k \right)$$

- **첫 번째 항**: 유니터리 코히런트 진화 (기존 CCE와 동일한 물리)
- **두 번째 항 (D[C_k]ρ)**: Lindblad dissipator — 배스 스핀의 비가역적 완화

### Jump Operator (이 코드에서 사용하는 것)

| 연산자 이름 | 행렬 | 물리적 의미 |
|------------|------|------------|
| `"+"` | √Γ · S⁺ | 배스 스핀 올림 (flip-up) |
| `"-"` | √Γ · S⁻ | 배스 스핀 내림 (flip-down) |
| `"z"` | √Γ · Sz | 위상 감쇄 (dephasing) |

> S⁺, S⁻를 동시에 같은 rate Γ로 주면: 배스 스핀이 열역학적으로 균형된 환경과 결합하는 것을 모사 (T₁ 완화)

---

## 4. Superoperator Formalism (C++ 코드 핵심 구현)

행렬 방정식을 벡터 방정식으로 변환하여 수치적으로 풀기 쉽게 만드는 기법.

### 벡터화 (row-major)

$$\text{vec}(\rho)_{i \cdot N + j} = \rho_{ij}$$

### 핵심 공식

| 연산 | Superoperator |
|------|--------------|
| ρ → A ρ B | S = A ⊗ B^T (= `op_to_supop(A, B)`) |
| ρ → -i[H, ρ] | L_coh = -i (H⊗I - I⊗H^T) |
| ρ → D[C]ρ | L_D = C⊗C* - ½(C†C)⊗I - ½I⊗(C†C)^T |
| ρ(t) | ρ(t) = exp(L·t) ρ(0) |

### 시간 전파

- **자유 진화**: `exp(L * t)` (행렬 지수함수)
- **Hahn echo (npulse=1)**: U_after · U_before
  - U_before = exp(L * t/2), U_after = exp(L_flipped * t/2)
  - L_flipped: α↔β 교환 (π 펄스 효과)

---

## 5. ME-CCE vs ME-gCCE 차이

### ME-CCE (projected / secular)
- 배스 공간에서만 density matrix ρ_bath (d_bath × d_bath)
- L_coh = projected_coherent_superoperator(H_α, H_β)
  → -i (H_α ⊗ I - I ⊗ H_β^T)
- 코히런스 = Tr(ρ_bath(t))
- **장점**: 작은 행렬, 빠른 계산
- **단점**: secular approximation이 성립해야 함

### ME-gCCE (full Hamiltonian)
- 전체 공간 density matrix ρ_total (d_q × d_bath × d_q × d_bath)
- L_coh = coherent_superoperator(H_total)
  → -i (H_total ⊗ I - I ⊗ H_total^T)
- 코히런스 = ⟨ψ_α|Tr_bath[ρ_q(t)]|ψ_β⟩ / L_0
- **장점**: 더 정확한 물리
- **단점**: 큰 행렬 (큐빗 × 배스 공간)

---

## 6. pCCE 클러스터링 (k-means 기반)

### Hash clustering (CCE, gCCE, ME-CCE, ME-gCCE)
- r_dipole 내에 있는 스핀 쌍만 클러스터로 묶음
- 디폴쌍 상호작용 강도 기준 정렬

### pCCE clustering (pCCE, pgCCE, pmeCCE, pmegCCE)
- 전체 배스 스핀을 **k-means**로 sK개씩 묶음
- r_dipole 제한 없이 **강제로 모든 스핀을 그룹화**
- `sK=1` → 기존 CCE/gCCE와 동일한 결과 (수치 검증 완료)
- `sK=2` → 2개씩 묶인 그룹으로 계산 (더 많은 상관관계 포함 가능)

#### 파라미터
- `sK`: 그룹당 배스 스핀 수
- `max_trial`: k-means 재시작 횟수
- `max_iter`: k-means 최대 반복 횟수
- `kmeans_pp`: k-means++ 초기화 (기본값: True)

---

## 7. 검증 결과 요약

| 검증 | 결과 | max|diff| |
|------|------|----------|
| ME-CCE(rate=0) = CCE | **PASS** | 5.2 × 10⁻⁹ |
| ME-gCCE(rate=0) = gCCE | **PASS** | 1.2 × 10⁻⁸ |
| pmeCCE(rate=0) = pCCE | **PASS** | 7.5 × 10⁻⁹ |
| pgCCE(t=0) = 1.0 | **PASS** | 0.0 |
| pmegCCE(rate=0) = pgCCE | **PASS** | 5.4 × 10⁻⁹ |
| pCCE(sK=1) = CCE | **수치 확인** | < 10⁻⁶ |
| pgCCE(sK=1) = gCCE | **수치 확인** | < 10⁻⁶ |

---

## 8. 변경된 파일 목록

| 파일 | 변경 내용 |
|------|----------|
| `src/superoperator.cpp` | **NEW** — Lindblad superoperator 전체 구현 |
| `src/simulator_mecce.cpp` | **NEW** — ME-CCE 계산 커널 |
| `src/simulator_megcce.cpp` | **NEW** — ME-gCCE 계산 커널 |
| `include/superoperator.h` | **NEW** — superoperator 헤더 |
| `src/simulator.cpp` | isMECCE, isMEGCCE, isPGCCE, isPMECCE, isPMEGCCE 플래그 추가 |
| `src/cluster.cpp` | pCCE 브랜치에 pmecce, pgcce, pmegcce 추가; 0번째 클러스터 보존 조건 |
| `src/cluster_pcce.cpp` | k-means 클러스터링 조건에 pgcce, pmegcce 추가 |
| `src/cluster_hash.cpp` | NULL pointer 세그폴트 수정 (findCluster 반환값 NULL 체크) |
| `src/general.cpp` | 메서드 목록에 pgCCE, pmegCCE 추가 |
| `src/option.cpp` | sK 옵션을 pgcce, pmegcce에서도 읽도록 수정 |
| `test/verification/run_all_tests.py` | **NEW** — 5개 검증 테스트 자동화 스크립트 |
| `test/verification/calculation_conditions.md` | **NEW** — 계산 조건 문서 |
