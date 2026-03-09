"""
ME-CCE Master Equation — C++ 코드와 1:1 대응하는 Python 구현
=============================================================

이 파일은 CCEX의 simulator_mecce.cpp + superoperator.cpp 를
Python으로 그대로 옮긴 것입니다.

간단한 모델: 큐빗(S=1/2) + 배스 스핀 1개(S=1/2)
  - FREE evolution (Ramsey, npulse=0): CCE vs ME-CCE(rate=0) 검증
  - Hahn echo (npulse=1) + rate>0: Lindblad dissipation의 물리적 효과 시연

※ 왜 Hahn echo + 단일 스핀으로는 CCE와 ME-CCE(rate=0)을 비교할 수 없나?
   Hahn echo는 단일 스핀의 dephasing을 완벽히 상쇄합니다.
   따라서 단일 스핀 클러스터에서 CCE = ME-CCE(rate=0) = 1.0 (에코 완벽 상쇄).
   Jump operator (rate>0)를 추가하면 에코가 불완전해져 L(t) < 1이 됩니다.
   → Lindblad 효과를 보려면 Hahn echo, 검증(CCE=ME-CCE)은 Ramsey를 사용.

목차:
  Part 1. 수학적 배경
  Part 2. C++ 함수들의 Python 구현 (1:1 대응)
  Part 3. Ramsey — CCE vs ME-CCE(rate=0) 수치 검증
  Part 4. Hahn echo — Lindblad의 물리적 효과
  Part 5. 결과 플롯
"""

import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# ======================================================================
# PART 1: 수학적 배경
# ======================================================================
"""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
【1. Lindblad Master Equation】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  dρ/dt = L[ρ] = -i[H, ρ] + Σ_k D[C_k]ρ

  Lindblad dissipator:
    D[C]ρ = C ρ C† - ½ C†C ρ - ½ ρ C†C

  ρ  : 배스 스핀의 밀도 행렬
  H  : 해밀토니안 (rad/ms 단위)
  C_k: jump operator (소산을 만드는 연산자)

  해석:
  - C ρ C†          : 환경에 의해 스핀이 새로운 상태로 점프
  - -½ C†C ρ        : 점프 전 상태가 줄어드는 효과 (확률 보존)
  - -½ ρ C†C        : (위와 동일, Hermitian 대칭성)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
【2. ME-CCE의 핵심: Projected Coherent Superoperator】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  CCE의 secular approximation에서:

  큐빗 상태로 projection된 해밀토니안:
    H_α = H_bath + ⟨α|H_qb|α⟩   (ms = +1/2)
    H_β = H_bath + ⟨β|H_qb|β⟩   (ms = -1/2)

  실제로 ME-CCE가 전파하는 방정식:

    dρ̃/dt = -i H_α ρ̃ + i ρ̃ H_β + Σ_k D[C_k]ρ̃

  여기서 ρ̃는 full density matrix의 off-diagonal 배스 블록:
    ρ̃ = ⟨α|ρ_total|β⟩

  코히런스: L(t) = Tr(ρ̃(t))

  → 큐빗 DOF를 명시적으로 추적하지 않아도 됨 (secular approx의 이점)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
【3. Superoperator Formalism】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  행렬 방정식을 벡터 방정식으로 변환 → 수치 계산에 편리.

  벡터화 (row-major):
    vec(ρ)[i*N + j] = ρ[i, j]

  핵심 공식:
    ρ → A ρ B   ⟺   vec(ρ) → (A ⊗ B^T) vec(ρ)
    이를 op_to_supop(A, B) = A ⊗ B^T 로 정의.

  코히런트 superoperator:
    L_coh = -i (H ⊗ I - I ⊗ H^T)
    ← dρ/dt = -i[H, ρ] 에 대응

  Projected coherent superoperator (ME-CCE):
    L_proj = -i (H_α ⊗ I - I ⊗ H_β^T)
    ← dρ̃/dt = -i H_α ρ̃ + i ρ̃ H_β 에 대응

  Lindblad dissipator superoperator:
    L_D[C] = (C ⊗ C*) - ½(C†C ⊗ I) - ½(I ⊗ (C†C)^T)
    ← D[C]ρ = C ρ C† - ½ C†C ρ - ½ ρ C†C 에 대응

  시간 전파:
    ρ(t) = vec_to_mat( exp(L·t) · vec(ρ(0)) )

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
【4. Hahn Echo에서 단일 스핀의 특별한 점】
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Hahn echo (π pulse at t/2):
    첫 τ = t/2: ρ̃(τ) = e^{-iH_α τ} ρ̃₀ e^{iH_β τ}
    π pulse: H_α ↔ H_β
    다음 τ: ρ̃(t) = e^{-iH_β τ} ρ̃(τ) e^{iH_α τ}
           = e^{-iH_β τ} e^{-iH_α τ} ρ̃₀ e^{iH_β τ} e^{iH_α τ}

  rate=0이면:
    L(t) = Tr(ρ̃(t)) = (1/N) Tr[e^{-iH_β τ} e^{-iH_α τ} e^{iH_β τ} e^{iH_α τ}]

  ★ 단일 배스 스핀, 대각 H_α, H_β의 경우:
    대각 원소들의 위상이 에코에서 정확히 상쇄 → L(t) = 1.0 for all t
    이것이 Hahn echo의 목적: 단일 스핀 dephasing 상쇄!

  rate>0이면:
    D[C_k]ρ̃ 항이 추가됨 → 에코가 불완전해짐 → L(t) < 1
    이것이 ME-CCE의 새로운 물리: Lindblad jump가 에코 상쇄를 방해

  ∴ CCE = ME-CCE(rate=0) 비교는 자유 진화(Ramsey)에서 해야 함
     Hahn echo에서는 rate>0 효과를 확인하는 용도로 사용.
"""

# ======================================================================
# PART 2: C++ 함수들의 Python 구현 (1:1 대응)
# ======================================================================

# ─── 벡터화 유틸리티 ───────────────────────────────────────────────────

def op_to_supop(A, B):
    """
    C++: op_to_supop(A, B) = kron(A, B^T)

    [superoperator.cpp line 64-68]
    ρ → A ρ B 연산의 superoperator 표현.
    Row-major 벡터화: kron(A, B^T) · vec(ρ) = vec(A ρ B)

    증명:
    [kron(A, B^T)]_{(i·N+j), (k·N+l)} = A_{ik} · (B^T)_{jl} = A_{ik} · B_{lj}
    [kron(A, B^T) · vec(ρ)]_{i·N+j} = Σ_{k,l} A_{ik} B_{lj} ρ_{kl} = [A ρ B]_{ij}  ✓
    """
    return np.kron(A, B.T)


def mat_to_vec(rho):
    """
    C++: mat_to_vec(mat) [superoperator.cpp line 70-80]

    행렬을 row-major 순서로 벡터화.
    vec[i*N + j] = rho[i, j]
    """
    return rho.flatten()   # numpy default: C-order (row-major)


def vec_to_mat(vec):
    """
    C++: vec_to_mat(vec) [superoperator.cpp line 82-92]

    벡터를 다시 행렬로 복원.
    """
    N = int(np.round(np.sqrt(len(vec))))
    return vec.reshape(N, N)


# ─── Superoperator 구성 ───────────────────────────────────────────────

def projected_coherent_superoperator(H_alpha, H_beta):
    """
    C++: projected_coherent_superoperator(Ha, Hb)
    [superoperator.cpp line 100-104]

    ME-CCE의 코히런트 부분:
      L_proj = -i (H_α ⊗ I - I ⊗ H_β^T)

    이 연산자를 ρ̃에 적용하면:
      dρ̃/dt = -i H_α ρ̃ + i ρ̃ H_β

    (op_to_supop(Ha, I) maps ρ → H_α ρ)
    (op_to_supop(I, Hb) maps ρ → ρ H_β)
    """
    N = H_alpha.shape[0]
    I = np.eye(N, dtype=complex)
    return -1j * (op_to_supop(H_alpha, I) - op_to_supop(I, H_beta))


def collapse_superoperator(C):
    """
    C++: collapse_superoperator_single(C)
    [superoperator.cpp line 106-119]

    Lindblad dissipator의 superoperator:
      D[C]ρ = C ρ C† - ½ C†C ρ - ½ ρ C†C

    각 항의 superoperator:
      C ρ C†    → op_to_supop(C, C†) = kron(C, (C†)^T) = kron(C, C*)
      ½ C†C ρ   → ½ op_to_supop(C†C, I)
      ½ ρ C†C   → ½ op_to_supop(I, C†C)

    왜 op_to_supop(C, C†) = kron(C, C*)?
      (C†)^T = (C^*)^T^T = C^*  (Hermitian conjugate의 transpose)
      kron(C, (C†)^T) = kron(C, C*)  ✓
    """
    N = C.shape[0]
    I = np.eye(N, dtype=complex)
    CdagC = C.conj().T @ C         # C† C

    term1 = op_to_supop(C, C.conj().T)    # C ⊗ C*   → C ρ C†
    term2 = 0.5 * op_to_supop(CdagC, I)  # ½ C†C ⊗ I → ½ C†C ρ
    term3 = 0.5 * op_to_supop(I, CdagC)  # ½ I ⊗ C†C → ½ ρ C†C

    return term1 - term2 - term3


# ─── Jump operator ───────────────────────────────────────────────────

def spin_operators(S=0.5):
    """
    스핀 S에 대한 Sx, Sy, Sz, S+, S- 행렬 생성.

    기저 순서: |m=S⟩, |m=S-1⟩, ..., |m=-S⟩  (위에서 아래)
    """
    dim = int(2*S + 1)
    Sz = np.zeros((dim, dim), dtype=complex)
    Sm = np.zeros((dim, dim), dtype=complex)   # S- (내림 연산자)

    for row in range(dim):
        m = S - row
        Sz[row, row] = m
        if row + 1 < dim:
            m_val = S - row   # |m⟩ 상태의 m 값
            # ⟨m-1|S-|m⟩ = sqrt(S(S+1) - m(m-1))
            Sm[row+1, row] = np.sqrt(S*(S+1) - m_val*(m_val-1))

    Sp = Sm.conj().T   # S+ = (S-)†
    Sx = (Sp + Sm) / 2
    Sy = (Sp - Sm) / (2j)

    return Sx, Sy, Sz, Sp, Sm


def get_jump_operator(op_type, rate, S=0.5):
    """
    C++: get_jump_operator_matrix(op_type, dim, rate)
    [superoperator.cpp line 121-187]

    op_type:
      "-+" → √Γ · S+  (올림 연산자: |↓⟩→|↑⟩ flip-up)
      "+-" → √Γ · S-  (내림 연산자: |↑⟩→|↓⟩ flip-down)
      "z"  → √Γ · Sz  (위상 감쇄, dephasing)

    rate는 rad/ms 단위. C++ 코드에서 rate는 radkHz = rad/ms.
    실제 rate의 단위는 [rate] = 1/s 이고, Γ = rate.
    Jump operator: C = √Γ · (S+, S-, 또는 Sz)
    Lindblad: D[C]ρ = Γ(SρS† - ½S†Sρ - ½ρS†S)
    """
    _, _, Sz, Sp, Sm = spin_operators(S)
    sqrt_rate = np.sqrt(rate)

    if op_type == "-+":
        return sqrt_rate * Sp   # √Γ · S+
    elif op_type == "+-":
        return sqrt_rate * Sm   # √Γ · S-
    elif op_type == "z":
        return sqrt_rate * Sz   # √Γ · Sz
    else:
        raise ValueError(f"Unknown op_type: {op_type}")


# ─── Lindblad superoperator 구성 ─────────────────────────────────────

def build_lindbladian(H_alpha, H_beta, jump_ops_list):
    """
    C++: simulator_mecce.cpp 의 L_coh + L_incoh 구성

    파라미터:
      H_alpha       : α 상태 projected 해밀토니안 (배스 스핀 공간)
      H_beta        : β 상태 projected 해밀토니안
      jump_ops_list : [C_1, C_2, ...] jump operator 행렬 리스트

    반환:
      L        : 전체 Lindblad superoperator
      L_flipped: α↔β 교환된 superoperator (Hahn echo π pulse 후)
    """
    N = H_alpha.shape[0]

    # ── 코히런트 부분 ────────────────────────────────────
    # C++: L_coh = projected_coherent_superoperator(Halpha, Hbeta)
    L_coh = projected_coherent_superoperator(H_alpha, H_beta)

    # ── 비코히런트 부분 (dissipation) ────────────────────
    # C++: L_incoh = incoherent_superoperator(ba, bsigmas, joa, bdim)
    L_incoh = np.zeros((N*N, N*N), dtype=complex)
    for C in jump_ops_list:
        # C++: collapse_superoperator_single(C_full)
        L_incoh += collapse_superoperator(C)

    L = L_coh + L_incoh

    # ── Flipped: π pulse 후 (H_α ↔ H_β 교환) ────────────
    # C++: L_flipped = projected_coherent_superoperator(Hbeta, Halpha) + L_incoh
    L_coh_flipped = projected_coherent_superoperator(H_beta, H_alpha)
    L_flipped = L_coh_flipped + L_incoh

    return L, L_flipped


# ─── 시간 전파 ─────────────────────────────────────────────────────────

def propagate_ramsey(L, t, rho0):
    """
    Ramsey (자유 진화, npulse=0): ρ̃(t) = exp(L·t) ρ̃₀

    C++: simple_incoherent_propagator(t, L) = expm(L * t)
    → ρ̃(t)의 trace = 코히런스 L(t)
    """
    if t == 0.0:
        return np.trace(rho0)
    prop = expm(L * t)               # C++: expm(L * t)
    rho_vec = prop @ mat_to_vec(rho0)
    rho_t = vec_to_mat(rho_vec)
    return np.trace(rho_t)


def propagate_hahn_echo(L, L_flipped, t, rho0):
    """
    Hahn Echo (npulse=1): π pulse at t/2

    시퀀스:
      1단계: 자유 진화 τ=t/2  → U_before = exp(L · τ)
      π pulse: H_α ↔ H_β     → L → L_flipped
      2단계: 자유 진화 τ=t/2  → U_after = exp(L_flipped · τ)

    전체 전파자: S_total = U_after · U_before

    C++: [simulator_mecce.cpp, npulse==1 분기]
      U_before = simple_incoherent_propagator(tau, L)
      U_after  = simple_incoherent_propagator(tau, L_flipped)
      super_prop = propagate_superpropagators(U_before, U_after, npulse=1)
               = U_after * U_before
    """
    if t == 0.0:
        return np.trace(rho0)

    tau = t / 2.0

    # C++: simple_incoherent_propagator(tau, L) = expm(L * tau)
    U_before = expm(L * tau)
    U_after  = expm(L_flipped * tau)

    # C++: propagate_superpropagators(U_before, U_after, 1) = U_after * U_before
    super_prop = U_after @ U_before

    rho_vec = super_prop @ mat_to_vec(rho0)
    rho_t = vec_to_mat(rho_vec)
    return np.trace(rho_t)


# ─── CCE 코히런스 (해석적 계산) ─────────────────────────────────────────

def cce_ramsey(H_alpha, H_beta, t, rho0):
    """
    CCE (Lindblad 없음) Ramsey 코히런스.

    L_CCE(t) = Tr[ρ₀ · e^{iH_β t} · e^{-iH_α t}]  ... *

    * 이것은 ME-CCE Ramsey와 동일한 결과:
      ME-CCE: L(t) = Tr[e^{-iH_α t} ρ₀ e^{iH_β t}]
                   = Tr[ρ₀ · e^{iH_β t} · e^{-iH_α t}]  (cyclic trace)  ✓
    """
    if t == 0.0:
        return 1.0 + 0j
    M = expm(1j * H_beta * t) @ expm(-1j * H_alpha * t)
    return np.trace(rho0 @ M)


def cce_hahn_echo(H_alpha, H_beta, t, rho0):
    """
    CCE (Lindblad 없음) Hahn echo 코히런스.

    단일 스핀, 대각 H_α, H_β의 경우 Hahn echo에서:
    L_CCE_echo(t) = 1.0  (에코가 단일 스핀 dephasing을 완벽히 상쇄)

    일반적으로:
    L_CCE_echo(t) = Tr[ρ₀ · e^{-iH_α τ}e^{-iH_β τ}e^{iH_α τ}e^{iH_β τ}]
    where τ = t/2

    ※ H_α와 H_β가 교환되지 않는 일반적인 경우에는 1 이 아닐 수 있음.
       이 모델에서는 H_α = diag(h_α, -h_α), H_β = diag(h_β, -h_β) 로
       둘 다 대각이므로 교환됨 → L=1.
    """
    if t == 0.0:
        return 1.0 + 0j
    tau = t / 2.0
    # Hahn echo CCE formula
    M = (expm(-1j * H_alpha * tau)
         @ expm(-1j * H_beta  * tau)
         @ expm( 1j * H_alpha * tau)
         @ expm( 1j * H_beta  * tau))
    return np.trace(rho0 @ M)


# ======================================================================
# PART 3: 간단한 모델 설정
# ======================================================================

def build_model(A_zz_radms, B_G=300.0, gyro_e=-17.608597):
    """
    배스 스핀 1개에 대한 해밀토니안을 구성합니다.

    물리:
      배스 Larmor 항: H_bath = ω_bath · Sz  (ω_bath = gyro_e · B)
      큐빗-배스 secular 결합: H_qb ~ A_zz · Sz_qubit ⊗ Sz_bath
        → 큐빗 α(ms=+1/2): H_qb → +A_zz/2 · Sz
        → 큐빗 β(ms=-1/2): H_qb → -A_zz/2 · Sz

    파라미터:
      A_zz_radms : A_zz [rad/ms]
      B_G        : 자기장 [G]
      gyro_e     : 전자 자이로자기비 [rad/(ms·G)]
                   (CCEX에서 -17608.597 rad/(ms·G) 이지만 ×10⁻³ 단위로 입력)

    반환:
      H_alpha, H_beta: 2×2 복소 행렬 [rad/ms 단위]
    """
    _, _, Sz, _, _ = spin_operators(S=0.5)

    omega_bath = gyro_e * B_G          # rad/ms
    H_bath = omega_bath * Sz

    H_qb_alpha = (A_zz_radms / 2.0) * Sz
    H_qb_beta  = -(A_zz_radms / 2.0) * Sz

    H_alpha = H_bath + H_qb_alpha
    H_beta  = H_bath + H_qb_beta

    return H_alpha, H_beta


# ======================================================================
# PART 4: 메인 계산
# ======================================================================

def run_all():

    print("=" * 62)
    print("  ME-CCE Python Tutorial — C++ superoperator.cpp 구현 확인")
    print("=" * 62)

    # ── 모델 파라미터 ──────────────────────────────────────────────────
    B_G     = 300.0         # 자기장 [G]
    gyro_e  = -17.608597    # 전자 자이로자기비 [rad/(ms·G)]  (≈ -17608.597 e-3)
    A_zz    = 0.05          # 하이퍼파인 결합 [rad/ms]
                            # (실제 전자-전자 쌍은 수 rad/ms 수준, 여기선 작은 값 사용)
    Gamma   = 2.0           # Jump rate [rad/ms]  (이 단위에서 1/(2π·1000) ≠ 1 이지만
                            #  물리 이해를 위해 직접 rad/ms 사용; T1 ≈ 1/Gamma)

    # 시간 그리드
    nstep  = 51
    deltat = 0.04           # [ms]
    times  = np.arange(1, nstep+1) * deltat

    print(f"\n  모델 파라미터")
    print(f"  {'자기장':20s} B     = {B_G} G")
    print(f"  {'Larmor 주파수':20s} ω     = {gyro_e*B_G:.4f} rad/ms")
    print(f"  {'HF 결합':20s} A_zz  = {A_zz} rad/ms")
    print(f"  {'Jump rate (S+,S-)':20s} Gamma = {Gamma} rad/ms")
    print(f"  {'시간 범위':20s} t     = 0 ~ {times[-1]:.2f} ms")

    # ── 해밀토니안 ──────────────────────────────────────────────────────
    H_alpha, H_beta = build_model(A_zz, B_G, gyro_e)

    print(f"\n  H_alpha (projected, α state) =")
    for row in H_alpha:
        print(f"    {np.array2string(row, precision=4)}")
    print(f"\n  H_beta  (projected, β state) =")
    for row in H_beta:
        print(f"    {np.array2string(row, precision=4)}")

    # ── Jump operators ──────────────────────────────────────────────────
    # C++ code: get_jump_operator_matrix("-+", dim=2, rate) → √Γ · S+
    C_plus  = get_jump_operator("-+", Gamma)   # √Γ · S+
    C_minus = get_jump_operator("+-", Gamma)   # √Γ · S-

    print(f"\n  Jump operator C+ = sqrt({Gamma}) · S+ =")
    print(f"    {C_plus}")
    print(f"  Jump operator C- = sqrt({Gamma}) · S- =")
    print(f"    {C_minus}")

    # ── Lindblad superoperator ──────────────────────────────────────────
    # rate=0 (검증용)
    L_zero,      L_flip_zero      = build_lindbladian(H_alpha, H_beta, [])
    # rate=Gamma
    L_with_rate, L_flip_with_rate = build_lindbladian(H_alpha, H_beta, [C_plus, C_minus])

    print(f"\n  Superoperator 크기: {L_zero.shape}  (N²×N², N=2)")
    print(f"\n  L_coh (rate=0) =")
    print(np.array2string(L_zero, precision=4, suppress_small=True))

    # ── 초기 상태 ──────────────────────────────────────────────────────
    # C++: rho0 = MatrixXcd::Identity(bdim, bdim) / doublec(bdim, 0.0)
    N = 2
    rho0 = np.eye(N, dtype=complex) / N
    print(f"\n  초기 상태 ρ̃₀ = I/2 (ensemble average):")
    print(f"    {rho0}")

    # ═══════════════════════════════════════════════════════════════════
    # 계산 1: RAMSEY (자유 진화) — CCE vs ME-CCE(rate=0) 검증
    # ═══════════════════════════════════════════════════════════════════
    print("\n" + "─" * 62)
    print("  계산 1: Ramsey (자유 진화) — CCE vs ME-CCE(rate=0)")
    print("─" * 62)
    print("  [이론] ME-CCE(rate=0) = CCE for free evolution")
    print("  L_CCE(t)   = Tr[ρ₀ · e^{iH_β t} · e^{-iH_α t}]")
    print("  L_MECCE(t) = Tr[ e^{-iH_α t} ρ̃₀ e^{iH_β t} ]")
    print("  Cyclic trace → 두 식은 동일 ✓")

    L_cce_ramsey      = np.zeros(len(times), dtype=complex)
    L_mecce_r0_ramsey = np.zeros(len(times), dtype=complex)
    L_mecce_ramsey    = np.zeros(len(times), dtype=complex)

    for i, t in enumerate(times):
        L_cce_ramsey[i]      = cce_ramsey(H_alpha, H_beta, t, rho0)
        L_mecce_r0_ramsey[i] = propagate_ramsey(L_zero,      t, rho0)
        L_mecce_ramsey[i]    = propagate_ramsey(L_with_rate, t, rho0)

    diff_ramsey = np.abs(L_cce_ramsey.real - L_mecce_r0_ramsey.real)
    print(f"\n  검증 결과: max|CCE - ME-CCE(rate=0)| = {diff_ramsey.max():.3e}")
    print(f"  → {'PASS ✓' if diff_ramsey.max() < 1e-10 else 'FAIL ✗'}  (기준: < 1e-10)")

    print(f"\n  {'t [ms]':>8}  {'CCE':>14}  {'MECCE(r=0)':>14}  {'MECCE(r>0)':>14}  {'|diff|':>12}")
    for i in range(0, len(times), 10):
        print(f"  {times[i]:8.2f}  "
              f"{L_cce_ramsey[i].real:+14.10f}  "
              f"{L_mecce_r0_ramsey[i].real:+14.10f}  "
              f"{L_mecce_ramsey[i].real:+14.10f}  "
              f"{diff_ramsey[i]:.6e}")

    # ═══════════════════════════════════════════════════════════════════
    # 계산 2: HAHN ECHO — Lindblad의 물리적 효과
    # ═══════════════════════════════════════════════════════════════════
    print("\n" + "─" * 62)
    print("  계산 2: Hahn Echo — Lindblad 효과")
    print("─" * 62)
    print("  [이론] 단일 스핀, rate=0 → L=1.0 (에코 완벽 상쇄)")
    print("         Jump operator (rate>0) → 에코가 불완전 → L<1")

    L_cce_echo      = np.zeros(len(times), dtype=complex)
    L_mecce_r0_echo = np.zeros(len(times), dtype=complex)
    L_mecce_echo    = np.zeros(len(times), dtype=complex)

    for i, t in enumerate(times):
        L_cce_echo[i]      = cce_hahn_echo(H_alpha, H_beta, t, rho0)
        L_mecce_r0_echo[i] = propagate_hahn_echo(L_zero,      L_flip_zero,      t, rho0)
        L_mecce_echo[i]    = propagate_hahn_echo(L_with_rate, L_flip_with_rate, t, rho0)

    print(f"\n  {'t [ms]':>8}  {'CCE echo':>12}  {'MECCE(r=0)':>12}  {'MECCE(r>0)':>12}")
    for i in range(0, len(times), 10):
        print(f"  {times[i]:8.2f}  "
              f"{L_cce_echo[i].real:12.8f}  "
              f"{L_mecce_r0_echo[i].real:12.8f}  "
              f"{L_mecce_echo[i].real:12.8f}")

    print(f"\n  CCE echo max deviation from 1.0: {np.max(np.abs(L_cce_echo.real - 1.0)):.2e}")
    print(f"  MECCE(r=0) echo max deviation from 1.0: {np.max(np.abs(L_mecce_r0_echo.real - 1.0)):.2e}")
    print(f"  MECCE(r>0) echo at t=2.0 ms: L = {L_mecce_echo[-1].real:.6f}")
    print(f"  → Jump operator가 에코를 방해: L < 1  (추가 디코히런스)")

    return times, (L_cce_ramsey, L_mecce_r0_ramsey, L_mecce_ramsey,
                   L_cce_echo, L_mecce_r0_echo, L_mecce_echo)


# ======================================================================
# PART 5: 수치 검증 — Superoperator 각 항 확인
# ======================================================================

def verify_superoperator_math():
    """
    Superoperator의 각 항이 올바르게 구현되었는지 수치 검증.
    C++ 코드와 1:1 대응 확인.
    """
    print("\n" + "=" * 62)
    print("  Superoperator 수학적 검증")
    print("=" * 62)

    # ── 검증 1: op_to_supop ─────────────────────────────────────────────
    print("\n[1] op_to_supop(A, B): ρ → A ρ B")
    A = np.array([[1, 2], [3, 4]], dtype=complex)
    B = np.array([[0, 1], [1, 0]], dtype=complex)
    rho = np.array([[0.5, 0.3], [0.3, 0.5]], dtype=complex)

    supop = op_to_supop(A, B)
    result_supop = vec_to_mat(supop @ mat_to_vec(rho))
    result_direct = A @ rho @ B

    print(f"   직접 계산 A·ρ·B:\n{result_direct}")
    print(f"   Superoperator 결과:\n{result_supop}")
    print(f"   최대 오차: {np.abs(result_supop - result_direct).max():.2e}  ✓" if
          np.abs(result_supop - result_direct).max() < 1e-12 else "   FAIL!")

    # ── 검증 2: coherent_superoperator ─────────────────────────────────
    print("\n[2] L_coh: dρ/dt = -i[H, ρ]")
    H = np.array([[1, 0], [0, -1]], dtype=complex)   # Sz
    rho = np.array([[0.5, 0.3+0.1j], [0.3-0.1j, 0.5]], dtype=complex)
    L_coh = projected_coherent_superoperator(H, H)   # H_α=H_β → same as full coherent

    drho_supop  = vec_to_mat(L_coh @ mat_to_vec(rho))
    drho_direct = -1j * (H @ rho - rho @ H)
    print(f"   최대 오차 (L_coh vs -i[H,ρ]): {np.abs(drho_supop - drho_direct).max():.2e}  ✓" if
          np.abs(drho_supop - drho_direct).max() < 1e-12 else "   FAIL!")

    # ── 검증 3: collapse_superoperator ─────────────────────────────────
    print("\n[3] D[C]ρ = C ρ C† - ½ C†C ρ - ½ ρ C†C")
    _, _, _, Sp, Sm = spin_operators(0.5)
    C = np.sqrt(0.5) * Sm    # √(0.5) · S-
    rho = np.array([[0.7, 0.2], [0.2, 0.3]], dtype=complex)

    L_D = collapse_superoperator(C)
    drho_supop  = vec_to_mat(L_D @ mat_to_vec(rho))
    drho_direct = C @ rho @ C.conj().T - 0.5*(C.conj().T @ C @ rho) - 0.5*(rho @ C.conj().T @ C)
    print(f"   최대 오차 (D[C]ρ): {np.abs(drho_supop - drho_direct).max():.2e}  ✓" if
          np.abs(drho_supop - drho_direct).max() < 1e-12 else "   FAIL!")

    # ── 검증 4: Trace 보존 ──────────────────────────────────────────────
    print("\n[4] Trace 보존 확인 (for Lindblad equation)")
    print("   Lindblad이 trace를 보존하려면 Tr[L[ρ]] = 0 이어야 함")
    _, _, Sz, Sp, Sm = spin_operators(0.5)
    H = 2.0 * Sz
    C1 = np.sqrt(0.3) * Sp
    C2 = np.sqrt(0.5) * Sm
    L_full_coh = projected_coherent_superoperator(H, H)
    L_full_D   = collapse_superoperator(C1) + collapse_superoperator(C2)
    L_full = L_full_coh + L_full_D

    rho_test = np.array([[0.6, 0.2+0.1j], [0.2-0.1j, 0.4]], dtype=complex)
    drho = vec_to_mat(L_full @ mat_to_vec(rho_test))
    print(f"   Tr[L[ρ]] = {np.trace(drho):.2e}  {'✓' if abs(np.trace(drho)) < 1e-12 else '✗'}")

    # ── 검증 5: exp(L*t) vs analytic solution ──────────────────────────
    print("\n[5] 해석 해와 비교 (단순 dephasing model)")
    print("   H=0, C = √Γ · Sz → dρ_{01}/dt = -Γ ρ_{01}")
    print("   → ρ_{01}(t) = ρ_{01}(0) · exp(-Γ t)")
    Gamma_test = 1.5   # rad/ms
    C_z = get_jump_operator("z", Gamma_test)
    # H=0이면 L_coh=0
    L_dephase, _ = build_lindbladian(
        np.zeros((2,2), dtype=complex),
        np.zeros((2,2), dtype=complex),
        [C_z]
    )
    t_test = 0.5   # ms
    prop = expm(L_dephase * t_test)
    rho_init = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)  # |+⟩⟨+|
    rho_t = vec_to_mat(prop @ mat_to_vec(rho_init))
    # C = √Γ Sz: D[C]ρ_01 = Γ(Sz_00·Sz_11 - ½(Sz²_00 + Sz²_11))·ρ_01
    #                      = Γ(1/2·(-1/2) - ½(1/4+1/4))·ρ_01 = -Γ/2 · ρ_01
    # → ρ_01(t) = ρ_01(0) · exp(-Γ/2 · t)
    analytic_offdiag = 0.5 * np.exp(-Gamma_test / 2.0 * t_test)
    print(f"   수치: ρ₁₂(t={t_test}) = {rho_t[0,1].real:.8f}")
    print(f"   해석: ρ₁₂(t={t_test}) = {analytic_offdiag:.8f}  [exp(-Γ/2 · t)]")
    print(f"   오차: {abs(rho_t[0,1].real - analytic_offdiag):.2e}  "
          f"{'✓' if abs(rho_t[0,1].real - analytic_offdiag) < 1e-6 else '✗'}")


# ======================================================================
# PART 6: 플롯
# ======================================================================

def plot_all(times, data):
    L_cce_ram, L_r0_ram, L_rm, L_cce_echo, L_r0_echo, L_echo = data

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle('ME-CCE Tutorial: 1 qubit + 1 bath spin (S=1/2)',
                 fontsize=13)

    # Panel 1: Ramsey — 검증
    ax = axes[0]
    ax.plot(times, L_cce_ram.real,  'b-',  lw=2.5, label='CCE (analytic)')
    ax.plot(times, L_r0_ram.real,   'r--', lw=2,   label='ME-CCE (rate=0)')
    ax.plot(times, L_rm.real,       'g-.', lw=2,   label=f'ME-CCE (rate=2.0 rad/ms)')
    ax.set_xlabel('Time (ms)'); ax.set_ylabel('Re[L(t)]')
    ax.set_title('Ramsey (free evolution)\nCCE vs ME-CCE verification')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_xlim([0, times[-1]])

    # Panel 2: Hahn Echo — 에코 상쇄와 Lindblad 효과
    ax = axes[1]
    ax.plot(times, np.ones(len(times)), 'b-',  lw=2.5, label='CCE = ME-CCE(rate=0) = 1.0')
    ax.plot(times, L_echo.real,          'g-.', lw=2,   label=f'ME-CCE (rate=2.0 rad/ms)')
    ax.set_xlabel('Time (ms)'); ax.set_ylabel('Re[L(t)]')
    ax.set_title('Hahn Echo (npulse=1)\nLindblad breaks echo refocusing')
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_xlim([0, times[-1]])
    ax.set_ylim([0, 1.05])

    # Panel 3: Lindblad 이론 다이어그램
    ax = axes[2]
    ax.axis('off')
    diagram = [
        "Lindblad Master Equation",
        "",
        "  drho/dt = L[rho]",
        "",
        "  L[rho] = -i[H, rho]     (coherent)",
        "         + D[C+]rho        (flip-up)",
        "         + D[C-]rho        (flip-down)",
        "",
        "  D[C]rho = C rho C† ",
        "          - 1/2 C†C rho",
        "          - 1/2 rho C†C",
        "",
        "Superoperator (vec form):",
        "  L_coh = -i(Ha x I - I x Hb^T)",
        "  L_D   = C x C* - 1/2 C†C x I",
        "                 - 1/2 I x (C†C)^T",
        "",
        "Propagator: exp(L * t)",
    ]
    ax.text(0.05, 0.95, '\n'.join(diagram),
            transform=ax.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    out = '/Users/istblue/Desktop/codes/CCEX/docs/lindblad_tutorial_result.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Plot saved: {out}")


# ======================================================================
# MAIN
# ======================================================================

if __name__ == "__main__":

    # 수학적 검증 먼저
    verify_superoperator_math()

    # 메인 계산
    times, data = run_all()

    # 플롯
    plot_all(times, data)

    # C++ ↔ Python 대응 요약
    print("\n" + "=" * 62)
    print("  C++ ↔ Python 함수 대응 요약")
    print("=" * 62)
    table = [
        ("C++ (superoperator.cpp)",             "Python",                       "기능"),
        ("─" * 32,                              "─" * 28,                       "─" * 22),
        ("op_to_supop(A, B)",                   "op_to_supop(A, B)",            "kron(A, B^T)"),
        ("mat_to_vec(mat)",                     "mat_to_vec(rho)",              "행렬→벡터"),
        ("vec_to_mat(vec)",                     "vec_to_mat(vec)",              "벡터→행렬"),
        ("projected_coherent_superoperator",    "projected_coherent_...(Ha,Hb)","ME-CCE 코히런트"),
        ("collapse_superoperator_single(C)",    "collapse_superoperator(C)",    "D[C]ρ dissipator"),
        ("get_jump_operator_matrix(t,d,r)",     "get_jump_operator(type, rate)","S+/S-/Sz 연산자"),
        ("simple_incoherent_propagator(t,L)",   "expm(L * t)",                  "exp(L·t) 전파자"),
        ("propagate_superpropagators(Ub,Ua,1)", "U_after @ U_before",           "Hahn echo 합성"),
    ]
    for c, p, f in table:
        print(f"  {c:<38}  {p:<30}  {f}")

    print("\n  Done!")
