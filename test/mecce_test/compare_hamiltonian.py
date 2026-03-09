"""
Compare single-cluster Hamiltonian between PyCCE and CCEX.
All computations done manually to avoid numba issues.
"""
import pycce
import numpy as np
from pycce.constants import HBAR_MU0_O4PI, PI2

np.set_printoptions(precision=10, linewidth=200, suppress=True)

# Setup
atoms = pycce.random_bath('e', [5e3,5e3,5e3], density=1e16, density_units='cm-3', seed=2)
cen = pycce.CenterArray(spin=0.5, alpha=[1, 0], beta=[0, 1])
calc = pycce.Simulator(spin=cen, bath=atoms, order=2, r_bath=1400, r_dipole=700, magnetic_field=300)
calc.generate_clusters()

bath = calc.bath
gyro_e = bath.types['e'].gyro  # -17608.59705
gyro_cen = float(cen.gyro.flatten()[0])  # scalar gyro for isotropic electron
Bz = 300  # Gauss

pair = calc.clusters[2][0]
i0, i1 = pair
xyz0 = bath[i0]['xyz']
xyz1 = bath[i1]['xyz']
print(f"Cluster: spins [{i0}, {i1}]")
print(f"  Spin {i0}: xyz = {xyz0}")
print(f"  Spin {i1}: xyz = {xyz1}")
print(f"  gyro_e = {gyro_e}, gyro_cen = {gyro_cen}")

# Spin-1/2 matrices
Sx = np.array([[0, 0.5], [0.5, 0]], dtype=np.complex128)
Sy = np.array([[0, -0.5j], [0.5j, 0]], dtype=np.complex128)
Sz = np.array([[0.5, 0], [0, -0.5]], dtype=np.complex128)
I2 = np.eye(2, dtype=np.complex128)

def kron(A, B):
    return np.kron(A, B)

# Expanded operators for 2-spin system (4x4)
S0x = kron(Sx, I2); S0y = kron(Sy, I2); S0z = kron(Sz, I2)
S1x = kron(I2, Sx); S1y = kron(I2, Sy); S1z = kron(I2, Sz)

def gen_pos_tensor(r1, r2):
    """PyCCE's position tensor: -(3*r⊗r - I*r²)/r⁵"""
    pos = r1 - r2
    r = np.linalg.norm(pos)
    return -(3 * np.outer(pos, pos) - np.eye(3) * r**2) / r**5

def vec_tensor_vec_manual(ivec0, tensor, ivec1):
    """H = Σ_ij tensor[i,j] * S0_i * S1_j"""
    ops = [ivec0, ivec1]  # Each is [Sx, Sy, Sz]
    H = np.zeros((4, 4), dtype=np.complex128)
    for i in range(3):
        for j in range(3):
            H += tensor[i, j] * ops[0][i] @ ops[1][j]
    return H

ivec0 = [S0x, S0y, S0z]
ivec1 = [S1x, S1y, S1z]

###############################################################################
# PyCCE-style computation (Hamiltonian in kHz)
###############################################################################
print("\n" + "="*60)
print("PyCCE-style: Hamiltonian in kHz")
print("="*60)

prefactor_pycce = HBAR_MU0_O4PI / PI2
print(f"HBAR_MU0_O4PI = {HBAR_MU0_O4PI}")
print(f"Prefactor = HBAR_MU0_O4PI / (2π) = {prefactor_pycce}")

# Bath-bath dipolar tensor (kHz)
dd_pycce = gen_pos_tensor(xyz0, xyz1) * prefactor_pycce * gyro_e * gyro_e
print(f"\nDipolar tensor dd (kHz):")
print(dd_pycce)

# Bath-bath Hamiltonian
H_dd_pycce = vec_tensor_vec_manual(ivec0, dd_pycce, ivec1)

# Zeeman
H_zeeman_pycce = gyro_e * Bz * (S0z + S1z)

# Total bath
H_bath_pycce = H_dd_pycce + H_zeeman_pycce

# Qubit-bath hyperfine (kHz)
A0_pycce = gen_pos_tensor(np.zeros(3), xyz0) * prefactor_pycce * gyro_cen * gyro_e
A1_pycce = gen_pos_tensor(np.zeros(3), xyz1) * prefactor_pycce * gyro_cen * gyro_e

# Project: H_qb = Σ_ij ⟨ψ|S_i|ψ⟩ * A_ij * I_j
# alpha: ⟨α|S_x|α⟩=0, ⟨α|S_y|α⟩=0, ⟨α|S_z|α⟩=+0.5
# beta:  ⟨β|S_x|β⟩=0, ⟨β|S_y|β⟩=0, ⟨β|S_z|β⟩=-0.5
proj_alpha = np.array([0., 0., 0.5])
proj_beta  = np.array([0., 0., -0.5])

def project_qb(A_tensor, proj, ivec):
    """H_qb = Σ_ij proj[i] * A[i,j] * I_j"""
    H = np.zeros((4, 4), dtype=np.complex128)
    for i in range(3):
        for j in range(3):
            H += proj[i] * A_tensor[i, j] * ivec[j]
    return H

H_qb_alpha_pycce = project_qb(A0_pycce, proj_alpha, ivec0) + project_qb(A1_pycce, proj_alpha, ivec1)
H_qb_beta_pycce  = project_qb(A0_pycce, proj_beta, ivec0)  + project_qb(A1_pycce, proj_beta, ivec1)

H_alpha_pycce = H_bath_pycce + H_qb_alpha_pycce
H_beta_pycce  = H_bath_pycce + H_qb_beta_pycce

###############################################################################
# CCEX-style computation (Hamiltonian in radkHz)
###############################################################################
print("\n" + "="*60)
print("CCEX-style: Hamiltonian in radkHz")
print("="*60)

H_BAR = 1.054571729
print(f"H_BAR = {H_BAR}")

# Bath-bath dipolar tensor (radkHz)
dd_ccex = gen_pos_tensor(xyz0, xyz1) * H_BAR * gyro_e * gyro_e
print(f"\nDipolar tensor dd (radkHz):")
print(dd_ccex)

# Bath-bath Hamiltonian
H_dd_ccex = vec_tensor_vec_manual(ivec0, dd_ccex, ivec1)

# Zeeman (radkHz) — CCEX converts Gauss to mT internally?
# Check: CCEX bfield handling
# gyro is in rad/(ms*mT), but bfield is in Gauss. 1 Gauss = 0.1 mT
# PyCCE: gyro is in rad/(ms*G) (since HBAR_MU0_O4PI is defined for Gauss)
# Wait — let me check this carefully
print(f"\n--- Zeeman term analysis ---")
print(f"gyro_e = {gyro_e}")
print(f"  PyCCE: gyro in rad/(ms·G), so Zeeman = gyro * B_Gauss * Sz")
print(f"  CCEX:  gyro in rad/(ms·mT)? or rad/(ms·G)?")
# From gyro file: "e  0.5  -17608.597050"
# PyCCE electron gyro = -17608.59705 — same value
# So both use same gyro units
# PyCCE: gamma_e = ge * MU_B / HBAR = 2.00231930436256 * 9.2740100783e-24 / 1.054571817e-34
#       = 1.760859644e11 rad/(s·T) = 17608.59644 rad/(ms·mT)
# 1 Gauss = 0.1 mT
# gamma_e * 1 Gauss = -17608.597 * 0.1 = -1760.8597 rad/ms = -1760.8597 radkHz
# But PyCCE says this is in kHz... hmm

# Let me check: does PyCCE use gyro in rad/(ms·G)?
# If gyro = -17608.59705 rad/(ms·G), then:
# Zeeman = -17608.59705 * 300 * Sz [kHz] = -5282579.115 * Sz [kHz]
# That's 5.28 MHz = huge, which makes sense for electron Larmor in 300G

# In CCEX: if gyro is also -17608.59705 and bfield=300:
# Does CCEX use bfield in Gauss directly?

H_zeeman_ccex = gyro_e * Bz * (S0z + S1z)  # Same formula, but interpreted as radkHz
# Wait - if PyCCE and CCEX use the same gyro and same B value,
# then Zeeman_PyCCE (kHz) == Zeeman_CCEX (radkHz) numerically
# But physically they can't be equal in different units!
# This means one of them has the wrong units...

print(f"\nZeeman = gyro * B = {gyro_e * Bz}")
print(f"  If this is in kHz:    {gyro_e * Bz} kHz = {gyro_e * Bz * 2 * np.pi} radkHz")
print(f"  If this is in radkHz: {gyro_e * Bz} radkHz = {gyro_e * Bz / (2*np.pi)} kHz")

###############################################################################
# KEY COMPARISON
###############################################################################
print("\n" + "="*60)
print("KEY COMPARISON")
print("="*60)

# Ratio of dipolar tensors
print("\n1. Dipolar tensor ratio (CCEX / PyCCE):")
ratio_dd = dd_ccex / dd_pycce
print(f"   All elements = {ratio_dd[0,0]:.10f}")
print(f"   2π = {2*np.pi:.10f}")
print(f"   Match: {np.allclose(ratio_dd, 2*np.pi)}")

# Zeeman is numerically identical
print(f"\n2. Zeeman (gyro*B): PyCCE and CCEX compute same number: {gyro_e * Bz}")
print(f"   PyCCE says this is kHz, CCEX says this is radkHz")
print(f"   → This is a UNIT INCONSISTENCY in one of the codes!")

# If Zeeman is correct in PyCCE (kHz), then in CCEX it should be 2π times larger (radkHz)
# But both use the same formula: gyro * B
# The question: is gyro in rad/(ms·G) or (cycles/ms·G)?
# If gyro = -17608.59705 is in rad/(ms·G):
#   Then gyro * B is in rad/ms = radkHz  → CCEX is correct
#   PyCCE divides by 2π somewhere to get kHz?

# In PyCCE propagator: exp(-i * H * t * 2π)
# So effectively: exp(-i * (H_kHz * 2π) * t) = exp(-i * H_radkHz * t)
# If PyCCE's H already contains gyro*B (which is in radkHz), and PyCCE multiplies by 2π:
#   exp(-i * gyro*B * Sz * 2π * t) → Wrong! Extra 2π
# Unless PyCCE's gyro is actually in kHz/(G), not rad/(ms·G)

# Let me check: gyro_e = -17608.59705
# True electron gyromagnetic ratio = 1.760859644e11 rad/(s·T)
# In rad/(ms·G): 1.760859644e11 * 1e-3 * 1e-4 = 1.760859644e4 = 17608.59644
# Yes! gyro_e = -17608.59705 IS in rad/(ms·G) ≈ rad/(ms·G)
# So gyro * B = -17608.59705 * 300 = -5282579.115 rad/ms = -5282579.115 radkHz

# But PyCCE treats this as kHz and then multiplies by 2π in propagator:
# exp(-i * H * t * 2π) where H in kHz = (gyro*B)/2π would make sense...
# But PyCCE doesn't divide gyro by 2π!

# Let me check PyCCE's Hamiltonian more carefully
print(f"\n3. PyCCE Hamiltonian unit check:")
print(f"   PyCCE H_zeeman = gyro*B*Sz = {gyro_e * Bz * 0.5:.4f}")
print(f"   If kHz: Larmor frequency = {gyro_e * Bz / (2*np.pi):.4f} kHz (cycles)")
print(f"   If radkHz: Larmor = {gyro_e * Bz:.4f} radkHz (angular)")
print(f"   Physical: e- Larmor at 300G = {abs(gyro_e * Bz / (2*np.pi)):.0f} kHz ≈ {abs(gyro_e * Bz / (2*np.pi))/1e3:.1f} MHz")
# Electron Larmor at 300G ≈ 2.8 MHz/G * 300G ≈ 840 MHz
# But gyro*B/2π = 17608*300/2π ≈ 840,764 kHz ≈ 840 MHz ✓
# So gyro*B is in radkHz, and gyro*B/(2π) is in kHz

# PyCCE: Hamiltonian is labeled as kHz but gyro*B gives radkHz??
# PyCCE then multiplies by 2π in propagator: exp(-i * H * 2π * t)
# If H = gyro*B*Sz (radkHz), exp(-i * radkHz * 2π * t) → extra 2π!

# Unless PyCCE's dipolar term (with /PI2 prefactor) compensates...
print(f"\n4. Dipolar vs Zeeman unit consistency:")
print(f"   PyCCE dipolar prefactor: HBAR_MU0_O4PI/PI2 = {prefactor_pycce:.10f}")
print(f"   PyCCE dipolar: {dd_pycce[2,2]:.10f} (labeled kHz)")
print(f"   CCEX dipolar:  {dd_ccex[2,2]:.10f} (labeled radkHz)")
print(f"   Ratio: {dd_ccex[2,2]/dd_pycce[2,2]:.10f} (= 2π)")

print(f"\n   PyCCE Zeeman:   {gyro_e * Bz:.4f} (labeled kHz)")
print(f"   CCEX Zeeman:    {gyro_e * Bz:.4f} (labeled radkHz)")
print(f"   Ratio: 1.0 (SAME NUMBER!)")

print(f"\n   *** INCONSISTENCY: Dipolar has 2π ratio but Zeeman has ratio 1 ***")
print(f"   This means in PyCCE: Zeeman is actually in radkHz, not kHz!")
print(f"   And PyCCE multiplies ALL of H by 2π in propagator:")
print(f"     exp(-i * H * 2π * t)")
print(f"   So Zeeman gets: exp(-i * gyro*B*Sz * 2π * t) → has extra 2π")
print(f"   And Dipolar gets: exp(-i * (dipolar/2π) * 2π * t) = exp(-i * dipolar * t) → correct")

print(f"\n5. In CCEX:")
print(f"   Zeeman: gyro*B (radkHz)")
print(f"   Dipolar: H_BAR * g1*g2 * geometric (radkHz)")
print(f"   Propagator: exp(-i * H * t) — no extra 2π")
print(f"   Zeeman: exp(-i * gyro*B*Sz * t) → correct")
print(f"   Dipolar: exp(-i * H_BAR*g1*g2*geo * t)")

# The question: does CCEX Zeeman use the right formula?
# In CCEX, how is bfield converted? Let me check if CCEX converts Gauss to mT
# From what we know: CCEX bfield is read in Gauss
# If CCEX uses gyro in rad/(ms·mT) and bfield in Gauss:
#   It should convert: B_mT = B_Gauss * 0.1
#   Then Zeeman = gyro * B_mT = gyro * B_Gauss * 0.1
# But if PyCCE uses gyro in rad/(ms·G):
#   Zeeman = gyro * B_Gauss

print(f"\n6. Checking CCEX bfield conversion:")
print(f"   Does CCEX convert Gauss to mT? If gyro is in rad/(ms·mT):")
print(f"   Zeeman_CCEX = gyro * B_Gauss * 0.1 = {gyro_e * Bz * 0.1:.4f} radkHz")
print(f"   Zeeman_PyCCE = gyro * B_Gauss = {gyro_e * Bz:.4f} 'kHz'")
print(f"   Ratio = 0.1 (factor of 10!)")
