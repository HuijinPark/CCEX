"""gCCE order 1 of an S = 1 qubit with a non-axial zero-field splitting (E != 0) and one 13C, against the exact answer (A91).
usage: gcce_zfs_exact.py write PHI IN.JSON BATHFILE     write the input (the system turned by PHI degrees about z)
       gcce_zfs_exact.py check PHI OUT_PREFIX            compare OUT_PREFIX_noDiv / _wiDiv with the exact coherence

The qubit at the origin, ZFS D = 2.87 GHz, E = 100 MHz (the tensor R diag(-D/3+E, -D/3-E, 2D/3) R^T, kHz), one 13C at
R (3, 1, 2) A (point dipole), B = 5 G along z, ensemble over the 13C, free induction, alpha_state ms [1], beta_state ms [0], 200 points to
0.1 ms. Order 1 with one bath spin is the whole system, so gCCE is exact; the exact coherence is written out here from
the full 6 x 6 Hamiltonian (S.T.S - gamma_e B S_z + S.A.I - gamma_n B I_z, A = H_BAR gamma_e gamma_n / r^3 (1 - 3 n n),
rad/ms), rho_0 = |psi0><psi0| x I/2 with psi0 = (|+1> + |0>)/sqrt 2, L(t) = <+1| Tr_n rho(t) |0> / L(0).
phi = 0 and 30: the turned system keeps the same m_s states, a different preparation phase relative to the ZFS axes
(check_rotation/check_zrot_gcce.py), and its exact |L| differs from phi = 0 by up to 3.6e-5; CCEX has to follow the exact
answer for each phi. The tolerance, 2e-8, is rounding amplified by the GHz qubit phase (measured 5e-10 and 7.6e-10).
Exit 1 above the tolerance, 2 when the output is missing or does not hold the grid."""
import json, sys
import numpy as np

H_BAR, GAMMA_E, GAMMA_N, B = 1.054571729, -17608.59705, 6.728284, 5.0
D, E, SPIN, NSTEP, DELTAT, TOL = 2870000.0, 100000.0, np.array([3.0, 1.0, 2.0]), 200, 0.0005, 2e-8
mode, phi = sys.argv[1], np.radians(float(sys.argv[2]))
R = np.array([[np.cos(phi), -np.sin(phi), 0], [np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
T, pos = R @ np.diag([-D/3 + E, -D/3 - E, 2*D/3]) @ R.T, R @ SPIN

if mode == 'write':
    inp = {'method': 'gCCE', 'quantity': 'coherence',
           'Qubit': {'qubit': [{'name': 'q', 'spin': 1.0, 'gyro': GAMMA_E, 'xyz': [0.0, 0.0, 0.0],
                                'detuning': 0.0}], 'alpha_state': {'ms': [1]}, 'beta_state': {'ms': [0]},
                     'intmap': [{'between': ['q', 'q'], 'tensor': T.tolist()}]},
           'configuration': {'gyrofile': 'gyro', 'bathfile': [sys.argv[4].split('/')[-1]]},
           'order': 1, 'bfield': B, 'rbath': 20, 'rdip': 10, 'deltat': DELTAT, 'nstep': NSTEP,
           'pulse': {'count': 0}, 'output': {'file': sys.argv[3].split('/')[-1][:-5]}}
    json.dump(inp, open(sys.argv[3], 'w'), indent=1)
    open(sys.argv[4], 'w').write(f'1.0\t0.0\t0.0\t0.0\n{pos[0]:.15f}\t{pos[1]:.15f}\t{pos[2]:.15f}\t13C\n')
    sys.exit(0)


def spin(S):
    m = np.arange(S, -S - 1, -1)
    Sp = np.diag(np.sqrt(S*(S + 1) - m[1:]*(m[1:] + 1)), 1).astype(complex)
    return [(Sp + Sp.conj().T)/2, (Sp - Sp.conj().T)/2j, np.diag(m).astype(complex)]


S, I = spin(1.0), spin(0.5)
n = pos/np.linalg.norm(pos)
A = H_BAR*GAMMA_E*GAMMA_N/np.linalg.norm(pos)**3*(np.eye(3) - 3*np.outer(n, n))
H = sum(2*np.pi*T[i, j]*np.kron(S[i] @ S[j], np.eye(2)) + A[i, j]*np.kron(S[i], I[j]) for i in range(3) for j in range(3)) \
    - GAMMA_E*B*np.kron(S[2], np.eye(2)) - GAMMA_N*B*np.kron(np.eye(3), I[2])
w, V = np.linalg.eigh(H)
psi = np.array([1, 1, 0])/np.sqrt(2)
rho0 = V.conj().T @ np.kron(np.outer(psi, psi), np.eye(2)/2) @ V
L = []
for k in range(NSTEP):
    u = np.exp(-1j*w*k*DELTAT)
    rho = V @ (u[:, None]*rho0*u.conj()[None, :]) @ V.conj().T
    L.append(np.trace(rho.reshape(3, 2, 3, 2), axis1=1, axis2=3)[0, 1])
L = np.array(L)/L[0]
worst = 0.0
for kind in ('noDiv', 'wiDiv'):
    try:
        rows = [l.split() for l in open(f'{sys.argv[3]}_{kind}') if l.strip()]
        t = np.array([float(r[0]) for r in rows]); x = np.array([complex(r[1]) for r in rows])
    except (OSError, ValueError, IndexError) as e:
        print(f'cannot read {sys.argv[3]}_{kind}: {e}'); sys.exit(2)
    if len(x) != NSTEP or np.max(np.abs(t - DELTAT*np.arange(NSTEP))) > 1e-9 or not np.all(np.isfinite(x)):
        print(f'{sys.argv[3]}_{kind} does not hold {NSTEP} finite points k * {DELTAT}'); sys.exit(2)
    worst = max(worst, np.max(np.abs(x - L)))
print(f'max |L - L_exact| = {worst:.2e} (tolerance {TOL:g}), |L(end)| exact {abs(L[-1]):.10f}')
sys.exit(0 if worst <= TOL else 1)
