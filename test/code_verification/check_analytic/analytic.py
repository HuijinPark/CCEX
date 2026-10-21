"""The closed-form coherence of one nuclear spin 1/2 beside the NV, against CCEX's output.
usage: analytic.py OUT_PREFIX --npulse 0|1 --nstep N --deltat MS [--tol 2e-10]
With one bath spin, order 1 is exact. For NV branch m (alpha_state ms [1], beta_state ms [0]) the spin precesses about
w_m = m A_z - gamma_n B z, A the point-dipole tensor H_BAR gamma_e gamma_n / r^3 (1 - 3 r r) written out here
independently of CCEX; n_m = w_m / |w_m|. Ensemble average over the nuclear state, t the total time:
  free induction  L(t) = cos(|w0| t/2) cos(|w1| t/2) + (n0 . n1) sin(|w0| t/2) sin(|w1| t/2)
  Hahn echo       L(t) = 1 - 2 |n0 x n1|^2 sin^2(|w0| t/4) sin^2(|w1| t/4)
The inputs (inputs/: qubit at (10,10,10), one 13C at (13,11.5,14) A, B = 500 G along z) are restated here; the
tolerance is the 10 decimals CCEX writes. The output has to hold the requested grid, t_k = k deltat for k < nstep,
every value finite. Exit 1 when the largest deviation is above the tolerance, 2 when the output is missing, cannot be
read or does not hold the grid (so a self-test can tell a caught error from a broken run)."""
import sys
import numpy as np

H_BAR, GAMMA_E, GAMMA_N, B = 1.054571729, -17608.59705, 6.728284, 500.0      # rad/ms/G, G
QUBIT, SPIN = np.array([10.0, 10.0, 10.0]), np.array([13.0, 11.5, 14.0])     # A


def opt(name, cast, default):
    return cast(sys.argv[sys.argv.index(name)+1]) if name in sys.argv else default


prefix, npulse, tol = sys.argv[1], opt('--npulse', int, None), opt('--tol', float, 2e-10)
nstep, deltat = opt('--nstep', int, None), opt('--deltat', float, None)


def broken(why):
    print(f'npulse {npulse}: {why}')
    sys.exit(2)


if npulse not in (0, 1) or nstep is None or deltat is None:
    broken('give --npulse 0|1, --nstep and --deltat')
r = SPIN-QUBIT; d = np.linalg.norm(r); u = r/d
A = H_BAR*GAMMA_E*GAMMA_N/d**3*(np.eye(3)-3*np.outer(u, u))
w0, w1 = (m*A[2]+np.array([0.0, 0.0, -GAMMA_N*B]) for m in (0, 1))
a0, a1 = np.linalg.norm(w0), np.linalg.norm(w1); n0, n1 = w0/a0, w1/a1
try:
    rows = [l.split() for l in open(prefix+'_noDiv') if l.strip()]
    t = np.array([float(x[0]) for x in rows]); L = np.array([complex(x[1]) for x in rows])
except (OSError, ValueError, IndexError) as e:
    broken(f'cannot read {prefix}_noDiv: {e}')
grid = np.arange(nstep)*deltat
if len(t) != nstep or not np.all(np.isfinite(t)) or np.max(np.abs(t-grid)) > 1e-12 or not np.all(np.isfinite(L)):
    broken(f'{len(t)} rows, expected the {nstep} points k*{deltat} ms, all finite')
if npulse == 0:
    exact = np.cos(a0*t/2)*np.cos(a1*t/2)+np.dot(n0, n1)*np.sin(a0*t/2)*np.sin(a1*t/2)
else:
    exact = 1-2*np.linalg.norm(np.cross(n0, n1))**2*np.sin(a0*t/4)**2*np.sin(a1*t/4)**2
dev = float(np.max(np.abs(L-exact)))
print(f'npulse {npulse}: {len(t)} points to t = {t[-1]:g} ms, max |L - closed form| = {dev:.3e} (tol {tol:g})')
sys.exit(0 if dev <= tol else 1)
