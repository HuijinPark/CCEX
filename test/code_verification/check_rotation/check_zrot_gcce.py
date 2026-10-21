"""gCCE with a non-axial zero-field splitting (E != 0) turned rigidly about the computational z axis (A91).
usage: check_zrot_gcce.py CCEX_BIN MPIRUN WORKDIR

The whole system -- the qubit positions, the ZFS tensors R T R^T, the 13C positions (point dipole, so the hyperfine
turns with them), the qubit-qubit point dipole (computed from the turned positions) -- is written out turned by phi
about z, by hand, with frame.rotation off; the field stays along z. The rotation is a symmetry of the Hamiltonian,
H' = U H U^dag with U = exp(-i phi J_z), but the run also needs the initial state and the observable, and CCEX keeps
those in the computational frame: |alpha>, |beta> are m_s states, psi0 = (|alpha> + |beta>)/norm. Turned with the
system they would be U|alpha> = exp(-i phi m_alpha)|alpha>, U|beta> = exp(-i phi m_beta)|beta>: the same states up to
a RELATIVE phase exp(-i phi (m_beta - m_alpha)) in psi0, i.e. another preparation phase. With E != 0 (or any coupling
that does not conserve J_z, the non-secular hyperfine too) the coherence depends on that phase, so with the m_s
states fixed the turned run is a different experiment; CCEX agrees with the exact answer for each of them
(check_analytic, gCCE + ZFS E). Here the turned state is given where it is real, phi (m_beta - m_alpha) a multiple of
180 degrees, as Qubit.alpha_state/beta_state vectors, and then |L| has to be the same at every point:

  fid_10_180   one qubit, +1 / 0, phi = 180 (ZFS unchanged, the bath turned), free induction, beta = -|0>
  fid_1m1_90   one qubit, +1 / -1, phi = 90 (E -> -E), free induction, beta = -|-1>
  hahn_10_180  as fid_10_180 with one pi pulse: the pulse axis x of the pair (alpha, -beta) is -x of (alpha, beta),
               a pi pulse about -x is the one about x times -1, so the echo is invariant too
  mq_fid_90    two qubits (E of opposite sign), (+1, +1) / (0, 0), phi = 90, free induction, beta = -|0 0>

The tolerance, 5e-8, is rounding: positions written to 15 digits and an eigendecomposition of a different matrix,
amplified by the GHz qubit phase (2 pi 2.87 GHz x 0.1 ms ~ 2e6 rad); measured 1e-10 to 8e-9, and turning by 360 degrees
(the same system up to the rounding of cos and sin) gives the same 6e-10 to 8e-9. That the comparison is not
empty is checked on fid_10_180 with the m_s states fixed (no turned state): that run has to differ by more than 1e-5
(measured 4.3e-4, the physical dependence on the preparation phase). For a general phi with pulses about X or Y the pulse
axis is a computational-frame direction like psi0: turning the system turns the pulse axes relative to it, so the same
kind of difference is expected and not checked. Exit 1 on a failed check, 2 when a run fails.
"""
import json, os, subprocess, sys
import numpy as np

BIN, MPIRUN, WORK = sys.argv[1], sys.argv[2], sys.argv[3]
TOL, MIN_FIXED = 5e-8, 1e-5
D = 2870000.0                                                # kHz
QUBITS = [([0.0, 0.0, 0.0], 100000.0), ([5.0, 2.0, 3.0], -60000.0)]   # position (A), E (kHz)
NUCLEI = [[2.0, 1.0, 1.5], [-1.5, 2.2, -1.0], [1.0, -2.5, 0.5], [-2.0, -1.0, 2.0]]


def run(name, phi, nq, ams, bms, turned, npulse):
    c, s = np.cos(np.radians(phi)), np.sin(np.radians(phi))
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    d = os.path.join(WORK, f'{name}_p{phi:g}{"_turned" if turned else ""}')
    os.makedirs(d, exist_ok=True)
    with open(f'{d}/bath', 'w') as f:
        f.write(f'{len(NUCLEI):12.5f}\t0.0\t0.0\t0.0\n')
        for p in NUCLEI:
            q = R @ np.array(p)
            f.write(f'{q[0]:.15f}\t{q[1]:.15f}\t{q[2]:.15f}\t13C\n')
    with open(f'{d}/gyro', 'w') as f:
        f.write('13C 0.5 6.728284\n')
    qubits, intmap = [], []
    for i in range(nq):
        pos, E = QUBITS[i]
        T = R @ np.diag([-D/3 + E, -D/3 - E, 2*D/3]) @ R.T
        qubits.append({'name': f'q{i}', 'spin': 1.0, 'gyro': -17608.59705, 'xyz': (R @ np.array(pos)).tolist(),
                       'detuning': 0.0})
        intmap.append({'between': [f'q{i}', f'q{i}'], 'tensor': T.tolist()})
    Q = {'qubit': qubits, 'intmap': intmap, 'alpha_state': {'ms': list(ams)}, 'beta_state': {'ms': list(bms)}}
    if turned:
        # U |m> = exp(-i phi m) |m>; the overall phase taken from alpha, beta keeps exp(-i phi (sum m_b - sum m_a))
        rel = np.exp(-1j * np.radians(phi) * (sum(bms) - sum(ams)))
        assert abs(rel.imag) < 1e-12, 'the turned state has to be real'
        def vec(ms, sign):
            v = np.ones(1)
            for m in ms:
                e = np.zeros(3); e[int(round(1 - m))] = 1.0; v = np.kron(v, e)
            return (sign * v).tolist()
        Q['alpha_state'] = {'vector': vec(ams, 1.0)}
        Q['beta_state'] = {'vector': vec(bms, round(rel.real))}
    inp = {'method': 'gCCE', 'quantity': 'coherence', 'Qubit': Q,
           'configuration': {'gyrofile': 'gyro', 'bathfile': ['bath']},
           'order': 1, 'bfield': 5, 'rbath': 50, 'rdip': 50, 'deltat': 0.0005, 'nstep': 200,
           'pulse': {'count': npulse}, 'output': {'file': 'out'}}
    with open(f'{d}/in.json', 'w') as f:
        json.dump(inp, f, indent=1)
    with open(f'{d}/run.log', 'w') as log:
        rc = subprocess.run([MPIRUN, '-n', '1', BIN, '-f', 'in.json'], cwd=d, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        print(f'  FAIL  {name} phi {phi:g}: CCEX exit {rc} (see {d}/run.log)'); sys.exit(2)
    out = {}
    for kind in ('noDiv', 'wiDiv'):
        with open(f'{d}/out_{kind}') as f:
            out[kind] = np.array([complex(l.split()[1]) for l in f if l.strip()])
        if len(out[kind]) != 200 or not np.all(np.isfinite(out[kind])):
            print(f'  FAIL  {name} phi {phi:g}: {d}/out_{kind} does not hold 200 finite values'); sys.exit(2)
    return out


def dev(a, b):
    return max(np.max(np.abs(np.abs(a[k]) - np.abs(b[k]))) for k in ('noDiv', 'wiDiv'))


nfail = 0
CASES = [('fid_10_180', 1, 180, [1], [0], 0), ('fid_1m1_90', 1, 90, [1], [-1], 0),
         ('hahn_10_180', 1, 180, [1], [0], 1), ('mq_fid_90', 2, 90, [1, 1], [0, 0], 0)]
for name, nq, phi, ams, bms, npulse in CASES:
    base = run(name, 0, nq, ams, bms, False, npulse)
    x = dev(run(name, phi, nq, ams, bms, True, npulse), base)
    ok = x <= TOL; nfail += not ok
    print(f'  {"ok  " if ok else "FAIL"}  {name}: turned by {phi} with the state turned too, max ||L| - |L_0|| {x:.2e} (<= {TOL:g})')
    if name == 'fid_10_180':
        y = dev(run(name, phi, nq, ams, bms, False, npulse), base)
        ok = y > MIN_FIXED; nfail += not ok
        print(f'  {"ok  " if ok else "FAIL"}  {name}: the m_s states kept (another preparation phase), max ||L| - |L_0|| '
              f'{y:.2e} (> {MIN_FIXED:g}: the comparison above is not empty)')
sys.exit(1 if nfail else 0)
