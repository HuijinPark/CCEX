"""
Compare CCEX ME-gCCE results against PyCCE ME-gCCE.
S=1 central spin, alpha=ms=-1, beta=ms=0 (NV center convention).
Electron spin bath with S+ and S- jump operators.
"""
import numpy as np
import json
import subprocess
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CCEX_BIN = "/Users/istblue/Desktop/codes/CCEX/bin/main.out"
TEST_DIR = "/Users/istblue/Desktop/codes/CCEX/test/mecce_test"
RESULT_DIR = os.path.join(TEST_DIR, "comparison_megcce_S1")
os.makedirs(RESULT_DIR, exist_ok=True)

###############################################################################
# Step 1: PyCCE reference calculation
###############################################################################
print("=" * 60)
print("Step 1: Running PyCCE reference calculation")
print("=" * 60)

import pycce

atoms = pycce.random_bath('e', [5e3, 5e3, 5e3],
                          density=1e16, density_units='cm-3', seed=2)
print(f"Total bath spins generated: {len(atoms)}")

# S=1 central spin, alpha=ms=-1, beta=ms=0
cen = pycce.CenterArray(spin=1, alpha=[0, 0, 1], beta=[0, 1, 0])

time_space = np.linspace(0, 2, 21)  # ms
magnetic_field = 300  # Gauss
order = 2
r_bath = 500
r_dipole = 300

# gCCE (no dissipation)
calc_gcce = pycce.Simulator(spin=cen, bath=atoms, order=order,
                            r_bath=r_bath, r_dipole=r_dipole,
                            magnetic_field=magnetic_field)
gcce_result = calc_gcce.compute(time_space, method='gcce', quantity='coherence',
                                nbstates=0, pulses=1)
print(f"PyCCE gCCE done. shape={gcce_result.shape}")

# ME-gCCE (with dissipation)
atoms_megcce = atoms.copy()
atoms_megcce.add_single_jump('p', rate=1.0, units='rad')  # S+
atoms_megcce.add_single_jump('m', rate=1.0, units='rad')  # S-

calc_megcce = pycce.Simulator(spin=cen, bath=atoms_megcce, order=order,
                              r_bath=r_bath, r_dipole=r_dipole,
                              magnetic_field=magnetic_field)
megcce_result = calc_megcce.compute(time_space, method='megcce', quantity='coherence',
                                    nbstates=0, pulses=1)
print(f"PyCCE ME-gCCE done. shape={megcce_result.shape}")

# Save PyCCE results
np.savetxt(os.path.join(RESULT_DIR, "pycce_gcce.dat"),
           np.column_stack([time_space, gcce_result.real, gcce_result.imag]),
           header="time(ms) Re(L) Im(L)", fmt="%.12e")
np.savetxt(os.path.join(RESULT_DIR, "pycce_megcce.dat"),
           np.column_stack([time_space, megcce_result.real, megcce_result.imag]),
           header="time(ms) Re(L) Im(L)", fmt="%.12e")
print("PyCCE results saved.")

###############################################################################
# Step 2: Export bath for CCEX
###############################################################################
print("\n" + "=" * 60)
print("Step 2: Exporting bath for CCEX")
print("=" * 60)

bath_xyz = atoms['xyz']
bath_names = atoms['N']

qubit_file = os.path.join(RESULT_DIR, "qubit_espin")
with open(qubit_file, 'w') as f:
    f.write("0.0\t0.0\t0.0\n")

gyro_file = os.path.join(RESULT_DIR, "gyro_espin")
with open(gyro_file, 'w') as f:
    f.write("e  0.5  -17608.597050\n")

bath_file = os.path.join(RESULT_DIR, "bath_espin")
nspin_total = len(bath_xyz)
with open(bath_file, 'w') as f:
    f.write(f"{nspin_total:.5f}\t0.00000\t0.00000\t0.00000\n")
    for i in range(nspin_total):
        f.write(f"{bath_xyz[i,0]:.5f}\t{bath_xyz[i,1]:.5f}\t{bath_xyz[i,2]:.5f}\t{bath_names[i]}\n")
print(f"Bath: {nspin_total} spins exported")

###############################################################################
# Step 3: Create CCEX configs and run
###############################################################################
print("\n" + "=" * 60)
print("Step 3: Running CCEX calculations")
print("=" * 60)

rate_pycce = 1.0  # rad/ms
rate_json = rate_pycce / (2.0 * np.pi * 1000.0)  # MHz (CCEX input unit)
deltat = time_space[1] - time_space[0]
nstep = len(time_space) - 1  # 50

common_config = {
    "quantity": "coherence",
    "qubitfile": qubit_file,
    "gyrofile": gyro_file,
    "bathfile": [bath_file],
    "order": order,
    "bfield": magnetic_field,
    "rbath": r_bath,
    "rdip": r_dipole,
    "deltat": deltat,
    "nstep": nstep,
    "nstate": 0,
    "qspin": 1.0,
    "qalphams": -1,
    "qbetams": 0,
    "addsubclus": True,
    "nk": [0, 0, 0],
    "npulse": 1,
    "savemode": "normal",
}

# gCCE
gcce_config = {**common_config, "method": "gcce",
               "outfile": os.path.join(RESULT_DIR, "ccex_gcce")}
gcce_config_file = os.path.join(RESULT_DIR, "ccex_gcce.json")
with open(gcce_config_file, 'w') as f:
    json.dump(gcce_config, f, indent=4)

print("Running CCEX gCCE...")
res = subprocess.run([CCEX_BIN, "-f", gcce_config_file],
                     capture_output=True, text=True, timeout=600)
print("gCCE done." if res.returncode == 0 else f"gCCE FAILED:\n{res.stderr[:500]}")
if res.returncode != 0:
    print(f"STDOUT:\n{res.stdout[-500:]}")

# ME-gCCE
megcce_config = {**common_config, "method": "megcce",
                 "outfile": os.path.join(RESULT_DIR, "ccex_megcce"),
                 "jump_operators": [
                     {"bath_name": "e", "operator": "-+", "rate": rate_json},
                     {"bath_name": "e", "operator": "+-", "rate": rate_json}
                 ]}
megcce_config_file = os.path.join(RESULT_DIR, "ccex_megcce.json")
with open(megcce_config_file, 'w') as f:
    json.dump(megcce_config, f, indent=4)

print("Running CCEX ME-gCCE...")
res = subprocess.run([CCEX_BIN, "-f", megcce_config_file],
                     capture_output=True, text=True, timeout=600)
print("ME-gCCE done." if res.returncode == 0 else f"ME-gCCE FAILED:\n{res.stderr[:500]}")
if res.returncode != 0:
    print(f"STDOUT:\n{res.stdout[-500:]}")

###############################################################################
# Step 4: Load and compare
###############################################################################
print("\n" + "=" * 60)
print("Step 4: Comparing results")
print("=" * 60)

def load_ccex(prefix):
    times, vals = [], []
    with open(prefix + "_wiDiv") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                times.append(float(parts[0]))
                vals.append(complex(parts[1]))
    return np.array(times), np.array(vals)

ccex_gcce_t, ccex_gcce_v = load_ccex(os.path.join(RESULT_DIR, "ccex_gcce"))
ccex_megcce_t, ccex_megcce_v = load_ccex(os.path.join(RESULT_DIR, "ccex_megcce"))

# Save CCEX results in comparable format
np.savetxt(os.path.join(RESULT_DIR, "ccex_gcce.dat"),
           np.column_stack([ccex_gcce_t, ccex_gcce_v.real, ccex_gcce_v.imag]),
           header="time(ms) Re(L) Im(L)", fmt="%.12e")
np.savetxt(os.path.join(RESULT_DIR, "ccex_megcce.dat"),
           np.column_stack([ccex_megcce_t, ccex_megcce_v.real, ccex_megcce_v.imag]),
           header="time(ms) Re(L) Im(L)", fmt="%.12e")

def match_and_compare(name, pycce_t, pycce_v, ccex_t, ccex_v):
    matched_t, matched_p, matched_c = [], [], []
    for i, tc in enumerate(ccex_t):
        idx = np.argmin(np.abs(pycce_t - tc))
        if np.abs(pycce_t[idx] - tc) < 1e-6:
            matched_t.append(tc)
            matched_p.append(pycce_v[idx])
            matched_c.append(ccex_v[i])
    matched_t = np.array(matched_t)
    matched_p = np.array(matched_p)
    matched_c = np.array(matched_c)
    diffs = np.abs(matched_p - matched_c)

    print(f"\n--- {name} ---")
    print(f"Matched {len(matched_t)} points")
    print(f"  Max |diff|:  {np.max(diffs):.6e}")
    print(f"  Mean |diff|: {np.mean(diffs):.6e}")

    print(f"  {'Time':>8}  {'PyCCE':>14}  {'CCEX':>14}  {'|Diff|':>12}")
    for k in [0, 1, 5, 12, 25, 37, 48, 49]:
        if k < len(matched_t):
            print(f"  {matched_t[k]:8.4f}  {matched_p[k].real:+14.10f}  {matched_c[k].real:+14.10f}  {diffs[k]:.6e}")

    return matched_t, matched_p, matched_c, diffs

mt_gcce, mp_gcce, mc_gcce, d_gcce = match_and_compare("gCCE", time_space, gcce_result, ccex_gcce_t, ccex_gcce_v)
mt_megcce, mp_megcce, mc_megcce, d_megcce = match_and_compare("ME-gCCE", time_space, megcce_result, ccex_megcce_t, ccex_megcce_v)

# Save diff data
np.savetxt(os.path.join(RESULT_DIR, "diff_gcce.dat"),
           np.column_stack([mt_gcce, d_gcce]),
           header="time(ms) |diff|", fmt="%.12e")
np.savetxt(os.path.join(RESULT_DIR, "diff_megcce.dat"),
           np.column_stack([mt_megcce, d_megcce]),
           header="time(ms) |diff|", fmt="%.12e")

###############################################################################
# Step 5: Generate comparison plots
###############################################################################
print("\n" + "=" * 60)
print("Step 5: Generating plots")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('PyCCE vs CCEX: gCCE and ME-gCCE Comparison\n(S=1, $m_s$=-1/0, e-spin bath, Hahn echo, B=300G)', fontsize=13)

# --- Plot 1: gCCE coherence ---
ax = axes[0, 0]
ax.plot(time_space, gcce_result.real, 'b-', linewidth=2, label='PyCCE gCCE')
ax.plot(ccex_gcce_t, ccex_gcce_v.real, 'r--', linewidth=2, label='CCEX gCCE')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Re[L(t)]')
ax.set_title('gCCE Coherence')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 2: ME-gCCE coherence ---
ax = axes[0, 1]
ax.plot(time_space, megcce_result.real, 'b-', linewidth=2, label='PyCCE ME-gCCE')
ax.plot(ccex_megcce_t, ccex_megcce_v.real, 'r--', linewidth=2, label='CCEX ME-gCCE')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Re[L(t)]')
ax.set_title('ME-gCCE Coherence (rate=1.0 rad/ms)')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 3: gCCE vs ME-gCCE overlay ---
ax = axes[1, 0]
ax.plot(time_space, gcce_result.real, 'b-', linewidth=2, label='PyCCE gCCE')
ax.plot(time_space, megcce_result.real, 'b--', linewidth=2, label='PyCCE ME-gCCE')
ax.plot(ccex_gcce_t, ccex_gcce_v.real, 'r-', linewidth=1.5, alpha=0.7, label='CCEX gCCE')
ax.plot(ccex_megcce_t, ccex_megcce_v.real, 'r--', linewidth=1.5, alpha=0.7, label='CCEX ME-gCCE')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Re[L(t)]')
ax.set_title('gCCE vs ME-gCCE (both codes)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# --- Plot 4: Absolute difference ---
ax = axes[1, 1]
ax.semilogy(mt_gcce, d_gcce, 'b-o', markersize=3, linewidth=1.5, label='gCCE |diff|')
ax.semilogy(mt_megcce, d_megcce, 'r-s', markersize=3, linewidth=1.5, label='ME-gCCE |diff|')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('|PyCCE - CCEX|')
ax.set_title('Absolute Difference')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=1e-14)

plt.tight_layout()
plot_file = os.path.join(RESULT_DIR, "comparison_plot.png")
plt.savefig(plot_file, dpi=150)
print(f"Plot saved: {plot_file}")

print(f"\nAll results saved in: {RESULT_DIR}")
print("=" * 60)
print("Done.")
