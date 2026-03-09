"""
Compare CCEX ME-CCE results against PyCCE ME-CCE.
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
RESULT_DIR = os.path.join(TEST_DIR, "comparison_S1_ms-1_0")
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

time_space = np.linspace(0, 2, 51)  # ms
magnetic_field = 300  # Gauss
order = 2
r_bath = 1400
r_dipole = 700

# CCE (no dissipation)
calc_cce = pycce.Simulator(spin=cen, bath=atoms, order=order,
                           r_bath=r_bath, r_dipole=r_dipole,
                           magnetic_field=magnetic_field)
cce_result = calc_cce.compute(time_space, method='cce', quantity='coherence',
                              nbstates=0, pulses=1)
print(f"PyCCE CCE done. shape={cce_result.shape}")

# ME-CCE (with dissipation)
atoms_mecce = atoms.copy()
atoms_mecce.add_single_jump('p', rate=1.0, units='rad')  # S+
atoms_mecce.add_single_jump('m', rate=1.0, units='rad')  # S-

calc_mecce = pycce.Simulator(spin=cen, bath=atoms_mecce, order=order,
                             r_bath=r_bath, r_dipole=r_dipole,
                             magnetic_field=magnetic_field)
mecce_result = calc_mecce.compute(time_space, method='mecce', quantity='coherence',
                                  nbstates=0, pulses=1)
print(f"PyCCE ME-CCE done. shape={mecce_result.shape}")

# Save PyCCE results
np.savetxt(os.path.join(RESULT_DIR, "pycce_cce.dat"),
           np.column_stack([time_space, cce_result.real, cce_result.imag]),
           header="time(ms) Re(L) Im(L)", fmt="%.12e")
np.savetxt(os.path.join(RESULT_DIR, "pycce_mecce.dat"),
           np.column_stack([time_space, mecce_result.real, mecce_result.imag]),
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

# CCE
cce_config = {**common_config, "method": "cce",
              "outfile": os.path.join(RESULT_DIR, "ccex_cce")}
cce_config_file = os.path.join(RESULT_DIR, "ccex_cce.json")
with open(cce_config_file, 'w') as f:
    json.dump(cce_config, f, indent=4)

print("Running CCEX CCE...")
res = subprocess.run([CCEX_BIN, "-f", cce_config_file],
                     capture_output=True, text=True, timeout=300)
print("CCE done." if res.returncode == 0 else f"CCE FAILED:\n{res.stderr[:300]}")

# ME-CCE
mecce_config = {**common_config, "method": "mecce",
                "outfile": os.path.join(RESULT_DIR, "ccex_mecce"),
                "jump_operators": [
                    {"bath_name": "e", "operator": "-+", "rate": rate_json},
                    {"bath_name": "e", "operator": "+-", "rate": rate_json}
                ]}
mecce_config_file = os.path.join(RESULT_DIR, "ccex_mecce.json")
with open(mecce_config_file, 'w') as f:
    json.dump(mecce_config, f, indent=4)

print("Running CCEX ME-CCE...")
res = subprocess.run([CCEX_BIN, "-f", mecce_config_file],
                     capture_output=True, text=True, timeout=300)
print("ME-CCE done." if res.returncode == 0 else f"ME-CCE FAILED:\n{res.stderr[:300]}")

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

ccex_cce_t, ccex_cce_v = load_ccex(os.path.join(RESULT_DIR, "ccex_cce"))
ccex_mecce_t, ccex_mecce_v = load_ccex(os.path.join(RESULT_DIR, "ccex_mecce"))

# Save CCEX results in comparable format
np.savetxt(os.path.join(RESULT_DIR, "ccex_cce.dat"),
           np.column_stack([ccex_cce_t, ccex_cce_v.real, ccex_cce_v.imag]),
           header="time(ms) Re(L) Im(L)", fmt="%.12e")
np.savetxt(os.path.join(RESULT_DIR, "ccex_mecce.dat"),
           np.column_stack([ccex_mecce_t, ccex_mecce_v.real, ccex_mecce_v.imag]),
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

    # Print selected points
    print(f"  {'Time':>8}  {'PyCCE':>14}  {'CCEX':>14}  {'|Diff|':>12}")
    for k in [0, 1, 5, 12, 25, 37, 48, 49]:
        if k < len(matched_t):
            print(f"  {matched_t[k]:8.4f}  {matched_p[k].real:+14.10f}  {matched_c[k].real:+14.10f}  {diffs[k]:.6e}")

    return matched_t, matched_p, matched_c, diffs

mt_cce, mp_cce, mc_cce, d_cce = match_and_compare("CCE", time_space, cce_result, ccex_cce_t, ccex_cce_v)
mt_mecce, mp_mecce, mc_mecce, d_mecce = match_and_compare("ME-CCE", time_space, mecce_result, ccex_mecce_t, ccex_mecce_v)

# Save diff data
np.savetxt(os.path.join(RESULT_DIR, "diff_cce.dat"),
           np.column_stack([mt_cce, d_cce]),
           header="time(ms) |diff|", fmt="%.12e")
np.savetxt(os.path.join(RESULT_DIR, "diff_mecce.dat"),
           np.column_stack([mt_mecce, d_mecce]),
           header="time(ms) |diff|", fmt="%.12e")

###############################################################################
# Step 5: Generate comparison plots
###############################################################################
print("\n" + "=" * 60)
print("Step 5: Generating plots")
print("=" * 60)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('PyCCE vs CCEX Comparison\n(S=1, $m_s$=-1/0, e-spin bath, Hahn echo, B=300G)', fontsize=13)

# --- Plot 1: CCE coherence ---
ax = axes[0, 0]
ax.plot(time_space, cce_result.real, 'b-', linewidth=2, label='PyCCE CCE')
ax.plot(ccex_cce_t, ccex_cce_v.real, 'r--', linewidth=2, label='CCEX CCE')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Re[L(t)]')
ax.set_title('CCE Coherence')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 2: ME-CCE coherence ---
ax = axes[0, 1]
ax.plot(time_space, mecce_result.real, 'b-', linewidth=2, label='PyCCE ME-CCE')
ax.plot(ccex_mecce_t, ccex_mecce_v.real, 'r--', linewidth=2, label='CCEX ME-CCE')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Re[L(t)]')
ax.set_title('ME-CCE Coherence (rate=1.0 rad/ms)')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Plot 3: CCE vs ME-CCE overlay ---
ax = axes[1, 0]
ax.plot(time_space, cce_result.real, 'b-', linewidth=2, label='PyCCE CCE')
ax.plot(time_space, mecce_result.real, 'b--', linewidth=2, label='PyCCE ME-CCE')
ax.plot(ccex_cce_t, ccex_cce_v.real, 'r-', linewidth=1.5, alpha=0.7, label='CCEX CCE')
ax.plot(ccex_mecce_t, ccex_mecce_v.real, 'r--', linewidth=1.5, alpha=0.7, label='CCEX ME-CCE')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Re[L(t)]')
ax.set_title('CCE vs ME-CCE (both codes)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# --- Plot 4: Absolute difference ---
ax = axes[1, 1]
ax.semilogy(mt_cce, d_cce, 'b-o', markersize=3, linewidth=1.5, label='CCE |diff|')
ax.semilogy(mt_mecce, d_mecce, 'r-s', markersize=3, linewidth=1.5, label='ME-CCE |diff|')
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

# --- Additional plot: zoomed difference ---
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
fig2.suptitle('Detailed Difference Analysis', fontsize=13)

ax = axes2[0]
ax.plot(mt_cce, d_cce, 'b-o', markersize=4, label='CCE')
ax.plot(mt_mecce, d_mecce, 'r-s', markersize=4, label='ME-CCE')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('|PyCCE - CCEX|')
ax.set_title('Absolute Difference (linear)')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes2[1]
# Relative difference (where |L| > 0.001)
rel_cce = np.full_like(d_cce, np.nan)
rel_mecce = np.full_like(d_mecce, np.nan)
mask_cce = np.abs(mp_cce) > 0.001
mask_mecce = np.abs(mp_mecce) > 0.001
rel_cce[mask_cce] = d_cce[mask_cce] / np.abs(mp_cce[mask_cce])
rel_mecce[mask_mecce] = d_mecce[mask_mecce] / np.abs(mp_mecce[mask_mecce])
ax.semilogy(mt_cce[mask_cce], rel_cce[mask_cce], 'b-o', markersize=4, label='CCE')
ax.semilogy(mt_mecce[mask_mecce], rel_mecce[mask_mecce], 'r-s', markersize=4, label='ME-CCE')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('|diff| / |L|')
ax.set_title('Relative Difference (where |L|>0.001)')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plot_file2 = os.path.join(RESULT_DIR, "difference_detail.png")
plt.savefig(plot_file2, dpi=150)
print(f"Plot saved: {plot_file2}")

print(f"\nAll results saved in: {RESULT_DIR}")
print("=" * 60)
print("Done.")
