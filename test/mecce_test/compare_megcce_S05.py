"""
Compare CCEX ME-gCCE results against PyCCE ME-gCCE.
S=0.5 central spin, alpha=ms=+0.5, beta=ms=-0.5.
Electron spin bath with S+ and S- jump operators.
Small bath for feasibility.
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
RESULT_DIR = os.path.join(TEST_DIR, "comparison_megcce_S05")
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

# S=0.5 central spin, alpha=ms=+0.5, beta=ms=-0.5
cen = pycce.CenterArray(spin=0.5, alpha=[1, 0], beta=[0, 1])

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

# Also compute CCE and ME-CCE for cross-check
cce_result = calc_gcce.compute(time_space, method='cce', quantity='coherence',
                               nbstates=0, pulses=1)
print(f"PyCCE CCE done. shape={cce_result.shape}")

atoms_mecce = atoms.copy()
atoms_mecce.add_single_jump('p', rate=1.0, units='rad')
atoms_mecce.add_single_jump('m', rate=1.0, units='rad')
calc_mecce = pycce.Simulator(spin=cen, bath=atoms_mecce, order=order,
                             r_bath=r_bath, r_dipole=r_dipole,
                             magnetic_field=magnetic_field)
mecce_result = calc_mecce.compute(time_space, method='mecce', quantity='coherence',
                                  nbstates=0, pulses=1)
print(f"PyCCE ME-CCE done. shape={mecce_result.shape}")

# Save PyCCE results
for name, data in [("pycce_gcce", gcce_result), ("pycce_megcce", megcce_result),
                   ("pycce_cce", cce_result), ("pycce_mecce", mecce_result)]:
    np.savetxt(os.path.join(RESULT_DIR, f"{name}.dat"),
               np.column_stack([time_space, data.real, data.imag]),
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
nstep = len(time_space) - 1  # 20

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
    "qspin": 0.5,
    "qalphams": 0.5,
    "qbetams": -0.5,
    "addsubclus": True,
    "nk": [0, 0, 0],
    "npulse": 1,
    "savemode": "normal",
}

methods = {
    "gcce": {"method": "gcce"},
    "megcce": {"method": "megcce",
               "jump_operators": [
                   {"bath_name": "e", "operator": "-+", "rate": rate_json},
                   {"bath_name": "e", "operator": "+-", "rate": rate_json}
               ]},
    "cce": {"method": "cce"},
    "mecce": {"method": "mecce",
              "jump_operators": [
                  {"bath_name": "e", "operator": "-+", "rate": rate_json},
                  {"bath_name": "e", "operator": "+-", "rate": rate_json}
              ]},
}

for mname, mconf in methods.items():
    config = {**common_config, **mconf,
              "outfile": os.path.join(RESULT_DIR, f"ccex_{mname}")}
    config_file = os.path.join(RESULT_DIR, f"ccex_{mname}.json")
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=4)

    print(f"Running CCEX {mname}...")
    res = subprocess.run([CCEX_BIN, "-f", config_file],
                         capture_output=True, text=True, timeout=600)
    if res.returncode == 0:
        print(f"  {mname} done.")
    else:
        print(f"  {mname} FAILED:\n{res.stderr[:300]}")
        print(f"  STDOUT (last 300):\n{res.stdout[-300:]}")

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

pycce_data = {}
ccex_data = {}
for mname in methods:
    ccex_data[mname] = load_ccex(os.path.join(RESULT_DIR, f"ccex_{mname}"))

pycce_results = {"gcce": gcce_result, "megcce": megcce_result,
                 "cce": cce_result, "mecce": mecce_result}

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
    if len(diffs) > 0:
        print(f"  Max |diff|:  {np.max(diffs):.6e}")
        print(f"  Mean |diff|: {np.mean(diffs):.6e}")
        print(f"  {'Time':>8}  {'PyCCE':>14}  {'CCEX':>14}  {'|Diff|':>12}")
        for k in [0, 1, 5, 10, 15, 19]:
            if k < len(matched_t):
                print(f"  {matched_t[k]:8.4f}  {matched_p[k].real:+14.10f}  {matched_c[k].real:+14.10f}  {diffs[k]:.6e}")

    return matched_t, matched_p, matched_c, diffs

results = {}
for mname in methods:
    ct, cv = ccex_data[mname]
    pv = pycce_results[mname]
    mt, mp, mc, d = match_and_compare(mname.upper(), time_space, pv, ct, cv)
    results[mname] = (mt, mp, mc, d)

###############################################################################
# Step 5: Generate comparison plots
###############################################################################
print("\n" + "=" * 60)
print("Step 5: Generating plots")
print("=" * 60)

fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('PyCCE vs CCEX: All Methods (S=0.5, e-spin bath, Hahn echo, B=300G, r_bath=500)', fontsize=13)

methods_plot = [("cce", "CCE"), ("gcce", "gCCE"), ("mecce", "ME-CCE"), ("megcce", "ME-gCCE")]

for idx, (mname, mlabel) in enumerate(methods_plot):
    row = idx // 3
    col = idx % 3
    ax = axes[row, col]

    pv = pycce_results[mname]
    ct, cv = ccex_data[mname]

    ax.plot(time_space, pv.real, 'b-', linewidth=2, label=f'PyCCE {mlabel}')
    ax.plot(ct, cv.real, 'r--', linewidth=2, label=f'CCEX {mlabel}')
    ax.set_xlabel('Time (ms)')
    ax.set_ylabel('Re[L(t)]')
    ax.set_title(f'{mlabel} Coherence')
    ax.legend()
    ax.grid(True, alpha=0.3)

# Plot 5: All methods overlay
ax = axes[1, 1]
for mname, mlabel in methods_plot:
    pv = pycce_results[mname]
    ct, cv = ccex_data[mname]
    ax.plot(time_space, pv.real, linewidth=2, label=f'PyCCE {mlabel}')
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Re[L(t)]')
ax.set_title('PyCCE: All Methods')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 6: Absolute differences
ax = axes[1, 2]
for mname, mlabel in methods_plot:
    mt, mp, mc, d = results[mname]
    if len(d) > 0:
        ax.semilogy(mt, d, '-o', markersize=3, linewidth=1.5, label=f'{mlabel} |diff|')
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
