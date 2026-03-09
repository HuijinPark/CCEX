"""
Verify pmecce implementation:
1. pmecce (rate=0) should match pcce
2. pmecce (with rate) should show additional decoherence vs pcce
"""
import numpy as np
import json
import subprocess
import os

CCEX_BIN = "/Users/istblue/Desktop/codes/CCEX/bin/main.out"
TEST_DIR = "/Users/istblue/Desktop/codes/CCEX/test/mecce_test"
RESULT_DIR = os.path.join(TEST_DIR, "verify_pmecce")
os.makedirs(RESULT_DIR, exist_ok=True)

# Use the bath files from the S=0.5 comparison
S05_DIR = os.path.join(TEST_DIR, "comparison_megcce_S05")

deltat = 0.04
nstep = 50

common_config = {
    "quantity": "coherence",
    "qubitfile": os.path.join(S05_DIR, "qubit_espin"),
    "gyrofile": os.path.join(S05_DIR, "gyro_espin"),
    "bathfile": [os.path.join(S05_DIR, "bath_espin")],
    "order": 2,
    "bfield": 300,
    "rbath": 500,
    "rdip": 300,
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
    "sK": 2,
    "max_trial": 100,
    "max_iter": 100,
}

rate_json = 1.0 / (2.0 * np.pi * 1000.0)  # 1 rad/ms → MHz

methods = {
    "pcce": {"method": "pcce"},
    "pmecce_rate0": {"method": "pmecce",
                     "jump_operators": [
                         {"bath_name": "e", "operator": "-+", "rate": 0.0},
                         {"bath_name": "e", "operator": "+-", "rate": 0.0}
                     ]},
    "pmecce": {"method": "pmecce",
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
    # Remove pCCE-specific keys for non-pCCE methods
    method_name = mconf.get("method", "")
    needs_pcce_params = method_name in ("pcce", "pmecce")
    if not needs_pcce_params:
        config.pop("sK", None)
        config.pop("max_trial", None)
        config.pop("max_iter", None)
    config_file = os.path.join(RESULT_DIR, f"ccex_{mname}.json")
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=4)

    print(f"Running CCEX {mname}...")
    res = subprocess.run(["mpirun", "-np", "1", CCEX_BIN, "-f", config_file],
                         capture_output=True, text=True, timeout=600)
    if res.returncode == 0:
        print(f"  {mname} done.")
    else:
        print(f"  {mname} FAILED:\n{res.stderr[:300]}")
        print(f"  STDOUT (last 300):\n{res.stdout[-300:]}")

# Load results
def load_ccex(prefix):
    times, vals = [], []
    with open(prefix + "_wiDiv") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                times.append(float(parts[0]))
                vals.append(complex(parts[1]))
    return np.array(times), np.array(vals)

data = {}
for mname in methods:
    data[mname] = load_ccex(os.path.join(RESULT_DIR, f"ccex_{mname}"))

# Verify: pmecce(rate=0) vs pcce
print("\n" + "=" * 60)
print("TEST 1: pmecce (rate=0) vs pcce — should match")
print("=" * 60)
t1, v1 = data["pcce"]
t2, v2 = data["pmecce_rate0"]
diffs = np.abs(v1 - v2)
print(f"Max |diff|:  {np.max(diffs):.6e}")
print(f"Mean |diff|: {np.mean(diffs):.6e}")
for i in [0, 5, 10, 25, 40, 49]:
    if i < len(t1):
        print(f"  t={t1[i]:6.3f}  pcce={v1[i].real:+14.10f}  pmecce(r=0)={v2[i].real:+14.10f}  |d|={diffs[i]:.3e}")

# Compare all methods
print("\n" + "=" * 60)
print("TEST 2: All methods comparison")
print("=" * 60)
for mname in methods:
    t, v = data[mname]
    print(f"\n{mname:>15s}:  t=0.0 → {v[0].real:.6f}  t=1.0 → {v[25].real:.6f}  t=2.0 → {v[-1].real:.6f}")

# Plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('pmecce Verification (S=0.5, e-spin bath, B=300G, sK=2)', fontsize=13)

ax = axes[0]
for mname in ["pcce", "pmecce_rate0", "pmecce"]:
    t, v = data[mname]
    label = {"pcce": "pCCE", "pmecce_rate0": "pmeCCE (rate=0)", "pmecce": "pmeCCE (rate=1rad/ms)"}[mname]
    ax.plot(t, v.real, linewidth=2, label=label)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Re[L(t)]')
ax.set_title('pCCE vs pmeCCE')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
for mname in ["cce", "mecce", "pcce", "pmecce"]:
    t, v = data[mname]
    label = {"cce": "CCE", "mecce": "ME-CCE", "pcce": "pCCE", "pmecce": "pmeCCE"}[mname]
    ax.plot(t, v.real, linewidth=2, label=label)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('Re[L(t)]')
ax.set_title('All Methods Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plot_file = os.path.join(RESULT_DIR, "verification_plot.png")
plt.savefig(plot_file, dpi=150)
print(f"\nPlot saved: {plot_file}")
print("Done.")
