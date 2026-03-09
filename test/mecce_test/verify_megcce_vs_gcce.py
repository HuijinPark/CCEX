"""
Verify ME-gCCE implementation by comparing with gCCE when rate=0.
If the Lindblad implementation is correct, ME-gCCE with zero dissipation
should give the same results as gCCE.
"""
import numpy as np
import json
import subprocess
import os

CCEX_BIN = "/Users/istblue/Desktop/codes/CCEX/bin/main.out"
TEST_DIR = "/Users/istblue/Desktop/codes/CCEX/test/mecce_test"
RESULT_DIR = os.path.join(TEST_DIR, "verify_megcce_rate0")
os.makedirs(RESULT_DIR, exist_ok=True)

# Use the bath files from the S=0.5 comparison
S05_DIR = os.path.join(TEST_DIR, "comparison_megcce_S05")

deltat = 0.1
nstep = 20

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
}

# Test 1: gCCE
gcce_config = {**common_config, "method": "gcce",
               "outfile": os.path.join(RESULT_DIR, "ccex_gcce")}
config_file = os.path.join(RESULT_DIR, "ccex_gcce.json")
with open(config_file, 'w') as f:
    json.dump(gcce_config, f, indent=4)
print("Running gCCE...")
res = subprocess.run([CCEX_BIN, "-f", config_file], capture_output=True, text=True, timeout=300)
print("gCCE:", "OK" if res.returncode == 0 else f"FAIL: {res.stderr[:200]}")

# Test 2: ME-gCCE with rate=0
megcce_config = {**common_config, "method": "megcce",
                 "outfile": os.path.join(RESULT_DIR, "ccex_megcce_rate0"),
                 "jump_operators": [
                     {"bath_name": "e", "operator": "-+", "rate": 0.0},
                     {"bath_name": "e", "operator": "+-", "rate": 0.0}
                 ]}
config_file = os.path.join(RESULT_DIR, "ccex_megcce_rate0.json")
with open(config_file, 'w') as f:
    json.dump(megcce_config, f, indent=4)
print("Running ME-gCCE (rate=0)...")
res = subprocess.run([CCEX_BIN, "-f", config_file], capture_output=True, text=True, timeout=600)
print("ME-gCCE rate=0:", "OK" if res.returncode == 0 else f"FAIL: {res.stderr[:200]}")

# Test 3: ME-gCCE with rate=0.001 (very small)
megcce_config2 = {**common_config, "method": "megcce",
                  "outfile": os.path.join(RESULT_DIR, "ccex_megcce_small_rate"),
                  "jump_operators": [
                      {"bath_name": "e", "operator": "-+", "rate": 1e-6},
                      {"bath_name": "e", "operator": "+-", "rate": 1e-6}
                  ]}
config_file = os.path.join(RESULT_DIR, "ccex_megcce_small_rate.json")
with open(config_file, 'w') as f:
    json.dump(megcce_config2, f, indent=4)
print("Running ME-gCCE (small rate)...")
res = subprocess.run([CCEX_BIN, "-f", config_file], capture_output=True, text=True, timeout=600)
print("ME-gCCE small rate:", "OK" if res.returncode == 0 else f"FAIL: {res.stderr[:200]}")

# Load and compare
def load_ccex(prefix):
    times, vals = [], []
    with open(prefix + "_wiDiv") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                times.append(float(parts[0]))
                vals.append(complex(parts[1]))
    return np.array(times), np.array(vals)

t_gcce, v_gcce = load_ccex(os.path.join(RESULT_DIR, "ccex_gcce"))
t_megcce0, v_megcce0 = load_ccex(os.path.join(RESULT_DIR, "ccex_megcce_rate0"))
t_megcce_s, v_megcce_s = load_ccex(os.path.join(RESULT_DIR, "ccex_megcce_small_rate"))

print("\n" + "=" * 60)
print("VERIFICATION: ME-gCCE (rate=0) vs gCCE")
print("=" * 60)
diffs = np.abs(v_gcce - v_megcce0)
print(f"Max |diff|:  {np.max(diffs):.6e}")
print(f"Mean |diff|: {np.mean(diffs):.6e}")
print(f"\n{'Time':>8}  {'gCCE':>18}  {'ME-gCCE(rate=0)':>18}  {'|Diff|':>12}")
for i in range(len(t_gcce)):
    print(f"{t_gcce[i]:8.2f}  {v_gcce[i].real:+18.12f}  {v_megcce0[i].real:+18.12f}  {diffs[i]:.6e}")

print("\n" + "=" * 60)
print("VERIFICATION: ME-gCCE (small rate) vs gCCE")
print("=" * 60)
diffs2 = np.abs(v_gcce - v_megcce_s)
print(f"Max |diff|:  {np.max(diffs2):.6e}")
print(f"Mean |diff|: {np.mean(diffs2):.6e}")
