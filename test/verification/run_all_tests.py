"""
CCEX Verification Test Suite
==============================
Tests all newly implemented methods: ME-CCE, ME-gCCE, pmeCCE
Each test runs CCEX and saves:
  - Numerical results (.dat)
  - Comparison plots (.png)
  - Summary text (summary.txt)

Usage:
    python run_all_tests.py            # CCEX-only tests (no PyCCE needed)
    python run_all_tests.py --pycce    # Also compare against PyCCE
"""

import numpy as np
import json
import subprocess
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────
CCEX_BIN  = "/Users/istblue/Desktop/codes/CCEX/bin/main.out"
BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
BATH_SRC  = "/Users/istblue/Desktop/codes/CCEX/test/mecce_test/comparison_megcce_S05"

USE_PYCCE = "--pycce" in sys.argv

results_summary = []   # (test_name, status, max_diff, note)

# ─────────────────────────────────────────────
def load_ccex(prefix):
    times, vals = [], []
    with open(prefix + "_wiDiv") as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                times.append(float(p[0]))
                vals.append(complex(p[1]))
    return np.array(times), np.array(vals)


def run_ccex(config, config_file, timeout=600):
    with open(config_file, 'w') as f:
        json.dump(config, f, indent=4)
    res = subprocess.run(["mpirun", "-np", "1", CCEX_BIN, "-f", config_file],
                         capture_output=True, text=True, timeout=timeout)
    return res.returncode == 0, res.stderr


def save_dat(path, t, v, header="time(ms) Re(L) Im(L)"):
    np.savetxt(path, np.column_stack([t, v.real, v.imag]),
               header=header, fmt="%.12e")


def bath_config(bath_dir):
    return {
        "qubitfile": os.path.join(bath_dir, "qubit_espin"),
        "gyrofile":  os.path.join(bath_dir, "gyro_espin"),
        "bathfile":  [os.path.join(bath_dir, "bath_espin")],
    }


print("=" * 60)
print("  CCEX Verification Test Suite")
print("=" * 60)
print(f"PyCCE mode: {USE_PYCCE}")
print()

# ╔══════════════════════════════════════════════════════════════╗
# ║  TEST 1 : ME-CCE  —  rate=0 must match CCE                  ║
# ╚══════════════════════════════════════════════════════════════╝
print("─" * 60)
print("TEST 1: ME-CCE verification  (rate=0 → CCE, rate>0 → decoherence)")
print("─" * 60)

T1_DIR = os.path.join(BASE_DIR, "1_verify_mecce")
os.makedirs(T1_DIR, exist_ok=True)

rate_json = 1.0 / (2.0 * np.pi * 1000.0)   # 1 rad/ms in MHz

base1 = {
    **bath_config(BATH_SRC),
    "quantity": "coherence",
    "order": 2, "bfield": 300,
    "rbath": 500, "rdip": 300,
    "deltat": 0.04, "nstep": 50,
    "nstate": 0,
    "qspin": 0.5, "qalphams": 0.5, "qbetams": -0.5,
    "addsubclus": True, "nk": [0, 0, 0],
    "npulse": 1, "savemode": "normal",
}

runs1 = {
    "cce"       : {"method": "cce"},
    "mecce_r0"  : {"method": "mecce",
                   "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": 0.0},
                                      {"bath_name": "e", "operator": "+-", "rate": 0.0}]},
    "mecce"     : {"method": "mecce",
                   "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": rate_json},
                                      {"bath_name": "e", "operator": "+-", "rate": rate_json}]},
}

data1 = {}
for mname, extra in runs1.items():
    cfg = {**base1, **extra, "outfile": os.path.join(T1_DIR, f"ccex_{mname}")}
    ok, err = run_ccex(cfg, os.path.join(T1_DIR, f"ccex_{mname}.json"))
    if ok:
        t, v = load_ccex(os.path.join(T1_DIR, f"ccex_{mname}"))
        data1[mname] = (t, v)
        save_dat(os.path.join(T1_DIR, f"{mname}.dat"), t, v)
        print(f"  {mname:12s}: OK")
    else:
        print(f"  {mname:12s}: FAILED  {err[:100]}")

# numerical comparison
t_cce, v_cce   = data1["cce"]
t_r0,  v_r0    = data1["mecce_r0"]
diffs1 = np.abs(v_cce - v_r0)
max_d1 = float(np.max(diffs1))
print(f"\n  mecce(rate=0) vs cce   max|diff| = {max_d1:.3e}")
results_summary.append(("ME-CCE rate=0 == CCE", "PASS" if max_d1 < 1e-6 else "FAIL", max_d1, ""))

# save numerical summary
with open(os.path.join(T1_DIR, "summary.txt"), "w") as f:
    f.write("TEST 1: ME-CCE Verification\n")
    f.write("=" * 50 + "\n")
    f.write(f"mecce(rate=0) vs cce   max|diff| = {max_d1:.3e}\n\n")
    f.write(f"{'Time':>8}  {'CCE':>16}  {'MECCE(r=0)':>16}  {'|diff|':>12}\n")
    for i in range(len(t_cce)):
        f.write(f"{t_cce[i]:8.3f}  {v_cce[i].real:+16.12f}  {v_r0[i].real:+16.12f}  {diffs1[i]:.6e}\n")

# plot
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle('ME-CCE Verification  (S=0.5, e-spin bath, B=300 G, Hahn echo)', fontsize=13)

ax = axes[0]
ax.plot(*data1["cce"],      lw=2,   label="CCE",                   color="steelblue")
ax.plot(*data1["mecce_r0"], lw=2,   label="ME-CCE (rate=0)",       color="tomato", ls="--")
ax.plot(*data1["mecce"],    lw=1.5, label="ME-CCE (rate=1 rad/ms)",color="seagreen",ls="-.")
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("Coherence"); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.semilogy(t_cce, diffs1, "r-o", ms=4, lw=1.5, label="ME-CCE(rate=0) vs CCE")
ax.set_xlabel("Time (ms)"); ax.set_ylabel("|CCE − ME-CCE(rate=0)|")
ax.set_title("Agreement Check"); ax.legend(); ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=1e-14)

plt.tight_layout()
plt.savefig(os.path.join(T1_DIR, "comparison_plot.png"), dpi=150)
plt.close()
print(f"  Plot  → {T1_DIR}/comparison_plot.png")


# ╔══════════════════════════════════════════════════════════════╗
# ║  TEST 2 : ME-gCCE  —  rate=0 must match gCCE               ║
# ╚══════════════════════════════════════════════════════════════╝
print("\n" + "─" * 60)
print("TEST 2: ME-gCCE verification  (rate=0 → gCCE, rate>0 → decoherence)")
print("─" * 60)

T2_DIR = os.path.join(BASE_DIR, "2_verify_megcce")
os.makedirs(T2_DIR, exist_ok=True)

base2 = {
    **bath_config(BATH_SRC),
    "quantity": "coherence",
    "order": 2, "bfield": 300,
    "rbath": 500, "rdip": 300,
    "deltat": 0.1, "nstep": 20,
    "nstate": 0,
    "qspin": 0.5, "qalphams": 0.5, "qbetams": -0.5,
    "addsubclus": True, "nk": [0, 0, 0],
    "npulse": 1, "savemode": "normal",
}

runs2 = {
    "gcce"       : {"method": "gcce"},
    "megcce_r0"  : {"method": "megcce",
                    "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": 0.0},
                                       {"bath_name": "e", "operator": "+-", "rate": 0.0}]},
    "megcce"     : {"method": "megcce",
                    "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": rate_json},
                                       {"bath_name": "e", "operator": "+-", "rate": rate_json}]},
}

data2 = {}
for mname, extra in runs2.items():
    cfg = {**base2, **extra, "outfile": os.path.join(T2_DIR, f"ccex_{mname}")}
    ok, err = run_ccex(cfg, os.path.join(T2_DIR, f"ccex_{mname}.json"))
    if ok:
        t, v = load_ccex(os.path.join(T2_DIR, f"ccex_{mname}"))
        data2[mname] = (t, v)
        save_dat(os.path.join(T2_DIR, f"{mname}.dat"), t, v)
        print(f"  {mname:12s}: OK")
    else:
        print(f"  {mname:12s}: FAILED  {err[:100]}")

t_gcce, v_gcce = data2["gcce"]
t_r0,   v_r0   = data2["megcce_r0"]
diffs2 = np.abs(v_gcce - v_r0)
max_d2 = float(np.max(diffs2))
print(f"\n  megcce(rate=0) vs gcce  max|diff| = {max_d2:.3e}")
results_summary.append(("ME-gCCE rate=0 == gCCE", "PASS" if max_d2 < 1e-6 else "FAIL", max_d2, ""))

with open(os.path.join(T2_DIR, "summary.txt"), "w") as f:
    f.write("TEST 2: ME-gCCE Verification\n")
    f.write("=" * 50 + "\n")
    f.write(f"megcce(rate=0) vs gcce   max|diff| = {max_d2:.3e}\n\n")
    f.write(f"{'Time':>8}  {'gCCE':>16}  {'MEGCCE(r=0)':>16}  {'|diff|':>12}\n")
    for i in range(len(t_gcce)):
        f.write(f"{t_gcce[i]:8.3f}  {v_gcce[i].real:+16.12f}  {v_r0[i].real:+16.12f}  {diffs2[i]:.6e}\n")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle('ME-gCCE Verification  (S=0.5, e-spin bath, B=300 G, Hahn echo)', fontsize=13)

ax = axes[0]
ax.plot(*data2["gcce"],      lw=2,   label="gCCE",                  color="steelblue")
ax.plot(*data2["megcce_r0"], lw=2,   label="ME-gCCE (rate=0)",      color="tomato",  ls="--")
ax.plot(*data2["megcce"],    lw=1.5, label="ME-gCCE (rate=1rad/ms)",color="seagreen",ls="-.")
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("Coherence"); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1]
ax.semilogy(t_gcce, diffs2, "r-o", ms=4, lw=1.5, label="ME-gCCE(rate=0) vs gCCE")
ax.set_xlabel("Time (ms)"); ax.set_ylabel("|gCCE − ME-gCCE(rate=0)|")
ax.set_title("Agreement Check"); ax.legend(); ax.grid(True, alpha=0.3)
ax.set_ylim(bottom=1e-14)

plt.tight_layout()
plt.savefig(os.path.join(T2_DIR, "comparison_plot.png"), dpi=150)
plt.close()
print(f"  Plot  → {T2_DIR}/comparison_plot.png")


# ╔══════════════════════════════════════════════════════════════╗
# ║  TEST 3 : pmeCCE  —  rate=0 must match pCCE                ║
# ╚══════════════════════════════════════════════════════════════╝
print("\n" + "─" * 60)
print("TEST 3: pmeCCE verification  (rate=0 → pCCE, rate>0 → decoherence)")
print("─" * 60)

T3_DIR = os.path.join(BASE_DIR, "3_verify_pmecce")
os.makedirs(T3_DIR, exist_ok=True)

base3 = {
    **bath_config(BATH_SRC),
    "quantity": "coherence",
    "order": 2, "bfield": 300,
    "rbath": 500, "rdip": 300,
    "deltat": 0.04, "nstep": 50,
    "nstate": 0,
    "qspin": 0.5, "qalphams": 0.5, "qbetams": -0.5,
    "addsubclus": True, "nk": [0, 0, 0],
    "npulse": 1, "savemode": "normal",
    "sK": 2, "max_trial": 5000, "max_iter": 5000,
}

base3_no_pcce = {k: v for k, v in base3.items()
                 if k not in ("sK", "max_trial", "max_iter")}

runs3 = {
    "cce"       : (base3_no_pcce, {"method": "cce"}),
    "mecce"     : (base3_no_pcce, {"method": "mecce",
                   "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": rate_json},
                                      {"bath_name": "e", "operator": "+-", "rate": rate_json}]}),
    "pcce"      : (base3,         {"method": "pcce"}),
    "pmecce_r0" : (base3,         {"method": "pmecce",
                   "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": 0.0},
                                      {"bath_name": "e", "operator": "+-", "rate": 0.0}]}),
    "pmecce"    : (base3,         {"method": "pmecce",
                   "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": rate_json},
                                      {"bath_name": "e", "operator": "+-", "rate": rate_json}]}),
}

data3 = {}
for mname, (bconf, extra) in runs3.items():
    cfg = {**bconf, **extra, "outfile": os.path.join(T3_DIR, f"ccex_{mname}")}
    ok, err = run_ccex(cfg, os.path.join(T3_DIR, f"ccex_{mname}.json"))
    if ok:
        t, v = load_ccex(os.path.join(T3_DIR, f"ccex_{mname}"))
        data3[mname] = (t, v)
        save_dat(os.path.join(T3_DIR, f"{mname}.dat"), t, v)
        print(f"  {mname:12s}: OK")
    else:
        print(f"  {mname:12s}: FAILED  {err[:100]}")

t_pcce,  v_pcce  = data3["pcce"]
t_pr0,   v_pr0   = data3["pmecce_r0"]
diffs3 = np.abs(v_pcce - v_pr0)
max_d3 = float(np.max(diffs3))
print(f"\n  pmecce(rate=0) vs pcce  max|diff| = {max_d3:.3e}")
results_summary.append(("pmeCCE rate=0 == pCCE", "PASS" if max_d3 < 1e-6 else "FAIL", max_d3, ""))

with open(os.path.join(T3_DIR, "summary.txt"), "w") as f:
    f.write("TEST 3: pmeCCE Verification\n")
    f.write("=" * 50 + "\n")
    f.write(f"pmecce(rate=0) vs pcce   max|diff| = {max_d3:.3e}\n\n")
    f.write(f"{'Time':>8}  {'pCCE':>16}  {'pmeCCE(r=0)':>16}  {'|diff|':>12}\n")
    for i in range(len(t_pcce)):
        f.write(f"{t_pcce[i]:8.3f}  {v_pcce[i].real:+16.12f}  {v_pr0[i].real:+16.12f}  {diffs3[i]:.6e}\n")
    f.write("\nAll Methods at selected time points:\n")
    f.write(f"{'Method':>14}  {'t=0.0':>12}  {'t=1.0':>12}  {'t=2.0':>12}\n")
    for mn in runs3:
        if mn in data3:
            t, v = data3[mn]
            i25 = 25 if 25 < len(t) else -1
            f.write(f"{mn:>14}  {v[0].real:12.8f}  {v[i25].real:12.8f}  {v[-1].real:12.8f}\n")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('pmeCCE Verification  (S=0.5, e-spin bath, B=300 G, sK=2, Hahn echo)', fontsize=13)

ax = axes[0]
colors = {"pcce": "steelblue", "pmecce_r0": "tomato", "pmecce": "seagreen"}
labels = {"pcce": "pCCE", "pmecce_r0": "pmeCCE (rate=0)", "pmecce": "pmeCCE (rate=1 rad/ms)"}
for mn in ["pcce", "pmecce_r0", "pmecce"]:
    if mn in data3:
        ax.plot(*data3[mn], lw=2, label=labels[mn], color=colors[mn],
                ls="--" if "r0" in mn else "-")
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("pCCE vs pmeCCE"); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1]
colors2 = {"cce": "gray", "mecce": "purple", "pcce": "steelblue", "pmecce": "seagreen"}
labels2 = {"cce": "CCE", "mecce": "ME-CCE", "pcce": "pCCE", "pmecce": "pmeCCE"}
for mn in ["cce", "mecce", "pcce", "pmecce"]:
    if mn in data3:
        ax.plot(*data3[mn], lw=2, label=labels2[mn], color=colors2[mn])
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("All Methods Comparison"); ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(T3_DIR, "comparison_plot.png"), dpi=150)
plt.close()
print(f"  Plot  → {T3_DIR}/comparison_plot.png")


# ╔══════════════════════════════════════════════════════════════╗
# ║  TEST 4 : pgCCE  —  rate=0 of pmeCCE and must match pCCE   ║
# ║           pgCCE(rate=0) == pCCE (via gCCE computation)      ║
# ╚══════════════════════════════════════════════════════════════╝
print("\n" + "─" * 60)
print("TEST 4: pgCCE verification  (pgCCE(rate=0) ≈ pCCE, pgCCE == pCCE+gCCE kernel)")
print("─" * 60)

T4_DIR = os.path.join(BASE_DIR, "4_verify_pgcce")
os.makedirs(T4_DIR, exist_ok=True)

base4_common = {
    **bath_config(BATH_SRC),
    "quantity": "coherence",
    "order": 2, "bfield": 300,
    "rbath": 500, "rdip": 300,
    "deltat": 0.04, "nstep": 50,
    "nstate": 0,
    "qspin": 0.5, "qalphams": 0.5, "qbetams": -0.5,
    "addsubclus": True, "nk": [0, 0, 0],
    "npulse": 1, "savemode": "normal",
}
base4_sk2 = {**base4_common, "sK": 2, "max_trial": 5000, "max_iter": 5000}
base4_sk1 = {**base4_common, "sK": 1, "max_trial": 5000, "max_iter": 5000}

runs4 = {
    "cce"       : (base4_common, {"method": "cce"}),
    "gcce"      : (base4_common, {"method": "gcce"}),
    "pcce_sK1"  : (base4_sk1,   {"method": "pcce"}),
    "pcce_sK2"  : (base4_sk2,   {"method": "pcce"}),
    "pgcce_sK1" : (base4_sk1,   {"method": "pgcce"}),
    "pgcce_sK2" : (base4_sk2,   {"method": "pgcce"}),
}

data4 = {}
for mname, (bconf, extra) in runs4.items():
    cfg = {**bconf, **extra, "outfile": os.path.join(T4_DIR, f"ccex_{mname}")}
    ok, err = run_ccex(cfg, os.path.join(T4_DIR, f"ccex_{mname}.json"))
    if ok:
        t, v = load_ccex(os.path.join(T4_DIR, f"ccex_{mname}"))
        data4[mname] = (t, v)
        save_dat(os.path.join(T4_DIR, f"{mname}.dat"), t, v)
        print(f"  {mname:12s}: OK")
    else:
        print(f"  {mname:12s}: FAILED  {err[:100]}")

t_pgcce2, v_pgcce2 = data4["pgcce_sK2"]
status4 = "PASS" if abs(v_pgcce2[0].real - 1.0) < 1e-6 else "FAIL"
results_summary.append(("pgCCE sanity check (t=0 → 1.0)", status4,
                        abs(v_pgcce2[0].real - 1.0), ""))

# numerical table
print(f"\n  {'Method':>12}  {'t=0.0':>10}  {'t=1.0':>10}  {'t=2.0':>10}")
for mn in runs4:
    if mn in data4:
        t, v = data4[mn]
        i25 = 25 if 25 < len(t) else -1
        print(f"  {mn:>12}  {v[0].real:10.6f}  {v[i25].real:10.6f}  {v[-1].real:10.6f}")

with open(os.path.join(T4_DIR, "summary.txt"), "w") as f:
    f.write("TEST 4: pgCCE Verification  (sK=1 and sK=2)\n")
    f.write("=" * 60 + "\n")
    f.write("pgCCE = pCCE clustering + gCCE computation kernel\n\n")
    f.write(f"{'Method':>12}  {'t=0.0':>12}  {'t=1.0':>12}  {'t=2.0':>12}\n")
    f.write("-" * 52 + "\n")
    for mn in runs4:
        if mn in data4:
            t, v = data4[mn]
            i25 = 25 if 25 < len(t) else -1
            f.write(f"{mn:>12}  {v[0].real:12.8f}  {v[i25].real:12.8f}  {v[-1].real:12.8f}\n")

fig, axes = plt.subplots(1, 3, figsize=(19, 5))
fig.suptitle('pgCCE Verification  (S=0.5, e-spin bath, B=300 G, Hahn echo, order=2)', fontsize=13)

colors_4 = {
    "cce":       ("gray",      "-",  "CCE"),
    "gcce":      ("steelblue", "-",  "gCCE"),
    "pcce_sK1":  ("tomato",    "--", "pCCE (sK=1)"),
    "pcce_sK2":  ("tomato",    "-.", "pCCE (sK=2)"),
    "pgcce_sK1": ("seagreen",  "--", "pgCCE (sK=1)"),
    "pgcce_sK2": ("seagreen",  "-.", "pgCCE (sK=2)"),
}

# Panel 1: CCE family (CCE, pCCE sK=1, pCCE sK=2)
ax = axes[0]
for mn in ["cce", "pcce_sK1", "pcce_sK2"]:
    if mn in data4:
        col, ls, lbl = colors_4[mn]
        ax.plot(*data4[mn], lw=2, color=col, ls=ls, label=lbl)
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("CCE kernel  (CCE vs pCCE)"); ax.legend(); ax.grid(True, alpha=0.3)

# Panel 2: gCCE family (gCCE, pgCCE sK=1, pgCCE sK=2)
ax = axes[1]
for mn in ["gcce", "pgcce_sK1", "pgcce_sK2"]:
    if mn in data4:
        col, ls, lbl = colors_4[mn]
        ax.plot(*data4[mn], lw=2, color=col, ls=ls, label=lbl)
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("gCCE kernel  (gCCE vs pgCCE)"); ax.legend(); ax.grid(True, alpha=0.3)

# Panel 3: sK=1 comparison (all methods)
ax = axes[2]
for mn in ["cce", "gcce", "pcce_sK1", "pgcce_sK1", "pcce_sK2", "pgcce_sK2"]:
    if mn in data4:
        col, ls, lbl = colors_4[mn]
        ax.plot(*data4[mn], lw=2, color=col, ls=ls, label=lbl)
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("All methods"); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(T4_DIR, "comparison_plot.png"), dpi=150)
plt.close()
print(f"  Plot  → {T4_DIR}/comparison_plot.png")


# ╔══════════════════════════════════════════════════════════════╗
# ║  TEST 5 : pmegCCE  —  rate=0 must match pgCCE               ║
# ║           pmegCCE = pCCE clustering + ME-gCCE kernel        ║
# ╚══════════════════════════════════════════════════════════════╝
print("\n" + "─" * 60)
print("TEST 5: pmegCCE verification  (pmegCCE(rate=0) ≈ pgCCE, rate>0 → decoherence)")
print("─" * 60)

T5_DIR = os.path.join(BASE_DIR, "5_verify_pmegcce")
os.makedirs(T5_DIR, exist_ok=True)

base5_common = {
    **bath_config(BATH_SRC),
    "quantity": "coherence",
    "order": 2, "bfield": 300,
    "rbath": 500, "rdip": 300,
    "deltat": 0.04, "nstep": 50,
    "nstate": 0,
    "qspin": 0.5, "qalphams": 0.5, "qbetams": -0.5,
    "addsubclus": True, "nk": [0, 0, 0],
    "npulse": 1, "savemode": "normal",
}
base5_sk2 = {**base5_common, "sK": 2, "max_trial": 5000, "max_iter": 5000}

runs5 = {
    "pgcce"      : (base5_sk2, {"method": "pgcce"}),
    "pmegcce_r0" : (base5_sk2, {"method": "pmegcce",
                                "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": 0.0},
                                                   {"bath_name": "e", "operator": "+-", "rate": 0.0}]}),
    "pmegcce"    : (base5_sk2, {"method": "pmegcce",
                                "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": rate_json},
                                                   {"bath_name": "e", "operator": "+-", "rate": rate_json}]}),
}

data5 = {}
for mname, (bconf, extra) in runs5.items():
    cfg = {**bconf, **extra, "outfile": os.path.join(T5_DIR, f"ccex_{mname}")}
    ok, err = run_ccex(cfg, os.path.join(T5_DIR, f"ccex_{mname}.json"))
    if ok:
        t, v = load_ccex(os.path.join(T5_DIR, f"ccex_{mname}"))
        data5[mname] = (t, v)
        save_dat(os.path.join(T5_DIR, f"{mname}.dat"), t, v)
        print(f"  {mname:14s}: OK")
    else:
        print(f"  {mname:14s}: FAILED  {err[:100]}")

if "pgcce" in data5 and "pmegcce_r0" in data5:
    t_pgcce, v_pgcce   = data5["pgcce"]
    t_pr0,   v_pr0     = data5["pmegcce_r0"]
    diffs5 = np.abs(v_pgcce.real - v_pr0.real)
    max_d5 = float(np.max(diffs5))
    status5 = "PASS" if max_d5 < 1e-6 else "FAIL"
else:
    max_d5, status5 = np.nan, "FAIL"
results_summary.append(("pmegCCE rate=0 == pgCCE", status5, max_d5, ""))

print(f"\n  pmegcce(rate=0) vs pgcce   max|diff| = {max_d5:.3e}  [{status5}]")

print(f"\n  {'Method':>14}  {'t=0.0':>10}  {'t=1.0':>10}  {'t=2.0':>10}")
for mn in runs5:
    if mn in data5:
        t, v = data5[mn]
        i25 = 25 if 25 < len(t) else -1
        print(f"  {mn:>14}  {v[0].real:10.6f}  {v[i25].real:10.6f}  {v[-1].real:10.6f}")

with open(os.path.join(T5_DIR, "summary.txt"), "w") as f:
    f.write("TEST 5: pmegCCE Verification  (sK=2)\n")
    f.write("=" * 55 + "\n")
    f.write("pmegCCE = pCCE clustering + ME-gCCE computation kernel\n\n")
    f.write(f"pmegcce(rate=0) vs pgcce  max|diff| = {max_d5:.3e}  [{status5}]\n\n")
    f.write(f"{'Method':>14}  {'t=0.0':>12}  {'t=1.0':>12}  {'t=2.0':>12}\n")
    f.write("-" * 56 + "\n")
    for mn in runs5:
        if mn in data5:
            t, v = data5[mn]
            i25 = 25 if 25 < len(t) else -1
            f.write(f"{mn:>14}  {v[0].real:12.8f}  {v[i25].real:12.8f}  {v[-1].real:12.8f}\n")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('pmegCCE Verification  (S=0.5, e-spin bath, B=300 G, sK=2, Hahn echo)', fontsize=13)

ax = axes[0]
colors_5 = {"pgcce": "steelblue", "pmegcce_r0": "tomato", "pmegcce": "seagreen"}
labels_5 = {"pgcce": "pgCCE", "pmegcce_r0": "pmegCCE (rate=0)", "pmegcce": "pmegCCE (rate=1 rad/ms)"}
for mn in ["pgcce", "pmegcce_r0", "pmegcce"]:
    if mn in data5:
        ax.plot(*data5[mn], lw=2, label=labels_5[mn], color=colors_5[mn],
                ls="--" if "r0" in mn else "-")
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("pgCCE vs pmegCCE"); ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1]
colors5b = {"pgcce": "steelblue", "pmegcce": "seagreen"}
labels5b = {"pgcce": "pgCCE (sK=2)", "pmegcce": "pmegCCE (sK=2)"}
# Also add pgcce_sK1 from TEST4 if available
all5b = {"pgcce": data5.get("pgcce"), "pmegcce": data5.get("pmegcce")}
for mn, d in all5b.items():
    if d is not None:
        ax.plot(*d, lw=2, label=labels5b[mn], color=colors5b[mn])
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("pgCCE vs pmegCCE (with dissipation)"); ax.legend(); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(T5_DIR, "comparison_plot.png"), dpi=150)
plt.close()
print(f"  Plot  → {T5_DIR}/comparison_plot.png")


# ╔══════════════════════════════════════════════════════════════╗
# ║  TEST 6 (optional): Compare all methods vs PyCCE            ║
# ╚══════════════════════════════════════════════════════════════╝
pycce_ok = False
if USE_PYCCE:
    print("\n" + "─" * 60)
    print("TEST 6: All methods vs PyCCE  (S=0.5)")
    print("─" * 60)

    T6_DIR = os.path.join(BASE_DIR, "6_pycce_comparison")
    os.makedirs(T6_DIR, exist_ok=True)

    try:
        import pycce

        atoms = pycce.random_bath('e', [5e3, 5e3, 5e3],
                                  density=1e16, density_units='cm-3', seed=2)
        cen = pycce.CenterArray(spin=0.5, alpha=[1, 0], beta=[0, 1])

        time_space = np.linspace(0, 2, 21)
        order, r_bath, r_dipole, B = 2, 500, 300, 300

        pycce_data = {}
        for pname, pmethod in [("cce", "cce"), ("gcce", "gcce"),
                                ("mecce", "mecce"), ("megcce", "megcce")]:
            atoms_run = atoms.copy()
            if pname in ("mecce", "megcce"):
                atoms_run.add_single_jump('p', rate=1.0, units='rad')
                atoms_run.add_single_jump('m', rate=1.0, units='rad')
            sim = pycce.Simulator(spin=cen, bath=atoms_run, order=order,
                                  r_bath=r_bath, r_dipole=r_dipole,
                                  magnetic_field=B)
            res = sim.compute(time_space, method=pmethod, quantity='coherence',
                              nbstates=0, pulses=1)
            pycce_data[pname] = (time_space, res)
            save_dat(os.path.join(T6_DIR, f"pycce_{pname}.dat"), time_space, res)
            print(f"  PyCCE {pname}: OK")

        # Export bath for CCEX
        bath_xyz = atoms['xyz']
        bath_names = atoms['N']
        nspin = len(bath_xyz)
        for fname, content in [
            ("qubit_espin", "0.0\t0.0\t0.0\n"),
            ("gyro_espin",  "e  0.5  -17608.597050\n"),
        ]:
            with open(os.path.join(T6_DIR, fname), 'w') as f:
                f.write(content)
        with open(os.path.join(T6_DIR, "bath_espin"), 'w') as f:
            f.write(f"{nspin:.5f}\t0.00000\t0.00000\t0.00000\n")
            for i in range(nspin):
                f.write(f"{bath_xyz[i,0]:.5f}\t{bath_xyz[i,1]:.5f}\t{bath_xyz[i,2]:.5f}\t{bath_names[i]}\n")

        # CCEX runs
        base6 = {
            **bath_config(T6_DIR),
            "quantity": "coherence",
            "order": order, "bfield": B,
            "rbath": r_bath, "rdip": r_dipole,
            "deltat": float(time_space[1] - time_space[0]),
            "nstep": len(time_space) - 1,
            "nstate": 0, "qspin": 0.5, "qalphams": 0.5, "qbetams": -0.5,
            "addsubclus": True, "nk": [0, 0, 0],
            "npulse": 1, "savemode": "normal",
        }
        runs6 = {
            "cce"   : {"method": "cce"},
            "gcce"  : {"method": "gcce"},
            "mecce" : {"method": "mecce",
                       "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": rate_json},
                                          {"bath_name": "e", "operator": "+-", "rate": rate_json}]},
            "megcce": {"method": "megcce",
                       "jump_operators": [{"bath_name": "e", "operator": "-+", "rate": rate_json},
                                          {"bath_name": "e", "operator": "+-", "rate": rate_json}]},
        }
        ccex_data6 = {}
        for mname, extra in runs6.items():
            cfg = {**base6, **extra, "outfile": os.path.join(T6_DIR, f"ccex_{mname}")}
            ok, err = run_ccex(cfg, os.path.join(T6_DIR, f"ccex_{mname}.json"))
            if ok:
                t, v = load_ccex(os.path.join(T6_DIR, f"ccex_{mname}"))
                ccex_data6[mname] = (t, v)
                save_dat(os.path.join(T6_DIR, f"ccex_{mname}.dat"), t, v)
                print(f"  CCEX  {mname}: OK")
            else:
                print(f"  CCEX  {mname}: FAILED  {err[:100]}")

        # Stats
        with open(os.path.join(T6_DIR, "summary.txt"), "w") as f:
            f.write("TEST 6: CCEX vs PyCCE Comparison  (S=0.5)\n")
            f.write("=" * 60 + "\n\n")
            for mn in ["cce", "gcce", "mecce", "megcce"]:
                if mn not in ccex_data6:
                    continue
                pt, pv = pycce_data[mn]
                ct, cv = ccex_data6[mn]
                diffs = []
                for i, tc in enumerate(ct):
                    idx = np.argmin(np.abs(pt - tc))
                    if np.abs(pt[idx] - tc) < 1e-6:
                        diffs.append(abs(pv[idx] - cv[i]))
                diffs = np.array(diffs)
                max_d = float(np.max(diffs)) if len(diffs) else np.nan
                status = "PASS" if max_d < 1e-3 else "WARN"
                f.write(f"{mn.upper():8s}  max|diff| = {max_d:.3e}  [{status}]\n")
                results_summary.append((f"CCEX {mn} vs PyCCE", status, max_d, ""))

        # Plot
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('CCEX vs PyCCE  (S=0.5, e-spin bath, B=300 G, Hahn echo)', fontsize=13)
        for idx, mn in enumerate(["cce", "gcce", "mecce", "megcce"]):
            ax = axes[idx // 2][idx % 2]
            pt, pv = pycce_data[mn]
            ax.plot(pt, pv.real, 'b-', lw=2, label=f'PyCCE {mn.upper()}')
            if mn in ccex_data6:
                ct, cv = ccex_data6[mn]
                ax.plot(ct, cv.real, 'r--', lw=2, label=f'CCEX {mn.upper()}')
            ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
            ax.set_title(mn.upper()); ax.legend(); ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(T6_DIR, "comparison_plot.png"), dpi=150)
        plt.close()
        print(f"  Plot  → {T6_DIR}/comparison_plot.png")
        pycce_ok = True

    except ImportError:
        print("  PyCCE not found — skipping Test 6")
        results_summary.append(("CCEX vs PyCCE", "SKIP", np.nan, "PyCCE not installed"))


# ╔══════════════════════════════════════════════════════════════╗
# ║  SUMMARY PLOT  (combines key results from all tests)        ║
# ╚══════════════════════════════════════════════════════════════╝
print("\n" + "─" * 60)
print("Generating overall summary plot ...")

fig, axes = plt.subplots(1, 5, figsize=(27, 5))
fig.suptitle('CCEX Method Verification Summary  (S=0.5, e-spin, B=300 G, Hahn echo)',
             fontsize=13)

# Panel 1 — ME-CCE
ax = axes[0]
if "cce" in data1 and "mecce_r0" in data1 and "mecce" in data1:
    ax.plot(*data1["cce"],      lw=2,   color="steelblue", label="CCE")
    ax.plot(*data1["mecce_r0"], lw=2,   color="tomato",  ls="--", label="ME-CCE (rate=0)")
    ax.plot(*data1["mecce"],    lw=1.5, color="seagreen", ls="-.", label="ME-CCE (rate=1 rad/ms)")
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("TEST 1: ME-CCE"); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# Panel 2 — ME-gCCE
ax = axes[1]
if "gcce" in data2 and "megcce_r0" in data2 and "megcce" in data2:
    ax.plot(*data2["gcce"],      lw=2,   color="steelblue", label="gCCE")
    ax.plot(*data2["megcce_r0"], lw=2,   color="tomato",  ls="--", label="ME-gCCE (rate=0)")
    ax.plot(*data2["megcce"],    lw=1.5, color="seagreen", ls="-.", label="ME-gCCE (rate=1 rad/ms)")
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("TEST 2: ME-gCCE"); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# Panel 3 — pmeCCE
ax = axes[2]
for mn, col, lbl, ls in [
    ("pcce",    "steelblue", "pCCE",              "-"),
    ("pmecce_r0","tomato",   "pmeCCE (rate=0)",   "--"),
    ("pmecce",  "seagreen",  "pmeCCE (rate=1 r/ms)", "-."),
]:
    if mn in data3:
        ax.plot(*data3[mn], lw=2, color=col, ls=ls, label=lbl)
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("TEST 3: pmeCCE"); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# Panel 4 — pgCCE
ax = axes[3]
for mn, col, lbl, ls in [
    ("cce",   "gray",       "CCE",  "-"),
    ("gcce",  "steelblue",  "gCCE", "--"),
    ("pcce_sK2",  "tomato",     "pCCE (sK=2)", "-."),
    ("pgcce_sK2", "seagreen",   "pgCCE (sK=2)","-"),
]:
    if mn in data4:
        ax.plot(*data4[mn], lw=2, color=col, ls=ls, label=lbl)
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("TEST 4: pgCCE"); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# Panel 5 — pmegCCE
ax = axes[4]
for mn, col, lbl, ls in [
    ("pgcce",      "steelblue", "pgCCE",             "-"),
    ("pmegcce_r0", "tomato",    "pmegCCE (rate=0)",  "--"),
    ("pmegcce",    "seagreen",  "pmegCCE (rate>0)", "-."),
]:
    if mn in data5:
        ax.plot(*data5[mn], lw=2, color=col, ls=ls, label=lbl)
ax.set_xlabel("Time (ms)"); ax.set_ylabel("Re[L(t)]")
ax.set_title("TEST 5: pmegCCE"); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(BASE_DIR, "summary_plot.png"), dpi=150)
plt.close()
print(f"  Summary plot → {BASE_DIR}/summary_plot.png")


# ╔══════════════════════════════════════════════════════════════╗
# ║  MASTER SUMMARY TEXT                                        ║
# ╚══════════════════════════════════════════════════════════════╝
with open(os.path.join(BASE_DIR, "summary.txt"), "w") as f:
    f.write("CCEX Verification Test Suite — Master Summary\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"{'Test':40s}  {'Status':6s}  {'Max|diff|':>12}  Note\n")
    f.write("-" * 70 + "\n")
    for name, status, maxd, note in results_summary:
        d_str = f"{maxd:.3e}" if not np.isnan(maxd) else "    N/A"
        f.write(f"{name:40s}  {status:6s}  {d_str:>12}  {note}\n")

print("\n" + "=" * 60)
print("  MASTER SUMMARY")
print("=" * 60)
print(f"  {'Test':40s}  {'Status':6s}  {'Max|diff|':>12}")
for name, status, maxd, note in results_summary:
    d_str = f"{maxd:.3e}" if not np.isnan(maxd) else "    N/A"
    print(f"  {name:40s}  {status:6s}  {d_str:>12}")
print()
print(f"  All results in: {BASE_DIR}/")
print("=" * 60)
print("Done.")
