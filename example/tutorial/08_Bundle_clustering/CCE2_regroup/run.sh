#!/bin/sh
# 08_Bundle_clustering/CCE2_regroup: the bundled CCE2 of the NV in six P1 spins (five sampled states) refined by regroup: where some state has |L| > 1.01, find the bundle pairs that cause it and regroup the bundles around them, until every curve is physical
#
# Added 2026-10-12 (no idm record). WSL, 4 ranks: 1 s (five CCEX steps, 49 probes).
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results

# exit status 0: physical (every |L| <= 1.01); 2: stopped first (no culprit, no candidate, budget or change limit)
mpirun -n $NP $CCEX -f ccein.json > results/process_CCE2_regroup_1 || echo "regroup stopped with status $?"

# --- post-processing
# the bundles before and after, the unphysical states (U) and the culprit pairs of each round, and the end
python3 -c 'import json; s = json.load(open("results/adaptive/summary.json")); g = s["regroup"]
for r in g["rounds"]:
    v = r.get("verification") or {}
    print("round %d: bundles %s, %d unphysical states (worst: state %d, time point %d, |L| - 1.01 = %.3g)" % (r["round"], r["units"], r["U"], r["s"], r["t"], r["excess"]))
    print("  culprit pairs (g > 0): %s" % "; ".join("%s-%s g=%.3g" % (a, b, f["g"]) for (a, b), f in zip(r["culprits"], r["factors"])))
    print("  accepted: %s -> rerun: %s unphysical states, excess %.3g" % (r.get("accepted"), v.get("U_actual"), v.get("excess_actual", float("nan"))))
print("status:", s["status"], "| changes:", g["changes"], "| final bundles:", s["final_bundles"])'
