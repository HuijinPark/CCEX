#!/bin/sh
# 08_Bundle_clustering/CCE2_adaptive: the bundled CCE2 of the NV of 02 (bath copied to ../bath) (rbath 20, 2 ms) refined by the merge policy: compute, measure the error of each bundle boundary, merge where it pays, repeat
#
# Added 2026-09-23 (no idm record). WSL, 2 ranks: 8 min, nearly all of it in the probes.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
mkdir -p results

# exit status 0: the curve settled within atol; 2: the policy stopped first (here: no merge left that improves it)
mpirun -n $NP $CCEX -f ccein.json > results/process_CCE2_adaptive_1 || echo "refinement stopped with status $?"

# --- post-processing
# the rounds: bundles before and after, the largest local error seen, and the reason for stopping
python3 -c 'import json; s = json.load(open("results/adaptive/summary.json"))
for r in s["rounds"]: print("round %d: %d bundles, max local scaled error %.3g" % (r["iteration"], r["bundles"], r["max_local_scaled_error"]))
print("status:", s["status"], "| final bundles:", len(s["final_bundles"]), "| steps:", len(s["runs"]))'
