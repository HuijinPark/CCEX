#!/bin/sh
# 12_rdip_ladder/ensemble_CCE2: CCE2 Hahn echo of the NV of 02 (rbath 50, rdip 12) with rdip_ladder 4, 6, 8, 10 - one run, five curves - checked against a separate run at rdip 6; then the same input with evolution_update exact and recursive, timed and compared
#
# Recorded run: 1 bath configuration (02's bath/), 4 MPI ranks on WSL (2026-10-12): ladder 6.8 s, rdip 6 1.0 s, exact 5.5 s, recursive 4.6 s (CCEX wall times; 21 s in all).
# Here NCONF configurations on NP ranks (02's bath/ holds configuration 1; more come from its bath/genbath/).
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
BATH=../../02_Diamond_NV_natab/bath/bath_Diamond_NV_natab
mkdir -p results

for i in $(seq 1 $NCONF); do
  for name in ladder rdip6 exact recursive; do
    mpirun -n $NP $CCEX -f ccein_$name.json -I ${BATH}_$i -o results/CCE2_${name}_$i > results/process_CCE2_${name}_$i
  done
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# the wall times CCEX prints
grep -H 'JOB DONE' results/process_CCE2_*_1
# every rung equals the separate run (here rdip 6), the main output equals the run without the ladder (both recursive),
# and exact against recursive: within 1e-9 at every point
python3 -c 'import sys, numpy as np; sys.path.insert(0, sys.argv[1]); from coherence import read_coherence
L = lambda n: read_coherence("results/" + n)[1]
for a, b in [("CCE2_ladder_1_rdip6_wiDiv", "CCE2_rdip6_1_wiDiv"), ("CCE2_ladder_1_wiDiv", "CCE2_recursive_1_wiDiv"),
             ("CCE2_exact_1_wiDiv", "CCE2_recursive_1_wiDiv"), ("CCE2_exact_1_noDiv", "CCE2_recursive_1_noDiv")]:
    print("max |%s - %s| = %.3e" % (a, b, abs(L(a) - L(b)).max()))' $(dirname $ANALYZER)
# plot: the curves of the one ladder run (top) and how far each rung is from the main output at rdip 12 (bottom);
# the fast ripple is the 13C Larmor precession (about 2 us at 500 G), which the 0.2 us steps resolve
python3 -c 'import sys, numpy as np; sys.path.insert(0, sys.argv[1]); from coherence import read_coherence
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, (ax, bx) = plt.subplots(2, 1, figsize=(6.4, 6.0), sharex=True)
t, ref = read_coherence("results/CCE2_ladder_1_wiDiv", time_scale=1e3)
for r in ["4", "6", "8", "10"]:
    L = read_coherence("results/CCE2_ladder_1_rdip%s_wiDiv" % r, time_scale=1e3)[1]
    ax.plot(t, np.abs(L), lw=0.6, label="rdip %s A" % r)
    bx.plot(t, np.abs(L - ref), lw=0.6, label="rdip %s A" % r)
ax.plot(t, np.abs(ref), "k", lw=0.6, label="rdip 12 A (main output)")
ax.set_ylabel("|L(t)|"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
bx.set_yscale("log"); bx.set_ylim(1e-6, 1); bx.set_xlabel("time (us)"); bx.set_ylabel("|L_rung - L_rdip12|"); bx.grid(alpha=0.3)
fig.tight_layout(); fig.savefig("results/CCE2_ladder_rungs.png", dpi=150)' $(dirname $ANALYZER)
