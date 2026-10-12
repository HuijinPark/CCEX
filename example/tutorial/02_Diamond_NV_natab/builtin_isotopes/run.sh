#!/bin/sh
# 02_Diamond_NV_natab/builtin_isotopes: ensemble_CCE's CCE2 Hahn echo with the gyro file (ccein_gyrofile.json) and without it (ccein_builtin.json: 13C from the built-in isotope table), compared and plotted together
#
# Added 2026-10-12 (no idm record). WSL, 4 ranks, 1 configuration: 4 s for the two runs and the plot.
# Here NCONF configurations on NP ranks.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  for name in gyrofile builtin; do
    mpirun -n $NP $CCEX -f ccein_$name.json -I ../bath/bath_Diamond_NV_natab_${i} -o results/CCE2_${name}_$i > results/process_CCE2_${name}_$i
  done
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# where the 13C spin and gamma came from (the run log says it per species)
grep -H -E 'built-in isotope table|Read GyroFile' results/process_CCE2_builtin_1 results/process_CCE2_gyrofile_1 || true
# the two gammas differ by 2e-7 (relative): over 2 ms at 500 G a phase of about 1e-3 rad on the 13C precession, so the curves differ by about 1e-4
python3 -c 'import sys, numpy as np; sys.path.insert(0, sys.argv[1]); from coherence import read_coherence
a = read_coherence("results/CCE2_gyrofile_1_wiDiv")[1]; b = read_coherence("results/CCE2_builtin_1_wiDiv")[1]
print("max |L_builtin - L_gyrofile| = %.3e" % abs(b - a).max())' $(dirname $ANALYZER)
# plot
python3 $ANALYZER -d results/ -fi CCE2_PREFIX_1_wiDiv -vn "PREFIX" -v gyrofile builtin -uc "ms" -ntn -fo results/CCE2_Comp_isotopes -pl_m
