#!/bin/sh
# 07_hBN_VB_natab/single-sample_gCCE: gCCE2 density matrix of one sampled bath state
#
# Recorded run: 1 bath configuration(s) x 1 sampled states, 10 MPI ranks on idm.
# Here NCONF configurations on NP ranks. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein.js -I ../bath/flat/bath_hBN_VB_flat_${i} -o results/gCCE2_hBN_VB_$i > results/process_gCCE2_hBN_VB_$i
done
