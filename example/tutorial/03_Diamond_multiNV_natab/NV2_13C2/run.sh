#!/bin/sh
# 03_Diamond_multiNV_natab/NV2_13C2: two NV centres with a 13C pair: gCCE density matrix, Ramsey (ccein_NV6_13C2_ramsey.json) and Hahn echo
#
# Recorded run: 1 bath configuration(s) x 1 sampled states, 1 MPI ranks on idm.
# Here NCONF configurations on NP ranks. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
mkdir -p results

mpirun -n $NP $CCEX -f ccein_NV6_13C2_hahnecho.json -I bath -o results/CCE2_Diamond_NV_natab_NV6_13C2_hahnecho_1 > results/process_CCE2_Diamond_NV_natab_NV6_13C2_hahnecho_1
mpirun -n $NP $CCEX -f ccein_NV6_13C2_ramsey.json -I bath -o results/CCE2_Diamond_NV_natab_NV6_13C2_ramsey_1 > results/process_CCE2_Diamond_NV_natab_NV6_13C2_ramsey_1
