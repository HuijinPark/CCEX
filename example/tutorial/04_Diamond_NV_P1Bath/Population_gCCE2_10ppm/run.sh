#!/bin/sh
# 04_Diamond_NV_P1Bath/Population_gCCE2_10ppm: gCCE2 population (density matrix, `pulse.count=0`) of the NV in a 10 ppm P1 bath
#
# Recorded run: 1 bath configuration(s), 4 MPI ranks on idm.
# Here NCONF configurations on NP ranks. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
mkdir -p results

mpirun -n $NP $CCEX -f ccein.json -o results/gCCE2_DiaP1_10ppm_1 > results/process_gCCE2_DiaP1_10ppm_1
