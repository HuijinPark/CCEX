#!/bin/sh
# 01_Arbitrary_system: gCCE order 1 on a hand-written three-file system; the density matrix (`quantity=dm`) of one bath state
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

mpirun -n $NP $CCEX -f ccein.json -I bath -o results/CCE2_Diamond_NV_natab_1 > results/process_CCE2_Diamond_NV_natab_1
