#!/bin/sh
# 01_ensemble_and_states/cCCE2_Hahn_states: 20 sampled bath states (nstate 20) against the ensemble (nstate 0) of the same bath
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
mpirun -n $NP $CCEX -f ccein.json -o results > results/ccein.log
# the same input with nstate 0 (the ensemble), written by set_input_keys.py as results/ensemble.json
python3 $ROOT/example/tools/InputKeys/set_input_keys.py ccein.json results/ensemble.json configuration.nstate=0
mpirun -n $NP $CCEX -f results/ensemble.json -o results > results/ensemble.log
python3 $PLOT results/coherence.png "ensemble=results/ensemble_wiDiv" "mean of 20 states=results/ccein_state*_wiDiv" \
  "state 1=results/ccein_state1_wiDiv"
