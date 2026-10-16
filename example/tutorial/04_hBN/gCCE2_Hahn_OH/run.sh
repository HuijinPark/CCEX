#!/bin/sh
# 04_hBN/gCCE2_Hahn_OH: gCCE2 with and without the Overhauser field of the bath on the qubit
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
PLOT=$ROOT/example/tutorial/plot_coherence.py
NP=${NP:-4}
mkdir -p results
for name in wi_OH wo_OH; do
  mpirun -n $NP $CCEX -f ccein_$name.json -o results > results/ccein_$name.log
done
python3 $PLOT results/coherence.png --unit ns "with OH=results/ccein_wi_OH_state1_wiDiv" "without=results/ccein_wo_OH_state1_wiDiv"
