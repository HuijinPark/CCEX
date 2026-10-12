#!/bin/sh
# 02_Diamond_NV_natab/ensemble_gCCE: gCCE2 on one configuration, single qubit (ccein_single.js) and two qubits (ccein_multi2.js), plotted together
#
# Recorded run: 1 bath configuration(s), 20 MPI ranks on idm.
# Here NCONF configurations on NP ranks. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_multi2.js -I ../bath/bath_Diamond_NV_natab_${i} -o results/CCE2_Diamond_NV_natab_multi2_$i > results/process_CCE2_Diamond_NV_natab_multi2_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_single.js -I ../bath/bath_Diamond_NV_natab_${i} -o results/CCE2_Diamond_NV_natab_single_$i > results/process_CCE2_Diamond_NV_natab_single_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# plot
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_PREFIX_1_wiDiv -vn "PREFIX" -v single multi2 -uc "us" -ntn -fo results/CCE2_Diamond_NV_natab_Comp_Multi -pl_m -ylim 0 1
