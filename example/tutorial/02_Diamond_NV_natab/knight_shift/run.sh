#!/bin/sh
# 02_Diamond_NV_natab/knight_shift: CCE2 with the Knight shift on (ccein_wi_knight.json) and off (ccein_wo_knight.json) in a tilted field, plotted together
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
  mpirun -n $NP $CCEX -f ccein_wi_knight.json -I ../bath/bath_Diamond_NV_natab_${i} -o results/CCE2_Diamond_NV_natab_wi_knight_$i > results/process_CCE2_Diamond_NV_natab_wi_knight_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_wo_knight.json -I ../bath/bath_Diamond_NV_natab_${i} -o results/CCE2_Diamond_NV_natab_wo_knight_$i > results/process_CCE2_Diamond_NV_natab_wo_knight_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# plot
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_PREFIX_1_wiDiv -vn "PREFIX" -v wi_knight wo_knight -uc "us" -ntn -fo results/CCE2_Diamond_NV_natab_Comp_Knight -pl_m -ylim 0 1
