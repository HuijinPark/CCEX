#!/bin/sh
# 07_hBN_VB_natab/single-sample_gCCE_DPS: gCCE2 with and without the `OH` term (DPS)
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
  mpirun -n $NP $CCEX -f ccein_wi_OH.json -I ../bath/flat/bath_hBN_VB_flat_${i} -o results/CCE2_hBN_VB_wi_OH_$i > results/process_CCE2_hBN_VB_wi_OH_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_wo_OH.json -I ../bath/flat/bath_hBN_VB_flat_${i} -o results/CCE2_hBN_VB_wo_OH_$i > results/process_CCE2_hBN_VB_wo_OH_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# plot
python3 $ANALYZER -d results/ -fi CCE2_hBN_VB_PREFIX_1_state1_wiDiv -vn "PREFIX" -v wi_OH wo_OH -uc "ns" -ntn -fo results/CCE2_hBN_VB_wi_wo_OH -pl_m
