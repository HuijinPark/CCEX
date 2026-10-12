#!/bin/sh
# 07_hBN_VB_natab/ensemble_CCE: CCE2 with the DFT hyperfine and quadrupole tensors: hyperfine order 3 with quadrupole order 0 and 2, and order 0 with 2
#
# Recorded run: 1 bath configuration(s), 10 MPI ranks on idm.
# Here NCONF configurations on NP ranks. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_hf0_qd2.json -I ../bath/flat/bath_hBN_VB_flat_${i} -o results/CCE2_hBN_VB_hf0_qd2_$i > results/process_CCE2_hBN_VB_hf0_qd2_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_hf3_qd0.json -I ../bath/flat/bath_hBN_VB_flat_${i} -o results/CCE2_hBN_VB_hf3_qd0_$i > results/process_CCE2_hBN_VB_hf3_qd0_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_hf3_qd2.json -I ../bath/flat/bath_hBN_VB_flat_${i} -o results/CCE2_hBN_VB_hf3_qd2_$i > results/process_CCE2_hBN_VB_hf3_qd2_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# plot
python3 $ANALYZER -d results/ -fi CCE2_hBN_VB_PREFIX_1_wiDiv -vn "PREFIX" -v hf3_qd2 hf0_qd2 -uc "us" -ntn -fo results/CCE2_hBN_VB_Comp_HF -pl_m
# plot: the two quadrupole orders at hyperfine order 3
python3 $ANALYZER -d results/ -fi CCE2_hBN_VB_PREFIX_1_wiDiv -vn "PREFIX" -v hf3_qd0 hf3_qd2 -uc "us" -ntn -fo results/CCE2_hBN_VB_Comp_Qd -pl_m
