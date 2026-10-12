#!/bin/sh
# 07_hBN_VB_natab/hyperfine_mediated_term: CCE2 with and without the hyperfine-mediated term
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
  mpirun -n $NP $CCEX -f ccein_wi_hfmedi.js -I ../bath/flat/bath_hBN_VB_flat_${i} -o results/CCE2_hBN_VB_wi_hfmedi_$i > results/process_CCE2_hBN_VB_wi_hfmedi_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_wo_hfmedi.js -I ../bath/flat/bath_hBN_VB_flat_${i} -o results/CCE2_hBN_VB_wo_hfmedi_$i > results/process_CCE2_hBN_VB_wo_hfmedi_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# plot
python3 $ANALYZER -d results/ -fi CCE2_hBN_VB_PREFIX_1_wiDiv -vn "PREFIX" -v wi_hfmedi wo_hfmedi -uc "us" -ntn -fo results/CCE2_hBN_VB_Comp_HFmedi -pl_m
