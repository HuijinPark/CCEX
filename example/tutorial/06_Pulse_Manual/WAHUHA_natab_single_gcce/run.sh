#!/bin/sh
# 06_Pulse_Manual/WAHUHA_natab_single_gcce: the four-pulse WAHUHA sequence, gCCE2
#
# Recorded run: 1 bath configuration(s) x 1 sampled states, 5 MPI ranks on idm.
# Here NCONF configurations on NP ranks. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein.json -I ../bath/bath_Diamond_NV_natab_${i} -o results/CCE2_Diamond_NV_natab_$i > results/process_CCE2_Diamond_NV_natab_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_Diamond_NV_natab_CONF_stateSTAT_noDiv -vn "CONF" "STAT" -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 1) -uc "ms" -err -ntn -fo results/CCE2_Diamond_NV_natab_CONF_stateSTAT_ErrCorr
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_CONF_stateSTAT_ErrCorr -vn "CONF" "STAT" -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 1) -uc "ms" -ea -ntn -fo results/CCE2_Diamond_NV_natab_EnsAvg -pl -fit -fit_initp 2 -fit_space 5
