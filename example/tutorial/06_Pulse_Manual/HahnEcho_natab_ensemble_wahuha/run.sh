#!/bin/sh
# 06_Pulse_Manual/HahnEcho_natab_ensemble_wahuha: the same sequence given in the WAHUHA style of the input
#
# Recorded run: 1 bath configuration(s), 5 MPI ranks on idm.
# Here NCONF configurations on NP ranks. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein.js -I ../bath/bath_Diamond_NV_natab_${i} -o results/CCE2_Diamond_NV_natab_$i > results/process_CCE2_Diamond_NV_natab_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_CONF_wiDiv -err_d results/ -err_fi CCE2_Diamond_NV_natab_CONF_noDiv -vn "CONF" -v $(seq -s ' ' 1 $NCONF) -uc "ms" -err -ntn -fo results/CCE2_Diamond_NV_natab_CONF_ErrCorr
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_CONF_ErrCorr -vn "CONF" -v $(seq -s ' ' 1 $NCONF) -uc "ms" -ea -ntn -fo results/CCE2_Diamond_NV_natab_EnsAvg -pl -fit -fit_initp 2 -fit_space 5
