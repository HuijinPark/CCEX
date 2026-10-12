#!/bin/sh
# 06_Pulse_Manual/SpinEcho_P1_1ppm: the spin echo on the 1 ppm P1 bath as an explicit sequence, original and WAHUHA-style
#
# Recorded run: 5 bath configuration(s) x 1 sampled states, 128 MPI ranks on idm.
# Here NCONF configurations on NP ranks; NSTATE sampled states. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}; NSTATE=${NSTATE:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_origin.js -I ../../04_Diamond_NV_P1Bath/bath/14N_1ppm/bath_DiaP1_1ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/CCE2_DiaP1_1ppm_Order_2_origin_$i > results/process_CCE2_DiaP1_1ppm_Order_2_origin_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_wahuha.js -I ../../04_Diamond_NV_P1Bath/bath/14N_1ppm/bath_DiaP1_1ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/CCE2_DiaP1_1ppm_Order_2_wahuha_$i > results/process_CCE2_DiaP1_1ppm_Order_2_wahuha_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_Order_2_origin_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_DiaP1_1ppm_Order_2_origin_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/CCE2_DiaP1_1ppm_Order_2_origin_CONF_stateSTAT_ErrCorr
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_Order_2_wahuha_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_DiaP1_1ppm_Order_2_wahuha_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/CCE2_DiaP1_1ppm_Order_2_wahuha_CONF_stateSTAT_ErrCorr
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_Order_2_origin_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/CCE2_DiaP1_1ppm_Order_2_origin_EnsAvg
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_Order_2_wahuha_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/CCE2_DiaP1_1ppm_Order_2_wahuha_EnsAvg
