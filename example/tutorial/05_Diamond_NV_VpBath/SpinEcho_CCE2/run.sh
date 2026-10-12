#!/bin/sh
# 05_Diamond_NV_VpBath/SpinEcho_CCE2: CCE2 spin echo of the NV in a Vp defect bath at 10, 12 and 50 ppm
#
# Recorded run: 1 bath configuration(s) x 1 sampled states, 128 MPI ranks on idm.
# Here NCONF configurations on NP ranks; NSTATE sampled states. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}; NSTATE=${NSTATE:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_10ppm.js -I ../bath/10ppm/bath_DiaVp_10ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/CCE2_DiaVp_10ppm_Order_2_$i > results/process_CCE2_DiaVp_10ppm_Order_2_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_12ppm.js -I ../bath/12ppm/bath_DiaVp_12ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/CCE2_DiaVp_12ppm_Order_2_$i > results/process_CCE2_DiaVp_12ppm_Order_2_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_50ppm.js -I ../bath/50ppm/bath_DiaVp_50ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/CCE2_DiaVp_50ppm_Order_2_$i > results/process_CCE2_DiaVp_50ppm_Order_2_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_DiaVp_10ppm_Order_2_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_DiaVp_10ppm_Order_2_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/CCE2_DiaVp_10ppm_Order_2_CONF_stateSTAT_ErrCorr
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_DiaVp_12ppm_Order_2_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_DiaVp_12ppm_Order_2_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/CCE2_DiaVp_12ppm_Order_2_CONF_stateSTAT_ErrCorr
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_DiaVp_50ppm_Order_2_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_DiaVp_50ppm_Order_2_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/CCE2_DiaVp_50ppm_Order_2_CONF_stateSTAT_ErrCorr
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_DiaVp_10ppm_Order_2_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/CCE2_DiaVp_10ppm_Order_2_EnsAvg
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_DiaVp_12ppm_Order_2_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/CCE2_DiaVp_12ppm_Order_2_EnsAvg
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_DiaVp_50ppm_Order_2_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/CCE2_DiaVp_50ppm_Order_2_EnsAvg
