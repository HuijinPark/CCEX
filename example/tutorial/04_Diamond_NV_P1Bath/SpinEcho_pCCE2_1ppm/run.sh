#!/bin/sh
# 04_Diamond_NV_P1Bath/SpinEcho_pCCE2_1ppm: pCCE2 with pcce.partition_size = 1..4 spins per centre, 1 ppm
#
# Recorded run: 50 bath configuration(s) x 50 sampled states, 128 MPI ranks on idm.
# Here NCONF configurations on NP ranks; NSTATE sampled states. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}; NSTATE=${NSTATE:-50}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_sK1.json -I ../bath/14N_1ppm/bath_DiaP1_1ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/CCE2_DiaP1_1ppm_sK_1_$i > results/process_CCE2_DiaP1_1ppm_sK_1_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_sK2.json -I ../bath/14N_1ppm/bath_DiaP1_1ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/CCE2_DiaP1_1ppm_sK_2_$i > results/process_CCE2_DiaP1_1ppm_sK_2_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_sK3.json -I ../bath/14N_1ppm/bath_DiaP1_1ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/CCE2_DiaP1_1ppm_sK_3_$i > results/process_CCE2_DiaP1_1ppm_sK_3_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_sK4.json -I ../bath/14N_1ppm/bath_DiaP1_1ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/CCE2_DiaP1_1ppm_sK_4_$i > results/process_CCE2_DiaP1_1ppm_sK_4_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_sK_1_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_DiaP1_1ppm_sK_1_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/CCE2_DiaP1_1ppm_sK_1_CONF_stateSTAT_ErrCorr
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_sK_2_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_DiaP1_1ppm_sK_2_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/CCE2_DiaP1_1ppm_sK_2_CONF_stateSTAT_ErrCorr
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_sK_3_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_DiaP1_1ppm_sK_3_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/CCE2_DiaP1_1ppm_sK_3_CONF_stateSTAT_ErrCorr
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_sK_4_CONF_stateSTAT_wiDiv -err_d results/ -err_fi CCE2_DiaP1_1ppm_sK_4_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/CCE2_DiaP1_1ppm_sK_4_CONF_stateSTAT_ErrCorr
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_sK_1_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/CCE2_DiaP1_1ppm_sK_1_EnsAvg
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_sK_2_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/CCE2_DiaP1_1ppm_sK_2_EnsAvg
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_sK_3_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/CCE2_DiaP1_1ppm_sK_3_EnsAvg
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_DiaP1_1ppm_sK_4_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/CCE2_DiaP1_1ppm_sK_4_EnsAvg
