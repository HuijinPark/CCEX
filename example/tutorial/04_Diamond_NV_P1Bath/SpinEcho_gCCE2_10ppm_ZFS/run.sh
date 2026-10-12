#!/bin/sh
# 04_Diamond_NV_P1Bath/SpinEcho_gCCE2_10ppm_ZFS: gCCE2 at 10 ppm and zero field, with and without the NV zero-field splitting (`job` is the field-sweep job the two came from)
#
# Recorded run: 1 bath configuration(s) x 1 sampled states, 4 MPI ranks on idm.
# Here NCONF configurations on NP ranks; NSTATE sampled states. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}; NSTATE=${NSTATE:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_with_ZFS.js -I ../bath/14N_10ppm/bath_DiaP1_10ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_$i > results/process_gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_with_ZFS_v2.js -I ../bath/14N_10ppm/bath_DiaP1_10ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/v2_gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_$i > results/process_v2_gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_$i
done
for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein_without_ZFS.js -I ../bath/14N_10ppm/bath_DiaP1_10ppm_${i} -a Random$i -s Random -S Random -N $NSTATE -o results/gCCE2_DiaP1_10ppm_B0_0G_without_ZFS_$i > results/process_gCCE2_DiaP1_10ppm_B0_0G_without_ZFS_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_CONF_stateSTAT_wiDiv -err_d results/ -err_fi gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_CONF_stateSTAT_ErrCorr
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi v2_gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_CONF_stateSTAT_wiDiv -err_d results/ -err_fi v2_gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/v2_gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_CONF_stateSTAT_ErrCorr
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi gCCE2_DiaP1_10ppm_B0_0G_without_ZFS_CONF_stateSTAT_wiDiv -err_d results/ -err_fi gCCE2_DiaP1_10ppm_B0_0G_without_ZFS_CONF_stateSTAT_noDiv -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -err -ntn -fo results/gCCE2_DiaP1_10ppm_B0_0G_without_ZFS_CONF_stateSTAT_ErrCorr
# ensemble average
python3 $ANALYZER -d results/ -fi gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_EnsAvg
# ensemble average
python3 $ANALYZER -d results/ -fi v2_gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/v2_gCCE2_DiaP1_10ppm_B0_0G_with_ZFS_EnsAvg
# ensemble average
python3 $ANALYZER -d results/ -fi gCCE2_DiaP1_10ppm_B0_0G_without_ZFS_CONF_stateSTAT_ErrCorr -vn CONF STAT -v $(seq -s ' ' 1 $NCONF) -v $(seq -s ' ' 1 $NSTATE) -uc "ms" -ea -ntn -fo results/gCCE2_DiaP1_10ppm_B0_0G_without_ZFS_EnsAvg
