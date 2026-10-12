#!/bin/sh
# 11_Auto_convergence/ensemble_CCE2: the automatic convergence of the CCE2 Hahn echo of the NV in natural-abundance diamond (02's system): the time window, rdip, rbath and the number of bath configurations found by main.out itself, then the production
#
# Recorded run: NBATH=40 bath configurations, 4 MPI ranks on WSL (2026-10-12): 3 min 24 s for the baths, 42 s for the convergence.
# Here NBATH configurations on NP ranks. Set them to repeat the recorded run (every config_<k> in results/bath/ is used).
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
NP=${NP:-4}
NBATH=${NBATH:-40}
mkdir -p results/bath

# the bath configurations: 13C at natural abundance in a 110 A cell around the NV (covers rbath up to 53 A), made with
# tools/BathGenerator from 02's genbath/ cell files, configuration k with the random seed k (the same files every time).
# The generator may abort at interpreter exit after writing ("double free"); the file is checked instead.
GEN=$ROOT/example/tutorial/02_Diamond_NV_natab/bath/genbath
for k in $(seq 1 $NBATH); do
  [ -s results/bath/config_$k ] && continue
  python3 -c 'import sys, random, runpy, numpy as np; k = int(sys.argv[1]); np.random.seed(k); random.seed(k); sys.argv = sys.argv[2:]; runpy.run_path(sys.argv[0], run_name="__main__")' \
    $k $ROOT/example/tools/BathGenerator/BathGenerator.py -o results/bath/config_$k -i $GEN/POSCAR_pure -I $GEN/POSCAR_NV -p $GEN/prob_natab \
    --size 110 --Defect V --Defect_center 0.5 0.5 0.5 > results/bath/process_$k 2>&1 || true
  [ -s results/bath/config_$k ] || { echo "bath $k was not made: see results/bath/process_$k"; exit 1; }
  mv results/bath/config_${k}_defect results/bath/config_defect
done

# the automatic convergence: rank 0 runs the runner (python3 with numpy, scipy, matplotlib), which starts every CCEX run
# with mpirun -n NP. Exit status 0: converged; 3: unresolved (the summary says what is missing); calling again continues.
mpirun -n $NP $CCEX -f ccein.json > results/process_auto || echo "automatic convergence ended with status $?"

# --- post-processing: what it decided (results/auto/result/summary.txt; figures 1_window.png ... 6_converged_coherence.png
# and the mean production curve coherence_mean.txt in the same folder)
cat results/auto/result/summary.txt
ls results/auto/result/*.png
