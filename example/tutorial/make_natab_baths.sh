#!/bin/sh
# make_natab_baths.sh N OUTDIR [SIZE]: N natural-abundance 13C bath configurations around an NV, OUTDIR/config_1 .. config_N
# and OUTDIR/config_defect (the NV position), made with tools/BathGenerator from 00_quickstart/bath/genbath/ (diamond
# cells, 1.1 % 13C) in a SIZE A cell (default 110: rbath up to 53 A). Configuration k uses the random seed k, so every
# call makes the same files; files already there are kept. About 5 s per configuration (SIZE 110).
set -e
N=$1; OUT=$2; SIZE=${3:-110}
[ -n "$OUT" ] || { echo "usage: make_natab_baths.sh N OUTDIR [SIZE]"; exit 2; }
ROOT=$(git rev-parse --show-toplevel)
GEN=$ROOT/example/tutorial/00_quickstart/bath/genbath
mkdir -p $OUT
for k in $(seq 1 $N); do
  [ -s $OUT/config_$k ] && continue
  # the generator may abort at interpreter exit after writing ("double free"); the file is checked instead
  python3 -c 'import sys, random, runpy, numpy as np; k = int(sys.argv[1]); np.random.seed(k); random.seed(k); sys.argv = sys.argv[2:]; runpy.run_path(sys.argv[0], run_name="__main__")' \
    $k $ROOT/example/tools/BathGenerator/BathGenerator.py -o $OUT/config_$k -i $GEN/POSCAR_pure -I $GEN/POSCAR_NV -p $GEN/prob_natab \
    --size $SIZE --Defect V --Defect_center 0.5 0.5 0.5 > $OUT/process_$k 2>&1 || true
  [ -s $OUT/config_$k ] || { echo "bath $k was not made: see $OUT/process_$k"; exit 1; }
  mv $OUT/config_${k}_defect $OUT/config_defect
done
