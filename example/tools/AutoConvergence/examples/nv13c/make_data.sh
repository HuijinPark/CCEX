#!/bin/sh
# NV-13C test data for the automatic convergence (WSL-sized): 70 diamond 13C baths (natural abundance, 110 A cell, rbath up to about 55 A)
# from BathGeneratorSpinCluster example 1_13C_only_diamond with a fixed seed, and 30 state files per bath.
#   sh make_data.sh [BGSC_DIR]        (default ~/bin/code_develop/BathGeneratorSpinCluster)
set -e
BGSC=${1:-$HOME/bin/code_develop/BathGeneratorSpinCluster}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$HERE/data" && cd "$HERE/data"
cp "$BGSC/examples/1_13C_only_diamond/POSCAR_diamond" "$BGSC/examples/1_13C_only_diamond/POSCAR_NV_defect" .
sed 's/^ cell_size \t = 30/ cell_size \t = 110/' "$BGSC/examples/1_13C_only_diamond/opt.txt" > opt.txt
# the generator writes every file and then may abort at interpreter exit ("double free"); the outputs are checked below
python "$BGSC/BathGeneratorSpinCluster.py" --opt opt.txt -i POSCAR_diamond -I POSCAR_NV_defect -o config --size 110 \
    --defect N --defect_center 0.0 0.0 0.0 --range 70 --mpi_seed_base 20261007 > generate.log 2>&1 || true
grep -q "Success creating the configure" generate.log && [ -s config_70 ] || { echo "generation failed: see data/generate.log"; exit 1; }
head -1 config_defect > qubit
printf '13C 0.5 6.728284\n' > gyro_13C
python "$HERE/make_states.py" . 30
