#!/bin/sh
# NV-13C test data for the automatic convergence (WSL-sized): 70 diamond 13C baths (natural abundance, 110 A cell, rbath up
# to about 53 A) made by the tutorials' example/tutorial/make_natab_baths.sh (example/tools/BathGenerator, seed k; with
# mpi4py over NP MPI ranks, default 4), and 30 state files per bath.
#   sh make_data.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE" && git rev-parse --show-toplevel)
sh "$ROOT/example/tutorial/make_natab_baths.sh" 70 "$HERE/data" 110
cd "$HERE/data"
cp config_defect qubit
printf '13C 0.5 6.728284\n' > gyro_13C
python "$HERE/make_states.py" . 30
