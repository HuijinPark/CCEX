#!/bin/sh
# DEER bath-pulse regression: six short runs, then check_deer.py verifies them.
#
# Twelve spins, cluster order 1, non-sampled (-N 0), 41 time steps - the whole thing takes a
# few seconds.  Order 1 is deliberate: the closed form is exact there, so the checks compare
# against an answer that is known rather than against a previous run.
#
#   sh run_test.sh                     # uses the binary at the repo root
#   sh run_test.sh /path/to/main.out   # or point it somewhere else
#
# Exit status is 0 only if every check passes.

set -e
HERE=$(cd "$(dirname "$0")" && pwd)
CODE=${1:-$HERE/../../../bin/main.out}
GYRO=${GYROFILE:-$HERE/../../Diamond_NV_P1Bath/bath/DiaP1_gyro}

if [ ! -x "$CODE" ]; then
    echo "binary not found: $CODE" >&2
    echo "usage: sh run_test.sh [/path/to/main.out]" >&2
    exit 1
fi
if [ ! -f "$GYRO" ]; then
    echo "gyro file not found: $GYRO   (override with GYROFILE=...)" >&2
    exit 1
fi

# numpy is needed by the checker; module purge can hide it, so look for an interpreter that has it
PY=${PYTHON:-}
if [ -z "$PY" ]; then
    for c in python3 /opt/anaconda3/2024.10-1/bin/python /usr/bin/python3; do
        if command -v "$c" > /dev/null 2>&1 && "$c" -c "import numpy" > /dev/null 2>&1; then
            PY=$c; break
        fi
    done
fi
if [ -z "$PY" ]; then
    echo "no python with numpy found (set PYTHON=/path/to/python)" >&2
    exit 1
fi

NP=${NP:-2}
for c in nopulse onres offres multi_hit multi_miss nearpulse; do
    d=$HERE/$c
    mkdir -p "$d/res_CCE" "$d/proc"
    rm -f "$d/res_CCE"/*
    sed "s|GYROFILE|$GYRO|" "$d/ccein.js" > "$d/ccein.run.js"
    ( cd "$d" && mpirun -np $NP "$CODE" -f ./ccein.run.js \
        -I "$HERE/bath_DEERtest_12spin.txt" -a "$HERE/avaax_DEERtest" \
        -s Random -S Random -o "./res_CCE/test" -N 0 \
        < /dev/null > proc/run.log 2>&1 )
    printf "  ran %-11s -> %s files\n" "$c" "$(ls "$d/res_CCE" | wc -l)"
done

"$PY" "$HERE/check_deer.py"
