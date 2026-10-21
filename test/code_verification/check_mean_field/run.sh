#!/usr/bin/env bash
#
# hamiltonian.common.bath_bath_mean_field "xyz" (vector_mean_field until A47, D45) against "z", where the answer is known.
# ---------------------------------------------------------------------------
# The p1_cce2e2n case (P1 with its 14N as a sub-spin, n_nuclear_spins 1, s_range [1, 1], state 1) with every P1 moved out of rbath but
# the nearest one whose axis is not the NV axis and whose 14N is not at m = 0: no bath spin is outside the cluster, so the mean field is zero and the vector mode must give the
# same curve as the z mode, byte for byte. Until review 158861 R1 the vector mode took the sub-spin's x, y field out of
# the main spin although it had never been added, and the curves differed.
#
# Usage:  bash test/code_verification/check_mean_field/run.sh      (bin/main.out; CCEX_BIN, CCEX_MPI as elsewhere)
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../../.." && pwd)"
export PYTHONPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)${PYTHONPATH:+:$PYTHONPATH}"   # input_layout.py (A64/A65)
# One run at a time per folder: a second run would delete the first one's work files (D29; review 20261008-003400 R3).
if [ "${CCEX_RUN_LOCK:-}" != "$HERE" ]; then
    exec 9>"$HERE/.run.lock"
    if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
        echo "ERROR: another check is running in $HERE (it would delete that run's work files); wait for it or use another checkout" >&2; exit 1
    fi
    CCEX_RUN_LOCK="$HERE"; export CCEX_RUN_LOCK
fi
BIN="${CCEX_BIN:-$REPO/bin/main.out}"; MPIRUN="${CCEX_MPI:-mpirun}"
CASE="$REPO/test/code_verification/check_numerics/cases/p1_cce2e2n"; INP="$REPO/test/code_verification/check_numerics/inputs/p1_1ppm"
WORK="$HERE/_work"; rm -rf "${WORK:?}"; mkdir -p "$WORK"
python3 - "$CASE" "$INP" "$WORK" <<'PY'
import json, math, sys; from input_layout import flat, nest
case, inp, work = sys.argv[1:]
q = [float(x) for x in open(f'{inp}/bath_DiaP1_1ppm_Defect').read().split()[:3]]
rows = open(f'{inp}/bath_DiaP1_1ppm_1').read().splitlines()
# a P1 whose 14N is not at m = 0 and whose Jahn-Teller axis is not the NV axis (avaax 1), so its hyperfine tensor has
# the x z, y z parts the bug took out of the main spin
axis = [l.strip() for l in open(f'{inp}/avaax_nav4').read().splitlines()]
ex = [l.strip() for l in open(f'{inp}/exstate_1').read().splitlines()]
ok = [i for i in range(1, len(rows)) if axis[i] != '1' and float(ex[i]) != 0.0]
near = min(ok, key=lambda i: math.dist([float(x) for x in rows[i].split()[:3]], q))
print(f'kept P1 line {near}: axis {axis[near]}, 14N m {ex[near]}')
out = [rows[0]]
for i in range(1, len(rows)):
    p = rows[i].split()
    out.append(rows[i] if i == near else ' %.5f\t %.5f\t %.5f\t        %s\t' % (float(p[0])+1e5, float(p[1])+1e5, float(p[2])+1e5, p[3]))
open(f'{work}/bath_one', 'w').write('\n'.join(out)+'\n')
c = flat(json.load(open(f'{case}/ccein.json')))
for block, k in ((c['Qubit']['qubit'][0], 'file'), (c['configuration'], 'gyrofile')):   # the qubit's position file (A94)
    block[k] = f'{case}/{block[k]}'
for block, k in ((c['configuration'], 'defect_axis_file'), (c['configuration'], 'bath_state_file'), (c['configuration'], 'defect_state_file')):
    block[k] = [f'{case}/{x}' for x in block[k]]       # one per bath file (A64)
c['configuration']['bathfile'] = [f'{work}/bath_one']; c['nstep'] = 40
for name, field in (('z', 'z'), ('vector', 'xyz')):
    json.dump(nest(dict(c, hamiltonian={'common': {'bath_bath_mean_field': field}}, outfile=f'{work}/{name}')),
              open(f'{work}/{name}.json', 'w'), indent=1)
PY
for m in z vector; do
    (cd "$WORK" && "$MPIRUN" -n 1 "$BIN" -f "$m.json" > "proc_$m" 2>&1 < /dev/null) || { echo "  FAIL  $m run: see $WORK/proc_$m"; exit 1; }
done
n=$(grep -c "Number of bath spins : 1 " "$WORK/proc_vector" || true)
if [ "$n" -lt 1 ]; then echo "  FAIL  the bath is not one P1: see $WORK/proc_vector"; exit 1; fi
if cmp -s "$WORK"/z_state1_wiDiv "$WORK"/vector_state1_wiDiv && cmp -s "$WORK"/z_state1_noDiv "$WORK"/vector_state1_noDiv; then
    echo "  PASS  one P1 with its 14N sub-spin: bath_bath_mean_field xyz gives the z curve byte for byte (no spin outside the cluster)"
else
    echo "  FAIL  one P1 with its 14N sub-spin: vector and z differ"; exit 1
fi
