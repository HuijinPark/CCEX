#!/usr/bin/env bash
#
# column_evolution (docs/user/InputFile_Backend.md#tag_column_evolution; on by default since 2026-10-15, A62):
# the check_numerics case p1_cce2e2n (one sampled state, NV + P1 with its 14N) at order 4, rbath 170, nstep 100 and one
# cluster of order 4 (max_clusters_per_order): one cluster of bath dimension 1296, the others 216 and below.
# C1  without the block the run is on with min_dim 256 and columns 64: the bytes of an explicit
#     {"enabled": true, "min_dim": 256, "columns": 64}, and the log prints both values;
# C2  {"enabled": false} is the path before: the bytes of a run whose min_dim no cluster reaches, and the log prints
#     neither value;
# C3  the default and the path before agree in the raw _noDiv to 1e-10 (the column path sums the same products in
#     another order; whether the written digits differ is only noted - on WSL, 2026-10-15, they do not, so the path is
#     judged by the log, the value cceEvolve dispatches on, as in evolution_update.sh).
# The numerics of the column path itself are test_cce_columns (make numeric-test); the refusals of the block are in
# ../check_refusals/.
#
# Usage:  bash test/code_verification/check_input_forms/column_evolution.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/column_evolution"
rm -rf "$W"; mkdir -p "$W"
CASE="$REPO/test/code_verification/check_numerics/cases/p1_cce2e2n"
BASE='{"order":4,"rbath":170,"nstep":100,"max_clusters_per_order":[0,0,0,0,1],"profile":true}'

# run <name> <column_evolution block or "">: the case with BASE and the block, 2 ranks, into $W/<name>/
run(){
    local o="$W/$1"; mkdir -p "$o"
    python3 -c 'import json,sys; from input_layout import flat, nest
j=flat(json.load(open(sys.argv[1]))); j.update(json.loads(sys.argv[3]))
if sys.argv[4]: j["column_evolution"]=json.loads(sys.argv[4])
j["outfile"]="out"; json.dump(nest(j),open(sys.argv[2],"w"),indent=1)' "$CASE/ccein.json" "$o.json" "$BASE" "$2"
    ( cd "$CASE" && "$MPIRUN" -n 2 "$BIN_DEV" -f "$o.json" -o "$o" > "$o.log" 2>&1 < /dev/null ) || { echo "    $1 exited $?"; return 1; }
}
same_bytes(){ [ -n "$(ls "$1")" ] && [ "$(ls "$1")" = "$(ls "$2")" ] && diff -rq "$1" "$2" >/dev/null; }
logs_on(){ grep -Eq "column_evolution\.min_dim *: +$2( |$)" "$1" && grep -Eq "column_evolution\.columns *: +$3( |$)" "$1"; }

ok=1
run default "" || ok=0
run explicit '{"enabled": true, "min_dim": 256, "columns": 64}' || ok=0
run off '{"enabled": false}' || ok=0
run unreached '{"min_dim": 100000}' || ok=0
grep -q "dim=1296 clusters=1 " "$W/default.log" || { echo "    the case has no cluster of dimension 1296 (see $W/default.log)"; ok=0; }

c1=$ok
same_bytes "$W/default" "$W/explicit" || { echo "    the default run is not the explicit block's bytes"; c1=0; }
logs_on "$W/default.log" 256 64 || { echo "    the default run's log does not print min_dim 256, columns 64"; c1=0; }
if [ $c1 -eq 1 ]; then echo "  PASS  C1 without the block: on, min_dim 256, columns 64"; PASS=$((PASS+1))
else echo "  FAIL  C1 (see $W)"; FAIL=$((FAIL+1)); fi

c2=$ok
same_bytes "$W/off" "$W/unreached" || { echo "    {\"enabled\": false} is not the bytes of the path before"; c2=0; }
grep -q "column_evolution\.min_dim" "$W/off.log" && { echo "    the log of {\"enabled\": false} prints column_evolution"; c2=0; }
if [ $c2 -eq 1 ]; then echo "  PASS  C2 {\"enabled\": false}: the path before"; PASS=$((PASS+1))
else echo "  FAIL  C2 (see $W)"; FAIL=$((FAIL+1)); fi

same_bytes "$W/default" "$W/off" && echo "    the default and the path before write the same bytes here"
if [ $ok -eq 1 ] && python3 - "$W/default/out_state1_noDiv" "$W/off/out_state1_noDiv" <<'PY'
import sys
def read(p): return [complex(l.split()[1].replace("+-", "-")) for l in open(p) if l.strip()]
a, b = read(sys.argv[1]), read(sys.argv[2])
assert len(a) == len(b) == 100, (len(a), len(b))
worst = max(abs(x - y) for x, y in zip(a, b))
print(f"    default - path before, raw _noDiv: max |difference| {worst:.3e}")
sys.exit(0 if worst <= 1e-10 else 1)
PY
then echo "  PASS  C3 the default agrees with the path before to 1e-10"; PASS=$((PASS+1))
else echo "  FAIL  C3 (see $W)"; FAIL=$((FAIL+1)); fi
echo "column_evolution: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
