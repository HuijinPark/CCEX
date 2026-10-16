#!/usr/bin/env bash
#
# A80 (2026-10-16): pCCE drops the spins left over when the bath does not divide into partitions of partition_size
# (shrink_restspins deletes the last ones), and the state files are read after that. Until A80 the reader set the state
# of the dropped spins too, in memory already freed (ASan: heap-buffer-overflow in setBathStates; a segfault at the end
# of a 64-rank run). P1 1 ppm cCCE2 pCCE (check_numerics/inputs/p1_1ppm, rbath 500), partition_size 3, two states read
# from bath_state_file, axes from defect_axis_file.
#   P1  the bath has spins left over (selected spins % 3 != 0), and the run finishes
#   P2  byte for byte the run whose bath, state and axis files lack the lines of the dropped spins
# A plain build need not show the old fault (it wrote into freed memory and went on); a build with CCEX_SANITIZE=1
# (CCEX_BIN=<that build>/bin/main.out) stops at it, and P1 fails.
#
# Usage:  bash test/code_verification/check_input_forms/pcce_rest.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/pcce_rest"
rm -rf "$W"; mkdir -p "$W/full" "$W/cut" "$W/out"
cd "$W"
python3 - "$REPO/test/code_verification/check_numerics" <<'PY'
import json, math, sys
num = sys.argv[1]; src = num + "/inputs/p1_1ppm/"
case = json.load(open(num + "/cases/p1_pcce2_sk2/ccein.json"))
rbath, sk = case["rbath"], 3
q = [float(x) for x in open(src + "bath_DiaP1_1ppm_Defect").read().split()[:3]]
bath = open(src + "bath_DiaP1_1ppm_1").read().splitlines()
# the spins the reader selects, as file line indices (line 0 is the header)
sel = [i for i, l in enumerate(bath) if i and len(l.split()) == 4
       and math.dist([float(x) for x in l.split()[:3]], q) <= rbath]
rest = len(sel) % sk
assert rest, f"{len(sel)} spins divide into partitions of {sk}: nothing is left over"
drop = set(sel[-rest:])                       # shrink_restspins deletes the last ones
open("rest", "w").write(f"{len(sel)} {rest}\n")
for name in ("bath_DiaP1_1ppm_1", "state_1", "state_2", "avaax_nav12"):
    lines = open(src + name).read().splitlines(True)
    open("full/" + name, "w").writelines(lines)
    open("cut/" + name, "w").writelines(l for i, l in enumerate(lines) if i not in drop)
for tag in ("full", "cut"):
    c = json.loads(json.dumps(case))
    c["cluster"]["pcce"].update(partition_size=sk, trials=100, max_iterations=100, print_trials=0)
    c["Qubit"]["qubit"] = [{"file": src + "bath_DiaP1_1ppm_Defect"}]
    c["configuration"].update(gyrofile=src + "DiaP1_gyro", bathfile=[tag + "/bath_DiaP1_1ppm_1"],
                              bath_state_file=[tag + "/state_"], defect_axis_file=[tag + "/avaax_nav12"], s_range=[1, 2])
    c["nstep"] = 50
    c["output"] = {"savemode": "coherence", "file": "./out/" + tag}
    json.dump(c, open(tag + ".json", "w"), indent=1)
PY
[ $? -eq 0 ] || { echo "  FAIL  inputs (see above)"; exit 1; }
read nsel rest < rest

ok=1
run_ccex "$BIN_DEV" full.json full.log || ok=0
grep -q "JOB DONE" full.log || ok=0
if [ $ok -eq 1 ]; then echo "  PASS  P1 $nsel spins, $rest left over by partition_size 3: the run with their states in the files finishes"; PASS=$((PASS+1))
else echo "  FAIL  P1 (see $W/full.log)"; tail -20 full.log | sed 's/^/      /'; FAIL=$((FAIL+1)); fi

# A89: the default log names them, once: "pCCE spins left out : <rest> of <selected> (partition_size 3), ..."
if [ "$(grep -c "^      pCCE spins left out : $rest of $nsel (partition_size 3), the last read" full.log)" = 1 ]; then
    echo "  PASS  P1 log: $(grep -m1 'pCCE spins left out' full.log | sed 's/^ *//')"; PASS=$((PASS+1))
else echo "  FAIL  P1 log has no single 'pCCE spins left out : $rest of $nsel' line (see $W/full.log)"; FAIL=$((FAIL+1)); fi

run_ccex "$BIN_DEV" cut.json cut.log
for s in 1 2; do cmp_raw "P2 state $s: == the bath without the $rest left-over spins" "out/full_state$s" "out/cut_state$s"; done
echo "pcce_rest: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
