#!/usr/bin/env bash
#
# gpu_options (docs/user/InputFile_Backend.md#tag_gpu_options, 2026-10-12, A26): the GPU settings in one block.
# G1  without the block, an empty block, the block with every default written out, and the block with every key changed
#     (min_dim, evolution_batch_size, evolution_batch_steps, cluster_batch_mb, cache_mb, packed_max_dim 0) give the same
#     output bytes on a CPU run (backend cpu: none of them is used) - single-state cCCE2 and gCCE2;
# G2  the log reports the block's values (gpu_options.min_dim ... gpu_options.packed_max_dim), so a run shows what it read;
# G3  without the key packed_max_dim is 0 in the CPU build, with backend cpu and auto (the default 4 is for the GPU build).
# The refused forms (the old flat keys, unknown keys, packed_max_dim other than 0, 2, 4, packed in bin/main.out) are in
# ../check_refusals/. The packed kernel itself needs the GPU: test/gpu_verification/check_gpu/test_small_cluster_cuda.
#
# Usage:  bash test/code_verification/check_input_forms/gpu_options.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/gpu_options"
rm -rf "$W"; mkdir -p "$W"
NUM="$REPO/test/code_verification/check_numerics"

run_case(){    # <case> <name> <json changes>
    local d="$NUM/cases/$1" o="$W/$2"; mkdir -p "$o"
    python3 -c 'import json,sys; j=json.load(open(sys.argv[1])); j.update(json.loads(sys.argv[3])); json.dump(j,open(sys.argv[2],"w"),indent=1)' \
        "$d/ccein.json" "$o.json" "$3"
    ( cd "$d" && "$MPIRUN" -n 2 "$BIN_DEV" -f "$o.json" -o "$o/out" > "$o.log" 2>&1 < /dev/null ) || { echo "    $2 exited $?"; return 1; }
}
same_bytes(){ [ -n "$(ls "$1")" ] && [ "$(ls "$1")" = "$(ls "$2")" ] && diff -rq "$1" "$2" >/dev/null; }

FORMS='none {}
empty {"gpu_options":{}}
defaults {"gpu_options":{"min_dim":256,"packed_max_dim":0,"evolution_batch_size":8,"evolution_batch_steps":64,"cluster_batch_mb":64,"cache_mb":64}}
changed {"gpu_options":{"min_dim":17,"packed_max_dim":0,"evolution_batch_size":3,"evolution_batch_steps":5,"cluster_batch_mb":7,"cache_mb":0}}'
for case in p1_cce2_single hbn_gcce2_oh; do
    ok=1
    while read -r form change; do
        run_case "$case" "$case/$form" "$(python3 -c 'import json,sys; print(json.dumps(dict(json.loads(sys.argv[1]), backend="cpu")))' "$change")" || ok=0
    done <<< "$FORMS"
    for form in empty defaults changed; do
        same_bytes "$W/$case/none" "$W/$case/$form" || { echo "    $case: the $form block changes the outputs"; ok=0; }
    done
    if [ $ok -eq 1 ]; then echo "  PASS  G1 $case: no block, an empty, a default and a changed block give the same bytes"; PASS=$((PASS+1))
    else echo "  FAIL  G1 $case (see $W/$case)"; FAIL=$((FAIL+1)); fi
    ok=1
    for want in "gpu_options.min_dim: +17" "gpu_options.evolution_batch_steps: +5" "gpu_options.cache_mb: +0" \
                "gpu_options.evolution_batch_size: +3" "gpu_options.cluster_batch_mb: +7" "gpu_options.packed_max_dim: +0"; do
        grep -Eq "^ +$want( |$)" "$W/$case/changed.log" || { echo "    $case: the log does not show $want"; ok=0; }
    done
    grep -q "input key \"gpu_options[^\"]*\" is not read" "$W/$case/changed.log" && { echo "    $case: a gpu_options key reported as not read"; ok=0; }
    if [ $ok -eq 1 ]; then echo "  PASS  G2 $case: the log reports the block's values"; PASS=$((PASS+1))
    else echo "  FAIL  G2 $case (see $W/$case/changed.log)"; FAIL=$((FAIL+1)); fi
done
# G3: without the key packed_max_dim is off in the CPU build, whatever the backend (default 4 only in bin/main_gpu.out with
# backend cuda or auto and method CCE; A26, user 2026-10-12), and backend auto runs (no refusal of the default)
ok=1
run_case p1_cce2_single auto_default '{"backend":"auto"}' || ok=0
run_case p1_cce2_single auto_exact '{"backend":"cpu","evolution_update":"exact"}' || ok=0
for log in "$W/p1_cce2_single/none.log" "$W/auto_default.log"; do
    grep -A1 -E "^ +gpu_options.packed_max_dim: +0( |$)" "$log" | grep -Eq "chosen by +: +default \(CPU build\)" ||
        { echo "    $log: packed_max_dim is not 0 by default (CPU build)"; ok=0; }
done
# backend auto runs exact (A29) and, in this build, everything on the CPU: the bytes of backend cpu with exact
same_bytes "$W/auto_exact" "$W/auto_default" || { echo "    backend auto (all on the CPU here) is not backend cpu exact"; ok=0; }
if [ $ok -eq 1 ]; then echo "  PASS  G3 without the key, packed_max_dim is 0 in the CPU build (backend cpu and auto)"; PASS=$((PASS+1))
else echo "  FAIL  G3 (see $W)"; FAIL=$((FAIL+1)); fi
echo "gpu_options: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
