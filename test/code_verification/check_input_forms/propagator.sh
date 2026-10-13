#!/usr/bin/env bash
#
# propagator (docs/user/InputFile_Backend.md; A2, 2026-10-12): the string and object forms, on inputs of
# ../check_numerics/cases/ (one sampled state, 1 MPI rank).
#   F1  "eigen" and {"method": "eigen"} give the default run's bytes; "sparse" and {"method": "sparse", "min_dim": 4096}
#       give those of evolution_update exact (every cluster here is below 4096, the default min_dim, so sparse leaves them
#       on the eigen path, and without the key sparse takes exact where the default is otherwise recursive, A29);
#   F2  {"method": "sparse", "min_dim": 1} (every cluster through the Krylov expansion) runs and is within 1e-8 of eigen
#       with exact in the raw _noDiv and in the raw _wiDiv at |L| <= 1, non-finite points equal on both sides: cCCE2 (P1, NV-13C) and gCCE2
#       evolution=vector (P1 with zero-field splitting). (The Manual-sequence gCCE2 and pCCE2 P1 inputs have diverged
#       _wiDiv points (|L| > 1) where eigen built with icpx and with g++ already differ by up to 2.8; see
#       (internal notes))
#   F3  malformed values are refused at the start with their message.
# The combinations sparse does not cover are refused in ../check_refusals/.
#
# Usage:  bash test/code_verification/check_input_forms/propagator.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/propagator"
rm -rf "$W"; mkdir -p "$W"
CASES="$REPO/test/code_verification/check_numerics/cases"

# variant <case> <name> <python edit of d>: the case's input with the edit, run in the case's folder (relative paths)
variant(){
    local case="$1" name="$2" edit="$3"
    mkdir -p "$W/$case"
    python3 - "$CASES/$case/ccein.json" "$W/$case/$name.json" <<PY
import json, sys
d = json.load(open(sys.argv[1]))
$edit
json.dump(d, open(sys.argv[2], "w"), indent=1)
PY
    ( cd "$CASES/$case" && run_ccex "$BIN_DEV" "$W/$case/$name.json" "$W/$case/$name.log" -o "$W/$case/$name" )
}

for case in p1_cce2_single p1_manual_hahn; do
    variant $case default "pass" || { echo "  FAIL  F1 $case: the default run exited $?"; FAIL=$((FAIL+1)); continue; }
    variant $case exact "d['evolution_update'] = 'exact'" || { echo "  FAIL  F1 $case: the exact run exited $?"; FAIL=$((FAIL+1)); continue; }
    for form in "eigen_string|d['propagator'] = 'eigen'" "eigen_object|d['propagator'] = {'method': 'eigen'}" \
                "sparse_string|d['propagator'] = 'sparse'" "sparse_4096|d['propagator'] = {'method': 'sparse', 'min_dim': 4096}"; do
        name="${form%%|*}"; edit="${form#*|}"
        # sparse takes evolution_update exact when the key is absent (the default is recursive elsewhere, A29)
        ref=default; [ "${name#sparse}" != "$name" ] && ref=exact
        if variant $case "$name" "$edit"; then cmp_raw "F1 $case: $name gives the bytes of the $ref run" "$W/$case/${ref}_state1" "$W/$case/${name}_state1"
        else echo "  FAIL  F1 $case $name exited"; FAIL=$((FAIL+1)); fi
    done
done

for case in p1_cce2_single natab_cce2_single p1_gcce2_zfs; do
    [ -f "$W/$case/exact.json" ] || variant $case exact "d['evolution_update'] = 'exact'"
    if variant $case sparse_all "d['propagator'] = {'method': 'sparse', 'min_dim': 1}" &&
       grep -Eq '^ +propagator +: +sparse' "$W/$case/sparse_all.log" &&
       python3 - "$W/$case" <<'PY'
import cmath, glob, math, sys
w = sys.argv[1]
def read(p): return [[complex(t) for t in l.split()[1:]] for l in open(p) if l.strip() and not l.startswith('#')]
def same(x, y): return all(a == b or (math.isnan(a) and math.isnan(b)) for a, b in ((x.real, y.real), (x.imag, y.imag)))
worst, bad, n = 0.0, 0, 0
# judged as the other checks: _noDiv at every finite point, _wiDiv where |L| <= 1 (diverged points apart)
for kind in ("noDiv", "wiDiv"):
    for f in sorted(glob.glob(f"{w}/exact_state*_{kind}")):
        for ra, rb in zip(read(f), read(f.replace("/exact_", "/sparse_all_"))):
            for x, y in zip(ra, rb):
                if cmath.isfinite(x) and cmath.isfinite(y):
                    if kind == "noDiv" or abs(x) <= 1: worst = max(worst, abs(x - y))
                elif not same(x, y): bad += 1
                n += 1
print(f"    sparse (min_dim 1) - eigen exact, raw _noDiv and _wiDiv at |L| <= 1: max |difference| {worst:.3e} over {n} values, {bad} non-finite mismatch(es)")
sys.exit(0 if n and worst <= 1e-8 and bad == 0 else 1)
PY
    then echo "  PASS  F2 $case: sparse through every cluster agrees with eigen (exact) to 1e-8"; PASS=$((PASS+1))
    else echo "  FAIL  F2 $case (see $W/$case/sparse_all.log)"; FAIL=$((FAIL+1)); fi
done

# refuse <name> <expected message fragment> <python edit>
refuse(){
    local name="$1" want="$2" edit="$3"
    variant p1_cce2_single "$name" "$edit"
    local rc=$?
    if [ $rc -eq 0 ]; then echo "  FAIL  F3 $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return; fi
    if ! grep -qF -- "$want" "$W/p1_cce2_single/$name.log"; then
        echo "  FAIL  F3 $name : refused (exit=$rc) but not with the expected message"; echo "        want: $want"
        grep -m1 -i "error" "$W/p1_cce2_single/$name.log" | sed 's/^/        got : /'; FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  F3 $name (exit=$rc)"; PASS=$((PASS+1))
}
refuse bad_value "current propagator options (fast) is not available" "d['propagator'] = 'fast'"
refuse bad_number "propagator must be a string (eigen, expm or sparse) or an object" "d['propagator'] = 1"
refuse bad_array "propagator must be a string (eigen, expm or sparse) or an object" "d['propagator'] = ['sparse']"
refuse no_method "propagator.method must be a string" "d['propagator'] = {'min_dim': 4096}"
refuse bad_method "current propagator options (krylov) is not available" "d['propagator'] = {'method': 'krylov'}"
refuse unknown_key "propagator has no key \"mindim\"" "d['propagator'] = {'method': 'sparse', 'mindim': 4096}"
refuse min_dim_eigen "propagator.min_dim and propagator.tol belong to method sparse" "d['propagator'] = {'method': 'eigen', 'min_dim': 4096}"
refuse tol_eigen "propagator.min_dim and propagator.tol belong to method sparse" "d['propagator'] = {'method': 'eigen', 'tol': 1e-12}"
refuse min_dim_zero "propagator.min_dim must be a positive 32-bit integer" "d['propagator'] = {'method': 'sparse', 'min_dim': 0}"
refuse min_dim_fraction "propagator.min_dim must be a positive 32-bit integer" "d['propagator'] = {'method': 'sparse', 'min_dim': 1.5}"
refuse min_dim_string "propagator.min_dim must be a finite JSON number" "d['propagator'] = {'method': 'sparse', 'min_dim': '4096'}"
refuse tol_zero "propagator.tol must be a number between 0 and 1" "d['propagator'] = {'method': 'sparse', 'tol': 0}"
refuse tol_one "propagator.tol must be a number between 0 and 1" "d['propagator'] = {'method': 'sparse', 'tol': 1}"
refuse tol_string "propagator.tol must be a finite JSON number" "d['propagator'] = {'method': 'sparse', 'tol': '1e-12'}"

echo "propagator: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
