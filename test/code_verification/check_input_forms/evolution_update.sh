#!/usr/bin/env bash
#
# evolution_update (docs/user/InputFile_Backend.md, 2026-10-11): an explicit "exact" gives the bytes of the default run;
# "recursive" (A16 b, the phases of the time grid carried over from the previous step) runs on the paths it covers
# and stays within 1e-8 of exact in the raw _noDiv. NV + P1 10 ppm (inputs/Diamond_NV_P1_10ppm_{cce,gcce}, nstep 40:
# the recurrence runs between steps 0, 16 and 32): cCCE2 and gCCE2 (evolution=vector) with one sampled state; since
# A18 (2026-10-11) also the ensemble cCCE2 (nstate 0: Hahn, CPMG-8, and a Manual sequence of unequal intervals) and
# gCCE2 evolution=matrix (the ensemble Hahn echo, and CPMG-8 with one state). The refusals ("recursive" elsewhere,
# other values) are in ../check_refusals/.
#
# Usage:  bash test/code_verification/check_input_forms/evolution_update.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/evolution_update"
rm -rf "$W"; mkdir -p "$W"
# name  input folder  input file  changes
VARIANTS="cce cce ccein_nstate1.json {}
gcce gcce ccein_nstate1.json {}
cce_ensemble cce ccein_nstate0.json {}
cce_ensemble_cpmg8 cce ccein_nstate0.json {\"npulse\":8}
cce_ensemble_manual cce ccein_nstate0.json {\"npulse\":2,\"pulsename\":\"Manual\",\"sequence\":[[0.2,\"X\",180],[0.5,\"X\",180],[0.3,\"I\",0]]}
gcce_matrix gcce ccein_nstate0.json {}
gcce_matrix_cpmg8 gcce ccein_nstate1.json {\"evolution\":\"matrix\",\"npulse\":8}"
NAMES=""
while read -r name sys file change; do
    NAMES="$NAMES $name"; mkdir -p "$W/$name/out"; cp -r "$REG/Diamond_NV_P1_10ppm_$sys/inputfiles" "$W/$name/"
    python3 - "$W/$name" "$file" "$change" <<'PY2'
import json, re, sys
d, file, change = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
j = json.loads(re.sub(r"(?m)^\s*#.*$", "", open(f"{d}/inputfiles/{file}").read()))
j.update(change)
json.dump(dict(j, outfile="./out/default"), open(f"{d}/default.json", "w"), indent=1)
for value in ("exact", "recursive"):
    json.dump(dict(j, outfile=f"./out/{value}", evolution_update=value), open(f"{d}/{value}.json", "w"), indent=1)
PY2
done <<< "$VARIANTS"
for sys in $NAMES; do
    cd "$W/$sys"; ok=1
    for v in default exact recursive; do run_ccex "$BIN_DEV" $v.json $v.log || { echo "    $sys $v exited $?"; ok=0; }; done
    ls out/default*_noDiv >/dev/null 2>&1 || ok=0
    for f in out/default*_noDiv out/default*_wiDiv; do
        cmp -s "$f" "${f/default/exact}" || { echo "    DIFF $sys $f ${f/default/exact}"; ok=0; }
    done
    if [ $ok -eq 1 ]; then echo "  PASS  E1 $sys: evolution_update exact gives the default run's bytes"; PASS=$((PASS+1))
    else echo "  FAIL  E1 $sys"; FAIL=$((FAIL+1)); fi
    if grep -Eq '^ +evolution_update +: +recursive' recursive.log && python3 - <<'PY'
import cmath, glob, math, sys
def read(p): return [complex(l.split()[1].replace("+-", "-")) for l in open(p) if l.strip()]
def same(x, y):    # nan in both is the same; nan on one side only, or nan vs inf, is not (|nan - y| would be skipped by max)
    return all(a == b or (math.isnan(a) and math.isnan(b)) for a, b in ((x.real, y.real), (x.imag, y.imag)))
worst, n, bad = 0.0, 0, 0
for f in sorted(glob.glob("out/exact*_noDiv")):
    a, b = read(f), read(f.replace("exact", "recursive"))
    assert len(a) == len(b) == 40, (f, len(a), len(b))
    for x, y in zip(a, b):
        if cmath.isfinite(x) and cmath.isfinite(y): worst = max(worst, abs(x - y))
        elif not same(x, y): bad += 1
    n += 1
print(f"    recursive - exact, raw _noDiv: max |difference| {worst:.3e} over {n} file(s), {bad} non-finite mismatch(es)")
sys.exit(0 if n and worst <= 1e-8 and bad == 0 else 1)
PY
    then echo "  PASS  E2 $sys: evolution_update recursive runs and agrees with exact to 1e-8"; PASS=$((PASS+1))
    else echo "  FAIL  E2 $sys (see $W/$sys/recursive.log)"; FAIL=$((FAIL+1)); fi
done
echo "evolution_update: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
