#!/usr/bin/env bash
#
# The built-in isotope table (src/input/isotope_table.cpp, docs/user/InputFile_General.md#gyrofile), on the natab
# cCCE2 system of inputs/Diamond_NV_natab_cce:
#   G1  no gyrofile vs a gyrofile holding the table's own 13C value        byte-identical
#   G2  a gyrofile with 13C only, some bath spins renamed 15N              runs, the log names 15N from the table
#   G3  a bath species in neither (XX)                                     refused, the log names it
#   G4  every entry of the generated table is an isotope with a positive spin
#
# Usage:  bash test/code_verification/check_input_forms/gyro.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/gyro"
rm -rf "$W"; mkdir -p "$W"
cp -r "$REG/Diamond_NV_natab_cce/inputfiles" "$W/"
cd "$W"
GAMMA=$(grep -o '"13C", 0.5f, [0-9.]*' "$REPO/src/input/isotope_table.cpp" | awk '{print $3}')
printf '13C 0.5 %s\n' "$GAMMA" > inputfiles/gyro_table
python3 - <<'PY'
import json, re; from input_layout import flat, nest
src = open("inputfiles/ccein_nstate0.json").read()
j = flat(json.loads(re.sub(r"(?m)^\s*#.*$", "", src)))
def write(name, **kw):
    k = dict(j, **kw)
    conf = dict(j["configuration"])                  # gyrofile and bathfile are the configuration block's (A64)
    for key in ("gyrofile", "bathfile"):
        if key in kw:
            conf[key] = k.pop(key)
    k["configuration"] = conf
    for key in [x for x, v in kw.items() if v is None]:
        (conf if key in ("gyrofile", "bathfile") else k).pop(key)
    open(name, "w").write(json.dumps(nest(k), indent=1))
lines = open("inputfiles/bath_Diamond_NV_natab_1").read().splitlines()
open("inputfiles/bath_15N", "w").write("\n".join([lines[0]] + [l.replace("13C", "15N") if i % 7 == 0 else l
                                                             for i, l in enumerate(lines[1:])]) + "\n")
open("inputfiles/bath_XX", "w").write("\n".join([lines[0], lines[1].replace("13C", "XX")] + lines[2:]) + "\n")
write("g1_file.json", gyrofile="./inputfiles/gyro_table", outfile="./out/g1_file")
write("g1_table.json", gyrofile=None, outfile="./out/g1_table")
write("g2.json", bathfile=["./inputfiles/bath_15N"], outfile="./out/g2")
write("g3.json", bathfile=["./inputfiles/bath_XX"], rbath=1000, outfile="./out/g3")     # XX is outside rbath 30
PY
mkdir -p out
run_ccex "$BIN_DEV" g1_file.json  g1_file.log  && run_ccex "$BIN_DEV" g1_table.json g1_table.log
cmp_raw "G1 no gyrofile == the table's value in a gyrofile" out/g1_file out/g1_table
if run_ccex "$BIN_DEV" g2.json g2.log && grep -q "15N from the built-in isotope table" g2.log; then
    echo "  PASS  G2 a species missing from the gyrofile comes from the table"; PASS=$((PASS+1))
else
    echo "  FAIL  G2 (see $W/g2.log)"; FAIL=$((FAIL+1))
fi
if ! run_ccex "$BIN_DEV" g3.json g3.log && grep -q 'species "XX" is neither in the gyrofile nor in the built-in' g3.log; then
    echo "  PASS  G3 a species in neither is refused"; PASS=$((PASS+1))
else
    echo "  FAIL  G3 (see $W/g3.log)"; FAIL=$((FAIL+1))
fi
# G4 every table entry is an isotope: mass number > 0, a positive integer or half-integer spin (review 20261008-003400 R1)
if python3 - "$REPO/src/input/isotope_table.cpp" <<'PY'
import re, sys
rows = re.findall(r'\{"(\d+)([A-Z][a-z]?)", (-?[0-9.]+)f, (-?[0-9.e+-]+)\}', open(sys.argv[1]).read())
bad = [r for r in rows if int(r[0]) <= 0 or float(r[2]) <= 0 or float(r[2]) * 2 != int(float(r[2]) * 2) or float(r[3]) == 0]
sys.exit(1 if bad or len(rows) < 100 else 0)
PY
then echo "  PASS  G4 every table entry is an isotope with a positive spin"; PASS=$((PASS+1))
else echo "  FAIL  G4 the table holds a non-isotope row"; FAIL=$((FAIL+1)); fi
echo "gyro: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
