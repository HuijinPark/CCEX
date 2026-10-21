#!/usr/bin/env bash
#
# Comments in the input (src/input/option.cpp cJSON_ReadFccein): '#' and '!' start a comment to the end of the line only
# outside a JSON string. Until 2026-10-08 a key such as "#note" or a value such as "run#1" was cut at the '#' and the input
# refused. The natab cCCE2 input of inputs/Diamond_NV_natab_cce, plain and with comments of both kinds, a "#note" key and an
# outfile holding '#' and '!', must give the same curves byte for byte.
#
# Usage:  bash test/code_verification/check_input_forms/comments.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/comments"
rm -rf "$W"; mkdir -p "$W/out"
cp -r "$REG/Diamond_NV_natab_cce/inputfiles" "$W/"
cd "$W"
python3 - <<'PY'
import json, re; from input_layout import flat, nest
j = flat(json.loads(re.sub(r"(?m)^\s*#.*$", "", open("inputfiles/ccein_nstate0.json").read())))
plain = dict(j, outfile="./out/plain")
open("plain.json", "w").write(json.dumps(nest(plain), indent=1))
lines = json.dumps(nest(dict(j, outfile="./out/run#1!x")), indent=1).splitlines()
lines.insert(1, '  "#note": "a key with # and ! inside a string",   # a comment after a value')
lines.insert(2, '  ! a comment line of the other kind')
open("commented.json", "w").write("# a comment line\n" + "\n".join(lines) + "\n")
PY
run_ccex "$BIN_DEV" plain.json plain.log && run_ccex "$BIN_DEV" commented.json commented.log
cmp_raw "C1 comments, a #note key and an outfile with # and ! give the plain input's curves" "out/plain" "out/run#1!x"
echo "comments: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
