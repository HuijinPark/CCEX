#!/usr/bin/env bash
#
# save_nodiv (docs/user/InputFile_Output.md, 2026-10-08): false writes _wiDiv alone, byte-identical to the _wiDiv of the
# default run; the default still writes both; a value that is not a JSON boolean is refused. natab cCCE2 input of
# inputs/Diamond_NV_natab_cce, ensemble and one sampled state.
#
# Usage:  bash test/code_verification/check_input_forms/nodiv.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/nodiv"
rm -rf "$W"; mkdir -p "$W/out"
cp -r "$REG/Diamond_NV_natab_cce/inputfiles" "$W/"
cd "$W"
python3 - <<'PY'
import json, re
for name in ("ccein_nstate0.json", "ccein_nstate1.json"):
    j = json.loads(re.sub(r"(?m)^\s*#.*$", "", open("inputfiles/" + name).read()))
    tag = name.split("_")[1].split(".")[0]
    json.dump(dict(j, outfile=f"./out/{tag}_both"), open(f"{tag}_both.json", "w"), indent=1)
    json.dump(dict(j, outfile=f"./out/{tag}_wi", save_nodiv=False), open(f"{tag}_wi.json", "w"), indent=1)
    json.dump(dict(j, outfile=f"./out/{tag}_bad", save_nodiv="false"), open(f"{tag}_bad.json", "w"), indent=1)
PY
for tag in nstate0 nstate1; do
    run_ccex "$BIN_DEV" ${tag}_both.json ${tag}_both.log
    run_ccex "$BIN_DEV" ${tag}_wi.json ${tag}_wi.log
    ok=1
    for f in out/${tag}_both*_wiDiv; do
        g=${f/_both/_wi}
        cmp -s "$f" "$g" || { echo "    DIFF $f $g"; ok=0; }
        [ -f "${f%_wiDiv}_noDiv" ] || { echo "    the default run wrote no ${f%_wiDiv}_noDiv"; ok=0; }
    done
    ls out/${tag}_wi*_noDiv >/dev/null 2>&1 && { echo "    save_nodiv false still wrote a _noDiv file"; ok=0; }
    ls out/${tag}_both*_wiDiv >/dev/null 2>&1 || ok=0
    if [ $ok -eq 1 ]; then echo "  PASS  N1 $tag: save_nodiv false gives the default _wiDiv and no _noDiv"; PASS=$((PASS+1))
    else echo "  FAIL  N1 $tag"; FAIL=$((FAIL+1)); fi
    if ! run_ccex "$BIN_DEV" ${tag}_bad.json ${tag}_bad.log && grep -q "save_nodiv must be a boolean" ${tag}_bad.log; then
        echo "  PASS  N2 $tag: save_nodiv \"false\" (a string) is refused"; PASS=$((PASS+1))
    else echo "  FAIL  N2 $tag (see $W/${tag}_bad.log)"; FAIL=$((FAIL+1)); fi
done
echo "nodiv: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
