#!/usr/bin/env bash
#
# auto_convergence.enabled (example/tools/AutoConvergence/README.md, 2026-10-08): main.out -f input.json runs the automatic
# convergence when the input holds an auto_convergence block whose "enabled" is true or absent; "enabled": false runs the
# ordinary calculation, byte-identical to the input without the block; under mpirun -n N rank 0 alone starts the runner;
# "enabled" that is not a JSON boolean is refused. natab cCCE2 input of inputs/Diamond_NV_natab_cce, one sampled state.
#
# Usage:  bash test/code_verification/check_input_forms/autoconv.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/autoconv"
rm -rf "$W"; mkdir -p "$W/out"
cp -r "$REG/Diamond_NV_natab_cce/inputfiles" "$W/"
cd "$W"
python3 - <<'PY'
import json, re
j = json.loads(re.sub(r"(?m)^\s*#.*$", "", open("inputfiles/ccein_nstate1.json").read()))
def out(name):                     # output.file (A65; outfile until then)
    return dict(j.get("output", {}), file=f"./out/{name}")
json.dump(dict(j, output=out("plain")), open("plain.json", "w"), indent=1)
json.dump(dict(j, output=out("off"), auto_convergence=dict(enabled=False, workdir="auto")), open("off.json", "w"), indent=1)
json.dump(dict(j, output=out("on"), auto_convergence=dict(workdir="auto")), open("on.json", "w"), indent=1)
json.dump(dict(j, output=out("bad"), auto_convergence=dict(enabled="yes", workdir="auto")), open("bad.json", "w"), indent=1)
PY
run_ccex "$BIN_DEV" plain.json plain.log
run_ccex "$BIN_DEV" off.json off.log
ok=1
for f in out/plain*_wiDiv out/plain*_noDiv; do
    [ -f "$f" ] || { ok=0; continue; }
    cmp -s "$f" "${f/plain/off}" || { echo "    DIFF $f ${f/plain/off}"; ok=0; }
done
ls out/plain*_wiDiv >/dev/null 2>&1 || ok=0
if [ $ok -eq 1 ]; then echo "  PASS  A1 enabled false: the ordinary calculation, the same bytes as without the block"; PASS=$((PASS+1))
else echo "  FAIL  A1 (see $W/off.log)"; FAIL=$((FAIL+1)); fi
# plan "mpirun으로 부르면 그 rank 수로" (2026-10-10): under mpirun -n 2 the runner starts once (rank 0), not refused
"$MPIRUN" -n 2 "$BIN_DEV" -f on.json > on_mpi.log 2>&1 < /dev/null
# (the input's bath file has no {k}: the runner starts and refuses it - A65, the configurations are the configuration block's)
if [ "$(grep -c "^\[auto\] runner" on_mpi.log)" = 1 ] && grep -q "configuration.bathfile must hold {k} once" on_mpi.log && ! grep -q "refused inside an MPI launch" on_mpi.log; then
    echo "  PASS  A2 an enabled block under mpirun -n 2: the runner starts once, from rank 0"; PASS=$((PASS+1))
else echo "  FAIL  A2 (see $W/on_mpi.log)"; FAIL=$((FAIL+1)); fi
if ! "$BIN_DEV" -f on.json > on.log 2>&1 < /dev/null && grep -q "^\[auto\] runner" on.log; then
    echo "  PASS  A3 main.out -f with an enabled block starts the runner"; PASS=$((PASS+1))
else echo "  FAIL  A3 (see $W/on.log)"; FAIL=$((FAIL+1)); fi
if ! "$BIN_DEV" -f bad.json > bad.log 2>&1 < /dev/null && grep -q "enabled must be true or false" bad.log; then
    echo "  PASS  A4 enabled \"yes\" (not a boolean) is refused"; PASS=$((PASS+1))
else echo "  FAIL  A4 (see $W/bad.log)"; FAIL=$((FAIL+1)); fi
# review 20261008-autoconv-autostart #1: another option would be dropped by the runner, so it is refused; -fFILE is -f FILE
# (-P until A63, 2026-10-15; -v is one of the options main.out keeps)
if ! "$BIN_DEV" -v -f on.json > opt.log 2>&1 < /dev/null && grep -q "takes only -f FILE" opt.log && ! grep -q "^\[auto\] runner" opt.log; then
    echo "  PASS  A5 -v with an enabled block is refused (the runner would drop it)"; PASS=$((PASS+1))
else echo "  FAIL  A5 (see $W/opt.log)"; FAIL=$((FAIL+1)); fi
if ! "$BIN_DEV" -fon.json > attached.log 2>&1 < /dev/null && grep -q "^\[auto\] runner" attached.log; then
    echo "  PASS  A6 -fFILE is read like -f FILE"; PASS=$((PASS+1))
else echo "  FAIL  A6 (see $W/attached.log)"; FAIL=$((FAIL+1)); fi
# plan "끊김 대비·진행 보기" R4: the status reads the workdir only (nothing has run here: not started); the tool's own command
# since A63 (main.out --auto_status was removed, 2026-10-15)
if PYTHONPATH="$REPO/example/tools/AutoConvergence" "${CCEX_PYTHON:-python3}" -m autoconv status on.json > status.log 2>&1 < /dev/null && grep -q "not started" status.log; then
    echo "  PASS  A7 python -m autoconv status tells a run that has not started"; PASS=$((PASS+1))
else echo "  FAIL  A7 (see $W/status.log)"; FAIL=$((FAIL+1)); fi
# 2026-10-09: a shared installation is allowed (the tool owned by this user, root or the binary's owner) but nothing in it may be
# writable by others - a .py file made writable by others in a copy of the tool is refused before Python starts
T="$W/tool"; cp -r "$REPO/example/tools/AutoConvergence" "$T"; rm -rf "$T/autoconv/__pycache__"; chmod -R go-w "$T"
chmod o+w "$T/autoconv/t2.py"
if ! CCEX_AUTOCONV_DIR="$T" CCEX_PYTHON="${CCEX_PYTHON:-python3}" "$BIN_DEV" -f on.json > perm.log 2>&1 < /dev/null \
   && grep -q "refuses $T/autoconv/t2.py" perm.log; then
    echo "  PASS  A8 a .py file of the tool writable by others is refused"; PASS=$((PASS+1))
else echo "  FAIL  A8 (see $W/perm.log)"; FAIL=$((FAIL+1)); fi
# a folder above the tool writable by others could have another tool put in its place after the check (TOCTOU, 2026-10-09)
O="$W/open"; mkdir -p "$O"; cp -r "$REPO/example/tools/AutoConvergence" "$O/tool"; rm -rf "$O/tool/autoconv/__pycache__"
chmod -R go-w "$O/tool"; chmod 777 "$O"
if ! CCEX_AUTOCONV_DIR="$O/tool" CCEX_PYTHON="${CCEX_PYTHON:-python3}" "$BIN_DEV" -f on.json > ancestor.log 2>&1 < /dev/null \
   && grep -q "refuses $O:" ancestor.log; then
    echo "  PASS  A9 a folder above the tool writable by others is refused"; PASS=$((PASS+1))
else echo "  FAIL  A9 (see $W/ancestor.log)"; FAIL=$((FAIL+1)); fi
chmod 755 "$O"
# the sticky bit excuses a folder above the tool (/tmp), never the tool itself (commit security review 2026-10-09)
chmod 1777 "$O/tool/autoconv"
if ! CCEX_AUTOCONV_DIR="$O/tool" CCEX_PYTHON="${CCEX_PYTHON:-python3}" "$BIN_DEV" -f on.json > sticky.log 2>&1 < /dev/null \
   && grep -q "refuses $O/tool/autoconv:" sticky.log; then
    echo "  PASS  A11 a sticky tool folder writable by others is refused"; PASS=$((PASS+1))
else echo "  FAIL  A11 (see $W/sticky.log)"; FAIL=$((FAIL+1)); fi
chmod 755 "$O/tool/autoconv"
# a symbolic link inside the tool points where the folders above are not checked (review 6548a9 R2)
L="$W/linked"; mkdir -p "$L/tool"; cp -r "$REPO/example/tools/AutoConvergence/autoconv" "$L/pkg"; rm -rf "$L/pkg/__pycache__"
ln -s "$L/pkg" "$L/tool/autoconv"; chmod -R go-w "$L"
if ! CCEX_AUTOCONV_DIR="$L/tool" CCEX_PYTHON="${CCEX_PYTHON:-python3}" "$BIN_DEV" -f on.json > link.log 2>&1 < /dev/null \
   && grep -q "symbolic link inside the tool" link.log; then
    echo "  PASS  A12 a symbolic link inside the tool is refused"; PASS=$((PASS+1))
else echo "  FAIL  A12 (see $W/link.log)"; FAIL=$((FAIL+1)); fi
# group-writable is allowed when the group is the owner's own private group (idm copies under umask 0002); checked only where
# this user's primary group has no other members; without getent (macOS, where the primary group is staff, every user's) it
# cannot tell and skips
G="$(id -gn)"; GID="$(id -g)"
if command -v getent >/dev/null 2>&1 && [ -z "$(getent group "$G" | cut -d: -f4)" ] && [ "$(getent passwd | awk -F: -v g="$GID" -v u="$(id -u)" '$4==g && $3!=u' | wc -l)" -eq 0 ]; then
    P="$W/private"; mkdir -p "$P"; cp -r "$REPO/example/tools/AutoConvergence" "$P/tool"; rm -rf "$P/tool/autoconv/__pycache__"
    chmod -R go-w "$P/tool"; chmod -R g+w "$P/tool"
    CCEX_AUTOCONV_DIR="$P/tool" CCEX_PYTHON="${CCEX_PYTHON:-python3}" "$BIN_DEV" -f on.json > private.log 2>&1 < /dev/null
    if grep -q "\[auto\] runner $P/tool/autoconv" private.log && ! grep -q "refuses" private.log; then
        echo "  PASS  A10 a tool writable by the owner's private group is accepted"; PASS=$((PASS+1))
    else echo "  FAIL  A10 (see $W/private.log)"; FAIL=$((FAIL+1)); fi
else echo "  skip  A10 (the primary group $G has other members or is another account's primary group)"; fi
echo "autoconv: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
