#!/usr/bin/env bash
# Shared helpers for the legacy <-> new input form equivalence suite (and ../check_refusals).
#
# Every check here is a COMPARISON, never an absolute expectation: no reference
# numbers are stored, so nothing can be "updated" into agreement by accident.

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
BIN_DEV="${CCEX_BIN:-$REPO/bin/main.out}"   # CCEX_BIN: the binary to test (a CMake build keeps it in its own folder)
# mpirun, as check_rotation finds it: CCEX_MPI when given (and then it has to exist), else the Intel MPI the
# loaded module points at, else whatever is on PATH, else the WSL location this suite was written against. It was
# that fixed WSL path alone until 2026-10-01, so the suite ran nowhere else.
if [ -n "${CCEX_MPI:-}" ]; then
    [ -x "$CCEX_MPI" ] || { echo "ERROR: CCEX_MPI is not an executable: $CCEX_MPI" >&2; exit 1; }
    MPIRUN="$CCEX_MPI"
else
    MPIRUN=""
    for __c in "${I_MPI_ROOT:-}/bin/mpirun" "$(command -v mpirun 2>/dev/null)" \
               /opt/intel/oneapi/mpi/2021.18/bin/mpirun; do
        [ -n "$__c" ] && [ -x "$__c" ] && { MPIRUN="$__c"; break; }
    done
    [ -n "$MPIRUN" ] || { echo "ERROR: no mpirun found (set CCEX_MPI)" >&2; exit 1; }
fi
echo "[check] mpirun: $MPIRUN; binary: $BIN_DEV" >&2
export I_MPI_FABRICS="${I_MPI_FABRICS:-shm}"

WORK="$REPO/test/code_verification/check_input_forms/_work"
# One run at a time per folder: a second run would delete the first one's work files (D29).
# A script this run starts (check_refusals runs legacy.sh) inherits the lock and does not take it again.
if [ "${CCEX_RUN_LOCK:-}" != "$REPO/test/code_verification/check_input_forms" ]; then
    exec 9>"$REPO/test/code_verification/check_input_forms/.run.lock"
    if command -v flock >/dev/null 2>&1 && ! flock -n 9; then
        echo "ERROR: another check is running in $REPO/test/code_verification/check_input_forms (it would delete that run's work files); wait for it or use another checkout" >&2; exit 1
    fi
    CCEX_RUN_LOCK="$REPO/test/code_verification/check_input_forms"; export CCEX_RUN_LOCK
fi
FIX="$REPO/test/code_verification/bath_data"
REG="$REPO/test/code_verification/check_input_forms/inputs"    # the four legacy-form inputs of the old regression suite (its v1.0.0 layout), the base of every variant here
EX="$REPO/example"

# run_ccex <binary> <ccein> <procfile> [extra args...]
# 1 MPI rank throughout: this is an equivalence check, not a performance test,
# and gather order is one less thing to have to argue about.
run_ccex(){
    local bin="$1" ccein="$2" proc="$3"; shift 3
    local t0 t1 rc
    t0=$(date +%s.%N)
    # stdin from /dev/null: mpirun drains it, and callers drive these runs from a
    # `while read` over a manifest, which would lose every line but the first.
    "$MPIRUN" -n 1 "$bin" -f "$ccein" "$@" > "$proc" 2>&1 < /dev/null
    rc=$?
    t1=$(date +%s.%N)
    CCEX_WALL=$(awk "BEGIN{printf \"%.2f\", $t1-$t0}")
    return $rc
}

# cmp_raw <label> <prefixA> <prefixB>
# Byte comparison of the raw CCEX coherence output. Deliberately NOT the
# error-corrected or ensemble-averaged file: those divide and threshold, which
# can hide a difference that is really there.
PASS=0; FAIL=0
cmp_raw(){
    local label="$1" a="$2" b="$3" ok=1 f
    for f in _noDiv _wiDiv; do
        if [ ! -f "$a$f" ]; then echo "    MISSING $a$f"; ok=0; continue; fi
        if [ ! -f "$b$f" ]; then echo "    MISSING $b$f"; ok=0; continue; fi
        if ! cmp -s "$a$f" "$b$f"; then
            echo "    DIFF $f"
            # first differing line, so a failure says what changed
            diff "$a$f" "$b$f" | head -4 | sed 's/^/      /'
            ok=0
        fi
    done
    if [ $ok -eq 1 ]; then echo "  PASS  $label"; PASS=$((PASS+1))
    else                   echo "  FAIL  $label"; FAIL=$((FAIL+1)); fi
}

# maxdiff <label> <prefixA> <prefixB>
# For pairs that reach the same numbers by different floating-point routes, where
# byte identity cannot be required. Reports the largest absolute deviation over
# every column; the caller decides what to make of it. No threshold is baked in.
# What is counted as a failure is only that the comparison did not happen.
maxdiff(){
    local label="$1" a="$2" b="$3" f
    for f in _noDiv _wiDiv; do
        [ -f "$a$f" ] && [ -f "$b$f" ] || { echo "  FAIL  $label ($f missing)"; FAIL=$((FAIL+1)); return; }
        python3 - "$a$f" "$b$f" "$label$f" <<'PYX'
import sys

def load(path):
    # CCEX writes complex columns as "+1.0000000000+0.0000000000j"; the time column
    # is a plain float. quantity="dm" rows are not all the same width, so keep the
    # rows ragged instead of forcing a rectangular array.
    rows = []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        rows.append([complex(tok) for tok in line.split()])
    return rows

A = load(sys.argv[1]); B = load(sys.argv[2])
if len(A) != len(B) or any(len(x) != len(y) for x, y in zip(A, B)):
    print(f"  ----  {sys.argv[3]}: layout differs ({len(A)} vs {len(B)} rows)"); sys.exit(1)
m = 0.0; scale = 0.0; nan_one = nan_both = 0
for ra, rb in zip(A, B):
    for va, vb in zip(ra, rb):
        na, nb = va != va, vb != vb          # max() skips NaN, so a NaN looked like no difference (D42)
        if na or nb:
            nan_both += na and nb; nan_one += na != nb
            continue
        m = max(m, abs(va - vb)); scale = max(scale, abs(va))
rel = m / scale if scale > 0 else m
tag = "byte-equal" if m == 0.0 and nan_one == nan_both == 0 else f"max|d|={m:.3e}  rel={rel:.3e}"
if nan_one: tag += f"  NaN in one file only at {nan_one} values"
if nan_both: tag += f"  NaN in both at {nan_both} values"
print(f"  ----  {sys.argv[3]}: {tag}")
PYX
        # a deviation is reported and read, not judged; a comparison that could not be
        # made (the script failed, or the two files are not laid out alike) is a failure
        [ $? -eq 0 ] || { echo "  FAIL  $label ($f not compared)"; FAIL=$((FAIL+1)); return; }
    done
}
