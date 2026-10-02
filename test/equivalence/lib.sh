#!/usr/bin/env bash
# Shared helpers for the v1.1.0 <-> axis-align equivalence suite.
#
# Every check here is a COMPARISON, never an absolute expectation: no reference
# numbers are stored, so nothing can be "updated" into agreement by accident.

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN_DEV="$REPO/bin/main.out"
BIN_V110="/home/jh0626/SourceCodes/CCEX-v1.1.0/bin/main.out.v110"
MPIRUN="${CCEX_MPI:-/opt/intel/oneapi/mpi/2021.18/bin/mpirun}"
export I_MPI_FABRICS="${I_MPI_FABRICS:-shm}"

WORK="$REPO/test/equivalence/_work"
FIX="$REPO/test/code_verification/CCE_Reprod/Bath_Data"
REG="$REPO/test/regression/v1.0.0"
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
    print(f"  ----  {sys.argv[3]}: layout differs ({len(A)} vs {len(B)} rows)"); sys.exit(0)
m = 0.0; scale = 0.0
for ra, rb in zip(A, B):
    for va, vb in zip(ra, rb):
        m = max(m, abs(va - vb)); scale = max(scale, abs(va))
rel = m / scale if scale > 0 else m
tag = "byte-equal" if m == 0.0 else f"max|d|={m:.3e}  rel={rel:.3e}"
print(f"  ----  {sys.argv[3]}: {tag}")
PYX
    done
}
