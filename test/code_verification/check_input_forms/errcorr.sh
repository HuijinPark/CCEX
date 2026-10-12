#!/usr/bin/env bash
#
# Error correction on every system's baseline run.
# ---------------------------------------------------------------------------
# CCE divides a super-cluster coherence by its sub-cluster coherences. Where an
# ESEEM modulation takes a denominator close to zero the quotient diverges, and
# the _wiDiv file alone will happily show |L| > 1. CCEX writes _wiDiv and _noDiv
# together precisely so the two can be compared; only the corrected result is
# physical.
#
# The equivalence suite compares the RAW files on purpose -- correction and
# averaging both smooth, and a difference that is really there could hide under
# them. This script is the other half: it says whether the systems those raw
# comparisons were made on are numerically sound in the first place. A system
# that diverges is not disqualified as an equivalence fixture, but it must not
# be read as physics.
#
# Usage:  bash test/code_verification/check_input_forms/errcorr.sh
# Needs:  numpy < 2 (the analyzer uses np.ComplexWarning, removed in 2.0).
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

PY="${CCEX_ANALYZER_PY:-python3}"   # any python3 with numpy: the analyzer runs on numpy 2 since 2026-09-23
RUN="$REPO/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py"
export MPLBACKEND=Agg
PY="$(command -v "$PY")" || { echo "ERROR: analyzer python not found: ${CCEX_ANALYZER_PY:-python3}" >&2; exit 1; }
[ -f "$RUN" ] || { echo "ERROR: analyzer not found: $RUN" >&2; exit 1; }

# err_tol is the analyzer's own default. Named here rather than left implicit so
# that a later change to it is visible in the log.
ERR_TOL=2.0

# A stage that fails has to fail the script (FAIL is lib.sh's counter, as in axis2.sh):
# the analyzer not running, a baseline or a corrected file that is not there, the
# reading of the corrected file. Until 2026-09-30 none of these did - the last echo
# made the exit status 0, and run.sh printed "all stages passed" over an analyzer
# that had not run.
DONE=0
# error_correct on the boundary inputs that raised IndexError until 2026-10-02 (D25): no points, and a
# last point exactly at the tolerance after a divergence; plus a curve with nothing to correct.
if "$PY" - "$REPO/example/analysis/Analyzer/CoherenceAnalyzer.v1" > "$WORK/errcorr_unit.log" 2>&1 <<'PYEOF'
import sys, numpy as np
sys.path.insert(0, sys.argv[1])
from coherence import error_correct
assert error_correct(np.array([], complex), np.array([], complex)).size == 0
out = error_correct(np.array([1, 5, 0.9, 2.5], complex), np.array([1, 1, 1, 0.5], complex), tolerance=2.0)
assert np.allclose(out, [1, 0.95, 0.9, 2.499999]), out
same = np.array([1, 0.8, 0.6], complex)
assert np.array_equal(error_correct(same.copy(), same.copy()), same)
PYEOF
then echo "  ok    error_correct boundary inputs"
else echo "  FAIL  error_correct boundary inputs"; sed 's/^/        /' "$WORK/errcorr_unit.log"; FAIL=$((FAIL+1)); fi
for sys in S1a S1b S2a S2b S2c S3 S4; do
    w="$WORK/$sys"
    # prefer the axis-2 baseline where it exists, else the axis-1 dev run
    base=$(ls "$w"/base*_wiDiv "$w"/dev*_wiDiv 2>/dev/null | head -1)
    if [ -z "$base" ]; then
        # skipped only when no run was made for the system; one the axis-2 stage ran
        # (it left a manifest) must have its baseline
        if [ -f "$w/manifest.tsv" ]; then echo "  FAIL  $sys : no baseline output"; FAIL=$((FAIL+1))
        else                              echo "  skip  $sys (no baseline output)"; fi
        continue
    fi
    name=$(basename "${base%_wiDiv}")
    echo "=== $sys : $name"
    [ -f "$w/${name}_noDiv" ] || { echo "  FAIL  $sys : no ${name}_noDiv beside the baseline"; FAIL=$((FAIL+1)); echo; continue; }
    # an _ErrCorr left by an earlier run must not pass for this run's
    rm -f "$w/${name}_ErrCorr"
    "$PY" "$RUN" -d "$w/" -fi "${name}_wiDiv" \
                 -err_d "$w/" -err_fi "${name}_noDiv" \
                 -vn CONF -v 1 \
                 -uc ms -err -err_tol "$ERR_TOL" -ntn \
                 -fo "$w/${name}_ErrCorr" > "$w/proc_errcorr" 2>&1
    rc=$?
    grep -E "Read file number|Done to correct|Write file" "$w/proc_errcorr" | sed 's/^/    /'
    if [ $rc -ne 0 ]; then
        echo "  FAIL  $sys : analyzer exit=$rc"
        tail -3 "$w/proc_errcorr" | sed 's/^/        /'
        FAIL=$((FAIL+1)); echo; continue
    fi
    [ -s "$w/${name}_ErrCorr" ] || { echo "  FAIL  $sys : the analyzer wrote no ${name}_ErrCorr"; FAIL=$((FAIL+1)); echo; continue; }

    "$PY" - "$w/${name}_wiDiv" "$w/${name}_ErrCorr" <<'PYX'
import sys
import numpy as np

def load(p):
    out = []
    for line in open(p):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        out.append([complex(t) for t in line.split()])
    return np.array(out, dtype=complex)

raw = load(sys.argv[1])
cor = load(sys.argv[2])
# column 0 is time; the coherence lives in the remaining columns
r = np.abs(raw[:, 1:]); c = np.abs(cor[:, 1:])
print(f"    max|L| raw = {r.max():.6f}   corrected = {c.max():.6f}   "
      f"points with |L|>1 : {int((r > 1.0 + 1e-9).sum())} -> {int((c > 1.0 + 1e-9).sum())}")
PYX
    rc=$?
    if [ $rc -ne 0 ]; then echo "  FAIL  $sys : reading the corrected file (exit=$rc)"; FAIL=$((FAIL+1))
    else                   DONE=$((DONE+1)); fi
    echo
done

# nothing to correct is not a pass either: there was then nothing this stage looked at
if [ "$DONE" -eq 0 ] && [ "$FAIL" -eq 0 ]; then
    echo "  FAIL  no system had a baseline to correct (run legacy.sh and axis2.sh first)"; FAIL=$((FAIL+1))
fi
echo "==========================================="
echo " error correction : corrected=$DONE FAIL=$FAIL"
echo "==========================================="
[ "$FAIL" -eq 0 ]
