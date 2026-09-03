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
# Usage:  bash test/equivalence/errcorr.sh
# Needs:  numpy < 2 (the analyzer uses np.ComplexWarning, removed in 2.0).
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

PY="${CCEX_ANALYZER_PY:-/home/jh0626/libs/venv_ccex_analyzer/bin/python}"
RUN="$REPO/example/z.Analyzer/CoherenceAnalyzer.v1/run.py"
export MPLBACKEND=Agg
[ -x "$PY"  ] || { echo "ERROR: analyzer python (numpy<2) not found: $PY" >&2; exit 1; }
[ -f "$RUN" ] || { echo "ERROR: analyzer not found: $RUN" >&2; exit 1; }

# err_tol is the analyzer's own default. Named here rather than left implicit so
# that a later change to it is visible in the log.
ERR_TOL=2.0

for sys in S1a S1b S2a S2b S2c S3 S4; do
    w="$WORK/$sys"
    # prefer the axis-2 baseline where it exists, else the axis-1 dev run
    base=$(ls "$w"/base*_wiDiv "$w"/dev*_wiDiv 2>/dev/null | head -1)
    [ -n "$base" ] || { echo "  skip  $sys (no baseline output)"; continue; }
    name=$(basename "${base%_wiDiv}")
    echo "=== $sys : $name"
    "$PY" "$RUN" -d "$w/" -fi "${name}_wiDiv" \
                 -err_d "$w/" -err_fi "${name}_noDiv" \
                 -vn CONF -v 1 \
                 -uc ms -err -err_tol "$ERR_TOL" -ntn \
                 -fo "$w/${name}_ErrCorr" 2>&1 \
        | grep -E "Read file number|Done to correct|Write file" | sed 's/^/    /'

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
try:
    cor = load(sys.argv[2])
except OSError:
    print("    ErrCorr file not written"); sys.exit(0)
# column 0 is time; the coherence lives in the remaining columns
r = np.abs(raw[:, 1:]); c = np.abs(cor[:, 1:])
print(f"    max|L| raw = {r.max():.6f}   corrected = {c.max():.6f}   "
      f"points with |L|>1 : {int((r > 1.0 + 1e-9).sum())} -> {int((c > 1.0 + 1e-9).sum())}")
PYX
    echo
done
