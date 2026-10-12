#!/usr/bin/env bash
#
# CCEX against a closed form: one 13C beside the NV, order 1 (exact for one spin),
# ensemble average, free induction and Hahn echo (analytic.py has the formulas).
# Everything else under code_verification compares CCEX with CCEX (another input
# form, another frame, a recorded run); this one does not, so a common error such as
# a wrong time scale in the free evolution is caught here (D22, 2026-10-02).
# A self-test runs the echo with a wrong gyro and requires the check to fail.
#
# Usage:  bash test/code_verification/check_analytic/run.sh      (CCEX_BIN, CCEX_MPI as for the others)
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../check_input_forms/lib.sh"
W="$HERE/_work"; rm -rf "${W:?}"; mkdir -p "$W"
cp "$HERE"/inputs/{qubit,bath,gyro} "$W/"
NSTEP=400; DELTAT=0.0005
# qubitfile form: the branch keys are qalphams/qbetams (alphams/betams belong to the Qubit section)
input(){ # input <npulse> <gyro file>
    cat <<JSON
{"method": "CCE", "quantity": "coherence", "qubitfile": "qubit", "gyrofile": "$2", "bathfile": ["bath"],
 "order": 1, "bfield": 500, "rbath": 20, "rdip": 10, "deltat": $DELTAT, "nstep": $NSTEP,
 "qalphams": 1, "qbetams": 0, "npulse": $1, "nstate": 0, "seed": 1}
JSON
}
PASS=0; FAIL=0
for np in 0 1; do
    input $np gyro > "$W/in$np.json"
    ( cd "$W" && run_ccex "$BIN_DEV" "in$np.json" "proc$np" -o "out$np" )
    if [ $? -ne 0 ]; then echo "  FAIL  npulse $np : CCEX exit"; FAIL=$((FAIL+1)); continue; fi
    if python3 "$HERE/analytic.py" "$W/out$np" --npulse $np --nstep $NSTEP --deltat $DELTAT; then PASS=$((PASS+1)); else echo "  FAIL  npulse $np"; FAIL=$((FAIL+1)); fi
done
# self-test: a 1 % wrong gyro has to be caught
printf '13C 0.5 6.795567\n' > "$W/gyro_wrong"
input 1 gyro_wrong > "$W/in_wrong.json"
# only a numerical mismatch (exit 1) counts: a run that failed, or an output that could not be read
# (exit 2), is a broken self-test, not a caught error
( cd "$W" && run_ccex "$BIN_DEV" "in_wrong.json" "proc_wrong" -o "out_wrong" )
rc=$?
python3 "$HERE/analytic.py" "$W/out_wrong" --npulse 1 --nstep $NSTEP --deltat $DELTAT > "$W/selftest.log"
ac=$?
if [ $rc -ne 0 ]; then echo "  FAIL  self-test : CCEX exit $rc"; FAIL=$((FAIL+1))
elif [ $ac -eq 1 ]; then echo "  PASS  self-test : a 1 % wrong gyro is caught ($(grep -o 'max .*= [0-9.e+-]*' "$W/selftest.log"))"; PASS=$((PASS+1))
elif [ $ac -eq 0 ]; then echo "  FAIL  self-test : a 1 % wrong gyro passed"; FAIL=$((FAIL+1))
else echo "  FAIL  self-test : $(cat "$W/selftest.log")"; FAIL=$((FAIL+1)); fi
echo " closed-form checks : PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
