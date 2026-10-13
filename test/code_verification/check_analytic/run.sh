#!/usr/bin/env bash
#
# CCEX against a closed form: one 13C beside the NV, order 1 (exact for one spin),
# ensemble average, free induction and Hahn echo (analytic.py has the formulas).
# Everything else under code_verification compares CCEX with CCEX (another input
# form, another frame, a recorded run); this one does not, so a common error such as
# a wrong time scale in the free evolution is caught here (D22, 2026-10-02).
# A self-test runs the echo with a wrong gyro and requires the check to fail.
# hfmedi runs at zero field (degenerate qubit levels) close it, see below.
#
# Usage:  bash test/code_verification/check_analytic/run.sh      (CCEX_BIN, CCEX_MPI as for the others)
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../check_input_forms/lib.sh"
W="$HERE/_work"; rm -rf "${W:?}"; mkdir -p "$W"
cp "$HERE"/inputs/{qubit,bath,gyro} "$W/"
NSTEP=400; DELTAT=0.0005
# Qubit.file form: the branch is Qubit.alphams/betams beside the file
input(){ # input <npulse> <gyro file>
    cat <<JSON
{"method": "cCCE", "quantity": "coherence", "Qubit": {"file": "qubit", "alphams": 1, "betams": 0}, "gyrofile": "$2",
 "bathfile": ["bath"], "order": 1, "bfield": 500, "rbath": 20, "rdip": 10, "deltat": $DELTAT, "nstep": $NSTEP,
 "pulse": {"count": $1}, "nstate": 0, "seed": 1}
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
# hfmedi with degenerate qubit levels (review 20261012-194049-af8e42 R1), not a closed form: at zero field an axial ZFS
# leaves +1 and -1 degenerate but uncoupled (the term of that pair is 0: finite output, exit 0, as it was NaN); without
# ZFS the levels 0 and +1 are degenerate and coupled, where the second-order term does not exist (refused). Review
# 20261012-233058-452278: P1-02, the same degenerate levels with a zero hyperfine (gyro 0) have no coupling and are
# accepted; P1-03, a ZFS shifted by 1e13 I (an energy shift) is decided as the unshifted one.
hfmedi(){ # hfmedi <name> <Qubit.zfs> [gyro file]
    cat > "$W/$1.json" <<JSON
{"method": "cCCE", "quantity": "coherence", "Qubit": {"file": "qubit", "alphams": 1, "betams": 0, "zfs": $2},
 "gyrofile": "${3:-gyro}", "bathfile": ["bath"], "order": 1, "bfield": 0, "rbath": 20, "rdip": 10, "deltat": 0.001,
 "nstep": 3, "hamiltonian": {"cCCE": {"hfmedi": true}}, "pulse": {"count": 1}, "nstate": 0, "seed": 1}
JSON
    ( cd "$W" && run_ccex "$BIN_DEV" "$1.json" "proc_$1" -o "out_$1" )
}
finite(){ # finite <name>: both outputs there, 3 points each, all finite
    python3 -c '
import math, sys
v = [complex(x) for f in sys.argv[1:] for l in open(f) for x in l.split()[1:]]
sys.exit(0 if len(v) == 6 and all(math.isfinite(z.real) and math.isfinite(z.imag) for z in v) else 1)' "$W/out_$1_noDiv" "$W/out_$1_wiDiv"
}
accepted(){ # accepted <name> <label> <exit>
    if [ "$3" -ne 0 ]; then echo "  FAIL  $2 : CCEX exit $3"; FAIL=$((FAIL+1))
    elif finite "$1"; then echo "  PASS  $2 : finite"; PASS=$((PASS+1))
    else echo "  FAIL  $2 : output missing or not finite"; FAIL=$((FAIL+1)); fi
}
hfmedi hfmedi_axial '[[-956666.6666666666,0,0],[0,-956666.6666666666,0],[0,0,1913333.3333333333]]'
accepted hfmedi_axial "hfmedi B=0 axial ZFS" $?
hfmedi hfmedi_nozfs '[[0,0,0],[0,0,0],[0,0,0]]'
rc=$?
if [ $rc -ne 0 ] && grep -q "hfmedi: qubit level .* is degenerate" "$W/proc_hfmedi_nozfs"; then
    echo "  PASS  hfmedi B=0 no ZFS : refused (exit=$rc)"; PASS=$((PASS+1))
else echo "  FAIL  hfmedi B=0 no ZFS : exit $rc, expected the degeneracy refusal"; FAIL=$((FAIL+1)); fi
printf '13C 0.5 0\n' > "$W/gyro_zero"
hfmedi hfmedi_zero_hf '[[0,0,0],[0,0,0],[0,0,0]]' gyro_zero
accepted hfmedi_zero_hf "hfmedi B=0 no ZFS, zero hyperfine" $?
hfmedi hfmedi_small '[[0,0,0],[0,0,0],[0,0,1]]'
rs=$?
hfmedi hfmedi_shift '[[1e13,0,0],[0,1e13,0],[0,0,10000000000001.0]]'
rh=$?
if [ $rs -eq 0 ] && [ $rh -eq 0 ] && finite hfmedi_small && finite hfmedi_shift; then
    echo "  PASS  hfmedi B=0 ZFS diag(0,0,1) and + 1e13 I : both accepted, finite"; PASS=$((PASS+1))
else echo "  FAIL  hfmedi B=0 ZFS diag(0,0,1) exit $rs, + 1e13 I exit $rh : expected the same decision (accepted, finite)"; FAIL=$((FAIL+1)); fi
# hamiltonian.cCCE.delta_A false (A47): the hyperfine tensors of a cluster's bath spins are replaced by their mean, so their
# Knight fields no longer differ (Zhao 2011 SI eq. 23, Maze 2008 eq. 8). Two 13C on the NV axis, 1.54 A apart at 14.5 A
# (inputs/bath_pair): their hyperfine tensors are axial (no x z, y z part), so one spin alone gives L = 1 in a Hahn echo and
# the pair decays only through the flip-flop the Knight-field difference (26 rad/ms, the size of the pair's coupling) lets
# through. delta_A false makes the cluster's two tensors equal: the flip-flop commutes with what is left and the pair gives
# L = 1 exactly (to the 10 decimals written); with delta_A true (the default) the same pair decays. Order 1 (a cluster of
# one spin keeps its own tensor) is byte for byte the run without the key.
pair(){ # pair <name> <order> <delta_A json value or "">
    local h=""; [ -n "$3" ] && h=", \"hamiltonian\": {\"cCCE\": {\"delta_A\": $3}}"
    cat > "$W/$1.json" <<JSON
{"method": "cCCE", "quantity": "coherence", "Qubit": {"file": "qubit", "alphams": 1, "betams": 0}, "gyrofile": "gyro",
 "bathfile": ["bath_pair"], "order": $2, "bfield": 500, "rbath": 20, "rdip": 10, "deltat": 0.005, "nstep": 400,
 "pulse": {"count": 1}, "nstate": 0, "seed": 1$h}
JSON
    ( cd "$W" && run_ccex "$BIN_DEV" "$1.json" "proc_$1" -o "out_$1" )
}
cp "$HERE"/inputs/bath_pair "$W/"
farthest(){ # farthest <name>: the largest |L - 1| of the noDiv output
    python3 -c '
import sys
print(max(abs(complex(l.split()[1].replace("+-", "-")) - 1) for l in open(sys.argv[1]) if l.strip()))' "$W/out_$1_noDiv"
}
pair dA_true 2 ""; r1=$?
pair dA_false 2 false; r2=$?
pair dA_true_o1 1 ""; r3=$?
pair dA_false_o1 1 false; r4=$?
if [ $r1 -ne 0 ] || [ $r2 -ne 0 ] || [ $r3 -ne 0 ] || [ $r4 -ne 0 ]; then
    echo "  FAIL  delta_A : CCEX exit ($r1 $r2 $r3 $r4)"; FAIL=$((FAIL+1))
else
    dt=$(farthest dA_true); df=$(farthest dA_false)
    if python3 -c "import sys; sys.exit(0 if float('$dt') > 1e-2 and float('$df') < 1e-9 else 1)"; then
        echo "  PASS  delta_A : the 13C pair on the axis decays with delta_A true (max |L-1| $dt), L = 1 with false (max |L-1| $df)"
        PASS=$((PASS+1))
    else echo "  FAIL  delta_A : max |L-1| $dt with delta_A true (expected > 1e-2), $df with false (expected < 1e-9)"; FAIL=$((FAIL+1)); fi
    if cmp -s "$W/out_dA_true_o1_noDiv" "$W/out_dA_false_o1_noDiv" && cmp -s "$W/out_dA_true_o1_wiDiv" "$W/out_dA_false_o1_wiDiv"; then
        echo "  PASS  delta_A : order 1 (clusters of one spin) is byte for byte the same with delta_A false"; PASS=$((PASS+1))
    else echo "  FAIL  delta_A : order 1 changes with delta_A false"; FAIL=$((FAIL+1)); fi
fi
# hamiltonian.common.bath_bath_mean_field "off" (A47): with nstate > 0 the static field of the bath spins outside a cluster
# is left out. Two electron-like spins (gyro e, inputs/bath_off, a fixed state each) off the NV axis, order 1: with "off"
# the coherence is the product of the two spins each alone in its bath (no other spin, so no field), point for point; with
# "z" (the default) the other spin's field moves it.
offrun(){ # offrun <name> <bath> <state file prefix> <bath_bath_mean_field or "">
    local h=""; [ -n "$4" ] && h=", \"hamiltonian\": {\"common\": {\"bath_bath_mean_field\": \"$4\"}}"
    cat > "$W/$1.json" <<JSON
{"method": "cCCE", "quantity": "coherence", "Qubit": {"file": "qubit", "alphams": 1, "betams": 0}, "gyrofile": "gyro_e",
 "bathfile": ["$2"], "order": 1, "bfield": 500, "rbath": 20, "rdip": 10, "deltat": 0.00002, "nstep": 200,
 "pulse": {"count": 1}, "nstate": 1, "bath_state_file": "$3", "seed": 1$h}
JSON
    ( cd "$W" && run_ccex "$BIN_DEV" "$1.json" "proc_$1" -o "out_$1" )
}
cp "$HERE"/inputs/{gyro_e,bath_off,bath_off_a,bath_off_b,st_ab1,st_a1,st_b1} "$W/"
offrun mf_off bath_off st_ab off; r1=$?
offrun mf_z bath_off st_ab ""; r2=$?
offrun mf_a bath_off_a st_a ""; r3=$?
offrun mf_b bath_off_b st_b ""; r4=$?
if [ $r1 -ne 0 ] || [ $r2 -ne 0 ] || [ $r3 -ne 0 ] || [ $r4 -ne 0 ]; then
    echo "  FAIL  bath_bath_mean_field off : CCEX exit ($r1 $r2 $r3 $r4)"; FAIL=$((FAIL+1))
else
    read -r doff dz < <(python3 -c '
import sys
def read(p): return [complex(l.split()[1].replace("+-", "-")) for l in open(p) if l.strip()]
off, z, a, b = (read(f"{sys.argv[1]}/out_{n}_state1_noDiv") for n in ("mf_off", "mf_z", "mf_a", "mf_b"))
prod = [x*y for x, y in zip(a, b)]
print(max(abs(x-y) for x, y in zip(off, prod)), max(abs(x-y) for x, y in zip(z, prod)))' "$W")
    if python3 -c "import sys; sys.exit(0 if float('$doff') < 1e-9 and float('$dz') > 1e-6 else 1)"; then
        echo "  PASS  bath_bath_mean_field off : the product of the spins alone (max difference $doff); z differs by $dz"
        PASS=$((PASS+1))
    else echo "  FAIL  bath_bath_mean_field off : off - product $doff (expected < 1e-9), z - product $dz (expected > 1e-6)"
         FAIL=$((FAIL+1)); fi
fi
echo " closed-form checks : PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
