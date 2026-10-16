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
# the qubit of a position file, a Qubit.qubit entry (A94; the Qubit.file form before): the branch is Qubit.alpha_state/beta_state {"ms": [m]}
input(){ # input <npulse> <gyro file> <outfile>
    cat <<JSON
{"method": "cCCE", "quantity": "coherence", "Qubit": {"qubit": [{"file": "qubit"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}},
 "configuration": {"gyrofile": "$2", "bathfile": ["bath"], "seed": 1},
 "order": 1, "bfield": 500, "rbath": 20, "rdip": 10, "deltat": $DELTAT, "nstep": $NSTEP,
 "pulse": {"count": $1}, "output": {"file": "$3"}}
JSON
}
PASS=0; FAIL=0
for np in 0 1; do
    input $np gyro "out$np" > "$W/in$np.json"
    ( cd "$W" && run_ccex "$BIN_DEV" "in$np.json" "proc$np" )
    if [ $? -ne 0 ]; then echo "  FAIL  npulse $np : CCEX exit"; FAIL=$((FAIL+1)); continue; fi
    if python3 "$HERE/analytic.py" "$W/out$np" --npulse $np --nstep $NSTEP --deltat $DELTAT; then PASS=$((PASS+1)); else echo "  FAIL  npulse $np"; FAIL=$((FAIL+1)); fi
done
# self-test: a 1 % wrong gyro has to be caught
printf '13C 0.5 6.795567\n' > "$W/gyro_wrong"
input 1 gyro_wrong out_wrong > "$W/in_wrong.json"
# only a numerical mismatch (exit 1) counts: a run that failed, or an output that could not be read
# (exit 2), is a broken self-test, not a caught error
( cd "$W" && run_ccex "$BIN_DEV" "in_wrong.json" "proc_wrong" )
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
hfmedi(){ # hfmedi <name> <ZFS tensor (kHz), the Qubit.intmap self entry> [gyro file]
    cat > "$W/$1.json" <<JSON
{"method": "cCCE", "quantity": "coherence", "Qubit": {"qubit": [{"file": "qubit"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}, "intmap": [{"between": ["q0", "q0"], "tensor": $2}]},
 "configuration": {"gyrofile": "${3:-gyro}", "bathfile": ["bath"], "seed": 1},
 "order": 1, "bfield": 0, "rbath": 20, "rdip": 10, "deltat": 0.001,
 "nstep": 3, "hamiltonian": {"cCCE": {"hfmedi": true}}, "pulse": {"count": 1}, "output": {"file": "out_$1"}}
JSON
    ( cd "$W" && run_ccex "$BIN_DEV" "$1.json" "proc_$1" )
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
# hamiltonian.cCCE.delta_A (A83): the second-order delta_A(m_s).I = -m_s/(2 Delta) (A_x x A_y).I of each bath spin. One 13C
# 2.55 A from the NV (inputs/bath_dA, |A| up to 9.4e3 rad/ms, |A_x x A_y|/(2 Delta) = 2.5 rad/ms), order 1, ensemble,
# 400 points to 1 ms: the cCCE against gCCE order 1, which keeps the whole S.A.I and is exact for one spin. Hahn echo:
# max |L - L_gCCE|; free induction: max ||L| - |L_gCCE|| (the second order also shifts the qubit frequency by a constant,
# which the cCCE leaves out, so the phase of L is not compared). B = 5 G along z (knight is 0 there): without delta_A
# 0.81 (echo) and 0.65 (free induction), with it 4.8e-3 and 5.3e-3. B = (2, 0, 5) G: knight alone 0.56 and 0.81, knight
# with delta_A 2.0e-2 and 2.9e-2. What is left is the qubit Zeeman energy b left out of the denominators (Delta, not
# Delta + m_s b, as for knight; b/Delta = 0.5 % at 5 G) and, across the axis, knight's own higher orders. The tolerances,
# 1e-2 and 5e-2, are about twice the measured residual; without delta_A the deviation has to be at least ten times the one with it
# (measured 27 to 170 times).
dArun(){ # dArun <name> <method> <npulse> <bfield> <hamiltonian.cCCE json or "">
    local h=""; [ -n "$5" ] && h=", \"hamiltonian\": {\"cCCE\": $5}"
    cat > "$W/$1.json" <<JSON
{"method": "$2", "quantity": "coherence",
 "Qubit": {"qubit": [{"file": "qubit"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}, "intmap": [{"between": ["q0", "q0"], "tensor": [[-956666.6666666666,0,0],[0,-956666.6666666666,0],[0,0,1913333.3333333333]]}]},
 "configuration": {"gyrofile": "gyro", "bathfile": ["bath_dA"], "seed": 1},
 "order": 1, "bfield": $4, "rbath": 20, "rdip": 10, "deltat": 0.0025, "nstep": 400,
 "pulse": {"count": $3}, "output": {"file": "out_$1"}$h}
JSON
    ( cd "$W" && run_ccex "$BIN_DEV" "$1.json" "proc_$1" )
}
cp "$HERE"/inputs/bath_dA "$W/"
dAcase(){ # dAcase <label> <bfield> <tolerance> <hamiltonian.cCCE without delta_A> <with delta_A>
    local rc=0 np
    for np in 0 1; do
        dArun "dA_${1}_g$np" gCCE $np "$2" "" || rc=1
        dArun "dA_${1}_off$np" cCCE $np "$2" "$4" || rc=1
        dArun "dA_${1}_on$np" cCCE $np "$2" "$5" || rc=1
    done
    if [ $rc -ne 0 ]; then echo "  FAIL  delta_A $1 : CCEX exit"; FAIL=$((FAIL+1)); return; fi
    for np in 0 1; do
        read -r doff don < <(python3 -c '
import sys
def read(p): return [complex(l.split()[1].replace("+-", "-")) for l in open(p) if l.strip()]
g, off, on = (read(f"{sys.argv[1]}/out_dA_{sys.argv[2]}_{n}{sys.argv[3]}_wiDiv") for n in ("g", "off", "on"))
d = (lambda x, y: abs(abs(x)-abs(y))) if sys.argv[3] == "0" else (lambda x, y: abs(x-y))
print(max(d(x, y) for x, y in zip(off, g)), max(d(x, y) for x, y in zip(on, g)))' "$W" "$1" $np)
        local what="Hahn echo max |L - L_gCCE|"; [ $np -eq 0 ] && what="free induction max ||L| - |L_gCCE||"
        if python3 -c "import sys; sys.exit(0 if float('$don') <= $3 and float('$doff') >= 10*float('$don') else 1)"; then
            echo "  PASS  delta_A $1 : $what $don with delta_A, $doff without"; PASS=$((PASS+1))
        else echo "  FAIL  delta_A $1 : $what $don with delta_A (expected <= $3), $doff without (expected >= 10 x the one with)"
             FAIL=$((FAIL+1)); fi
    done
}
dAcase "B=5G_z" "[0, 0, 5]" 1e-2 "" '{"delta_A": true}'
dAcase "B=(2,0,5)G" "[2, 0, 5]" 5e-2 '{"knight": true}' '{"knight": true, "delta_A": true}'
# hamiltonian.common.bath_bath_mean_field "off" (A47): with sampled states (s_range) the static field of the bath spins outside a cluster
# is left out. Two electron-like spins (gyro e, inputs/bath_off, a fixed state each) off the NV axis, order 1: with "off"
# the coherence is the product of the two spins each alone in its bath (no other spin, so no field), point for point; with
# "z" (the default) the other spin's field moves it.
offrun(){ # offrun <name> <bath> <state file prefix> <bath_bath_mean_field or "">
    local h=""; [ -n "$4" ] && h=", \"hamiltonian\": {\"common\": {\"bath_bath_mean_field\": \"$4\"}}"
    cat > "$W/$1.json" <<JSON
{"method": "cCCE", "quantity": "coherence", "Qubit": {"qubit": [{"file": "qubit"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}},
 "configuration": {"gyrofile": "gyro_e", "bathfile": ["$2"], "bath_state_file": ["$3"], "s_range": [1, 1], "seed": 1},
 "order": 1, "bfield": 500, "rbath": 20, "rdip": 10, "deltat": 0.00002, "nstep": 200,
 "pulse": {"count": 1}, "output": {"file": "out_$1"}$h}
JSON
    ( cd "$W" && run_ccex "$BIN_DEV" "$1.json" "proc_$1" )
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
# gCCE with a non-axial ZFS (E != 0) against the exact 6 x 6 answer (A91), the system as written and turned by 30 degrees
# about z (gcce_zfs_exact.py): the turned run keeps the m_s states, so its coherence differs (up to 3.6e-5), and CCEX has to
# follow the exact answer in both.
for phi in 0 30; do
    python3 "$HERE/gcce_zfs_exact.py" write $phi "$W/zfs$phi.json" "$W/bath_zfs$phi"
    ( cd "$W" && run_ccex "$BIN_DEV" "zfs$phi.json" "proc_zfs$phi" )
    if [ $? -ne 0 ]; then echo "  FAIL  gCCE ZFS E phi $phi : CCEX exit"; FAIL=$((FAIL+1)); continue; fi
    if r=$(python3 "$HERE/gcce_zfs_exact.py" check $phi "$W/zfs$phi"); then echo "  PASS  gCCE ZFS E phi $phi : $r"; PASS=$((PASS+1))
    else echo "  FAIL  gCCE ZFS E phi $phi : $r"; FAIL=$((FAIL+1)); fi
done
echo " closed-form checks : PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
