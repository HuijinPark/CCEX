#!/usr/bin/env bash
#
# A67: a key in another unit, {"value": ..., "unit": ...} (docs/user/InputFile.md#Units), gives the bytes of the plain
# number, on check_numerics cases. The values are chosen so that the one exact operation (x 10^k, or / 10^k) gives the
# double the parser gives for the plain number: an integer divided by a power of ten is the correctly rounded quotient,
# which is what strtod makes of the decimal, and the products below were checked to land on the plain number.
#   U1  natab_cce2_hahn: deltat 1 us, bfield 50 mT, rbath 3 nm, rdip 1 nm          == 0.001 ms, 500 G, 30 A, 10 A
#   U2  natab_cce2_hahn: every key in its own unit, written out ({"value": 0.001, "unit": "ms"}, ...)  == plain
#   U3  arbitrary_gcce1_dm: Qubit.qubit[0].xyz [0.15, 0, 0] nm, detuning 2000 Hz   == [1.5, 0, 0] A, 2 kHz
#   U4  natab_gcce2_bathadjust (two bath files): bath_shift [[0.1, 0.05, 0], [0, 0, 0]] nm == [[1, 0.5, 0], [0, 0, 0]] A
#   U5  natab_cce2_hahn as meCCE: jump rate 1 kHz                                   == 0.001 MHz
#   U6  arbitrary_gcce1_dm with a finite pulse (detuning_factor 0.98): duration 0.01 us == 10 ns
# pulse.bath.energy_shift and bandwidth are checked in ../check_deer; malformed objects and unknown units are refused in
# ../check_refusals.
#
# Usage:  bash test/code_verification/check_input_forms/units.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/units"
rm -rf "$W"; mkdir -p "$W"

# same <case> <label> <python edit of the unit run> <python edit of the plain run>: the two runs' outputs, byte for byte
same(){
    local case="$1" label="$2" c="$REPO/test/code_verification/check_numerics/cases/$1" ok=1 n=0 f
    python3 - "$c/ccein.json" "$W/$label" "$3" "$4" <<'PY'
import json, sys; from input_layout import flat, nest
src, out, edit, ref = sys.argv[1:5]
for name, code in ((out, edit), (out + "_ref", ref)):
    d = flat(json.load(open(src))); exec(code); d["outfile"] = name
    json.dump(nest(d), open(name + ".json", "w"), indent=1)
PY
    ( cd "$c" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$W/${label}_ref.json" > "$W/proc_${label}_ref" 2>&1 < /dev/null ) || ok=0
    ( cd "$c" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$W/$label.json" > "$W/proc_$label" 2>&1 < /dev/null ) || ok=0
    for f in "$W/${label}_ref"_*; do
        n=$((n+1)); cmp -s "$f" "$W/$label${f#$W/${label}_ref}" || { echo "    DIFF ${f#$W/}"; ok=0; }
    done
    [ $n -gt 0 ] || ok=0
    if [ $ok -eq 1 ]; then echo "  PASS  $label ($n files)"; PASS=$((PASS+1))
    else echo "  FAIL  $label : see $W/proc_$label"; FAIL=$((FAIL+1)); fi
}

same natab_cce2_hahn U1_other_units \
  "d['deltat'] = {'value': 1, 'unit': 'us'}; d['bfield'] = {'value': 50, 'unit': 'mT'}; d['rbath'] = {'value': 3, 'unit': 'nm'}; d['rdip'] = {'value': 1, 'unit': 'nm'}" \
  "d['deltat'] = 0.001; d['bfield'] = 500; d['rbath'] = 30; d['rdip'] = 10"
same natab_cce2_hahn U2_own_units_written \
  "d['deltat'] = {'value': d['deltat'], 'unit': 'ms'}; d['bfield'] = {'value': d['bfield'], 'unit': 'G'}; d['rbath'] = {'value': d['rbath'], 'unit': 'A'}; d['rdip'] = {'value': d['rdip'], 'unit': 'A'}" \
  "pass"
same arbitrary_gcce1_dm U3_qubit_xyz_detuning \
  "q = d['Qubit']['qubit'][0]; q['xyz'] = {'value': [0.15, 0, 0], 'unit': 'nm'}; q['detuning'] = {'value': 2000, 'unit': 'Hz'}" \
  "q = d['Qubit']['qubit'][0]; q['xyz'] = [1.5, 0, 0]; q['detuning'] = 2.0"
same natab_gcce2_bathadjust U4_bath_shift \
  "d['bath_shift'] = {'value': [[0.1, 0.05, 0], [0, 0, 0]], 'unit': 'nm'}" \
  "d['bath_shift'] = [[1.0, 0.5, 0], [0, 0, 0]]"
same natab_cce2_hahn U5_jump_rate \
  "d['method'] = 'meCCE'; d['nstep'] = 50; d['hamiltonian'] = {'meCCE': {'jump_operators': [{'bath_name': '13C', 'operator': 'z', 'rate': {'value': 1, 'unit': 'kHz'}}]}}" \
  "d['method'] = 'meCCE'; d['nstep'] = 50; d['hamiltonian'] = {'meCCE': {'jump_operators': [{'bath_name': '13C', 'operator': 'z', 'rate': 0.001}]}}"
same arbitrary_gcce1_dm U6_pulse_duration \
  "d['pulse'] = {'count': 1, 'detuning_factor': 0.98, 'duration': {'value': 0.01, 'unit': 'us'}}" \
  "d['pulse'] = {'count': 1, 'detuning_factor': 0.98, 'duration': 10.0}"

echo
echo "==========================================="
echo " units : PASS=$PASS FAIL=$FAIL"
echo "==========================================="
[ "$FAIL" -eq 0 ]
