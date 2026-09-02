#!/usr/bin/env bash
#
# Refused combinations, on the PHYSICAL fixtures.
# ---------------------------------------------------------------------------
# test/rotation/run.sh [10] already covers the refusals on its own toy systems.
# What is checked here is that the same guards fire on the real inputs -- a DFT
# tensor file, two qubits with explicit intmap tensors -- and that they fire
# BEFORE any tensor file is opened, which is the part that would regress quietly
# if the check ever drifted later into the run.
#
# A guard that passes here is one that cannot be reached by accident: the run
# stops, says why, and has not read anything it would have had to interpret in
# the wrong frame.
#
# Usage:  bash test/equivalence/reject.sh
# Needs:  axis1.sh already run (it writes the legacy cceins).
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

# reject <sys> <name> <expected message fragment> <python edit>
reject(){
    local sys="$1" name="$2" want="$3" edit="$4"
    local w="$WORK/$sys"
    [ -f "$w/legacy.json" ] || { echo "  FAIL  $name : no legacy.json"; FAIL=$((FAIL+1)); return; }
    python3 - "$w" "$name" <<PY
import json, sys
w, name = sys.argv[1:3]
d = json.load(open(f"{w}/legacy.json"))
$edit
d["outfile"] = f"./{name}"
json.dump(d, open(f"{w}/{name}.json","w"), indent=4)
PY
    local extra=""
    [ "$sys" = "S3" ] && extra="-s $FIX/1.single_data/state_DiaP1_1ppm_ -S $FIX/1.single_data/stateEx_DiaP1_1ppm_ -N 1"
    [ "$sys" = "S5" ] && extra="-s $EX/Diamond_multiNV_natab/NV2_13C2/state -N 1"
    ( cd "$w" && run_ccex "$BIN_DEV" "$name.json" "proc_$name" $extra )
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "  FAIL  $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return
    fi
    if ! grep -qF "$want" "$w/proc_$name"; then
        echo "  FAIL  $name : refused (exit=$rc) but not with the expected message"
        echo "        want: $want"
        grep -m1 -i "error" "$w/proc_$name" | sed 's/^/        got : /'
        FAIL=$((FAIL+1)); return
    fi
    # the guard must run before either tensor reader announces itself
    if grep -qE "Read the (Hyperfine|Quadrupole) interaction from DFT inputfile" "$w/proc_$name"; then
        echo "  FAIL  $name : refused, but a DFT tensor file was read first"; FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}

R111='{"enabled": True, "bath_axis": [0.0,0.0,1.0], "qubit_axis": [1.0,1.0,1.0]}'

echo "=== S4 : DFT tensor file + rotation"
reject S4 rj_hf_frame_missing "hf_readmode = 3 with the coordinate frame rotation on" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['qd_tensor_frame'] = 'bath'"
reject S4 rj_qd_frame_missing "qd_readmode" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['hf_tensor_frame'] = 'bath'"
reject S4 rj_tilted_bfield "bfield" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['hf_tensor_frame'] = 'bath'
d['coordinate_frame_rotation']['qd_tensor_frame'] = 'bath'
d['bfield'] = [100.0, 0.0, 30000.0]"
reject S4 rj_pos_frame_qubit "qubit_position_frame" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['hf_tensor_frame'] = 'bath'
d['coordinate_frame_rotation']['qd_tensor_frame'] = 'bath'
d['coordinate_frame_rotation']['qubit_position_frame'] = 'qubit'"

echo
echo "=== S5 : two qubits + rotation"
reject S5 rj_refq_missing "reference_qubit is required" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['qubit_position_frame'] = 'bath'
for e in d['Qubit']['intmap']: e['tensor_frame'] = 'bath'"
reject S5 rj_intmap_frame_missing "intmap tensor" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['reference_qubit'] = 'q1'
d['coordinate_frame_rotation']['qubit_position_frame'] = 'bath'"
reject S5 rj_refq_unknown "reference_qubit" \
  "d['coordinate_frame_rotation'] = $R111
d['coordinate_frame_rotation']['reference_qubit'] = 'q9'
d['coordinate_frame_rotation']['qubit_position_frame'] = 'bath'
for e in d['Qubit']['intmap']: e['tensor_frame'] = 'bath'"

echo
echo "=== S3 : removed key, and the Defect object forms"
reject S3 rj_magnetic_field_axis "magnetic_field_axis" \
  "d['magnetic_field_axis'] = [1.0,1.0,1.0]"
reject S3 rj_defaxis_without_rot "defect_axis_reference" \
  "d['defect_axis_reference'] = 'qubit_axis'"
reject S3 rj_zfs_de_no_axis "axis is missing" \
  "df = d['Defect'][0]
df['coordinate_frame'] = 'qubit'
df['zfs'] = {'D': 2870.0, 'unit': 'MHz'}"
reject S3 rj_zfs_E_nonzero "zfs.E is non-zero" \
  "df = d['Defect'][0]
df['coordinate_frame'] = 'qubit'
df['axis'] = [1.0,1.0,1.0]
df['zfs'] = {'D': 2870.0, 'E': 5.0, 'unit': 'MHz'}"
# axis is validated before zfs, so these are two different guards, not one.
reject S3 rj_axis_no_frame "axis is present but coordinate_frame is missing" \
  "df = d['Defect'][0]
df['axis'] = [1.0,1.0,1.0]"
reject S3 rj_zfs_object_no_frame "object-form zfs but has no coordinate_frame" \
  "df = d['Defect'][0]
df['zfs'] = {'D': 2870.0, 'unit': 'MHz'}"
reject S3 rj_bad_unit "unit" \
  "df = d['Defect'][0]
df['coordinate_frame'] = 'qubit'
df['eqs'] = {'values': [2.044], 'unit': 'furlongs'}"
reject S3 rj_rot_alias_both "both" \
  "d['coordinate_frame_rotation'] = $R111
d['bath_coordinate_rotation'] = $R111"

echo
echo "==========================================="
echo " refused combinations : PASS=$PASS FAIL=$FAIL"
echo "==========================================="
[ "$FAIL" -eq 0 ]
