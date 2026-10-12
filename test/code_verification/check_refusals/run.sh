#!/usr/bin/env bash
#
# Refused combinations, on the PHYSICAL fixtures.
# ---------------------------------------------------------------------------
# test/code_verification/check_rotation/run.sh [10] already covers the refusals on its own toy systems.
# What is checked here is that the same guards fire on the real inputs -- a DFT
# tensor file, two qubits with explicit intmap tensors -- and that they fire
# BEFORE any tensor file is opened, which is the part that would regress quietly
# if the check ever drifted later into the run.
#
# A guard that passes here is one that cannot be reached by accident: the run
# stops, says why, and has not read anything it would have had to interpret in
# the wrong frame.
#
# Usage:  bash test/code_verification/check_refusals/run.sh
# Needs:  bin/main.out; ../check_input_forms/legacy.sh is run first (it writes the legacy cceins).
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/../check_input_forms/lib.sh"
# the legacy-form inputs the refused variants are edited from
bash "$HERE/../check_input_forms/legacy.sh" > /dev/null || { echo "legacy inputs could not be written"; exit 1; }

# reject <sys> <name> <expected message fragment> <python edit> [read]
# "read": the refusal is about a tensor file, so it comes after a reader has started
reject(){
    local sys="$1" name="$2" want="$3" edit="$4" when="${5:-}"
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
    [ "$sys" = "S3" ] && extra="-s $FIX/single_data/state_DiaP1_1ppm_ -S $FIX/single_data/stateEx_DiaP1_1ppm_ -N 1"
    [ "$sys" = "S5" ] && extra="-s $EX/tutorial/03_Diamond_multiNV_natab/NV2_13C2/state -N 1"
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
    if [ "$when" != read ] && grep -qE "Read the (Hyperfine|Quadrupole) interaction from DFT inputfile" "$w/proc_$name"; then
        echo "  FAIL  $name : refused, but a DFT tensor file was read first"; FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}

R111='{"enabled": True, "bath_axis": [0.0,0.0,1.0], "qubit_axis": [1.0,1.0,1.0]}'

# A valid input with a command-line option that is not a number (review 78d583 R5): -N and -B were read by
# atoi/atof, so "-N nope" ran an ensemble (nstate 0) and "-N 1.5" one state.
reject_cli(){
    local sys="$1" name="$2" want="$3"; shift 3
    local w="$WORK/$sys"
    ( cd "$w" && run_ccex "$BIN_DEV" legacy.json "proc_$name" "$@" )
    local rc=$?
    if [ $rc -eq 0 ]; then echo "  FAIL  $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return; fi
    if ! grep -qF -- "$want" "$w/proc_$name"; then
        echo "  FAIL  $name : refused (exit=$rc) but not with the expected message"; FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}
echo "=== S1a : command-line numbers"
reject_cli S1a rj_cli_nstate_word "-N must be an integer" -N nope
reject_cli S1a rj_cli_nstate_fraction "-N must be an integer" -N 1.5
reject_cli S1a rj_cli_bfield_word "-B must be a finite number" -B nope
reject_cli S1a rj_cli_bfield_trailing "-B must be a finite number" -B 500G
# The command line wins over the input file whatever the order (review b9eb09 R1): -N before -f took nstate from the file.
order_check(){
    local w="$WORK/S1a" a b
    a="$(cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -N 3 -f legacy.json -P -o ./order_a < /dev/null 2>&1 | grep -E '^ +nstate +:')"
    b="$(cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -f legacy.json -N 3 -P -o ./order_b < /dev/null 2>&1 | grep -E '^ +nstate +:')"
    if [ -n "$a" ] && [ "$a" = "$b" ] && grep -Eq ':[[:space:]]+3([[:space:]]|$)' <<< "$a"; then echo "  PASS  cli_order_nstate"; PASS=$((PASS+1))
    else echo "  FAIL  cli_order_nstate : '-N 3 -f' gave [$a], '-f -N 3' gave [$b]"; FAIL=$((FAIL+1)); fi
}
order_check
# a qubitfile the input names but that is not there ran the qubit at (0,0,0) (review 4dd5fd R1)
reject S1a rj_qubitfile_missing "cannot open the qubit file" "d['qubitfile'] = d['qubitfile'] + '.missing'"
# nstep 0 left the result arrays empty and the output read result[0] (review 7eed38 R1)
reject S1a rj_nstep_zero "nstep (0 < 1)" "d['nstep'] = 0"
# A Manual sequence (review 5f7834 R2): a sum of fractions short of 1, a fraction that is not a/b or a number all
# through, a row that is not [fraction, axis, angle] used to pass (or skip the row); a number fraction used to throw
# and an angle of 90.5 became 90 - sequence_accepted checks those two now go through as given.
MANUAL="d['method'] = 'gCCE'; d['npulse'] = 2; d['pulsename'] = 'Manual'"   # gCCE: CCE refuses an I pulse on its own
reject S1a rj_seq_sum_short "sum to 0.75" "$MANUAL; d['sequence'] = [['1/4','I',0],['1/4','X',180],['1/4','I',0]]"
reject S1a rj_seq_fraction_text "must be a/b" "$MANUAL; d['sequence'] = [['1/2x','I',0],['1/4','X',180],['1/4','I',0]]"
reject S1a rj_seq_row_short "must be [fraction, axis, angle]" "$MANUAL; d['sequence'] = [['1/2','I'],['1/4','X',180],['1/4','I',0]]"
sequence_accepted(){
    local w="$WORK/S1a"
    python3 - "$w" <<'PY'
import json, sys
w = sys.argv[1]; d = json.load(open(f"{w}/legacy.json"))
d.update(method='gCCE', npulse=2, pulsename='Manual', sequence=[[0.25,'I',0],['1/2','X',90.5],[0.25,'I',0]], outfile='./seq_ok')
json.dump(d, open(f"{w}/seq_ok.json", "w"))
PY
    if (cd "$w" && "$MPIRUN" -n 1 "$BIN_DEV" -f seq_ok.json -P -o ./seq_ok < /dev/null > proc_seq_ok 2>&1) &&
       grep -Eq '0\.250 +0\.750 +0\.500 +X +90\.50' "$w/proc_seq_ok"; then
        echo "  PASS  sequence_accepted (a number fraction, an angle of 90.5)"; PASS=$((PASS+1))
    else echo "  FAIL  sequence_accepted : see $w/proc_seq_ok"; FAIL=$((FAIL+1)); fi
}
sequence_accepted
# a finite double outside the float it is stored in became inf (review fa7f39 R3)
reject S1a rj_bfield_float_overflow "within the float range" "d['bfield'] = 1e40"

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
# D43: alphaidx/betaidx count the qubit eigenstates by descending energy and override alphams/betams; written together and
# naming other states, the run stops (S5 writes alphams 1 / betams 0 for both NVs; index 2 is not |0,0>, the lowest level)
reject S5 rj_alphaidx_vs_alphams "Keep one of the two" \
  "d['Qubit']['alphaidx'] = 1; d['Qubit']['betaidx'] = 2" read
# only what is written constrains (review b467fe R1): q1's betams 0 alone against betaidx 0, the highest level |+1,+1>
reject S5 rj_alphaidx_one_side "Keep one of the two" \
  "d['Qubit']['alphaidx'] = 1; d['Qubit']['betaidx'] = 0
for q in d['Qubit']['qubit']: q.pop('alphams', None); q.pop('betams', None)
d['Qubit']['qubit'][0]['betams'] = 0" read
# and the same kind of input that agrees runs: q1's alphams 1 alone with alphaidx 0 (|+1,+1>); betaidx 1 is not written
# against, so it is free. A short run (nstep 2), since the check sits where the run picks psia/psib.
idx_ms_accepted(){
    local w="$WORK/S5"
    python3 - "$w" <<'PY'
import json, sys
w = sys.argv[1]; d = json.load(open(f"{w}/legacy.json"))
for q in d['Qubit']['qubit']: q.pop('alphams', None); q.pop('betams', None)
d['Qubit']['qubit'][0]['alphams'] = 1
d['Qubit'].update(alphaidx=0, betaidx=1); d.update(nstep=2, outfile='./idx_ms_ok')
json.dump(d, open(f"{w}/idx_ms_ok.json", "w"))
PY
    if (cd "$w" && run_ccex "$BIN_DEV" idx_ms_ok.json proc_idx_ms_ok -s "$EX/tutorial/03_Diamond_multiNV_natab/NV2_13C2/state" -N 1) &&
       grep -Eq 'psia is mostly +: +m_s \|\+1, \+1>' "$w/proc_idx_ms_ok"; then
        echo "  PASS  idx_ms_accepted (alphaidx with a partly written alphams that agrees)"; PASS=$((PASS+1))
    else echo "  FAIL  idx_ms_accepted : see $w/proc_idx_ms_ok"; FAIL=$((FAIL+1)); fi
}
idx_ms_accepted
# D28: Qubit.intmap[].tensor shorter than 3x3 was padded with zeros
reject S5 rj_intmap_tensor_short "3x3" \
  "d['Qubit']['intmap'][0]['tensor'] = [[1.0]]"
reject S5 rj_intmap_tensor_string "finite JSON number" \
  "d['Qubit']['intmap'][0]['tensor'] = [[1.0,0,0],[0,'a',0],[0,0,1.0]]"
# review 20261008-autoconv-d42 #5: a qubit array shorter than nqubit left the last qubit allocated but never set
reject S5 rj_qubit_array_short "Qubit.qubit must be an array of nqubit (2) entries, 1 given" \
  "d['Qubit']['qubit'] = d['Qubit']['qubit'][:1]"
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
# The one behaviour change that can stop a v1.1.0 input which used to complete. The
# arrays hold navaax+1 entries (index 0 is the main spin), so navaax+1 is one past the
# end; it used to be accepted and written to. S3's own fixture uses 1..navaax, and
# axis 1 shows that side is untouched.
reject S3 rj_iax_out_of_range "out of range" \
  "df = d['Defect'][0]
df['hypf'] = df['hypf'] + [[df['navaax'] + 1, '14N', [0.0]*9]]"

# D28 (2026-10-02): the entries of the indexed form are checked as the object form checks them. A string
# or fractional index, a non-string spin name, a non-number value and a field that is not an array were
# read as 0, truncated, crashed or skipped.
reject S3 rj_iax_string "out of range" \
  "df = d['Defect'][0]
df['hypf'] = df['hypf'] + [['1', '14N', [0.0]*9]]"
reject S3 rj_iax_fraction "out of range" \
  "df = d['Defect'][0]
df['hypf'] = df['hypf'] + [[1.5, '14N', [0.0]*9]]"
reject S3 rj_spin_not_string "spin name must be a string" \
  "df = d['Defect'][0]
df['hypf'] = df['hypf'] + [[1, 7, [0.0]*9]]"
reject S3 rj_value_string "must be a finite JSON number" \
  "df = d['Defect'][0]
df['hypf'] = df['hypf'] + [[1, '14N', ['x'] + [0.0]*8]]"
reject S3 rj_rxyzs_not_array "must be an array of [index, name, values] entries" \
  "df = d['Defect'][0]
df['rxyzs'] = 5"

# D42 (2026-10-08): an array shorter than the input needs was filled with empty strings or 0 and used. The bathadjust cases
# ran to the end before (zero offsets); the types case was refused later anyway (the hypf entries name 14N, which the empty
# type then lacked), so it checks the earlier, clearer refusal.
reject S3 rj_types_short "types has 0 entries" \
  "d['Defect'][0]['types'] = []"
reject S3 rj_bathadjust_row_short "bathadjust[0] has 2 numbers" \
  "d['bathadjust'] = [[0.0, 0.0]]"
reject S3 rj_bathadjust_rows_short "bathadjust has 0 rows" \
  "d['bathadjust'] = []"
# navaax 4 means axes 1-4: an avaax file with 5 was read one past the tables (2026-10-09); refused while the file is read
reject S3 rj_avaax_past_navaax "The principal axis is out of range" \
  "lines = open(d['avaaxfile']).read().split('\n')
open(f'{w}/avaax_axis5','w').write('\n'.join([lines[0]] + ['5' if l.strip() else l for l in lines[1:]]))
d['avaaxfile'] = f'{w}/avaax_axis5'" read

reject S3 rj_bad_unit "unit" \
  "df = d['Defect'][0]
df['coordinate_frame'] = 'qubit'
df['eqs'] = {'values': [2.044], 'unit': 'furlongs'}"
reject S3 rj_rot_alias_both "both" \
  "d['coordinate_frame_rotation'] = $R111
d['bath_coordinate_rotation'] = $R111"

echo
echo "=== S4 : a tensor file without a condition the reader needs"
# READ_Tensor_const / READ_Tensor_etc said so and then ended with exit 0 until 2026-10-01 (D12)
TENSOR_EDIT='
import re
key, cond = ARGS
src = d[key]
dst = f"{w}/{name}_" + src.rsplit("/", 1)[1]
open(dst, "w").write("".join(l for l in open(src) if cond not in l))
d[key] = dst
if key == "qd_tensorfile": d["qd_tensorfile_woqubit"] = dst'
reject S4 rj_hf_no_gfactor "Tensor Condition (g-factor___) can not read" \
  "ARGS = ('hf_tensorfile', 'g-factor___')
$TENSOR_EDIT" read
reject S4 rj_qd_no_etc "Tensor Condition (etc____) can not read" \
  "ARGS = ('qd_tensorfile', 'etc____')
$TENSOR_EDIT" read

echo
echo "=== a present key of the wrong type"
# The defaulting readers took a present key of another type for a missing one and ran
# with the default until 2026-10-02 (D18): "knight": "true" was false, a quoted gyro the
# default gyro, nstep 10.5 was 10. Absent keys keep their defaults; a boolean may still
# be 0 or 1 ("addsubclus": 0 in 29 order-2 inputs) and a qubit's detuning [x].
reject S1a rj_type_knight "knight must be true or false" "d['knight'] = 'true'"
reject S1a rj_type_hfmedi "hfmedi must be true or false" "d['hfmedi'] = 'true'"
reject S1a rj_type_bool_number "addsubclus must be true or false" "d['addsubclus'] = 2"
reject S5 rj_type_qubit_gyro "gyro must be a number" "d['Qubit']['qubit'][0]['gyro'] = str(d['Qubit']['qubit'][0]['gyro'])"
reject S1a rj_type_nstep_fraction "nstep must be an integer" "d['nstep'] = d['nstep'] + 0.5"
reject S1a rj_type_bfield_string "bfield must be a number" "d['bfield'] = str(d['bfield'])"
reject S1a rj_type_nk_element "nk[1] must be an integer" "d['nk'] = [0, '1', 0]"
reject S1a rj_type_bathadjust_element "bathadjust[0][2] must be a number" "d['bathadjust'] = [[0.0, 0.0, '1.0']]"
reject S1a rj_type_method "method must be a string" "d['method'] = 1"

echo
echo "=== a key nobody reads, close to one that took its default"
# Until 2026-10-02 (D24) such a key was ignored without a word: "nstates" left nstate at 0, and
# the top-level "alphams" of the qubitfile form (which reads qalphams) left the NV branch at +1.
reject S1a rj_typo_nstates "did you mean \"nstate\"" "d['nstates'] = d.pop('nstate')"
reject S1a rj_typo_qubitfile_alphams "did you mean \"qalphams\"" "d['alphams'] = d.pop('qalphams')"

echo
echo "=== a result that cannot be written in full"
# Output_save and Output_save_info checked only that the files opened (and
# Output_save_info only the first of the two): with a full disk, a quota or a
# file-size limit the run left a truncated result and ended with exit 0 until
# 2026-10-01 (D10). On a small case whose result is 18 KB (47 KB in info mode):
#   full     the result files are links to /dev/full, so every write fails with
#            ENOSPC, as on a full disk. Until 2026-10-02 a 4 KB file-size limit
#            (ulimit -f) did this, which Open MPI's mpirun did not survive: its own
#            shared-memory files met the limit and MPI_Init failed (idm, job 16853)
#   nodiv    a directory where the _noDiv file has to go, so only the second file
#            cannot be opened
CASE="$REPO/test/code_verification/check_numerics/cases/natab_cce2_hahn"
write_fail(){
    local name="$1" savemode="$2" how="$3" want="$4"
    local outw="$WORK/write_failure/$name"
    rm -rf "${outw:?}"; mkdir -p "$outw"
    python3 - "$CASE/ccein.json" "$outw/ccein.json" "$savemode" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); d["savemode"] = sys.argv[3]
json.dump(d, open(sys.argv[2], "w"), indent=1)
PY
    [ "$how" = nodiv ] && mkdir "$outw/out_noDiv"
    [ "$how" = full ] && ln -s /dev/full "$outw/out_noDiv" && ln -s /dev/full "$outw/out_wiDiv"
    # relative paths in the case's input are read from the case's folder
    ( cd "$CASE" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$outw/ccein.json" -o "$outw/out" 2>&1 < /dev/null ) > "$outw/proc"
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "  FAIL  $name : ended with exit 0"; FAIL=$((FAIL+1))
    elif ! grep -qF "$want $outw/out_" "$outw/proc"; then
        echo "  FAIL  $name : exit=$rc but no message naming the file"; FAIL=$((FAIL+1))
    else
        echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
    fi
}
write_fail write_full_normal  normal full  "cannot write file"
write_fail write_full_info    info   full  "cannot write file"
write_fail open_nodiv_normal  normal nodiv "cannot open file"
write_fail open_nodiv_info    info   nodiv "cannot open file"

echo
echo "=== the factors and except work modes of method D's regroup"
# adaptive-regroup-stage3.md 4절, -stage3b.md 1절: a request CCEX cannot answer for these modes is refused at the
# request, before any state is calculated. The driver checks the same before it asks (an in-process refusal would end
# the job), so these requests are written by hand: a P1 input of the regroup regression (six sites, five states).
P1IN="$REPO/test/cluster_verification/check_adaptive/regression/inputs/regroup-p1-0.json"
reject_work(){
    local name="$1" want="$2" edit="$3" request="$4"; shift 4
    local w="$WORK/adaptive_work"; mkdir -p "$w"
    python3 - "$P1IN" "$w" "$name" "$request" <<PY
import json, sys
src, w, name, request = sys.argv[1:5]
d = json.load(open(src))
$edit
r = json.loads(request); r["response"] = f"{w}/{name}.response"
json.dump(r, open(f"{w}/{name}.request.json", "w"))
d["adaptive_work_file"] = f"{w}/{name}.request.json"
d["outfile"] = f"{w}/{name}"
json.dump(d, open(f"{w}/{name}.json", "w"), indent=1)
PY
    # relative paths in the input are the repository's
    ( cd "$REPO" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/$name.json" "$@" > "$w/proc_$name" 2>&1 < /dev/null )
    local rc=$?
    if [ $rc -eq 0 ]; then echo "  FAIL  $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return; fi
    if ! grep -qF -- "$want" "$w/proc_$name"; then
        echo "  FAIL  $name : refused (exit=$rc) but not with the expected message"
        echo "        want: $want"
        grep -m1 -i "error" "$w/proc_$name" | sed 's/^/        got : /'
        FAIL=$((FAIL+1)); return
    fi
    if ls "$w/$name.response".*state-* > /dev/null 2>&1; then
        echo "  FAIL  $name : refused, but after a state was calculated"; FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}
B='"bundles": [[0, 3], [1, 5], [2, 4]]'
F="{\"mode\": \"factors\", $B, \"state\": 1, \"step\": 3}"
E="{\"mode\": \"except\", $B, \"except\": [[0, 1, 3, 5]]}"
reject_work rj_factors_state_zero "factors state is not a state of this run" "" "${F/\"state\": 1/\"state\": 0}"
reject_work rj_factors_state_past "factors state is not a state of this run" "" "${F/\"state\": 1/\"state\": 6}"
reject_work rj_factors_step "factors step is outside 0..nstep-1" "" "${F/\"step\": 3/\"step\": 30}"
reject_work rj_factors_step_text "factors work requires integer state and step" "" "${F/\"step\": 3/\"step\": \"3\"}"
reject_work rj_factors_gcce "factors/except work requires method cce" "d['method'] = 'gCCE'" "$F"
reject_work rj_factors_order1 "factors/except work requires bundle order 2" "d['order'] = 1" "$F"
reject_work rj_factors_probes "probes require probe mode" "" "${F%\}}, \"probes\": [[0]]}"
reject_work rj_factors_contributions "cannot be combined with cluster_contributions" "d['cluster_contributions'] = True" "$F"
reject_work rj_factors_plan_only "probe, factors and except work require dynamics" "" "$F" -P
reject_work rj_state_without_factors "state and step require factors mode" "" "{\"mode\": \"run\", $B, \"state\": 1, \"step\": 3}"
reject_work rj_except_gcce "factors/except work requires method cce" "d['method'] = 'gCCE'" "$E"
reject_work rj_except_order1 "factors/except work requires bundle order 2" "d['order'] = 1" "$E"
reject_work rj_except_probes "probes require probe mode" "" "${E%\}}, \"probes\": [[0]]}"
reject_work rj_except_empty "except mode requires a nonempty except list" "" "{\"mode\": \"except\", $B, \"except\": []}"
reject_work rj_except_twice "except lists a term twice" "" "{\"mode\": \"except\", $B, \"except\": [[0, 3], [0, 3]]}"
reject_work rj_except_not_a_term "an except entry is not a term of this partition's expansion" "" "{\"mode\": \"except\", $B, \"except\": [[0, 1]]}"
reject_work rj_except_plan_only "probe, factors and except work require dynamics" "" "$E" -P
reject_work rj_except_without_mode "except requires except mode" "" "{\"mode\": \"run\", $B, \"except\": [[0, 3]]}"

echo
echo "=== evolution_update"
# exact | recursive (2026-10-11): recursive (A16 b, A18) exists on the CPU for cCCE (single state and ensemble) and for
# gCCE with propagator=eigen (evolution vector and matrix); on the GPU the run would compute as exact and gCCE expm has no
# phases to carry over, so a given recursive there is refused (without the key those run exact, A29). On check_numerics
# cases, refused at the start. (The ensemble cCCE and gCCE matrix were refused until A18; they run now,
# ../check_input_forms/evolution_update.sh.) Since A36 (2026-10-12) backend auto in this CPU-only binary computes every
# cluster on the CPU and a given recursive runs (check_input_forms E5); backend cuda stays refused here, and recursive
# with backend auto or cuda in bin/main_gpu.out is refused in test/gpu_verification/check_gpu/check_options.py.
reject_case(){
    local case="$1" name="$2" want="$3" edit="$4"
    local w="$WORK/evolution_update_refusals"; mkdir -p "$w"
    python3 - "$REPO/test/code_verification/check_numerics/cases/$case/ccein.json" "$w/$name.json" <<PY
import json, sys
d = json.load(open(sys.argv[1]))
$edit
json.dump(d, open(sys.argv[2], "w"), indent=1)
PY
    # relative paths in the case's input are read from the case's folder
    ( cd "$REPO/test/code_verification/check_numerics/cases/$case" &&
      "$MPIRUN" -n 1 "$BIN_DEV" -f "$w/$name.json" -o "$w/$name" > "$w/proc_$name" 2>&1 < /dev/null )
    local rc=$?
    if [ $rc -eq 0 ]; then echo "  FAIL  $name : accepted (exit 0), expected refusal"; FAIL=$((FAIL+1)); return; fi
    if ! grep -qF -- "$want" "$w/proc_$name"; then
        echo "  FAIL  $name : refused (exit=$rc) but not with the expected message"
        echo "        want: $want"
        grep -m1 -i "error" "$w/proc_$name" | sed 's/^/        got : /'
        FAIL=$((FAIL+1)); return
    fi
    echo "  PASS  $name (exit=$rc)"; PASS=$((PASS+1))
}
reject_case p1_cce2_single rj_update_value "evolution_update must be exact or recursive" "d['evolution_update'] = 'fast'"
reject_case p1_cce2_single rj_update_type "evolution_update must be exact or recursive" "d['evolution_update'] = 1"
reject_case p1_manual_hahn rj_update_gcce_expm "this run has method=gCCE with propagator=expm" \
  "d['evolution_update'] = 'recursive'; d['evolution'] = 'matrix'; d['propagator'] = 'expm'"
reject_case p1_cce2_single rj_update_backend "this run has backend=cuda" \
  "d['evolution_update'] = 'recursive'; d['backend'] = 'cuda'"

echo "=== gpu_options"
# The GPU settings are one block since 2026-10-12 (A26, docs/user/InputFile_Backend.md#tag_gpu_options): the flat keys
# before it are refused with the new form named, and so are an unknown key of the block, a block that is not an object,
# and a packed_max_dim other than 0, 2 or 4. packed_max_dim needs method=CCE and the GPU binary (this check runs
# bin/main.out), both refused before the bath is read. The accepted forms are in ../check_input_forms/gpu_options.sh.
reject_case p1_cce2_single rj_gpu_min_dim "gpu_min_dim is now gpu_options.min_dim" "d['gpu_min_dim'] = 256"
reject_case p1_cce2_single rj_gpu_batch_steps "gpu_batch_steps is now gpu_options.evolution_batch_steps" "d['gpu_batch_steps'] = 64"
reject_case p1_cce2_single rj_gpu_batch_clusters "gpu_batch_clusters is now gpu_options.evolution_batch_size" "d['gpu_batch_clusters'] = 8"
reject_case p1_cce2_single rj_gpu_cluster_batch_mb "gpu_cluster_batch_mb is now gpu_options.cluster_batch_mb" "d['gpu_cluster_batch_mb'] = 64"
reject_case p1_cce2_single rj_gpu_cache_mb "gpu_cache_mb is now gpu_options.cache_mb" "d['gpu_cache_mb'] = 64"
reject_case p1_cce2_single rj_gpu_options_type "gpu_options must be an object" "d['gpu_options'] = 256"
reject_case p1_cce2_single rj_gpu_options_key "unknown key gpu_options.packed_dim" "d['gpu_options'] = {'packed_dim': 4}"
reject_case p1_cce2_single rj_packed_3 "gpu_options.packed_max_dim must be 0 (off), 2 or 4" "d['gpu_options'] = {'packed_max_dim': 3}"
reject_case p1_cce2_single rj_packed_8 "gpu_options.packed_max_dim must be 0 (off), 2 or 4" "d['gpu_options'] = {'packed_max_dim': 8}"
reject_case p1_cce2_single rj_packed_string "finite JSON number" "d['gpu_options'] = {'packed_max_dim': '4'}"
reject_case p1_cce2_single rj_min_dim_zero "gpu_options.min_dim must be a positive 32-bit integer" "d['gpu_options'] = {'min_dim': 0}"
reject_case p1_cce2_single rj_packed_cpu_binary "gpu_options.packed_max_dim needs the GPU binary" "d['gpu_options'] = {'packed_max_dim': 4}"
reject_case p1_manual_hahn rj_packed_gcce "gpu_options.packed_max_dim works with method=CCE" \
  "d['gpu_options'] = {'packed_max_dim': 4}; d['method'] = 'gCCE'"

echo
echo "==========================================="
echo " refused combinations : PASS=$PASS FAIL=$FAIL"
echo "==========================================="
[ "$FAIL" -eq 0 ]
