#!/usr/bin/env bash
# DEER bath pulses: the sequence merge (unit), then the closed form and the evolution paths (Task A3).
# Usage: bash test/code_verification/check_deer/run.sh      (CCEX_BIN, CCEX_MPI as for the others)
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../../.." && pwd)"
source "$HERE/../check_input_forms/lib.sh"
W="$HERE/_work"; rm -rf "${W:?}"; mkdir -p "$W"
PASS=0; FAIL=0
UNIT_ENV_FILE="$(mktemp)"
bash "$REPO/test/unit_build_env.sh" > "$UNIT_ENV_FILE" || { rm -f "$UNIT_ENV_FILE"; exit 1; }
unit_value(){ sed -n "s/^$1='\(.*\)'\$/\1/p" "$UNIT_ENV_FILE"; }
CXX="$(unit_value CCEX_UNIT_CXX)"; CXXFLAGS="$(unit_value CCEX_UNIT_CXXFLAGS) $(unit_value CCEX_UNIT_INCLUDES)"
UNIT_LIBS="$(unit_value CCEX_UNIT_LIBS)"; UNIT_OBJS="$(unit_value CCEX_UNIT_OBJS)"
for kv in $(unit_value CCEX_UNIT_ENV); do export "$kv"; done; rm -f "$UNIT_ENV_FILE"
"$CXX" "$HERE/test_merge_unit.cpp" $UNIT_OBJS -o "$W/test_merge_unit.out" $CXXFLAGS -w \
    -I"$REPO/external/eigen" -I"$REPO/external/uthash/include/" -I"$REPO/include/" $UNIT_LIBS \
    && "$W/test_merge_unit.out" && PASS=$((PASS+1)) || { echo "  FAIL  merge unit"; FAIL=$((FAIL+1)); }
# explicit "pulse.bath": {"count": 0} is the same run as no bath block, byte for byte (outputs and nothing else)
CASE="$REPO/test/code_verification/check_numerics/cases/p1_manual_hahn"
python3 - "$CASE/ccein.json" "$W" <<'PY'
import json, sys; from input_layout import flat, nest
w = sys.argv[2]
d = flat(json.load(open(sys.argv[1]))); d["outfile"] = f"{w}/bn_absent"; json.dump(nest(d), open(f"{w}/bn_absent.json", "w"), indent=1)
d["pulse"]["bath"] = {"count": 0}; d["outfile"] = f"{w}/bn_zero"; json.dump(nest(d), open(f"{w}/bn0.json", "w"), indent=1)
PY
# s_range [1, 1] (nstate 1 until A94) with bath_state_file "../../inputs/p1_1ppm/state_" (until A63 given as -s as well), so the outputs are <o>_state1_*
( cd "$CASE" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$W/bn_absent.json" > "$W/proc_bn_absent" 2>&1 < /dev/null )
( cd "$CASE" && "$MPIRUN" -n 1 "$BIN_DEV" -f "$W/bn0.json" > "$W/proc_bn_zero" 2>&1 < /dev/null )
if [ -s "$W/bn_absent_state1_noDiv" ] && cmp -s "$W/bn_absent_state1_noDiv" "$W/bn_zero_state1_noDiv" \
   && cmp -s "$W/bn_absent_state1_wiDiv" "$W/bn_zero_state1_wiDiv"; then
    echo "  PASS  pulse.bath.count 0 == absent (bytes)"; PASS=$((PASS+1))
else echo "  FAIL  pulse.bath.count 0 differs from absent (or no output)"; FAIL=$((FAIL+1)); fi
# --- Task A3 appends the run checks here ---
# six branch cases (adding_MECCE example/z.Tests/DEER): gCCE order 1 ensemble, closed form at order 1 (check_deer.py)
python3 - "$HERE/inputs" "$REPO/test/code_verification/bath_data/single_data/DiaP1_gyro" "$W" <<'PY'
import json, sys, copy; from input_layout import flat, nest
inp, gyro, w = sys.argv[1:4]

def merge(d, extra):   # nested blocks (pulse, Qubit, hamiltonian) merge key by key
    for k, v in extra.items():
        if isinstance(v, dict) and isinstance(d.get(k), dict):
            merge(d[k], v)
        else:
            d[k] = copy.deepcopy(v)
    return d

# the branch's ccein.js in the A47 keys: avaaxfile -> Defect.axis_file (configuration.defect_axis_file since A64), Defect list -> Defect.species
# (dfname -> name, naddspin -> n_nuclear_spins, navaax -> n_configurations; apprx was removed, "it changed nothing"),
# overhaus -> hamiltonian.common.qubit_bath_mean_field, npulse -> pulse.count, addsubclus -> add_subclusters;
# nk [0, 0], hf_readmode 0 and savemode normal are the defaults and are left out
base = {
 "method": "gCCE", "quantity": "coherence",
 "configuration": {"gyrofile": gyro, "bathfile": [f"{inp}/bath_DEERtest_12spin.txt"], "seed": 3528,
                   "defect_axis_file": [f"{inp}/avaax_DEERtest"]},
 "Qubit": {"qubit": [{"name": "qubit1", "spin": 1.0, "gyro": -17608.597050, "xyz": [0.0, 0.0, 0.0],
            "detuning": 0.0}],
           "alpha_state": {"ms": [-1]}, "beta_state": {"ms": [0]},
           "intmap": [{"between": ["qubit1", "qubit1"],
                       "tensor": [[-960000.0, 0, 0], [0, -960000.0, 0], [0, 0, 1920000.0]]}]},
 "hamiltonian": {"common": {"qubit_bath_mean_field": True}},
 "order": 1, "bfield": 180, "rbath": 200, "rdip": 200, "deltat": 2.5e-04, "nstep": 41,
 "add_subclusters": False, "pulse": {"count": 1},
 "Defect": {"species": [{"name": "P1", "n_nuclear_spins": 0, "n_configurations": 8, "detuning": [
     [1, "e", 79.975539], [2, "e", 59.579552], [3, "e", 59.579552], [4, "e", 59.579552],
     [5, "e", -79.975539], [6, "e", -59.579552], [7, "e", -59.579552], [8, "e", -59.579552]]}]}}

def pump(tones, seq):
    return {"pulse": {"bath": {"count": 1, "defect": "P1", "levels": [0.5, -0.5],   # S is the P1's own (A94)
                               "energy_shift": tones, "sequence": seq}}}
HAHN = [["0.50", "X", 180.0], ["0.50", "I", 0.0]]
cases = {
 "nopulse":    {},
 "onres":      pump([79.975539], HAHN),
 "offres":     pump([84.98], HAHN),
 "multi_hit":  pump([79.975539, 200], HAHN),   # the branch value: 79.98 is 28 rad*kHz off, outside the window
 "multi_miss": pump([200, 300], HAHN),
 "nearpulse":  pump([79.975539], [["0.505", "X", 180.0], ["0.495", "I", 0.0]]),
}
# the evolution paths on a pure state (s_range [1, 1]): matrix, vector, exact and recursive must agree with bath pulses.
# evolution_update is written out: its default is recursive on these paths since 15a136555 (general.cpp:539-548)
st = {"configuration": {"s_range": [1, 1], "bath_state_file": [f"{inp}/state_"], "seed": 3528}}
paths = {"path_matrix":           {"evolution": "matrix", "evolution_update": "exact"},
         "path_matrix_recursive": {"evolution": "matrix", "evolution_update": "recursive"},
         "path_vector":           {"evolution": "vector", "evolution_update": "exact"},
         "path_recursive":        {"evolution": "vector", "evolution_update": "recursive"},
         "path_nopump":           {"evolution": "vector", "evolution_update": "exact"}}
for path, extra in paths.items():
    cases[path] = merge(merge({}, cases["onres"] if path != "path_nopump" else {}), {**st, **extra})
# levels of a spin 1 on the P1 (S = 1/2): pulse.bath.spin was written beside them until A94 and refused when it was not the P1's
cases["bad_bspin"] = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"levels": [1.0, 0.0]}}})
# a defect name no spin has, and levels that are not two levels of S
cases["bad_defect"]   = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"defect": "P2"}}})
cases["bad_ms_equal"] = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"levels": [0.5, 0.5]}}})
cases["bad_ms_range"] = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"levels": [1.5, -0.5]}}})
# --- megCCE (rate 0: no jump operators; evolution_update exact by default for ME) and cluster_mode pcce
cases["meg_onres"] = merge(merge({}, cases["onres"]), {"method": "megCCE"})
# the bath pi on an edge of its own (qubit identity there): megCCE follows the merged sequence, each pulse once
cases["meg_nearpulse"] = merge(merge({}, cases["nearpulse"]), {"method": "megCCE"})
# the resonant tone second: every tone is tried, not the first one only (A3 review)
cases["multi_hit_rev"] = pump([200, 79.975539], HAHN)
# pcce with partition_size 1 is the standard single-spin set (base seed 3528 fixes the partition search)
PC1 = {"cluster_mode": "pcce", "pcce": {"partition_size": 1}}
cases["pg_onres"]   = merge(merge({}, cases["onres"]), PC1)
cases["pmeg_onres"] = merge(merge({}, cases["meg_onres"]), PC1)
# order 2 on the pure state, and the same with the bath file's spins in reverse order (the axis file and the state
# file reordered with them, line for line): the same _noDiv (A3 review: kron order and identity placement). Any order
# of two spins either keeps a pair or swaps it, and a U_b built in reverse kron order swaps the two spins' operators the
# same way in both runs; so gCCE also runs order 3 with the spins rotated by 5 (a 3-spin cluster that wraps around is
# rotated, not reversed, which such a U_b does not survive)
def reorder(src, dst, perm):   # the header line, then the spin lines in the order perm gives
    lines = open(src).read().splitlines()
    open(dst, "w").write("\n".join([lines[0]] + perm(lines[1:])) + "\n")
for tag, perm in (("rev", lambda l: l[::-1]), ("rot", lambda l: l[-5:] + l[:-5])):
    reorder(f"{inp}/bath_DEERtest_12spin.txt", f"{w}/{tag}_bath.txt", perm)
    reorder(f"{inp}/avaax_DEERtest", f"{w}/{tag}_avaax", perm)
    reorder(f"{inp}/state_1", f"{w}/{tag}_state_1", perm)
O2 = {"order": 2, "configuration": {"s_range": [1, 1], "bath_state_file": [f"{inp}/state_"]}}
def files(tag):
    return {"configuration": {"bathfile": [f"{w}/{tag}_bath.txt"], "bath_state_file": [f"{w}/{tag}_state_"],
                              "defect_axis_file": [f"{w}/{tag}_avaax"]}}
cases["o2"]         = merge(merge({}, cases["onres"]), O2)
cases["o2_rev"]     = merge(merge({}, cases["o2"]), files("rev"))
cases["meg_o2"]     = merge(merge({}, cases["meg_onres"]), O2)
cases["meg_o2_rev"] = merge(merge({}, cases["meg_o2"]), files("rev"))
cases["o3"]         = merge(merge({}, cases["o2"]), {"order": 3})
cases["o3_rot"]     = merge(merge({}, cases["o3"]), files("rot"))
# the gCCE expm propagator with the pump (README: propagator eigen/expm), and cluster_mode pcce at order 2
# (multi-spin pcce clusters with bath pulses; partition_size 1, the base seed 3528 fixes the partition search)
cases["onres_expm"] = merge(merge({}, cases["onres"]), {"propagator": "expm"})
cases["pg_o2"]      = merge(merge({}, cases["o2"]), PC1)
# A67: the tones in kHz and the default bandwidth written out (10 rad*kHz = 10/(2 pi 1000) MHz) are the onres run; a tone
# 28 rad*kHz (4.5 kHz) off the line drives it with a 5 kHz bandwidth and not with the default; a bandwidth <= 0 is refused
cases["units_onres"] = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"energy_shift": {"value": [79975.539], "unit": "kHz"},
                                                                          "bandwidth": 0.0015915494309189536}}})
cases["bw_wide"]     = merge(pump([79.98], HAHN), {"pulse": {"bath": {"bandwidth": {"value": 5, "unit": "kHz"}}}})
cases["bw_default"]  = pump([79.98], HAHN)
cases["bad_bw"]      = merge(merge({}, cases["onres"]), {"pulse": {"bath": {"bandwidth": 0}}})
for name, extra in cases.items():
    d = merge(copy.deepcopy(base), extra); d["outfile"] = f"{w}/{name}"
    json.dump(nest(d), open(f"{w}/{name}.json", "w"), indent=1)
PY
for c in nopulse onres offres multi_hit multi_miss nearpulse; do
    run_ccex "$BIN_DEV" "$W/$c.json" "$W/proc_$c" || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
done
for c in path_matrix path_matrix_recursive path_vector path_recursive path_nopump; do   # bath_state_file and s_range [1, 1] in the input (also -s, -N 1 until A63)
    run_ccex "$BIN_DEV" "$W/$c.json" "$W/proc_$c" || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
done
# megCCE and cluster_mode pcce with bath pulses; the pure-state runs at order 2 and 3 and their reordered
# baths (the bath_state_file key, also -s until A63, and the log must name the state file read)
for c in meg_onres meg_nearpulse multi_hit_rev pg_onres pmeg_onres onres_expm; do
    run_ccex "$BIN_DEV" "$W/$c.json" "$W/proc_$c" || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
done
for c in o2 o2_rev meg_o2 meg_o2_rev o3 o3_rot pg_o2; do
    st="$HERE/inputs/state_"
    case "$c" in *_rev) st="$W/rev_state_" ;; *_rot) st="$W/rot_state_" ;; esac
    run_ccex "$BIN_DEV" "$W/$c.json" "$W/proc_$c" || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
    grep -qF "Read StateFile : ${st}1" "$W/proc_$c" || { echo "  FAIL  $c : did not read ${st}1"; FAIL=$((FAIL+1)); }
done
if python3 "$HERE/check_deer.py" "$W" "$HERE/inputs"; then PASS=$((PASS+1)); else FAIL=$((FAIL+1)); fi
# a pump whose tones reach no spin of pulse.bath.defect runs (a control), with one warning
for c in offres multi_miss; do
    if grep -qF "Warning: pulse.bath drives no spin" "$W/proc_$c"; then echo "  PASS  $c : off-resonant pump warned"; PASS=$((PASS+1))
    else echo "  FAIL  $c : no off-resonant warning"; FAIL=$((FAIL+1)); fi
done
# refused <case> <message fragment> <label>: the run stops with that message, not later with a kron size error
refused(){
    run_ccex "$BIN_DEV" "$W/$1.json" "$W/proc_$1"; local rc=$?
    if [ $rc -ne 0 ] && grep -qF "$2" "$W/proc_$1"; then echo "  PASS  $3"; PASS=$((PASS+1))
    else echo "  FAIL  $3 (exit $rc)"; grep -m1 -i error "$W/proc_$1"; FAIL=$((FAIL+1)); fi
}
# levels that are not two levels of the driven spins' S (A94: S is theirs, pulse.bath.spin was removed)
refused bad_bspin    "pulse.bath.levels [1, 0] must be two different levels of spin 0.5" "pulse.bath.levels of another spin refused"
refused bad_defect   'pulse.bath.defect "P2" names no bath spin' "pulse.bath.defect naming no spin refused"
refused bad_ms_equal "pulse.bath.levels [0.5, 0.5] must be two different levels of spin 0.5" "pulse.bath.levels equal refused"
refused bad_ms_range "pulse.bath.levels [1.5, -0.5] must be two different levels of spin 0.5" "pulse.bath.levels |m| > spin refused"
# A67 pulse.bath.bandwidth and the unit form of energy_shift (bytes against onres)
for c in units_onres bw_wide bw_default; do
    run_ccex "$BIN_DEV" "$W/$c.json" "$W/proc_$c" || { echo "  FAIL  $c : CCEX exit"; FAIL=$((FAIL+1)); }
done
for c in units_onres bw_wide; do
    if [ -s "$W/onres_noDiv" ] && cmp -s "$W/onres_noDiv" "$W/${c}_noDiv" && cmp -s "$W/onres_wiDiv" "$W/${c}_wiDiv"; then
        echo "  PASS  $c == onres (bytes)"; PASS=$((PASS+1))
    else echo "  FAIL  $c differs from onres (or no output)"; FAIL=$((FAIL+1)); fi
done
if grep -qF "Warning: pulse.bath drives no spin" "$W/proc_bw_default"; then echo "  PASS  bw_default : 4.5 kHz off the line, not driven"; PASS=$((PASS+1))
else echo "  FAIL  bw_default : no off-resonant warning"; FAIL=$((FAIL+1)); fi
refused bad_bw "pulse.bath.bandwidth must be a positive number" "pulse.bath.bandwidth 0 refused"
echo " DEER checks : PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
