#!/usr/bin/env bash
#
# The configuration block of A64 (2026-10-15): path, k_range and {k}, and one state or axis file per bath file, on the
# P1 10 ppm cCCE2 system of inputs/Diamond_NV_P1_10ppm_cce (one bath state, random defect axes and nuclear states unless
# files are given). Two bath configurations: cfg_1 (the bath file) and cfg_2 (every spin moved by 1 A in x).
#   K1  k_range [1, 2] (configuration.path data, bathfile ["cfg_{k}"]): <outfile>_k1 and _k2, each byte for byte the run
#       of that file alone - random axes and states drawn as in a run of its own
#   K2  k_range [2, 2]: byte for byte the run of cfg_2 alone
#   K3  k_range [1, 2] with bath_state_file, defect_axis_file and defect_state_file per k: each k byte for byte its run
#       alone, and the two differ
#   K4  two bath files with a state, axis and defect state file each == the two merged into one file (the lists are read
#       in the order of the bath files)
#   K4  the log names the merge: "Bath = data/part_a (n spins) + data/part_b (m spins)"
#   K8  "" in the lists: that bath file has no state, axis or defect state file - the first file's are read, the second
#       file's spins drawn from the seed; "" for every bath file is the key absent (byte for byte)
#   K9  k_range [1, 2] without a seed (A73): each configuration has its own generated seed, printed per k in the log, and
#       configuration k is byte for byte the run of its file alone with that seed (with a seed, K1: the same seed for all)
#   K10 a file of a per-bath-file list that does not open (A73): a warning naming it, and only that bath file's spins drawn
#       from the seed - byte for byte K8, where that entry is ""
#   K5  {{k}} is a literal {k} with k_range: a bath file named lit_{k} is read
#   K6  without k_range {k} is literal: a bath file named cfg_{k} is read as it is
#   K7  cluster_mode pcce with the partition cache, cache.file "cache_{k}": two cache files, each k byte for byte its run
#       alone with a fresh cache
#   R1-R9 refused before anything is calculated: k_range without {k}, a missing file of one k, k_range [2, 1], a list of
#       the wrong length, a string where a list is, an enabled cache without {k}, the old top-level bathfile and nstate,
#       Defect.axis_file (configuration.defect_axis_file since A64)
#   R10 two bath files with spins at the same position (the same file twice): refused, naming both files (A73: infinite
#       coupling, nan)
#   W10 two bath files with spins 0.3 A apart (cfg_1 and cfg_1 shifted by 0.3 A in x; A73, user 2026-10-16): calculated,
#       with a warning naming the files and positions
#
# Usage:  bash test/code_verification/check_input_forms/configuration.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/configuration"
rm -rf "$W"; mkdir -p "$W/data" "$W/out"
cp -r "$REG/Diamond_NV_P1_10ppm_cce/inputfiles" "$W/"
cd "$W"
ok(){ echo "  PASS  $1"; PASS=$((PASS+1)); }
no(){ echo "  FAIL  $1"; FAIL=$((FAIL+1)); }
python3 - "$REPO" <<'PY'
import json, random, sys; from input_layout import flat, nest
sys.path.insert(0, sys.argv[1] + "/example/tools/InputKeys")
from convert_input_keys import strip_comments
j = flat(json.loads(strip_comments(open("inputfiles/ccein_nstate1.json").read())))
j.pop("outfile", None)
j["nstep"] = 20
lines = open("inputfiles/bath_DiaP1_10ppm_1").read().splitlines()
head, body = lines[0], [l for l in lines[1:] if l.strip()]
def shifted(l):
    x, y, z, name = l.split()
    return f"{float(x) + 1.0:.5f} {y} {z} {name}"
cfg = {1: body, 2: [shifted(l) for l in body]}
def write(path, rows, bath=False):
    # a bath file starts with its line "N 0 0 0"; a state or axis file with one token, the count
    open(path, "w").write("\n".join([head if bath else str(len(rows))] + rows) + "\n")
rng = random.Random(7)
for k, rows in cfg.items():
    write(f"data/cfg_{k}", rows, True)
    write(f"data/st{k}_1", [rng.choice(["0.5", "-0.5"]) for _ in rows])        # bath state 1, one line per bath line
    write(f"data/ax_{k}", [str(rng.randint(1, 3)) for _ in rows])             # P1 axes 1 .. 3 (n_configurations 4)
    write(f"data/ds{k}_1", [rng.choice(["0.5", "-0.5"]) for _ in rows])       # the 15N of each P1, state 1
# K4: the first half and the second half of cfg_1 as two bath files, and the same merged
half = len(body) // 2
for name, rows in (("a", body[:half]), ("b", body[half:])):
    write(f"data/part_{name}", rows, True)
    write(f"data/stp_{name}_1", [rng.choice(["0.5", "-0.5"]) for _ in rows])
    write(f"data/axp_{name}", [str(rng.randint(1, 3)) for _ in rows])
    write(f"data/dsp_{name}_1", [rng.choice(["0.5", "-0.5"]) for _ in rows])
for kind, suffix in (("part", ""), ("stp", "_1"), ("axp", ""), ("dsp", "_1")):
    rows = []
    for name in ("a", "b"):
        rows += open(f"data/{kind}_{name}{suffix}").read().splitlines()[1:]
    write(f"data/{kind}_merged{suffix}", rows, kind == "part")
open("data/lit_{k}", "w").write(open("data/cfg_2").read())
open("data/cfg_{k}", "w").write(open("data/cfg_2").read())

def inp(name, conf, defect=None, extra=None):
    d = json.loads(json.dumps(j))
    c = d["configuration"]
    for key in ("bathfile",):
        c.pop(key)
    c.update(conf)
    if "path" in conf:                     # the gyrofile is in the block too: relative to configuration.path
        c["gyrofile"] = "../inputfiles/gyro_P1"
    if defect:
        c.update(defect)                   # the Defect files are the configuration's (A64)
    d.update(extra or {})
    d["outfile"] = f"out/{name}"
    open(f"{name}.json", "w").write(json.dumps(nest(d), indent=1))

inp("k1", {"path": "data", "k_range": [1, 2], "bathfile": ["cfg_{k}"]})
inp("plain1", {"bathfile": ["data/cfg_1"]})
inp("plain2", {"path": "data", "bathfile": ["cfg_2"]})
inp("k2", {"path": "data", "k_range": [2, 2], "bathfile": ["cfg_{k}"]})
inp("k3", {"path": "data", "k_range": [1, 2], "bathfile": ["cfg_{k}"], "bath_state_file": ["st{k}_"]},
    {"defect_axis_file": ["ax_{k}"], "defect_state_file": ["ds{k}_"]})
for k in (1, 2):
    inp(f"k3plain{k}", {"bathfile": [f"data/cfg_{k}"], "bath_state_file": [f"data/st{k}_"]},
        {"defect_axis_file": [f"data/ax_{k}"], "defect_state_file": [f"data/ds{k}_"]})
inp("k4two", {"path": "data", "bathfile": ["part_a", "part_b"], "bath_state_file": ["stp_a_", "stp_b_"],
              "bath_shift": [[0, 0, 0], [0, 0, 0]]},
    {"defect_axis_file": ["axp_a", "axp_b"], "defect_state_file": ["dsp_a_", "dsp_b_"]})
inp("k4one", {"path": "data", "bathfile": ["part_merged"], "bath_state_file": ["stp_merged_"]},
    {"defect_axis_file": ["axp_merged"], "defect_state_file": ["dsp_merged_"]})
inp("k8", {"path": "data", "bathfile": ["part_a", "part_b"], "bath_state_file": ["stp_a_", ""],
           "bath_shift": [[0, 0, 0], [0, 0, 0]]}, {"defect_axis_file": ["axp_a", ""], "defect_state_file": ["dsp_a_", ""]})
inp("k8empty", {"path": "data", "bathfile": ["part_a", "part_b"], "bath_state_file": ["", ""],
                "bath_shift": [[0, 0, 0], [0, 0, 0]]}, {"defect_axis_file": ["", ""], "defect_state_file": ["", ""]})
inp("k8absent", {"path": "data", "bathfile": ["part_a", "part_b"], "bath_shift": [[0, 0, 0], [0, 0, 0]]})
inp("k10", {"path": "data", "bathfile": ["part_a", "part_b"], "bath_state_file": ["stp_a_", "missing_st_"],
            "bath_shift": [[0, 0, 0], [0, 0, 0]]},
    {"defect_axis_file": ["axp_a", "missing_ax"], "defect_state_file": ["dsp_a_", "missing_ds_"]})
inp("k9", {"path": "data", "k_range": [1, 2], "bathfile": ["cfg_{k}"]})
d = json.load(open("k9.json")); d["configuration"].pop("seed"); json.dump(d, open("k9.json", "w"), indent=1)
inp("r10", {"path": "data", "bathfile": ["cfg_1", "cfg_1"], "bath_shift": [[0, 0, 0], [0, 0, 0]]})
inp("w10", {"path": "data", "bathfile": ["cfg_1", "cfg_1"], "bath_shift": [[0, 0, 0], [0.3, 0, 0]]})
inp("k5", {"path": "data", "k_range": [1, 1], "bathfile": ["lit_{{k}}", "cfg_{k}"], "bath_shift": [[0, 0, 0], [0, 0, 0]]})
inp("k6", {"path": "data", "bathfile": ["cfg_{k}"]})
pcce = {"method": "cCCE", "cluster_mode": "pcce"}
inp("k7", {"path": "data", "k_range": [1, 2], "bathfile": ["cfg_{k}"]}, None,
    dict(pcce, pcce={"partition_size": 2, "cache": {"mode": "readwrite", "file": "out/cache_{k}"}}))
for k in (1, 2):
    inp(f"k7plain{k}", {"bathfile": [f"data/cfg_{k}"]}, None,
        dict(pcce, pcce={"partition_size": 2, "cache": {"mode": "readwrite", "file": f"out/fresh_cache_{k}"}}))
# refusals
inp("r1", {"path": "data", "k_range": [1, 2], "bathfile": ["cfg_1"]})
inp("r2", {"path": "data", "k_range": [1, 3], "bathfile": ["cfg_{k}"]})
inp("r3", {"path": "data", "k_range": [2, 1], "bathfile": ["cfg_{k}"]})
inp("r4", {"path": "data", "bathfile": ["cfg_1"], "bath_state_file": ["st1_", "st2_"]})
inp("r5", {"path": "data", "bathfile": ["cfg_1"], "defect_axis_file": "ax_1"}, None, {"__keep_old__": ["defect_axis_file"]})
inp("r6", {"path": "data", "k_range": [1, 2], "bathfile": ["cfg_{k}"]}, None,
    dict(pcce, pcce={"partition_size": 2, "cache": {"mode": "readwrite", "file": "out/cache"}}))
d = json.loads(json.dumps(j)); d["bathfile"] = d["configuration"].pop("bathfile"); d["outfile"] = "out/r7"; d["__keep_old__"] = ["bathfile"]
open("r7.json", "w").write(json.dumps(nest(d), indent=1))
d = json.loads(json.dumps(j)); d["nstate"] = d["configuration"].pop("nstate"); d["outfile"] = "out/r8"; d["__keep_old__"] = ["nstate"]
open("r8.json", "w").write(json.dumps(nest(d), indent=1))
d = json.loads(json.dumps(j)); d["Defect"]["axis_file"] = ["data/ax_1"]; d["outfile"] = "out/r9"; d["__keep_old__"] = ["axis_file"]
open("r9.json", "w").write(json.dumps(nest(d), indent=1))
PY

for n in k1 plain1 plain2 k2 k3 k3plain1 k3plain2 k4two k4one k8 k8empty k8absent k10 k9 w10 k5 k6 k7 k7plain1 k7plain2; do
    run_ccex "$BIN_DEV" $n.json $n.log || no "$n exited $? (see $W/$n.log)"
done
cmp_raw "K1 k_range: configuration 1 == its run alone" out/k1_k1_state1 out/plain1_state1
cmp_raw "K1 k_range: configuration 2 == its run alone" out/k1_k2_state1 out/plain2_state1
cmp_raw "K2 k_range [2, 2] == the run of cfg_2 alone" out/k2_k2_state1 out/plain2_state1
cmp_raw "K3 state, axis and defect state files per k: configuration 1 == its run alone" out/k3_k1_state1 out/k3plain1_state1
cmp_raw "K3 state, axis and defect state files per k: configuration 2 == its run alone" out/k3_k2_state1 out/k3plain2_state1
if cmp -s out/k3_k1_state1_wiDiv out/k3_k2_state1_wiDiv; then no "K3 the two configurations gave the same curve"; else ok "K3 the two configurations differ"; fi
# the files were read, not drawn at random (a random draw with the same seed would agree too)
read_all(){ local log=$1; shift; for f in "$@"; do grep -qF "$f" $log || return 1; done; }
if read_all k3.log "Read StateFile : data/st2_1" "Read AvaaxFile : data/ax_2" "Read ExstateFile : data/ds2_1" \
   && read_all k4two.log "Read StateFile : data/stp_b_1" "Read AvaaxFile : data/axp_b" "Read ExstateFile : data/dsp_b_1" \
   && ! grep -q "randomly generated" k3.log k4two.log
then ok "K3/K4 every state, axis and defect state file was read"; else no "K3/K4 a file was not read (see k3.log, k4two.log)"; fi
cmp_raw "K4 two bath files with their state, axis and defect state lists == the files merged" out/k4two_state1 out/k4one_state1
if grep -q "Bath = data/part_a ([0-9]* spins) + data/part_b ([0-9]* spins)" k4two.log
then ok "K4 the log names the merged bath files"; else no "K4 the merge line (see $W/k4two.log)"; fi
if read_all k8.log "Read StateFile : data/stp_a_1" "Read AvaaxFile : data/axp_a" "Read ExstateFile : data/dsp_a_1" \
   "bath states of the bath files without a state file" "paxes of the bath files without an axis file" \
   "subbath states of the bath files without a state file" && [ -s out/k8_state1_wiDiv ]
then ok "K8 \"\" in a list: the first file's files read, the second file's spins drawn"; else no "K8 (see $W/k8.log)"; fi
cmp_raw "K8 \"\" for every bath file == the keys absent" out/k8empty_state1 out/k8absent_state1
if grep -qF "Cannot open the state file (data/missing_st_1)" k10.log && grep -qF "Cannot open the avaax file (data/missing_ax)" k10.log \
   && grep -qF "Cannot open the exstate file (data/missing_ds_1)" k10.log && grep -qF "Read StateFile : data/stp_a_1" k10.log
then ok "K10 a missing file of one bath file: warned by name, the other bath file's files read"; else no "K10 (see $W/k10.log)"; fi
cmp_raw "K10 a missing file of one bath file == \"\" there (K8)" out/k10_state1 out/k8_state1
# K9: the seeds of the two configurations, from the log; then each configuration alone with its seed
s1=$(sed -n 's/^\[seed\] configuration k=1: seed \(-\{0,1\}[0-9]*\) .*/\1/p' k9.log)
s2=$(sed -n 's/^\[seed\] configuration k=2: seed \(-\{0,1\}[0-9]*\) .*/\1/p' k9.log)
if [ -n "$s1" ] && [ -n "$s2" ] && [ "$s1" != "$s2" ]; then ok "K9 k_range without a seed: a seed per k in the log ($s1, $s2)"
else no "K9 the per-k seed lines (see $W/k9.log)"; fi
for k in 1 2; do
    seed=$([ $k = 1 ] && echo "$s1" || echo "$s2")
    python3 - $k "${seed:-1}" <<'PY'
import json, sys
k, seed = sys.argv[1], int(sys.argv[2])
d = json.load(open("k9.json")); c = d["configuration"]
c.pop("k_range"); c["bathfile"] = [f"cfg_{k}"]; c["seed"] = seed; d["output"]["file"] = f"out/k9plain{k}"
json.dump(d, open(f"k9plain{k}.json", "w"), indent=1)
PY
    run_ccex "$BIN_DEV" k9plain$k.json k9plain$k.log || no "k9plain$k exited $? (see $W/k9plain$k.log)"
    cmp_raw "K9 configuration $k without a seed == its run alone with the seed of the log" out/k9_k${k}_state1 out/k9plain${k}_state1
done
if cmp -s out/k9_k1_state1_wiDiv out/k1_k1_state1_wiDiv; then no "K9 configuration 1 drew the same as seed 1 (seed not generated?)"; fi
if grep -qF "Warning(readBathFile): two bath spins of different bath files are 0.3000 Angstrom apart" w10.log \
   && grep -qF "in data/cfg_1 and" w10.log && grep -q "two spins at one site" w10.log && [ -s out/w10_state1_wiDiv ] \
   && ! grep -q nan out/w10_state1_wiDiv
then ok "W10 spins of two bath files 0.3 A apart: a warning naming files and positions, and calculated"; else no "W10 (see $W/w10.log)"; fi
if grep -q "lit_{k}" k5.log && [ -s out/k5_k1_state1_wiDiv ]; then ok "K5 {{k}} is a literal {k} with k_range"; else no "K5 (see $W/k5.log)"; fi
if grep -q "cfg_{k}" k6.log && [ -s out/k6_state1_wiDiv ]; then ok "K6 without k_range {k} is read as written"; else no "K6 (see $W/k6.log)"; fi
if [ -s out/cache_1 ] && [ -s out/cache_2 ] && ! cmp -s out/cache_1 out/cache_2; then ok "K7 two configurations, two cache files"
else no "K7 the cache files cache_1 and cache_2 (see $W/k7.log)"; fi
cmp_raw "K7 pcce with a cache per k: configuration 1 == its run alone" out/k7_k1_state1 out/k7plain1_state1
cmp_raw "K7 pcce with a cache per k: configuration 2 == its run alone" out/k7_k2_state1 out/k7plain2_state1

refused(){   # refused <name> <label> <message fragment>
    if run_ccex "$BIN_DEV" $1.json $1.log; then no "$2: accepted"; return; fi
    if ! grep -qF -- "$3" $1.log; then no "$2: refused, but not with \"$3\" (see $W/$1.log)"; return; fi
    if ls out/$1* > /dev/null 2>&1; then no "$2: refused after writing outputs"; return; fi
    if grep -q "Clusterize" $1.log; then no "$2: refused after the calculation started"; return; fi
    ok "$2"
}
refused r1 "R1 k_range without {k}" "configuration.k_range is given, but no file name has {k}"
refused r2 "R2 a missing file of one k, before any configuration runs" "configuration k = 3: cannot read data/cfg_3"
refused r3 "R3 k_range [2, 1]" "configuration.k_range must be [first, last]"
refused r4 "R4 a list of the wrong length" "configuration.bath_state_file must be a list of one name per bath file (1"
refused r5 "R5 a string where the list is" "configuration.defect_axis_file must be a list of one name per bath file"
refused r6 "R6 an enabled pcce cache without {k}" "pcce.cache.file \"out/cache\" needs {k}"
refused r7 "R7 the top-level bathfile" "input key \"bathfile\" is now \"configuration.bathfile\""
refused r8 "R8 the top-level nstate" "input key \"nstate\" is now \"configuration.nstate\""
refused r10 "R10 two bath files with spins at the same position" "two bath spins of different bath files are at the same position"
refused r9 "R9 Defect.axis_file" "input key \"Defect.axis_file\" is now \"configuration.defect_axis_file\""
echo "configuration: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
