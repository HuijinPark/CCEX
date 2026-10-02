#!/usr/bin/env bash
#
# Axis (1) : binary equivalence.
# ---------------------------------------------------------------------------
# The SAME legacy input, with no new key anywhere in it, run under the v1.1.0
# binary and under the axis-align binary. The two must produce byte-identical
# raw coherence output.
#
# This is the check that makes the rest of the suite mean anything. The new
# input forms cannot be run under v1.1.0 at all -- v1.1.0's parser rejects them
# -- so their agreement with v1.1.0 is established transitively:
#
#     legacy(v1.1.0) == legacy(dev)      <- here
#     legacy(dev)    == new-form(dev)    <- axis2.sh
#     => new-form(dev) == legacy(v1.1.0)
#
# so a failure here invalidates every later conclusion, and the suite stops.
#
# Usage:  bash test/equivalence/axis1.sh [system ...]
# Needs:  bin/main.out built, and the preserved v1.1.0 binary (see lib.sh).
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

[ -x "$BIN_DEV"  ] || { echo "ERROR: dev binary missing: $BIN_DEV" >&2; exit 1; }
[ -x "$BIN_V110" ] || { echo "ERROR: v1.1.0 binary missing: $BIN_V110" >&2; exit 1; }

echo "dev   : $BIN_DEV"
echo "        $(md5sum "$BIN_DEV" | cut -c1-32)"
echo "v1.1.0: $BIN_V110"
echo "        $(md5sum "$BIN_V110" | cut -c1-32)"
echo

# --- ccein generation -------------------------------------------------------
# Each system writes ONE legacy ccein; the two runs differ only in the binary and
# in "outfile", which is patched per run so the outputs sit side by side.
gen(){  # gen <sysdir> <python body on stdin>
    local sys="$1"; mkdir -p "$WORK/$sys"; python3 - "$WORK/$sys" "$REPO" "$FIX" "$REG" "$EX"
}

# regression case -> absolute-path ccein, config 1 only (byte identity does not
# need the whole ensemble)
gen_reg(){  # gen_reg <sysname> <case> <tag>
    local sys="$1" case="$2" tag="$3"
    mkdir -p "$WORK/$sys"
    python3 - "$WORK/$sys" "$REG/$case" "$tag" <<'PY'
import json, sys, os, glob
work, casedir, tag = sys.argv[1:4]
d = json.load(open(f"{casedir}/inputfiles/ccein_{tag}.json"))
for k in ("qubitfile","gyrofile","avaaxfile","statefile","exstatefile"):
    if k in d and isinstance(d[k], str):
        d[k] = os.path.abspath(os.path.join(casedir, d[k]))
baths = sorted(glob.glob(f"{casedir}/inputfiles/bath_*_1"))
d["bathfile"] = [os.path.abspath(baths[0])]
d["outfile"] = "./OUT"
json.dump(d, open(f"{work}/legacy.json","w"), indent=4)
print(f"    {os.path.basename(casedir)}/{tag}: bath={os.path.basename(baths[0])} "
      f"order={d.get('order')} rbath={d.get('rbath')} nstep={d.get('nstep')} nstate={d.get('nstate')}")
PY
}

run_pair(){  # run_pair <sysname> [extra ccex args...]
    local sys="$1"; shift
    local w="$WORK/$sys"
    ( cd "$w"
      python3 -c "import json;d=json.load(open('legacy.json'));d['outfile']='./v110';json.dump(d,open('r_v110.json','w'),indent=4)"
      python3 -c "import json;d=json.load(open('legacy.json'));d['outfile']='./dev' ;json.dump(d,open('r_dev.json' ,'w'),indent=4)"
      run_ccex "$BIN_V110" r_v110.json proc_v110 "$@"; rc1=$?; w1=$CCEX_WALL
      run_ccex "$BIN_DEV"  r_dev.json  proc_dev  "$@"; rc2=$?; w2=$CCEX_WALL
      echo "    v1.1.0 exit=$rc1 wall=${w1}s | dev exit=$rc2 wall=${w2}s"
      [ $rc1 -eq 0 ] && [ $rc2 -eq 0 ] || { echo "  FAIL  $sys (non-zero exit)"; exit 9; }
    )
    [ $? -eq 9 ] && { FAIL=$((FAIL+1)); return; }
    # nstate>0 appends _state<n> to the name; take whichever pair exists
    local pv pd
    pv=$(ls "$WORK/$sys"/v110*_noDiv 2>/dev/null | head -1); pv=${pv%_noDiv}
    pd=$(ls "$WORK/$sys"/dev*_noDiv  2>/dev/null | head -1); pd=${pd%_noDiv}
    cmp_raw "$sys" "$pv" "$pd"
}

# --- the systems ------------------------------------------------------------
declare -A SYSTEMS=(
  [S1a]="NV natab, CCE-2, ensemble bath"
  [S1b]="NV natab, CCE-2, sampled bath"
  [S2a]="NV natab, gCCE-2"
  [S2b]="NV + P1 10ppm, CCE-2, Defect naddspin=0"
  [S2c]="NV + P1 10ppm, gCCE-2, Defect naddspin=0"
  [S3]="P1 bath, Defect naddspin=1 (14N), navaax=4, avaax+state files"
  [S4]="hBN VB, hf_readmode 3 + qd_readmode 2 (DFT tensor files, vertex boundary)"
  [S4b]="hBN VB, the same but with the MinDif/MaxDif axis-aligned boundary"
  [S5]="2 NV, gCCE, explicit intmap self/ZFS tensors"
)
ORDER=(S1a S1b S2a S2b S2c S3 S4 S4b S5)
WANT=("${@:-${ORDER[@]}}")

for sys in "${WANT[@]}"; do
    echo "=== $sys : ${SYSTEMS[$sys]:-?}"
    case "$sys" in
      S1a) gen_reg S1a Diamond_NV_natab_cce   nstate0 ; run_pair S1a ;;
      S1b) gen_reg S1b Diamond_NV_natab_cce   nstate1 ; run_pair S1b ;;
      S2a) gen_reg S2a Diamond_NV_natab_gcce  nstate0 ; run_pair S2a ;;
      S2b) gen_reg S2b Diamond_NV_P1_10ppm_cce  nstate0 ; run_pair S2b ;;
      S2c) gen_reg S2c Diamond_NV_P1_10ppm_gcce nstate0 ; run_pair S2c ;;

      S3)  mkdir -p "$WORK/S3"
           python3 - "$WORK/S3" "$REPO" "$FIX" <<'PY'
import json, re, sys
work, repo, fix = sys.argv[1:4]
src = open(f"{repo}/test/code_verification/1.3.FullP1SpinBath/cce.json").read()
src = re.sub(r'!\S*','',src); src = re.sub(r'#[^\n]*','',src)   # CCEX json dialect
d = json.loads(src)
F = f"{fix}/1.single_data"
d["qubitfile"]=f"{F}/bath_DiaP1_1ppm_Defect"; d["gyrofile"]=f"{F}/DiaP1_gyro"
d["bathfile"]=[f"{F}/bath_DiaP1_1ppm_1"];     d["avaaxfile"]=f"{F}/bathJT_DiaP1_1ppm_1"
d["statefile"]=f"{F}/state_DiaP1_1ppm_";      d["exstatefile"]=f"{F}/stateEx_DiaP1_1ppm_"
# Reduced from the original (rbath 850 / rdip 620 / nstep 300, 469 spins, sized for
# 24 ranks). Equivalence needs reproducibility, not convergence; rbath 360 keeps
# 45 spins and 518 pairs, which still exercises every Defect code path.
d["rbath"]=360; d["rdip"]=360; d["nstep"]=20; d["nstate"]=1
d["outfile"]="./OUT"
json.dump(d, open(f"{work}/legacy.json","w"), indent=4)
print(f"    P1 Defect: rbath={d['rbath']} nstep={d['nstep']} navaax={d['Defect'][0]['navaax']} naddspin={d['Defect'][0]['naddspin']}")
PY
           # statefile passed BOTH in the ccein and on the command line, on purpose:
           # a missing -s falls back to randomized states that still look plausible.
           run_pair S3 -s "$FIX/1.single_data/state_DiaP1_1ppm_" \
                       -S "$FIX/1.single_data/stateEx_DiaP1_1ppm_" -N 1 ;;

      S4|S4b)
           # S4  : Afile_vertex / Qfile_vertex -> CheckBD_vertex
           # S4b : Afile        / Qfile        -> CheckBD_Range, the axis-aligned
           #       MinDif/MaxDif box. Two different boundary tests over the same
           #       system, and only one of them is reachable per file.
           SUF=""; [ "$sys" = "S4" ] && SUF="_vertex"
           mkdir -p "$WORK/$sys"
           python3 - "$WORK/$sys" "$FIX" "$SUF" <<'PY'
import json, sys
work, fix, suf = sys.argv[1:4]
H=f"{fix}/4.hexagonal_data"
# rbath 20 -> 6 : 3897 spins -> 101, and 121 s -> 0.4 s on one core. The A/Q files
# cover this region comfortably, so both readers stay on their file path.
d = {"method":"CCE","quantity":"coherence",
     "qubitfile":f"{H}/defect","gyrofile":f"{fix}/2.ensemble_data/Gyro_h_BN",
     "bathfile":[f"{H}/bath_1"],
     "order":2,"bfield":30000,"rbath":6,"rdip":4,"deltat":0.002,"nstep":20,
     "nstate":0,"alphams":1,"betams":0,"npulse":1,
     # hf_ignore_oor=1: the MinDif/MaxDif box in Afile is wider than the set of atoms
     # the file lists, so a spin can sit inside the box with no row of its own. Without
     # this the range-boundary run stops there and never reaches a coherence to compare.
     "hf_readmode":3,"hf_tensorfile":f"{H}/Afile{suf}","hf_cutoff":0,
     "hf_ignore_oor":(0 if suf else 1),
     # The QD reader has no counterpart to hf_ignore_oor, and Qfile's box is wide in
     # the same way, so the range-boundary case runs without the quadrupole. S4 covers
     # the QD reader over the vertex boundary; the QD range boundary stays uncovered.
     "qd_readmode":(2 if suf else 0),
     "qd_tensorfile":f"{H}/Qfile{suf}",
     "qd_tensorfile_woqubit":f"{H}/Qfile{suf}",
     "outfile":"./OUT"}
json.dump(d, open(f"{work}/legacy.json","w"), indent=4)
print(f"    hBN VB: hf_readmode={d['hf_readmode']} qd_readmode={d['qd_readmode']} "
      f"rbath={d['rbath']} boundary={'vertex' if suf else 'MinDif/MaxDif'}")
PY
           run_pair "$sys" ;;

      S5)  mkdir -p "$WORK/S5"
           python3 - "$WORK/S5" "$EX" <<'PY'
import json, re, sys
work, ex = sys.argv[1:3]
E=f"{ex}/Diamond_multiNV_natab/NV2_13C2"
s=open(f"{E}/ccein_NV6_13C2_hahnecho.js").read()
s=re.sub(r'#[^\n]*','',s); s=re.sub(r'!\S*','',s)
d=json.loads(s)
# bathfile as a LIST. Given as a bare string it is read as zero bath files and the
# run segfaults with "Total 1-Cluster # : 0" -- v1.1.0 behaves the same way, so this
# is the example being stale, not anything this branch changed.
d["bathfile"]=[f"{E}/bath"]; d["statefile"]=f"{E}/state"; d["gyrofile"]=f"{E}/gyro_13C"
d["nstep"]=20; d["outfile"]="./OUT"
json.dump(d, open(f"{work}/legacy.json","w"), indent=4)
print(f"    2 NV: nqubit={d['Qubit']['nqubit']} intmap={len(d['Qubit']['intmap'])} method={d['method']}")
PY
           run_pair S5 -s "$EX/Diamond_multiNV_natab/NV2_13C2/state" -N 1 ;;
      *) echo "  unknown system: $sys"; FAIL=$((FAIL+1)) ;;
    esac
    echo
done

echo "==========================================="
echo " axis (1) binary equivalence : PASS=$PASS FAIL=$FAIL"
echo "==========================================="
[ "$FAIL" -eq 0 ]
