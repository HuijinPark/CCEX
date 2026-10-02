#!/usr/bin/env bash
#
# DFT tensor-path legacy baseline
# ---------------------------------------------------------------------------
# The four regression tests under test/regression/ all run hf_readmode = 0 and
# qd_readmode = 0, so the DFT tensor readers -- the code the rotation feature is
# about to start transforming -- have NO byte-level baseline anywhere in the
# repository. The committed references under test/code_verification/ do not
# serve either: they already disagree with what the current build produces.
#
# So the baseline has to be made, not found. This script runs the same four
# inputs through two binaries and compares them byte for byte:
#
#     hf_readmode 1   Fermi contact + point dipole
#     hf_readmode 2   DFT dipolar tensor
#     hf_readmode 3   Fermi contact + DFT dipole
#     qd_readmode 2   DFT quadrupole (QuadFileWiDefect)
#
# None of them uses any of the new options, so every one must stay identical
# for as long as the rotation work goes on.
#
# The inputs are the h-BN V_B set already committed under
# test/code_verification/CCE_Reprod/Bath_Data/4.hexagonal_data/ -- real DFT
# tensors, and genuinely non-diagonal, which is what makes them worth using.
# Nothing there is written to; only read.
#
# Usage:
#   dft_legacy_check.sh <ref-binary>                 make the baseline only
#   dft_legacy_check.sh <ref-binary> <test-binary>   make it and compare
#
# Results and checksums go to _work/dft/, which run.sh wipes.
# ---------------------------------------------------------------------------

set -uo pipefail

# Resolve the binary arguments BEFORE cd'ing: they are usually given relative to
# the repository root, which is not where this script runs from.
abspath(){ [ -z "${1:-}" ] && return 0; ( cd "$(dirname "$1")" 2>/dev/null && printf '%s/%s\n' "$(pwd)" "$(basename "$1")" ) || printf '%s\n' "$1"; }
set -- "$(abspath "${1:-}")" "$(abspath "${2:-}")"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
REPO="$(cd "$HERE/../.." && pwd)"
DATA="$REPO/test/code_verification/CCE_Reprod/Bath_Data"
HEX="$DATA/4.hexagonal_data"

MPIRUN="${CCEX_MPI:-/opt/intel/oneapi/mpi/2021.18/bin/mpirun}"
export I_MPI_FABRICS="${I_MPI_FABRICS:-shm}"

REF_BIN="${1:-}"
TEST_BIN="${2:-}"

[ -n "$REF_BIN" ] || { echo "usage: $0 <ref-binary> [<test-binary>]" >&2; exit 1; }
[ -x "$REF_BIN" ] || { echo "ERROR: not executable: $REF_BIN" >&2; exit 1; }
[ -d "$HEX" ]     || { echo "ERROR: h-BN tensor data not found: $HEX" >&2; exit 1; }

OUT="$HERE/_work/dft"
rm -rf "$OUT"; mkdir -p "$OUT"

# --- the four cases ----------------------------------------------------------
# order 1 and nstep 20 on purpose: this is a check on the TENSOR READERS, and a
# deeper cluster expansion would only add runtime without touching more of them.
write_ccein(){
    local tag="$1" hf="$2" qd="$3" dest="$4" outfile="$5"
    python3 - "$dest" "$hf" "$qd" "$outfile" "$HEX" <<'PY'
import json, sys
dest, hf, qd, outfile, hex_dir = sys.argv[1:6]
cfg = {
    "method"    : "CCE",
    "quantity"  : "coherence",
    "qubitfile" : hex_dir + "/defect",
    "gyrofile"  : hex_dir + "/../2.ensemble_data/Gyro_h_BN",
    "bathfile"  : [hex_dir + "/bath_1"],
    "order"     : 1,
    "bfield"    : 30000,
    "rbath"     : 15,
    "rdip"      : 10,
    "deltat"    : 0.002,
    "nstep"     : 20,
    "npulse"    : 1,
    "nstate"    : 0,
    "seed"      : 1,
    "qalphams"  : 1,
    "qbetams"   : 0,
    "hf_readmode": int(hf),
    "qd_readmode": int(qd),
    "outfile"   : outfile,
}
if int(hf) != 0:
    cfg["hf_tensorfile"] = hex_dir + "/Afile_vertex"
    cfg["hf_cutoff"] = 0
    cfg["hf_ignore_oor"] = 0
if int(qd) != 0:
    cfg["qd_tensorfile"] = hex_dir + "/Qfile_vertex"
    cfg["qd_tensorfile_woqubit"] = hex_dir + "/Qfile_vertex"
json.dump(cfg, open(dest, "w"), indent=4)
PY
}

CASES="hf1:1:0 hf2:2:0 hf3:3:0 qd2:0:2"

run_all(){
    local bin="$1" label="$2"
    mkdir -p "$OUT/$label"
    for c in $CASES; do
        local tag=${c%%:*} rest=${c#*:}
        local hf=${rest%%:*} qd=${rest##*:}
        write_ccein "$tag" "$hf" "$qd" "$OUT/$label/$tag.json" "$OUT/$label/$tag"
        "$MPIRUN" -n 1 "$bin" -f "$OUT/$label/$tag.json" > "$OUT/$label/$tag.log" 2>&1
        printf "  %-4s hf=%s qd=%s  rc=%d\n" "$tag" "$hf" "$qd" $?
    done
}

echo "=== reference binary : $REF_BIN"
echo "    md5 $(md5sum "$REF_BIN" | cut -d' ' -f1)"
run_all "$REF_BIN" ref
( cd "$OUT/ref" && md5sum *_noDiv *_wiDiv 2>/dev/null | sort ) > "$OUT/baseline.md5"
echo
echo "  baseline checksums -> $OUT/baseline.md5"
sed 's/^/    /' "$OUT/baseline.md5"

[ -n "$TEST_BIN" ] || exit 0

echo
echo "=== test binary : $TEST_BIN"
echo "    md5 $(md5sum "$TEST_BIN" | cut -d' ' -f1)"
run_all "$TEST_BIN" new

echo
echo "=== byte comparison ==="
nfail=0
for c in $CASES; do
    tag=${c%%:*}
    for suf in _noDiv _wiDiv; do
        if cmp -s "$OUT/ref/$tag$suf" "$OUT/new/$tag$suf"; then
            echo "  IDENTICAL  $tag$suf"
        else
            echo "  DIFFERS    $tag$suf"
            nfail=$((nfail+1))
        fi
    done
done

echo
[ $nfail -eq 0 ] && echo "DFT legacy baseline holds." || echo "DFT legacy baseline BROKEN ($nfail file(s))."
exit $nfail
