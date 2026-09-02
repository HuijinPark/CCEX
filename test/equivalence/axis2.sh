#!/usr/bin/env bash
#
# Axis (2) : input-form equivalence.
# ---------------------------------------------------------------------------
# For every new input form, the legacy way of saying the same thing is run beside
# it under the SAME binary. Axis (1) already established that the legacy side is
# byte-identical to v1.1.0, so agreement here carries the new form back to v1.1.0
# without ever having to hand v1.1.0 a key it cannot parse.
#
# gen_variants.py writes the variants and a manifest saying, per variant, what it
# is compared against and whether byte identity is required or a numeric
# comparison is the most that can be asked (a unit conversion in between).
#
# Usage:  bash test/equivalence/axis2.sh [system ...]
# Needs:  bin/main.out, and axis1.sh already run (it writes the legacy cceins).
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

ORDER=(S1a S3 S4 S5)
WANT=("${@:-${ORDER[@]}}")

declare -A EXTRA=(
  [S3]="-s $FIX/1.single_data/state_DiaP1_1ppm_ -S $FIX/1.single_data/stateEx_DiaP1_1ppm_ -N 1"
  [S5]="-s $EX/Diamond_multiNV_natab/NV2_13C2/state -N 1"
)

# out_prefix <workdir> <name> : nstate>0 appends _state<n>, so glob for it
out_prefix(){ local p; p=$(ls "$1/$2"*_noDiv 2>/dev/null | head -1); echo "${p%_noDiv}"; }

for sys in "${WANT[@]}"; do
    w="$WORK/$sys"
    echo "=== $sys"
    [ -f "$w/legacy.json" ] || { echo "  FAIL  $sys : no legacy.json (run axis1.sh first)"; FAIL=$((FAIL+1)); echo; continue; }

    python3 "$HERE/gen_variants.py" "$sys" "$w" || { FAIL=$((FAIL+1)); continue; }

    # run every variant once
    while IFS=$'\t' read -r name baseline mode desc; do
        ( cd "$w" && run_ccex "$BIN_DEV" "$name.json" "proc_$name" ${EXTRA[$sys]:-} )
        rc=$?
        if [ $rc -ne 0 ]; then
            echo "  FAIL  $name : exit=$rc"
            sed -n '/Error/{p;q}' "$w/proc_$name" | sed 's/^/        /'
            FAIL=$((FAIL+1))
        fi
    done < "$w/manifest.tsv"

    # compare
    while IFS=$'\t' read -r name baseline mode desc; do
        [ "$name" = "base" ] && continue
        [ "$name" = "$baseline" ] && continue
        a=$(out_prefix "$w" "$baseline"); b=$(out_prefix "$w" "$name")
        if [ -z "$a" ] || [ -z "$b" ]; then
            echo "  FAIL  $name : output missing"; FAIL=$((FAIL+1)); continue
        fi
        echo "  -- $name  [$mode vs $baseline]  $desc"
        if [ "$mode" = "byte" ]; then cmp_raw "$name" "$a" "$b"
        else                          maxdiff "$name" "$a" "$b"; fi
    done < "$w/manifest.tsv"
    echo
done

echo "==========================================="
echo " axis (2) input-form equivalence : PASS=$PASS FAIL=$FAIL"
echo " (num-mode rows report a deviation and are read, not auto-judged)"
echo "==========================================="
[ "$FAIL" -eq 0 ]
