#!/usr/bin/env bash
#
# savemode cluster_density_matrix (docs/user/InputFile_Output.md, fixed 2026-10-13; "allfull" until the names of 2026-10-13,
# and the files and messages below keep the short tags: all = cluster_coherence, allfull = cluster_density_matrix): it
# writes every file savemode cluster_coherence writes, byte for
# byte (the result of the expansion and each cluster's file), and with gCCE also one full density-matrix file per cluster,
# named after the cluster (<outfile>_<bath indices>_full_state1_*, the zeroth cluster <outfile>_full_state1_*), no two
# alike; with cCCE nothing more. 1 and 3 ranks write the same files (no collision across ranks). natab gCCE2 and cCCE2,
# one sampled state, inputs/Diamond_NV_natab_{gcce,cce}.
#
# Usage:  bash test/code_verification/check_input_forms/allfull.sh
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib.sh"

W="$WORK/allfull"
rm -rf "$W"; mkdir -p "$W"
cd "$W"
for m in gcce cce; do
    mkdir -p $m
    cp -r "$REG/Diamond_NV_natab_$m/inputfiles" $m/
done
python3 - <<'PY'
import json, re; from input_layout import flat, nest
for m in ("gcce", "cce"):
    j = flat(json.loads(re.sub(r"(?m)^\s*#.*$", "", open(f"{m}/inputfiles/ccein_nstate1.json").read())))
    for sm, mode in (("all", "cluster_coherence"), ("allfull", "cluster_density_matrix")):
        json.dump(nest(dict(j, savemode=mode, outfile=f"./out_{sm}/r")), open(f"{m}/{sm}.json", "w"), indent=1)
PY
for m in gcce cce; do
    for np in 1 3; do
        for sm in all allfull; do
            mkdir -p $m/np$np && cp -r $m/inputfiles $m/np$np/ && cp $m/$sm.json $m/np$np/
            mkdir -p $m/np$np/out_$sm
            (cd $m/np$np && "$MPIRUN" -n $np "$BIN_DEV" -f $sm.json > $sm.log 2>&1 < /dev/null) \
                || { echo "  FAIL  $m np$np $sm: exit status (see $W/$m/np$np/$sm.log)"; FAIL=$((FAIL+1)); }
        done
        d=$m/np$np
        ok=1
        [ -f $d/out_all/r_state1_wiDiv ] || { echo "    savemode all wrote no result"; ok=0; }
        for f in $d/out_all/*; do
            cmp -s "$f" "$d/out_allfull/$(basename "$f")" || { echo "    allfull differs from all or lacks $(basename "$f")"; ok=0; }
        done
        nall=$(ls $d/out_all | wc -l)
        nfull=$(ls $d/out_allfull | grep -c '_full_' || true)
        nrest=$(ls $d/out_allfull | grep -vc '_full_' || true)
        [ "$nrest" -eq "$nall" ] || { echo "    allfull wrote $nrest files that are not _full_, all $nall"; ok=0; }
        if [ $m = gcce ]; then
            # one per cluster of all (its files less the two of the result, wiDiv and noDiv) and the zeroth cluster
            [ "$nfull" -eq "$nall" ] || { echo "    $nfull full files, expected $nall"; ok=0; }
            for f in $d/out_all/r_*_state1_wiDiv; do
                g=${f#$d/out_all/r}; g=$d/out_allfull/r${g%_state1_wiDiv}_full_state1_wiDiv
                [ -f "$g" ] || { echo "    no $g"; ok=0; break; }
            done
            ndist=$(md5sum $d/out_allfull/*_full_state1_wiDiv | awk '{print $1}' | sort -u | wc -l)
            [ "$ndist" -eq $((nall/2)) ] || { echo "    $ndist distinct full files of $((nall/2))"; ok=0; }
        else
            [ "$nfull" -eq 0 ] || { echo "    cCCE wrote $nfull full files"; ok=0; }
        fi
        if [ $ok -eq 1 ]; then echo "  PASS  F1 $m np$np: allfull = all ($nall files, the result among them) + $nfull full files, one per cluster"; PASS=$((PASS+1))
        else echo "  FAIL  F1 $m np$np"; FAIL=$((FAIL+1)); fi
    done
    if diff -rq $m/np1/out_allfull $m/np3/out_allfull > /dev/null; then
        echo "  PASS  F2 $m: 1 and 3 ranks write the same allfull files"; PASS=$((PASS+1))
    else echo "  FAIL  F2 $m: 1 and 3 ranks differ"; FAIL=$((FAIL+1)); fi
done
echo "allfull: $PASS passed, $FAIL failed"
[ $FAIL -eq 0 ]
