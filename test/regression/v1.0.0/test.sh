#!/usr/bin/env bash
#
# CCEX v1.0.0 regression suite
# ---------------------------------------------------------------------------
# Runs every test in this directory. A "test" is any subdirectory holding an
# executable `job`, so adding a case needs no edit here.
#
# The jobs stay completely independent: this script only cd's into each one and
# calls ./job, exactly as you would by hand. Nothing is passed in, nothing is
# patched, no environment is exported -- so
#
#     cd Diamond_NV_natab_cce && ./job
#
# still works on its own and produces the same output as running it from here.
#
# Usage:
#   ./test.sh                       run every test, in the order listed by -l
#   ./test.sh <name> [<name> ...]   run only those (directory names)
#   ./test.sh -l | --list           list the tests and exit
#   ./test.sh -h | --help           this text
#
# Each test's console output is streamed AND saved to <test>/job.log.
# Exit status is 0 only if every test that ran returned 0.
#
# NOTE ON COST: the tests are not the same size. Diamond_NV_P1_10ppm_gcce is a
# 10 ppm P1 ELECTRON bath -- 201 sites, 1823 pairs, and gCCE propagates 48x48
# blocks for a pair (3 x 4 x 4) against 12x12 for the 13C gCCE test. Run it on
# its own first if you have not timed it before.
# ---------------------------------------------------------------------------

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

# --- discover tests: any subdir with an executable job -----------------------
mapfile -t ALL < <(for d in */; do [ -x "${d}job" ] && basename "$d"; done | sort)

usage() { sed -n '3,29p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

case "${1-}" in
    -h|--help) usage; exit 0 ;;
    -l|--list)
        echo "tests in $(basename "$HERE"):"
        for t in "${ALL[@]}"; do
            m=$(sed -n 's/^[[:space:]]*"method"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' \
                    "$t"/inputfiles/ccein_*.json 2>/dev/null | head -1)
            n=$(ls -1 "$t"/inputfiles/ccein_*.json 2>/dev/null | wc -l)
            printf "    %-28s method=%-5s ccein=%s\n" "$t" "${m:-?}" "$n"
        done
        exit 0 ;;
esac

if [ $# -gt 0 ]; then
    TESTS=("$@")
    for t in "${TESTS[@]}"; do
        [ -x "$t/job" ] || { echo "ERROR: no executable job in ./$t" >&2
                             echo "       available: ${ALL[*]}" >&2; exit 1; }
    done
else
    TESTS=("${ALL[@]}")
fi
[ ${#TESTS[@]} -gt 0 ] || { echo "ERROR: no test found in $HERE" >&2; exit 1; }

# --- run ---------------------------------------------------------------------
echo "############################################################"
echo "# CCEX v1.0.0 regression suite -- ${#TESTS[@]} test(s)"
echo "############################################################"
printf '  %s\n' "${TESTS[@]}"
echo

names=(); codes=(); walls=()
rc_all=0

for t in "${TESTS[@]}"; do
    echo "============================================================"
    echo "== $t"
    echo "============================================================"

    t0=$(date +%s.%N)
    ( cd "$t" && ./job ) 2>&1 | tee "$t/job.log"
    rc=${PIPESTATUS[0]}
    t1=$(date +%s.%N)
    wall=$(awk "BEGIN{printf \"%.2f\", $t1-$t0}")

    names+=("$t"); codes+=("$rc"); walls+=("$wall")
    [ "$rc" -ne 0 ] && rc_all=1

    echo
    printf ">> %s : exit=%s, wall=%ss  (log: %s/job.log)\n\n" "$t" "$rc" "$wall" "$t"
done

# --- summary -----------------------------------------------------------------
echo "############################################################"
echo "# summary"
echo "############################################################"
printf "  %-30s %6s %10s   %s\n" "test" "exit" "wall[s]" "result"
for i in "${!names[@]}"; do
    if [ "${codes[$i]}" -eq 0 ]; then res="PASS"; else res="FAIL"; fi
    printf "  %-30s %6s %10s   %s\n" "${names[$i]}" "${codes[$i]}" "${walls[$i]}" "$res"
done
echo
if [ $rc_all -eq 0 ]; then
    echo "  all ${#names[@]} test(s) returned 0"
else
    echo "  at least one test failed"
fi

exit $rc_all
