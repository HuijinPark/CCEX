#!/usr/bin/env bash
#
# v1.1.0 <-> axis-align equivalence suite.
#
#     bash test/equivalence/run.sh          run everything, in order
#     bash test/equivalence/run.sh axis1    just one stage
#
# Needs bin/main.out built (bash do_compile.sh wsl) and the preserved v1.1.0
# binary named in lib.sh. Everything written lands in _work/, which is gitignored.
#
# The stages are ordered by what depends on what, and the suite stops at the first
# failure: axis (2) is only meaningful once axis (1) has tied the legacy inputs
# back to v1.1.0, and nothing is meaningful if the binary does not build.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAGES=(axis1 rotation axis2 reject errcorr)
WANT=("${@:-${STAGES[@]}}")
rc_all=0

for stage in "${WANT[@]}"; do
    echo
    echo "###########################################################"
    echo "# $stage"
    echo "###########################################################"
    case "$stage" in
        axis1)    bash "$HERE/axis1.sh" ;;
        # the upstream suite: the geometry unit test, the run-level cases, and the
        # A/B check against a hand-transformed input. It is dev-internal by design,
        # which is why axis1 exists beside it.
        rotation) bash "$HERE/../rotation/run.sh" ;;
        axis2)    bash "$HERE/axis2.sh" ;;
        reject)   bash "$HERE/reject.sh" ;;
        errcorr)  bash "$HERE/errcorr.sh" ;;
        *) echo "unknown stage: $stage"; rc_all=1; continue ;;
    esac
    rc=$?
    [ $rc -ne 0 ] && { echo "STAGE FAILED: $stage (exit=$rc)"; rc_all=$rc; break; }
done

echo
echo "==========================================================="
[ $rc_all -eq 0 ] && echo " all stages passed" || echo " FAILED (exit=$rc_all)"
echo "==========================================================="
exit $rc_all
