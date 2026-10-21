#!/usr/bin/env bash
# mps.sh without a GPU (review 20261012-233058-452278 R1/P3-07): a stand-in nvidia-cuda-mps-control on PATH, mps.sh
# sourced by a shell under set -euo pipefail as the jobs do. Every case must remove the private /tmp/ccex-mps-* directory
# and end with the job's own status: a quit that fails changes neither.
#   bash test/gpu_verification/server/test_mps.sh
set -u
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
work="$(mktemp -d)"
trap 'rm -rf -- "$work"' EXIT
mkdir "$work/bin"
cat > "$work/bin/nvidia-cuda-mps-control" << 'EOF'
#!/bin/bash
if [ "${1-}" = -d ]; then printf '%s\n' "$CUDA_MPS_PIPE_DIRECTORY" > "$MOCK_RECORD"; exit "${MOCK_START_RC:-0}"; fi
read -r command
exit "${MOCK_QUIT_RC:-0}"
EOF
chmod +x "$work/bin/nvidia-cuda-mps-control"

failed=0
# case NAME START_RC QUIT_RC BODY EXPECTED_STATUS ("nonzero" for any)
case_() {
    local name=$1 start=$2 quit=$3 body=$4 want=$5 got dir
    rm -f "$work/record"
    PATH="$work/bin:$PATH" MOCK_RECORD="$work/record" MOCK_START_RC=$start MOCK_QUIT_RC=$quit SLURM_JOB_ID="test$$" \
        bash -euo pipefail -c 'source "$1"; '"$body" bash "$here/mps.sh" > "$work/out" 2>&1
    got=$?
    dir="$(cat "$work/record" 2>/dev/null)"; dir="${dir%/pipe}"
    local ok=1
    case "$dir" in /tmp/ccex-mps-test$$.*) ;; *) ok=0;; esac
    [ -n "$dir" ] && [ -e "$dir" ] && ok=0
    if [ "$want" = nonzero ]; then [ "$got" -ne 0 ] || ok=0; else [ "$got" -eq "$want" ] || ok=0; fi
    if [ "$ok" = 1 ]; then
        echo "PASS mps $name: status $got, $dir removed"
    else
        echo "FAIL mps $name: status $got (want $want), directory '$dir' $( [ -e "$dir" ] && echo left || echo removed )"
        sed 's/^/    /' "$work/out"
        [ -n "$dir" ] && case "$dir" in /tmp/ccex-mps-test$$.*) rm -rf -- "$dir";; esac
        failed=1
    fi
}
# the body also checks that the directory is private (0700) while the job runs
private='[ "$(stat -c %a "${CUDA_MPS_PIPE_DIRECTORY%/pipe}")" = 700 ]'
case_ "quit ok"                 0 0 "$private; true"    0
case_ "quit fails"              0 1 "$private; true"    0
case_ "quit fails, job fails"   0 1 "$private; exit 3"  3
case_ "job fails"               0 0 "exit 4"            4
case_ "start fails"             1 0 'true'              0
case_ "SIGTERM"                 0 1 'kill -TERM $$; sleep 5' nonzero
[ "$failed" = 0 ] && echo "PASS test_mps" || echo "FAIL test_mps"
exit "$failed"
