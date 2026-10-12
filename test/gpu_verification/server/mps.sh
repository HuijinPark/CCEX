#!/usr/bin/env bash
# CUDA MPS for one Slurm job on idm4-gpu: source after environment.sh. With MPS the ranks of a job share the GPU at
# once instead of taking turns (job 698: GPU runs of 32 ranks took 0.49-0.62 of the time, curves byte-identical).
# The daemon runs as the user, its pipe and log in a private mktemp directory on the node's /tmp, and stops (and the
# directory goes) when the job ends or is cancelled.
# CCEX_MPS=0 leaves it off.
if [ "${CCEX_MPS:-1}" != 0 ] && command -v nvidia-cuda-mps-control > /dev/null 2>&1; then
    # a fresh directory only this user can open (mktemp: unpredictable name, mode 0700), never a fixed path in /tmp
    ccex_mps_dir="$(mktemp -d "/tmp/ccex-mps-${SLURM_JOB_ID:-$$}.XXXXXX")" || { echo "[mps] no temporary directory; running without it"; return 0 2>/dev/null || exit 0; }
    export CUDA_MPS_PIPE_DIRECTORY="$ccex_mps_dir/pipe"
    export CUDA_MPS_LOG_DIRECTORY="$ccex_mps_dir/log"
    mkdir -m 700 "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
    if nvidia-cuda-mps-control -d; then
        trap 'echo quit | nvidia-cuda-mps-control > /dev/null 2>&1; sleep 1; case "$ccex_mps_dir" in /tmp/ccex-mps-*) rm -rf -- "$ccex_mps_dir";; esac' EXIT
        echo "[mps] on (pipe $CUDA_MPS_PIPE_DIRECTORY)"
    else
        unset CUDA_MPS_PIPE_DIRECTORY CUDA_MPS_LOG_DIRECTORY
        case "$ccex_mps_dir" in /tmp/ccex-mps-*) rm -rf -- "$ccex_mps_dir";; esac
        echo "[mps] could not start; running without it"
    fi
else
    echo "[mps] off"
fi
