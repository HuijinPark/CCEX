#!/usr/bin/env bash
#SBATCH --job-name=ccex-gpu-sharing
#SBATCH --partition=gpu
#SBATCH --qos=gpu-short
#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=1-00:00:00
#SBATCH --output=server-results/sharing-%j.log

# A3: gpu_options.min_dim against the number of ranks sharing the GPU (1, 4, 8, 16, 32); benchmark_gpu_sharing.py.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
source test/gpu_verification/server/environment.sh
run_dir="server-results/sharing-${SLURM_JOB_ID}"
mkdir "$run_dir"
bash do_compile.sh gpu > "$run_dir/build.log" 2>&1
python3 test/gpu_verification/server/record_environment.py --gpu --output "$run_dir/environment.json"
python3 -u test/gpu_verification/gpu_benchmarks/benchmark_gpu_sharing.py --binary bin/main_gpu.out \
    --mpirun "$(command -v mpirun)" --output "$run_dir/conditions" "$@"
