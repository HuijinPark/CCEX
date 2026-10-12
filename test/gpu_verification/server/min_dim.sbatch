#!/usr/bin/env bash
#SBATCH --job-name=ccex-gpu-min-dim
#SBATCH --partition=gpu
#SBATCH --qos=gpu-long
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=2-00:00:00
#SBATCH --output=server-results/min-dim-%j.log

# A3: per-dimension CPU vs CUDA time for gpu_options.min_dim; both backends from bin/main_gpu.out.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
source test/gpu_verification/server/environment.sh
run_dir="server-results/min-dim-${SLURM_JOB_ID}"
mkdir "$run_dir"
bash do_compile.sh gpu > "$run_dir/build.log" 2>&1
python3 test/gpu_verification/server/record_environment.py --gpu --output "$run_dir/environment.json"
python3 -u test/gpu_verification/gpu_benchmarks/benchmark_min_dim.py --cpu bin/main_gpu.out --cuda bin/main_gpu.out \
    --mpirun "$(command -v mpirun)" --output "$run_dir/conditions" "$@"
