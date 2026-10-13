#!/usr/bin/env bash
#SBATCH --job-name=ccex-cce2-measure
#SBATCH --partition=gpu
#SBATCH --qos=gpu-short
#SBATCH --export=NONE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=server-results/measure-%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
source test/gpu_verification/server/environment.sh
source test/gpu_verification/server/mps.sh
run_dir="server-results/measure-${SLURM_JOB_ID}"
mkdir "$run_dir"
python3 test/gpu_verification/server/validate_inputs.py > "$run_dir/inputs.json"
python3 test/gpu_verification/server/record_environment.py --gpu --output "$run_dir/environment.json"

./bin/test_cce_cuda.out
./bin/test_gcce_cuda.out

# Repeat CPU references within this GPU allocation for both comparisons.
python3 test/gpu_verification/check_gpu/benchmark.py --binary bin/main_gpu.out --reference bin/main.out \
    --backend cuda --method all --repeats 3 --output "$run_dir/cuda"
python3 test/gpu_verification/check_gpu/benchmark.py --binary bin/main_gpu.out --reference bin/main.out \
    --backend auto --method all --repeats 3 --output "$run_dir/auto"

