#!/usr/bin/env bash
#SBATCH --job-name=ccex-cce2-build
#SBATCH --partition=gpu
#SBATCH --qos=gpu-short
#SBATCH --export=NONE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=server-results/build-%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
source test/gpu_verification/server/environment.sh
python3 test/gpu_verification/server/validate_inputs.py

# bin/main.out (CPU) and bin/main_gpu.out (the same sources with CUDA), then the CUDA unit
# tests that measure.sbatch runs on the GPU
CCEX_JOBS=2 bash do_compile.sh gpu
make bin/test_cce_cuda.out bin/test_gcce_cuda.out

python3 test/gpu_verification/server/record_environment.py \
    --output "server-results/build-${SLURM_JOB_ID}-environment.json"

# Both executables use CPU here. GPU computation requires the subsequent allocation.
python3 test/gpu_verification/check_gpu/benchmark.py --binary bin/main_gpu.out --reference bin/main.out \
    --backend cpu --method all --output "server-results/build-${SLURM_JOB_ID}-cpu-checks"

