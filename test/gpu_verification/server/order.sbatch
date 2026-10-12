#!/usr/bin/env bash
#SBATCH --job-name=ccex-high-order
#SBATCH --partition=gpu
#SBATCH --qos=gpu-long
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=3-00:00:00
#SBATCH --output=server-results/order-%j.log

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
source test/gpu_verification/server/environment.sh
source test/gpu_verification/server/mps.sh
order="${1:?Supply order 3, 4 or 5}"
case "$order" in 3|4|5) ;; *) exit 2 ;; esac
run_dir="server-results/order${order}-${SLURM_JOB_ID}"
mkdir "$run_dir"
python3 test/gpu_verification/server/validate_inputs.py > "$run_dir/base-inputs.json"
python3 test/gpu_verification/server/record_environment.py --gpu --output "$run_dir/environment.json"
python3 test/gpu_verification/server/run_orders.py --order "$order" --ranks 32 --batch 4 \
    --output "$run_dir/comparisons"
