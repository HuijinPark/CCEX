#!/bin/sh
# 09_GPU_backend/CCE5_P1: CCE5 Hahn echo of the NV in a five-spin P1 bath, once per backend: cpu, auto, and cuda when bin/main_gpu.out exists
#
# Added 2026-09-23 (no idm record). On a CPU-only binary "auto" falls back to the CPU and the
# curves agree to rounding; on idm4-gpu, "bash do_compile.sh gpu" also makes bin/main_gpu.out.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
CCEX_GPU=${CCEX_GPU_BIN:-$ROOT/bin/main_gpu.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
mkdir -p results

mpirun -n $NP $CCEX -f ccein_cpu.js -o results/CCE5_P1_cpu_1 > results/process_CCE5_P1_cpu_1
mpirun -n $NP $CCEX -f ccein_auto.js -o results/CCE5_P1_auto_1 > results/process_CCE5_P1_auto_1
BACKENDS="cpu auto"
if [ -x $CCEX_GPU ]; then
  mpirun -n $NP $CCEX_GPU -f ccein_cuda.js -o results/CCE5_P1_cuda_1 > results/process_CCE5_P1_cuda_1
  BACKENDS="cpu auto cuda"
fi

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# the backends must agree: largest difference from the cpu curve
python3 -c 'import sys, numpy as np; sys.path.insert(0, sys.argv[1]); from coherence import read_coherence
ref = read_coherence("results/CCE5_P1_cpu_1_wiDiv")[1]
[print("max |L_%s - L_cpu| = %.3e" % (b, abs(read_coherence("results/CCE5_P1_%s_1_wiDiv" % b)[1] - ref).max())) for b in sys.argv[3:]]' $(dirname $ANALYZER) $BACKENDS
# timings per backend and dimension: the [profile] lines of each process file
grep -h '\[profile\]' results/process_CCE5_P1_* | head -20 || true
# plot
python3 $ANALYZER -d results/ -fi CCE5_P1_PREFIX_1_wiDiv -vn "PREFIX" -v $BACKENDS -uc "us" -ntn -fo results/CCE5_P1_Comp_backend -pl_m
