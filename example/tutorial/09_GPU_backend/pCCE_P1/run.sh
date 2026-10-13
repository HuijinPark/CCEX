#!/bin/sh
# 09_GPU_backend/pCCE_P1: pCCE2 (groups of pcce.partition_size = 2 P1) Hahn echo of the NV in the five-spin P1 bath of CCE5_P1, once per backend: cpu, auto, and cuda when bin/main_gpu.out exists
#
# Added 2026-10-12. idm4-gpu job 695 (32 ranks, RTX PRO 6000): cpu 0.43 s, auto 0.44 s, cuda 4.2 s, curves identical.
# WSL, 4 ranks: cpu 0.4 s, auto 0.4 s (CPU-only binary or RTX 3050 build), cuda 0.4 s (RTX 3050).
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
CCEX_GPU=${CCEX_GPU_BIN:-$ROOT/bin/main_gpu.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
mkdir -p results

mpirun -n $NP $CCEX -f ccein_cpu.json -o results/pCCE2_P1_cpu_1 > results/process_pCCE2_P1_cpu_1
BACKENDS="cpu"
if [ -x $CCEX_GPU ]; then
  mpirun -n $NP $CCEX_GPU -f ccein_auto.json -o results/pCCE2_P1_auto_1 > results/process_pCCE2_P1_auto_1
  mpirun -n $NP $CCEX_GPU -f ccein_cuda.json -o results/pCCE2_P1_cuda_1 > results/process_pCCE2_P1_cuda_1
  BACKENDS="cpu auto cuda"
else
  mpirun -n $NP $CCEX -f ccein_auto.json -o results/pCCE2_P1_auto_1 > results/process_pCCE2_P1_auto_1
  BACKENDS="cpu auto"
fi

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# the backends must agree: largest difference from the cpu curve (the partition is the same: deterministic, from the seed)
python3 -c 'import sys, numpy as np; sys.path.insert(0, sys.argv[1]); from coherence import read_coherence
ref = read_coherence("results/pCCE2_P1_cpu_1_wiDiv")[1]
[print("max |L_%s - L_cpu| = %.3e" % (b, abs(read_coherence("results/pCCE2_P1_%s_1_wiDiv" % b)[1] - ref).max())) for b in sys.argv[3:]]' $(dirname $ANALYZER) $BACKENDS
# wall time and the timings per backend and dimension: the [profile] lines of each process file
grep -H 'JOB DONE' results/process_pCCE2_P1_*
grep -h '\[profile\]' results/process_pCCE2_P1_* | head -30 || true
# plot
python3 $ANALYZER -d results/ -fi pCCE2_P1_PREFIX_1_wiDiv -vn "PREFIX" -v $BACKENDS -uc "us" -ntn -fo results/pCCE2_P1_Comp_backend -pl_m
