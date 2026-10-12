#!/bin/sh
# 09_GPU_backend/gCCE_P1: gCCE4 Hahn echo of the NV (spin 1 with its zero-field splitting) in four P1 spins (one cluster of dimension 768), once per backend: cpu, auto, and cuda when bin/main_gpu.out exists
#
# Added 2026-10-12. idm4-gpu job 695 (32 ranks, RTX PRO 6000): cpu 35.5 s, auto 2.1 s, cuda 6.2 s, GPU curves within 1e-10 of cpu.
# WSL, 4 ranks: CPU-only binary cpu 40 s, auto 41 s (falls back to the CPU); with an RTX 3050 build auto and cuda 16.8 s.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
CCEX_GPU=${CCEX_GPU_BIN:-$ROOT/bin/main_gpu.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
mkdir -p results

mpirun -n $NP $CCEX -f ccein_cpu.json -o results/gCCE4_P1_cpu_1 > results/process_gCCE4_P1_cpu_1
BACKENDS="cpu"
if [ -x $CCEX_GPU ]; then
  mpirun -n $NP $CCEX_GPU -f ccein_auto.json -o results/gCCE4_P1_auto_1 > results/process_gCCE4_P1_auto_1
  mpirun -n $NP $CCEX_GPU -f ccein_cuda.json -o results/gCCE4_P1_cuda_1 > results/process_gCCE4_P1_cuda_1
  BACKENDS="cpu auto cuda"
else
  mpirun -n $NP $CCEX -f ccein_auto.json -o results/gCCE4_P1_auto_1 > results/process_gCCE4_P1_auto_1
  BACKENDS="cpu auto"
fi

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# the backends must agree: largest difference from the cpu curve
python3 -c 'import sys, numpy as np; sys.path.insert(0, sys.argv[1]); from coherence import read_coherence
ref = read_coherence("results/gCCE4_P1_cpu_1_wiDiv")[1]
[print("max |L_%s - L_cpu| = %.3e" % (b, abs(read_coherence("results/gCCE4_P1_%s_1_wiDiv" % b)[1] - ref).max())) for b in sys.argv[3:]]' $(dirname $ANALYZER) $BACKENDS
# wall time and the timings per backend and dimension: the [profile] lines of each process file
grep -H 'JOB DONE' results/process_gCCE4_P1_*
grep -h '\[profile\]' results/process_gCCE4_P1_* | head -30 || true
# plot
python3 $ANALYZER -d results/ -fi gCCE4_P1_PREFIX_1_wiDiv -vn "PREFIX" -v $BACKENDS -uc "us" -ntn -fo results/gCCE4_P1_Comp_backend -pl_m
