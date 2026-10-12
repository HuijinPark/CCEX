#!/bin/sh
# 10_Coordinate_rotation/HahnEcho_CCE2: CCE2 Hahn echo of the NV (02's system, bath copied to ../bath) three ways: the bath as generated (NV axis along z), the same bath written in the cubic frame with coordinate_frame_rotation on, and that cubic-frame bath with the rotation off
#
# Added 2026-09-23 (no idm record). aligned and cubic_rotated must agree to rounding;
# cubic_norotation puts the field 54.7 degrees off the NV axis and differs.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
mkdir -p results

mpirun -n $NP $CCEX -f ccein_aligned.js -o results/CCE2_aligned_1 > results/process_CCE2_aligned_1
mpirun -n $NP $CCEX -f ccein_cubic_rotated.js -o results/CCE2_cubic_rotated_1 > results/process_CCE2_cubic_rotated_1
mpirun -n $NP $CCEX -f ccein_cubic_norotation.js -o results/CCE2_cubic_norotation_1 > results/process_CCE2_cubic_norotation_1

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# largest difference from the aligned curve: ~1e-12 with the rotation on, of order 0.1 with it off
python3 -c 'import sys, numpy as np; sys.path.insert(0, sys.argv[1]); from coherence import read_coherence
ref = read_coherence("results/CCE2_aligned_1_wiDiv")[1]
[print("max |L_%s - L_aligned| = %.3e" % (b, abs(read_coherence("results/CCE2_%s_1_wiDiv" % b)[1] - ref).max())) for b in sys.argv[2:]]' $(dirname $ANALYZER) cubic_rotated cubic_norotation
# plot
python3 $ANALYZER -d results/ -fi CCE2_PREFIX_1_wiDiv -vn "PREFIX" -v aligned cubic_rotated cubic_norotation -uc "us" -ntn -fo results/CCE2_Comp_frame -pl_m
