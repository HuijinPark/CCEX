#!/bin/sh
# 08_Bundle_clustering/CCE2_bundle: CCE2 Hahn echo of the NV of 02 (bath copied to ../bath) (rbath 30) with ordinary clusters and with bundles of up to three 13C, plotted together
#
# Added 2026-09-23 (no idm record). WSL, 2 ranks: 1 s standard, 21 s bundled.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=${CCEX_BIN:-$ROOT/bin/main.out}
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
mkdir -p results

mpirun -n $NP $CCEX -f ccein_standard.json -o results/CCE2_standard_1 > results/process_CCE2_standard_1
mpirun -n $NP $CCEX -f ccein_bundle.json -o results/CCE2_bundle_1 > results/process_CCE2_bundle_1

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# how the bath was bundled, and how far the two curves are apart
grep -h '^\[adaptive\]' results/process_CCE2_bundle_1
python3 -c 'import sys, numpy as np; sys.path.insert(0, sys.argv[1]); from coherence import read_coherence
a = read_coherence("results/CCE2_standard_1_wiDiv")[1]; b = read_coherence("results/CCE2_bundle_1_wiDiv")[1]
print("max |L_bundle - L_standard| = %.3e" % abs(b - a).max())' $(dirname $ANALYZER)
# plot
python3 $ANALYZER -d results/ -fi CCE2_PREFIX_1_wiDiv -vn "PREFIX" -v standard bundle -uc "us" -ntn -fo results/CCE2_Comp_bundle -pl_m
