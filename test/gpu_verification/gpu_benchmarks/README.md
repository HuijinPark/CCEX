# gpu_benchmarks/ — how much faster, not whether right

| script | times |
|---|---|
| `benchmark_acceleration.py` | gCCE4 clusters by dimension, CPU against CUDA |
| `benchmark_pcce_pl6.py`, `benchmark_pl6_acceleration.py` | the PL6 pCCE case, where the GPU helps most |
| `benchmark_rank_scaling.py` | wall time against the number of MPI ranks |

All want a CUDA build and an idle GPU (`../server/` has the jobs that ran them on
idm4-gpu). They import `../check_gpu/benchmark.py`, the helper every check shares
(it runs a build on the four standard calculations and reads what it wrote);
`benchmark_pcce_pl6.py` reads coherence files through
`example/analysis/Analyzer/CoherenceAnalyzer.v1/coherence.py`.
