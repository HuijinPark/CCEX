# 09_GPU_backend

The CUDA backend. `backend` in the input chooses it: `cpu` (the default), `auto`
(CUDA for every cluster whose Hilbert dimension is at least `gpu_min_dim`, the CPU
for the rest, and the CPU everywhere in a binary built without CUDA) or `cuda`
(required; refused when it is not available). `docs/methods/gpu.md` has the details and
`test/gpu_verification/check_gpu/` the checks and benchmarks.

| Calculation | Inputs | What it shows | Recorded run |
|---|---|---|---|
| `CCE5_P1/` | ccein_cpu.js, ccein_auto.js, ccein_cuda.js | CCE5 Hahn echo of the NV in a five-spin P1 bath (electron + ¹⁵N each, so one cluster of dimension 1024): the same run per backend, the largest difference between the curves, the `[profile]` timings, and a plot | added 2026-09-23; WSL, 2 ranks, CPU only |

`bath/` is the five-spin bath of `test/gpu_verification/check_gpu/examples/` with the NV position and the
P1 gyromagnetic ratio. The single 1024-dimensional cluster is the case the GPU is
for: the CPU spends about a minute on its 100 time steps.

`run.sh` runs `cpu` and `auto` with `bin/main.out`, and `cuda` with `bin/main_gpu.out`
when that exists (`bash do_compile.sh gpu` on idm4-gpu makes both). On a machine
without CUDA, `auto` falls back to the CPU and the two curves are identical.

    cd CCE5_P1 && NP=4 bash run.sh

`example/tutorial/README.md` explains the layout.
