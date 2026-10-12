# 09_GPU_backend

The CUDA backend. `backend` in the input chooses it: `cpu` (the default), `auto`
(CUDA for every cluster whose Hilbert dimension is at least `gpu_options.min_dim`, the CPU
for the rest, and the CPU everywhere in a binary built without CUDA) or `cuda`
(required; refused when it is not available). In the GPU build, `auto` and `cuda` with method CCE (coherence) also
send cCCE clusters of bath dimension 2 and 4 to a packed kernel by default (`gpu_options.packed_max_dim` 4). `docs/methods/gpu.md` has the details and
`test/gpu_verification/check_gpu/` the checks and benchmarks. With many ranks on one GPU, start CUDA MPS in the job
(`test/gpu_verification/server/mps.sh`): the ranks then share the GPU at once (GPU runs 0.49-0.62 of the time on idm4-gpu).

| Calculation | Inputs | What it shows | Recorded run |
|---|---|---|---|
| `CCE5_P1/` | ccein_cpu.json, ccein_auto.json, ccein_cuda.json | CCE5 Hahn echo of the NV in a five-spin P1 bath (electron + ¹⁵N each, so one cluster of dimension 1024): the same run per backend, the largest difference between the curves, the `[profile]` timings, and a plot | added 2026-09-23; WSL, 2 ranks, CPU only |
| `gCCE_P1/` | ccein_cpu.json, ccein_auto.json, ccein_cuda.json | gCCE4 Hahn echo of the NV (spin 1 with its zero-field splitting, the `Qubit` block) in four of those P1 (`bath/bath_p1_4`): the cluster of all four has dimension 768 with the qubit, the case the GPU is for; per backend, as `CCE5_P1`. The gCCE example of `test/gpu_verification/check_gpu/examples/p1_gcce4.json` with 100 time points | added 2026-10-12; idm4-gpu job 695, 32 ranks + RTX PRO 6000: cpu 35.5 s, auto 2.1 s, cuda 6.2 s, both GPU curves within 1.0e-10 of cpu; WSL, 4 ranks + RTX 3050: cpu 40.9 s, auto 16.8 s, cuda 16.8 s, 1.0e-10 |
| `pCCE_P1/` | ccein_cpu.json, ccein_auto.json, ccein_cuda.json | pCCE2 (`sK` 2) on `CCE5_P1`'s five P1: two groups, whose pair is one cluster of dimension 256 (on the GPU under `auto` too, `min_dim` 256); per backend, as `CCE5_P1` | added 2026-10-12; idm4-gpu job 695, 32 ranks: cpu 0.43 s, auto 0.44 s, cuda 4.2 s, difference 0; WSL, 4 ranks + RTX 3050: 0.38, 0.40, 0.39 s, 0 |

`bath/` is the five-spin bath of `test/gpu_verification/check_gpu/examples/` with the NV position and the
P1 gyromagnetic ratio, and the four-spin bath `bath_p1_4` of the same folder (the last four of the five). The single 1024-dimensional cluster is the case the GPU is
for: the CPU spends about a minute on its 100 time steps.

Each `run.sh` runs `cpu` with `bin/main.out`, and `auto` and `cuda` with `bin/main_gpu.out`
when that exists (`bash do_compile.sh gpu` on idm4-gpu makes both), so that `auto` actually
uses the GPU. Without it, `auto` runs with `bin/main.out`, falls back to the CPU and its curve
is identical to the `cpu` one, and `cuda` is skipped. `CCEX_GPU_BIN` names another GPU binary (a laptop card, e.g.
`CUDA_HOME=/usr/local/cuda-12.8 CUDA_ARCH=86 bash do_compile.sh generic` for an RTX 3050).

What the 2026-10-12 runs say: the one 768-dimensional gCCE cluster is where the GPU pays
(its evolution took 1.6 s on the RTX PRO 6000; the whole CPU run 35.5 s, where that cluster sits on one rank). `cuda` sends
every cluster to the GPU, the qubit alone (dimension 3) included, and was slower than `auto`
there (6.2 s against 2.1 s) and for the 0.4 s pCCE run (4.2 s); the extra time is presumably
the CUDA start-up of those small clusters, not measured separately. On idm4-gpu the
analyzer's plot at the end of `run.sh` stopped with `matplotlib.cm has no attribute
get_cmap` (the python/3.13 module's newer matplotlib); the CCEX runs and the printed
comparisons before it were complete. The analyzer now uses `pyplot.get_cmap`, which old and new matplotlib both have.

    cd CCE5_P1 && NP=4 bash run.sh

`example/tutorial/README.md` explains the layout.
