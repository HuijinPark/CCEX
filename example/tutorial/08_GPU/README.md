# 08_GPU

The CUDA backend. `execution.backend` chooses it: `cpu` (the default), `auto` (CUDA for every cluster whose Hilbert
dimension is at least `execution.gpu_options.min_dim`, the CPU for the rest, and the CPU everywhere in a binary built
without CUDA) or `cuda` (every cluster; refused when CUDA is not available). In the GPU build, `auto` and `cuda` with
cCCE also send clusters of bath dimension 2 and 4 to a packed kernel (`gpu_options.packed_max_dim` 4).
`docs/methods/gpu.md` has the details and `test/gpu_verification/check_gpu/` the checks. Each folder has one input with
`backend` cpu; `run.sh` runs it with `bin/main.out`, then with `backend` auto and cuda (`set_input_keys.py`) with
`bin/main_gpu.out` when that exists (`bash do_compile.sh gpu` on idm4-gpu makes both; `CCEX_GPU_BIN` names another GPU
binary, e.g. a laptop card: `CUDA_HOME=/usr/local/cuda-12.8 CUDA_ARCH=86 bash do_compile.sh generic`). Without it
`auto` runs on the CPU and `cuda` is skipped. With many ranks on one GPU, start CUDA MPS in the job
(`test/gpu_verification/server/mps.sh`). Every key not written in an input has its default (`main.out -e`).

`bath/` is the P1 bath of `test/gpu_verification/check_gpu/examples/`: five P1 (`bath_p1_5`, electron + ¹⁵N each), the
last four of them (`bath_p1_4`), the NV position (`qubit_p1`) and the P1 gyromagnetic ratio (`gyro_P1`).

## cCCE5_Hahn_backends

- **What it shows:** cCCE5 Hahn echo of the NV in the five P1: one cluster of dimension 1024, the case the GPU is for
  (about a minute of CPU for its 100 time steps); the same run per backend, the `[profile]` timings and the curves.
- **Keys that matter:** `execution.backend`, `execution.gpu_options.min_dim` 256, `execution.profile` true, `order` 5,
  the `Defect` block (¹⁵N, four Jahn-Teller axes).
- **How to run:** on a GPU node, `NP=32 bash run.sh` (CPU part alone: 24 s on 4 ranks, WSL).
- **Expected result:** pending (server): (internal notes) makes
  `expected.png`; the backends agree to about 1e-10.
- **What to change next:** `gpu_options.min_dim` (which clusters go to the GPU).

## gCCE4_Hahn_backends

- **What it shows:** gCCE4 Hahn echo of the NV (spin 1 with its zero-field splitting) in four P1: the cluster of all
  four has dimension 768 with the qubit. Recorded 2026-10-12 (idm4-gpu job 695, 32 ranks + RTX PRO 6000): cpu 35.5 s,
  auto 2.1 s, cuda 6.2 s, both GPU curves within 1.0e-10 of cpu; WSL, 4 ranks + RTX 3050: 40.9 / 16.8 / 16.8 s.
- **Keys that matter:** `method` gCCE, `execution.backend`, `gpu_options.min_dim` 256, `evolution_batch_steps` 64.
- **How to run:** on a GPU node, `NP=32 bash run.sh` (CPU part alone: 85 s on 4 ranks, WSL).
- **Expected result:** pending (server), as above.
- **What to change next:** `cuda` sends every cluster, the qubit alone included, to the GPU and was slower than `auto`
  here (6.2 s against 2.1 s).

## pCCE2_Hahn_backends

- **What it shows:** pCCE2 (partitions of two P1) on the five P1: the pair of partitions is one cluster of dimension
  256, on the GPU under `auto` too (`min_dim` 256). Recorded 2026-10-12 (job 695, 32 ranks): cpu 0.43 s, auto 0.44 s,
  cuda 4.2 s, difference 0.
- **Keys that matter:** `cluster.mode` pcce, `cluster.pcce.partition_size` 2, `execution.backend`.
- **How to run:** on a GPU node, `NP=32 bash run.sh` (CPU part alone: 1 s on 4 ranks, WSL).
- **Expected result:** pending (server), as above.
- **What to change next:** `partition_size` 3 (larger clusters, more of them on the GPU).
