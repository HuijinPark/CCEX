# check_gpu/ — the CUDA backend gives the CPU's numbers

| script | checks | run |
|---|---|---|
| `test_*.cpp` (14) | unit tests of the kernels, the batching, the scaled product, the Hamiltonian check before the CUDA cCCE path, the seed sharing and the host-side lifetimes | `make numeric-test`, `cpu-test`, `gpu-test` (not `bash do_compile.sh gpu`, which builds `bin/main.out` and `bin/main_gpu.out` only): see "Unit tests" below |
| `check_cce_pulses.py --cpu BIN --cuda BIN --output DIR` | cCCE with every pulse count, CPU against CUDA | CUDA build |
| `check_acceleration.py`, `check_pcce.py`, `check_options.py` | gCCE4 accuracy and batching, pCCE, `backend` and the `gpu_options` block, CPU against CUDA | CUDA build |
| `check_seed_states.py --cpu BIN --output DIR` | seeds and sampled bath states are reproduced across ranks | CPU build is enough |
| `check_review_controls.py --cpu BIN --output DIR` | the controls the 2026-09-12 review asked for | CPU build is enough |
| `check_review_supplement.py`, `check_review_baseline.py`, `check_host_leaks.py` | the review's supplement, its baseline, host memory over repeated runs | CUDA build, and a `--baseline` from an earlier run |

## Unit tests: which command builds and runs which

After `bash do_compile.sh <preset>` has written the Makefile (or a CMake configure):

| tests | built and run by | needs |
|---|---|---|
| `test_scaled_product`, `test_hamiltonian_check` | `make numeric-test` (CMake: `cmake --build build --target numeric-test`) | a CPU |
| `test_cluster_log`, `test_review_controls`, `test_mpi_seed` (1 and 4 ranks), `test_host_lifetimes` (1 and 4 ranks) | `make cpu-test` (CMake: `--target cpu-test`; also in the CI workflow) | a CPU and `mpirun` (`make cpu-test MPIRUN=...` picks another; CMake uses the launcher FindMPI found) |
| the numeric two and `test_gcce_cuda`, `test_cce_cuda`, `test_gcce_accuracy` (twice: built-in case and `data/gcce_scale_separation_three_spins.txt`), `test_gcce_refinement` (the parallel refinement at dim 192, 15, 2 and 1), `test_gcce_batch`, `test_review_failure`, `test_review_memory`, `test_small_cluster_cuda` (the packed kernel of `gpu_options.packed_max_dim`, A26: the same arithmetic on the host and on the card against the CPU coherences of the 11 files in `data/small_cluster/`, 1e-8 + 2e-6 |CPU| per point, both lanes the same bits; `data/small_cluster/make_dumps.py` writes the files again with a CPU build, byte for byte; `small_cluster_host_emulation.cpp` stands in for the card so the node feeder and shared-memory windows can be run end to end without CUDA, as its header says) | `make gpu-test` (CMake: `--target gpu-test`) | a CUDA build and the GPU, so inside a job; `server/build.sbatch` builds `test_cce_cuda` and `test_gcce_cuda`, and `server/measure.sbatch` runs those two only |

`cpu-test` compiles the four with the build's own flags and `-UNDEBUG` for their asserts;
`test_host_lifetimes` links the objects of `bin/main.out` but `main.o` (CMake compiles those
sources again for it). Until 2026-10-02 these four and the two CUDA review tests had no target
and were built by hand (D17); built without `-march=native` against objects built with it,
`test_host_lifetimes` crashed (exit 139) on 2026-10-01.

On 2026-09-24 `check_review_controls.py --cpu` passed on the CPU build;
`check_seed_states.py` wrote `adaptive_cluster_max_dim` in its adaptive cases until
2026-09-24; it writes `bundle_max_dim` now and its nine cases pass on the CPU build.
One run in two failed in `benchmark.parse_profiles` because MPI interleaved two
ranks' `[profile]` lines in the p1 case; that is the output forwarding, not the
code, and a rerun passed.

`benchmark.py` is the helper the checks and `../gpu_benchmarks/` import: it runs a
build on the four standard calculations of `code_verification/check_input_forms/inputs`
and reads what it wrote. `examples/` holds the P1 and pCCE inputs the checks use, `data/` the gCCE accuracy
test's files and, in `data/small_cluster/`, the 11 `.bin` dumps of `test_small_cluster_cuda` with `make_dumps.py`. The seven `validation-*.json` files are the
results recorded on idm4-gpu and the RTX 3050 on 2026-09-08..11 that (internal notes)
cites; the other records of 2026-09-12 and the probe that led to `scaled_product.h`
were removed on 2026-09-24 (git history up to 8f2bcbd75).
