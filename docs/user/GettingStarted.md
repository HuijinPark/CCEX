# Getting started
<a id="GettingStarted"></a>

How to build CCEX, run it, and what a run does with its input.

<a id="Building"></a>
## Building

```
bash do_compile.sh <preset>                     # wsl-gnu, wsl-icpx, idm4, gpu, nurion, mac, or generic
cmake -B build && cmake --build build -j 8      # any machine: finds the MPI wrapper and, if there is one, nvcc
```

**What is built.** `do_compile.sh` writes a `Makefile` at the top (generated: change `do_compile.sh`, never the
`Makefile`), `obj/` and `bin/main.out`, and, where there is CUDA, `bin/main_gpu.out` for `execution.backend` `cuda` or
`auto`. CMake gives the same two in `build/bin/` and keeps its objects in `build/` (`-DCCEX_BIN_DIR=$PWD/bin` puts the
binaries in `bin/`). `bin/`, `obj/`, `build/` and the `Makefile` are build output and not in git. Needed: a C++17
compiler behind an MPI wrapper (`mpicxx`, `mpic++`, `mpiicpc`) and, for the GPU binary, the CUDA toolkit (nvcc,
cuSOLVER, cuBLAS, cudart; [GPU execution](../methods/gpu.md)); Eigen, uthash and SQLite are bundled in `external/`.

**Presets.** Each preset loads the modules of its machine itself; the header of `do_compile.sh` describes each.

| preset | machine |
|---|---|
| `wsl-gnu` (old name `wsl`) | WSL/Ubuntu with Intel oneAPI: Intel MPI's `mpicxx` with g++, `I_MPI_FABRICS=shm` |
| `wsl-icpx` | the same with icpx behind Intel MPI (`apt install intel-oneapi-compiler-dpcpp-cpp`) |
| `idm4` | idm: module intel-23.2, `mpiicpc` with icpx |
| `gpu` | idm4-gpu: OpenMPI + CUDA 12.8 (icpx when oneAPI is there, else gcc); also `bin/main_gpu.out`, compute capability 12.0 by default |
| `nurion` | KISTI Nurion: intel/19.0.5 + impi/19.0.5 |
| `mac` | macOS with an MPI wrapper (Homebrew open-mpi or mpich) |
| `generic` | any machine with an MPI C++ wrapper on `PATH`; `CXX` (default `mpicxx`), `CC`, `CXXFLAGS`, and for the GPU binary `CUDA_HOME` (or nvcc on `PATH`) and `CUDA_ARCH` from the environment |

`do_compile.sh` also reads `CCEX_JOBS` (parallel compile jobs, default 5), `CCEX_BLAS` (below), `CXXFLAGS_EXTRA`,
`LIBRARY_EXTRA`, and `CCEX_SANITIZE=1` (address and undefined-behaviour checks).

**CMake options.** `-DCCEX_CUDA=OFF`, `-DCMAKE_CUDA_ARCHITECTURES=<cc>`, `-DMPI_CXX_COMPILER=<wrapper>`,
`-DCCEX_NATIVE=OFF`, `-DCCEX_SANITIZE=ON`, `-DCCEX_BLAS=OFF|MKL`, `-DCCEX_BIN_DIR=<dir>` and
`-DCCEX_STRICT_WARNINGS=ON` (warnings on, format, sequence-point and misleading-indentation ones as errors; the CI
builds with it).

**Intel MKL.** Linked when the build finds it (`$MKLROOT` from setvars.sh or the intel module, with its libraries
there); the build prints `[build] MKL on (...)` or `MKL off (why)`, and every run log's header says which. It takes the
large cCCE/gCCE clusters: zgemm in the recursive path of the cCCE ensemble Hahn echo from bath dimension 21, zheevd from
dimension 63 (the bath dimension for cCCE, the whole qubit × bath dimension for gCCE). On idm it ran h-BN in 0.72–0.88
of the time and left P1 unchanged; the curves move in the last digits (up to 9e-8 at late times of h-BN gCCE).
`CCEX_BLAS=off bash do_compile.sh <preset>` (CMake `-DCCEX_BLAS=OFF`) builds without it, `CCEX_BLAS=mkl`
(`-DCCEX_BLAS=MKL`) stops when it is missing.

**The binary of a CMake build** is in `build/bin/`, so tell the tutorials (`export CCEX_BIN=<repository>/build/bin/main.out`,
and `CCEX_GPU_BIN` for `08_GPU`) and the checks (`CCEX_BIN=$PWD/build/bin/main.out CCEX_BUILD_DIR=$PWD/build`; the
unit-test drivers take the compiler, flags and objects from the `unit-build.env` the build writes, any generator,
through `test/unit_build_env.sh`).

<hr>

<a id="HowToRun"></a>
## How to run

Run the binary under MPI with the input file:

```
mpirun -n 4 bin/main.out -f ccein.json
```

- `-n 4`: the number of MPI processes.
- `bin/main.out`: the binary (`build/bin/main.out` from CMake, `bin/main_gpu.out` for the CUDA backend).
- `ccein.json`: the [input file](InputFile.md), JSON with `//`, `#` or `!` comments (`.json` or `.jsonc`).

**Command-line options**
- `-h`, `--help` : Print a short usage and exit
- -f FILE : [Input file](InputFile.md) option (Required to calculate) (See [Input file](InputFile.md))
- -o DIR : Output folder, made when it is missing. The outputs are `DIR/<output.file>`, or `DIR/<input file name
  without .json or .jsonc>` when the input has no [output.file](InputFile_Output.md#tag_outfile); an absolute `output.file` together
  with `-o` is refused. Without `-o` the outputs are in [output.path](InputFile_Output.md#tag_output_path), else in the
  working directory (`-o DIR` > `output.path` > the working directory).
- -v : Verbose output (longer log)
- -e : Write the input template `ccex_example.jsonc` in the working directory (in `DIR` with `-o DIR`) and stop. Its
  uncommented part is the setting of `example/tutorial/00_quickstart` (NV-<sup>13</sup>C cCCE2 Hahn echo, ensemble,
  rbath 50 Å, 2 ms) with a `bathfile` for you to fill in (the template does not run until it names a bath file, a
  ready one is `example/tutorial/00_quickstart/bath/bath_Diamond_NV_natab_1`); `gyrofile` is commented (the built-in
  isotope table). Every other input key is a commented line with its allowed values and its default, the comments are
  `//`. An existing file is not overwritten; `-e` with `-f` is refused.

Every other setting of a run is a key of the input file. The options below, which older versions took, are refused
with the key that replaces them:

| removed option | input key |
|---|---|
| `-I FILE` | [configuration.bathfile](InputFile_General.md#tag_bathfile): `["FILE"]` |
| `-N N` | [configuration.nstate](InputFile_General.md#tag_nstate) |
| `-m METHOD` | [method](InputFile_General.md#tag_method) (with [cluster.mode](InputFile_Cluster.md#tag_cluster_mode)) |
| `-q QUANTITY` | [quantity](InputFile_General.md#tag_quantity) |
| `-B BZ` | [bfield](InputFile_General.md#tag_bfield) (a number is Bz) |
| `-s FILE` | [configuration.bath_state_file](InputFile_General.md#tag_bath_state_file): `["FILE"]` |
| `-S FILE` | [configuration.defect_state_file](InputFile_General.md#tag_defect_state_file): `["FILE"]` |
| `-a FILE` | [configuration.defect_axis_file](InputFile_General.md#tag_defect_axis_file): `["FILE"]` |
| `-P` | [output.cluster_only](InputFile_Output.md#tag_cluster_only): `true` |
| `-o FILE` (the output name) | [output.file](InputFile_Output.md#tag_outfile); `-o` names the output folder |
| `--auto_conv` | an [auto_convergence](InputFile_AutoConvergence.md) block (`main.out -f FILE` starts it) |
| `--auto_status` | `PYTHONPATH=example/tools/AutoConvergence python3 -m autoconv status FILE` (the tool's own command) |
| `--adaptive ...` | a [cluster.adaptive](InputFile_Cluster.md#tag_adaptive) block |

Two run modes take the input file and replace the ordinary calculation: an `auto_convergence` block in the input runs
the automatic convergence (See [Automatic convergence tags](InputFile_AutoConvergence.md); it takes no option but
`-f FILE`), and an `adaptive` block the adaptive refinement (See [cluster.adaptive](InputFile_Cluster.md#tag_adaptive);
`-o` does not apply to it, its steps go to the block's `output`).

<a id="EnvironmentVariables"></a>
**Environment** **variables**

None is needed for an ordinary run. They switch an optimisation off for a comparison, point the automatic convergence
at its tool, or serve the checks; set them in front of the command (`CCEX_SMALL_EIG=0 mpirun -n 32 main.out -f ccein.json`)
or export them in a job script.

| variable | when it is read | effect |
|---|---|---|
| `CCEX_CPU_BATCH=0` | every run | the CPU batch of small cCCE clusters (bath dimension 2 and 4) off: each cluster alone. The batch gives the same bits as that path; where it would not (another compiler or its options), it checks itself at start and stays off by itself, with a line in the log |
| `CCEX_SMALL_EIG=0` | every run | the fixed-sweep Jacobi eigensolver of the same small clusters off: Eigen solves them as well (the reference for a comparison). Matrices outside a safe range or not converged already go to Eigen one by one |
| `CUDA_VISIBLE_DEVICES` | GPU build (`main_gpu.out`) | which card the run may use; empty or not matching any card hides them all, and `execution.backend` `auto` then stays on the CPU |
| `CCEX_PYTHON` | a run with an `auto_convergence` block | the Python that runs the automatic-convergence tool (default `python3`; it needs numpy, scipy, matplotlib) |
| `CCEX_AUTOCONV_DIR` | a run with an `auto_convergence` block | the tool's folder when it is not found next to the binary (`example/tools/AutoConvergence`); a folder others can write to is refused |
| `CCEX_CHECK_DISORDER=1` | every run with `nstate` > 0 | the checks only: also builds the static mean field of the bath on one rank the old way and stops on any bit that differs from the split one |
| `CCEX_SMALL_CLUSTER_DUMP=FILE`, `CCEX_SMALL_CLUSTER_DUMP_N=N` | cCCE runs | the checks only: write the first N (400) small clusters with the CPU result to FILE (the CPU batch is off meanwhile) |
| `CCEX_WARP_CLUSTER_DUMP=FILE`, `CCEX_WARP_CLUSTER_DUMP_N=N` | gCCE runs | the checks only: the same for the clusters of the GPU warp kernel |
| `CCEX_TEST_GPU_OOM_BYTES=N` | GPU build | the checks only: every GPU allocation from N bytes fails as if the card were full, to exercise the fall-back to the CPU |

The variables `do_compile.sh` reads at build time are in [Building](#Building).

<hr>
<a id="Workflow"></a>
## Workflow

A run goes through these steps:

1. The command line: the input file (`-f`), the output folder (`-o`) and the log length (`-v`).
2. Every key of the input file, grouped as in [Input file](InputFile.md); a likely typo of a key or a combination
   that is not supported stops the run here, before any file is read, and any other key the run does not read is
   reported in the log and ignored.
3. The files the input names: the qubit file, the gyromagnetic ratios (or the built-in isotope table), the bath files,
   hyperfine and quadrupole tensor files, and the bath states and defect axes when they are given.
4. The clusters of the bath (`cluster.mode`, `order`), their propagation (CPU or `execution.backend` `cuda`) and the
   product of the cluster contributions, for each bath configuration.
5. The output files (`output.path`, `output.file`): the coherence or the density matrix as a function of time.

<img src="CCEX_workflow.png" alt="CCEX workflow" width="600">
