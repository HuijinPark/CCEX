# Overview

CCEX computes the coherence of spin qubits in a bath of nuclear and electron
spins by the cluster-correlation expansion (CCE). It is a C++ and MPI code of the
CQML group at Sungkyunkwan University: <https://github.com/CQML-SKKU/CCEX>.

## What it computes

The input names a qubit (or several), the bath around it and a pulse sequence.
From the positions, gyromagnetic ratios and interaction tensors CCEX builds the
spin Hamiltonian, divides the bath into clusters, solves the dynamics of every
cluster exactly and combines them into the quantity asked for as a function of
time:

- the coherence function `L(t)` of the qubit under the pulse sequence, or with
  gCCE its density matrix (`quantity`);
- by conventional CCE, generalized CCE (gCCE) or partition CCE (pCCE), with the
  disjoint-cluster and inter-bath variants dsj, itb and dsjitb (`method`);
- for a Hahn echo, a Ramsey experiment, CPMG (`npulse`) or any sequence written
  out pulse by pulse (`sequence`);
- for one bath configuration or an ensemble of them, and for electron-spin baths
  with sampled bath states (`nstate`).

Baths can be nuclear spins with point-dipole or DFT tensors (hyperfine and
quadrupole), electron-spin defects such as P1 centres with their own on-site
spins and axes, or both. Beyond the standard clustering by order and distance
CCEX has bundle clustering (`cluster_mode`, `bundle_max_dim`) and an adaptive
refinement of the clusters (the `Adaptive` block), a CUDA backend for the
cluster propagation (`backend`), and a coordinate frame rotation that lets the
bath file stay in the crystal frame (`coordinate_frame_rotation`).

## Building

```
cmake -B build && cmake --build build -j 8      # any machine: finds the MPI wrapper and, if there is one, nvcc
bash do_compile.sh <preset>                     # the group's machines: wsl-gnu, wsl-icpx, idm4, gpu, nurion, mac (or generic)
```

Both write `bin/main.out` and, where there is CUDA, `bin/main_gpu.out` for
`backend` cuda or auto. The build needs a C++17 compiler behind an MPI wrapper
and, for the GPU binary, the CUDA toolkit. Eigen, uthash and SQLite are
bundled in `external/`, and nothing else is linked. `CMakeLists.txt` takes
`-DCCEX_CUDA=OFF`, `-DCMAKE_CUDA_ARCHITECTURES=<cc>` and
`-DMPI_CXX_COMPILER=<wrapper>`; `do_compile.sh` loads the modules of the
group's machines itself, and its `generic` preset takes `CXX`, `CUDA_HOME` and
`CUDA_ARCH` from the environment.

## Running

```
mpirun -n 4 bin/main.out -f ccein.json
```

The input file is JSON with `#` and `!` comments; its keys are the reference in
[Input file](InputFile.md), and the command-line options that override them are
in [Getting started](GettingStarted.md). A run writes the coherence (or density
matrix) per bath configuration into text files, which the analyzer in
`example/analysis/` error-corrects, averages over the ensemble and plots.

## Where to go next

- [Getting started](GettingStarted.md): the command-line options, the input
  files (bath, gyromagnetic ratios, qubit, states, axes) and the order in which
  a run reads them.
- [Input file](InputFile.md): every key, by group.
- [Examples](Examples.md): ten worked examples in `example/tutorial/`, from a
  hand-written bath to the GPU backend and the frame rotation.
- [Analysis](Analysis.md) and [Tools](Tools.md): what runs on the output, and what
  makes the inputs (`example/analysis/`, `example/tools/`).
- Data structures and Source files: the code itself, for developers.
- `docs/methods/` in the repository: bundle and adaptive clustering, the GPU
  backend and the frame rotation in detail.

## Citing

The authors and the version are in `CITATION.cff` at the repository root.
