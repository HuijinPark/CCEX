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
- by conventional CCE (cCCE) or generalized CCE (gCCE), or their master-equation forms meCCE and megCCE
  (`method`), with partition CCE (pCCE), bundle
  clustering and the disjoint-cluster and inter-bath variants dsj, itb and dsjitb (`cluster.mode`);
- for a Hahn echo, a Ramsey experiment, CPMG (`pulse.count`) or any sequence written
  out pulse by pulse (`pulse.sequence`), with DEER pulses on the bath (`pulse.bath`);
- for one bath configuration or an ensemble of them, and for electron-spin baths
  with sampled bath states (`configuration.nstate`).

Baths can be nuclear spins with point-dipole or DFT tensors (hyperfine and
quadrupole), electron-spin defects such as P1 centres with their own on-site
spins and axes, or both. Beyond the standard clustering by order and distance
CCEX has bundle clustering (`cluster.mode` bundle, `cluster.bundle.max_dim`) and an adaptive
refinement of the clusters (the `cluster.adaptive` block), a CUDA backend for the
cluster propagation (`execution.backend`), and a coordinate frame rotation that lets the
bath file stay in the crystal frame (`frame.rotation`).

## Building

`bash do_compile.sh <preset>` or `cmake -B build && cmake --build build -j 8`; the presets, what is needed and the
options are in [Getting started](GettingStarted.md#Building).

## Running

```
mpirun -n 4 bin/main.out -f ccein.json
```

The input file is JSON (`.json` or `.jsonc`) with `//`, `#` and `!` comments; its keys are
the reference in [Input file](InputFile.md). Every setting of a run is a key of the input
file: the command-line options only name the input (`-f`), the output folder (`-o`, in
place of `output.path`), a longer log (`-v`), the input template (`-e`) and the usage (`-h`);
they are in [Getting started](GettingStarted.md#HowToRun). A run writes the coherence (or density
matrix) per bath configuration into text files, which the analyzer in
`example/analysis/` error-corrects, averages over the ensemble and plots.

## Where to go next

- [Getting started](GettingStarted.md): building, the command-line options, the
  environment variables and the order in which a run reads its input.
- [Input file](InputFile.md): every key, by group.
- [Examples](Examples.md): eleven worked examples in `example/tutorial/`, in learning
  order from a quickstart of seconds to the GPU backend, the automatic convergence and the
  frame rotation.
- [Analysis](Analysis.md) and [Tools](Tools.md): what runs on the output, and what
  makes the inputs (`example/analysis/`, `example/tools/`).
- Data structures and Source files: the code itself, for developers.
- Methods: [the coordinate frame rotation](../methods/coordinate-frame-rotation.md),
  [DEER bath pulses](../methods/deer.md), [master-equation CCE](../methods/master-equation.md) and
  [GPU execution](../methods/gpu.md), in detail (`docs/methods/` in the repository). Bundle clustering and
  the adaptive refinement are described with their keys in [Cluster tags](InputFile_Cluster.md).

## Citing

The authors and the version are in `CITATION.cff` at the repository root.
