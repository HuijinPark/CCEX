# example/

`tutorial/` is what you run, `analysis/` is what you run on the output, `tools/` is what makes the input.

## Start here

    bash do_compile.sh wsl-gnu                     # 1. build bin/main.out (presets: wsl-gnu, wsl-icpx, idm4, gpu, nurion, mac, generic; or CMake)
    cd example/tutorial/00_quickstart && NP=4 bash run.sh     # 2. run (about 1 s): results/ccein_wiDiv
    #                                                3. run.sh plots |L(t)| to results/coherence.png and prints T2
    #                                                4. compare with expected.png in the same folder (T2 1.068 ms)

`bin/main.out -e` writes an input template, `ccex_example.jsonc`: 00_quickstart's settings with every other key
commented and its default beside it, to start an input of your own (it runs once its `bathfile` names a bath file).

Then follow `tutorial/` in its order (`tutorial/README.md` has the list); every topic's `README.md` says what each
calculation shows, which keys matter, how long it runs and what to change next.

## What should I use for my case?

| Situation | Method and keys | Example |
|---|---|---|
| a defect in a nuclear-spin bath (¹³C, ²⁹Si, ...), high field along the defect axis | cCCE (the default), `order` 2, ensemble | `00_quickstart`, `01_ensemble_and_states/cCCE2_Hahn_ensemble` |
| several bath configurations from one input | `configuration.bathfile` with `{k}` and `configuration.k_range` | `01_ensemble_and_states/cCCE2_Hahn_ensemble` |
| the bath in sampled pure states (electron baths, polarised or single-shot) | `configuration.nstate` > 0 and `seed` (or `bath_state_file`) | `01_ensemble_and_states/cCCE2_Hahn_states`, `03_P1` |
| dynamical decoupling | `pulse.count` (1 Hahn, N CPMG-N), `pulse.name` `Manual` with `pulse.sequence` for any other π sequence | `02_pulses` |
| an electron-spin bath (P1, other paramagnetic defects) | the `Defect` block (axes, detunings, on-site nuclei), `nstate` > 0, check `order` 4 | `03_P1` |
| DFT hyperfine and quadrupole tensors (2D materials, defects with strong contact terms) | `hamiltonian.common.hyperfine` / `quadrupole` with the files of `tools/SpinTensorReader` | `04_hBN` |
| low or tilted field, a spin-1 qubit whose zero-field splitting matters, the qubit density matrix | gCCE (`method` gCCE, `Qubit.zfs` or `intmap`), `quantity` dm | `06_cCCE_vs_gCCE/gCCE2_Hahn_vsCCE`, `03_P1/gCCE2_Hahn_zfs` |
| cCCE in a tilted field | `hamiltonian.cCCE.knight` | `06_cCCE_vs_gCCE/cCCE2_Hahn_knightShift` |
| several qubits | `Qubit.qubit` list with `intmap`, gCCE | `05_multiNV` |
| strongly coupled groups of bath spins (dense electron baths) | pCCE (`cluster.mode` pcce, `partition_size`) or bundles (`cluster.mode` bundle, `bundle.max_dim`, optionally `cluster.adaptive`) | `07_pCCE_and_bundle` |
| large clusters (Hilbert dimension of hundreds or more) and a GPU | `execution.backend` auto with `bin/main_gpu.out` | `08_GPU` |
| not sure which rbath, rdip, time window or number of configurations | the `auto_convergence` block | `09_auto_convergence` |
| a bath file written in the crystal frame (defect axis not along z) | the `frame.rotation` block, `qubit_axis` | `10_coordinate_rotation` |

## analysis/ and tools/

| Folder | Does |
|---|---|
| `analysis/Analyzer/` | `CoherenceAnalyzer.v1/run.py` reads the `_noDiv` and `_wiDiv` files, corrects the division errors (`_ErrCorr`), averages over configurations (`EnsAvg`), fits and plots; `coherence.py` is the same as a library. `test/code_verification/check_input_forms` (`errcorr.sh`) calls `run.py` |
| `tools/AutoConvergence/` | the automatic convergence: `mpirun -n N main.out -f input.json` with an `auto_convergence` block runs it (time window, rdip, rbath, configurations, then the production); its `README.md` has the keys and outputs |
| `tools/BathGenerator/` | `BathGenerator.py`: a bath file from a POSCAR unit cell (`tutorial/make_natab_baths.sh` calls it) |
| `tools/SpinTensorReader/` | `main.py` with the `spHamilReader/` package: hyperfine and quadrupole tensors from VASP or QE output |
| `tools/IsotopeTable/` | `make_isotope_table.py`: writes `src/input/isotope_table.cpp` (the isotopes CCEX knows without a gyrofile) |
| `tools/InputKeys/` | `convert_input_keys.py` rewrites inputs written with older keys; `set_input_keys.py` writes a copy of an input with keys set |

Removed on 2026-10-16 (git history has them): `tools/Runner/` (job folders per variable value) and
`analysis/ConvTest/` (convergence over one variable); `configuration.k_range` and the automatic convergence took their
place. The adaptive-refinement study scripts moved from `analysis/AdaptiveStudy/` to (internal notes).
