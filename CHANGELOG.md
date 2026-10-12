# Changelog

## 2.0.0 — unreleased (changes since v1.2.2, 2026-09-04)

### Added
- **GPU backends** (`backend`: `cpu`, `cuda`, `auto`; `gpu_min_dim`, `gpu_batch_clusters`) for cCCE, gCCE and pCCE, with
  `propagator` (`eigen`, `expm`), `evolution` (`vector`, `matrix`) and `gcce_refinement` for gCCE. See `docs/methods/gpu.md`.
- **Adaptive clustering** (`cluster_mode: adaptive` with an `Adaptive` input block): bundles of strongly coupled spins with
  the merge policy, and `"regroup": true` in that block (command line `--regroup`) to regroup where a bundle expansion
  gives |L| > 1.
  See `docs/methods/adaptive-cluster-mode.md` and `docs/methods/adaptive-input-block.md`.
- **pCCE**: a deterministic partition search, a partition cache (`pcce_cache_file`, `pcce_cache_mode`) and
  `pcce_connections` (`centers`, `spins`).
- **Automatic convergence**: an `auto_convergence` block in the input (`main.out -f input.json`; `"enabled": false` runs
  the ordinary calculation; command line `--auto_conv`) chooses rbath, rdip, the time window and the number of bath
  configurations and states (`example/tools/AutoConvergence/README.md`).
- **`rdip_ladder`**: one order-2 run (cCCE, gCCE or pCCE) also writes the curves of smaller rdip values.
- **`evolution_update`**: `exact` (default) or `recursive` (the phases by a recurrence; 1.0–1.6x faster on the tested systems,
  where the curves stay within 1e-9 of `exact`; not bit-identical).
- **Built-in isotope table**: `gyrofile` is optional; species it does not name come from the table.
- **`save_nodiv`**: `false` writes `_wiDiv` alone.
- **Builds**: `do_compile.sh` presets (`wsl-gnu`, `wsl-icpx`, `idm4`, `gpu`, `nurion`, `mac`, `generic`) and a CMake build.

### Changed
- **Faster CPU paths with unchanged results** (byte-identical output): the cCCE and gCCE product assembly, the state preparation,
  fixed-size paths for small clusters, the order-2 cluster list built without duplicated work, the small-cluster
  Hamiltonians, and the pCCE partition search (up to 9x faster on large baths). See Performance below.
- **Repository layout**: `src/` and `include/` in five parts (input, physical, clustering, calculation, common);
  `example/{tutorial,analysis,tools}`; `external/`; `test/{code_verification,cluster_verification,gpu_verification}`.
  The tutorials were rewritten as ten numbered folders with `run.sh`.
- In input files, `#` and `!` start a comment only outside JSON strings.

### Performance (measured during the October 2026 work)

Each row compares the code before and after that change, on the same inputs and server (idm, 64 ranks unless noted);
the results are byte-identical, or within the stated bound for `recursive`. A comparison with v1.2.2 is not measured yet.

| What | Case | Before -> after | Speed-up | Baseline |
|---|---|---|---|---|
| P1 production cCCE2, 30 states | CR60 sigma 4 | 229 s -> 85.9 s | 2.7x | code of 2026-10-09 |
| | CR0 | 21.3 s -> 16.5 s | 1.3x | code of 2026-10-11 |
| P1 production gCCE2, 30 states | CR60 sigma 4 | 313.3 s -> 224.7 s | 1.4x | code of 2026-10-11 |
| | CR0 | 56.1 s -> 41.2 s | 1.4x | code of 2026-10-11 |
| pCCE partition (P1, 30,000 trials) | rbath 483 A | 1,616 s -> 173 s | 9.3x | code of 2026-10-11 |
| Automatic convergence, h-BN gCCE2 ensemble | whole run | 14,406 s -> 4,281 s (CCEX runs 596 -> 103) | 3.4x | code of 2026-10-11 |
| Automatic convergence, pCCE runs differing only in rdip | NV-13C, 4 ranks (local) | 8.2 s -> 0.3 s | 27x | code of 2026-10-11 |
| `evolution_update: recursive` (option) | 16 inputs (NV-13C, h-BN, P1) | | 1.0-1.6x | `exact` |

### Fixed
- Ensemble cCCE: the conjugation of the coherence and the half interval computed in double precision.
- Density-matrix output: a component that divides noise by zero is set to 0 and reported, instead of `nan`.
- `seed` had no effect on sampled states when `nstate` > 0.
- Builds with Intel `icpx`: the adaptive code is compiled without fast floating point, which had folded `isnan`/`isinf`.

### Refused now (inputs that used to be accepted)
- `alphams` at the top of a `qubitfile` (it had been ignored since v1.1.0).
- Arrays shorter than needed (Defect types shorter than `naddspin`, `bathadjust` rows with fewer than 3 numbers).
- An `avaax` axis beyond the defect's `navaax`.
- Combinations a method does not support; `test/code_verification/check_refusals` lists them.

### Removed
- The MKL dependency, `hf_cutoff`, the bundled `zlib/`, the top-level `main.cpp` and `config/` (moved into `src/` and the tools).
