# Changelog

## 2.0.0 — unreleased (changes since v1.2.2, 2026-09-04)

### Added
- **GPU backends** (`backend`: `cpu`, `cuda`, `auto`; settings in the `gpu_options` block) for cCCE, gCCE and pCCE, with
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
- **`evolution_update`**: `exact` or `recursive` (the phases by a recurrence; 1.0–1.6x faster on the tested systems,
  where the curves stay within 1e-9 of `exact`; not bit-identical). The default is `recursive` where it exists (see Changed).
- **Built-in isotope table**: `gyrofile` is optional; species it does not name come from the table.
- **`save_nodiv`**: `false` writes `_wiDiv` alone.
- **Packed GPU kernel for small cCCE clusters** (`gpu_options.packed_max_dim`: 4 by default in `bin/main_gpu.out` with
  `backend` `cuda` or `auto` and method CCE, 0 with `backend` `cpu` and in `bin/main.out`; 2 or 4 may be given): every cCCE
  cluster of bath dimension 2 and 4 is computed by one GPU thread, thousands per kernel, for the whole node through one CUDA
  context (node rank 0 feeds the GPU from node shared memory), the rest unchanged; single state and ensemble, Hahn, CPMG,
  Ramsey, `exact` and `recursive`; `bin/main_gpu.out` and `method` CCE. Within 1e-8 + 2e-6 |CPU| of the CPU, not
  bit-identical. On idm4-gpu (32 ranks, one GPU) the production P1 cCCE2 runs took 0.59–0.84 of the CPU's time; a run of a
  few seconds gains nothing (CR0 ensemble, 1.8 s: 0.95). See `docs/user/InputFile_Backend.md#tag_gpu_options` and
  (internal notes).
- **Builds**: `do_compile.sh` presets (`wsl-gnu`, `wsl-icpx`, `idm4`, `gpu`, `nurion`, `mac`, `generic`) and a CMake build.

### Changed
- **Faster CPU paths with unchanged results** (byte-identical output): the cCCE and gCCE product assembly, the state preparation,
  fixed-size paths for small clusters, the order-2 cluster list built without duplicated work, the small-cluster
  Hamiltonians, and the pCCE partition search (up to 9x faster on large baths). See Performance below.
- **Repository layout**: `src/` and `include/` in five parts (input, physical, clustering, calculation, common);
  `example/{tutorial,analysis,tools}`; `external/`; `test/{code_verification,cluster_verification,gpu_verification}`.
  The tutorials were rewritten as ten numbered folders with `run.sh`.
- In input files, `#` and `!` start a comment only outside JSON strings.
- **GPU settings moved into one block, `gpu_options`**: `min_dim` (was `gpu_min_dim`), `evolution_batch_size`
  (`gpu_batch_clusters`), `evolution_batch_steps` (`gpu_batch_steps`), `cluster_batch_mb` (`gpu_cluster_batch_mb`),
  `cache_mb` (`gpu_cache_mb`), and the new `packed_max_dim`. The old flat keys are refused with a message naming the new
  form (e.g. `gpu_min_dim is now gpu_options.min_dim`); `backend` keeps `cpu`, `cuda`, `auto`.
- **`bin/main_gpu.out` with `backend` `cuda` or `auto` runs the cCCE clusters of bath dimension 2 and 4 in the packed
  kernel by default** (`gpu_options.packed_max_dim` 4; `"packed_max_dim": 0` gives the earlier path). Its results agree with
  the CPU within 1e-8 + 2e-6 |CPU| and are not bit-identical to those of the earlier GPU path. `backend` `cpu` and
  `bin/main.out` are unchanged.
- **The default of `evolution_update` is `recursive` where it exists** (the CPU cCCE, single state and ensemble, and the CPU
  gCCE with `propagator` = `eigen`, `evolution` `vector` or `matrix`), and `exact` elsewhere (GPU backends, gCCE `expm`),
  without a message; the run log says which and why. The outputs are therefore not bit-identical to earlier runs (within
  1e-9 on the tested curves at |L| ≤ 1). `"evolution_update": "exact"` gives the earlier bits.
- **Clusters of an order are split over the ranks by cost when their sizes differ** (dimension D, cost D²(D + nstep); each
  rank takes a contiguous run of the list): mixed baths such as h-BN run 1.6–3x faster on 64 ranks (cCCE2 1.6x, gCCE2 2.3x,
  cCCE3 3.1x). Orders whose clusters all have one size (P1, ¹³C) keep the equal-count split, byte for byte; in a mixed bath
  the product of each rank's clusters is grouped differently, so the last printed digit can change (h-BN gCCE2 dm: 5 of
  36,000 numbers, 1e-10). Data (internal notes).

### Compared with v1.2.2 (same inputs and node, idm Xeon Gold 6530, 64 ranks, median of 3)

| Case | v1.2.2 | 2.0.0 (`exact`) | 2.0.0 (`recursive`) |
|---|---|---|---|
| P1 production cCCE2, CR0 | 73.2 s | 16.7 s (4.4x) | 10.7 s (6.8x) |
| P1 production cCCE2, CR60 sigma 4 | 452.1 s | 84.7 s (5.3x) | 53.3 s (8.5x) |
| P1 production gCCE2, CR0 | 91.0 s | 41.4 s (2.2x) | 31.1 s (2.9x) |
| P1 production gCCE2, CR60 sigma 4 | 546.4 s | 224.9 s (2.4x) | 168.0 s (3.3x) |

The coherence differs from v1.2.2 only where a fix in this release changes it (deltat read as a float in v1.2.2,
the sampled states and `seed`, the ensemble half interval in double, the density-matrix `nan`, the pCCE partition).

### Performance (measured during the October 2026 work)

Each row compares the code before and after that change, on the same inputs and server (idm, 64 ranks unless noted);
the results are byte-identical, or within the stated bound for `recursive`. For v1.2.2 see the table above.

| What | Case | Before -> after | Speed-up | Baseline |
|---|---|---|---|---|
| P1 production cCCE2, 30 states | CR60 sigma 4 | 229 s -> 85.9 s | 2.7x | code of 2026-10-09 |
| | CR0 | 21.3 s -> 16.5 s | 1.3x | code of 2026-10-11 |
| P1 production gCCE2, 30 states | CR60 sigma 4 | 313.3 s -> 224.7 s | 1.4x | code of 2026-10-11 |
| | CR0 | 56.1 s -> 41.2 s | 1.4x | code of 2026-10-11 |
| pCCE partition (P1, 30,000 trials) | rbath 483 A | 1,616 s -> 173 s | 9.3x | code of 2026-10-11 |
| Automatic convergence, h-BN gCCE2 ensemble | whole run | 14,406 s -> 4,281 s (CCEX runs 596 -> 103) | 3.4x | code of 2026-10-11 |
| Automatic convergence, pCCE runs differing only in rdip | NV-13C, 4 ranks (local) | 8.2 s -> 0.3 s | 27x | code of 2026-10-11 |
| Product assembly, P1 production cCCE2, 30 states, `recursive` | CR60 sigma 4 | 55.7 s -> 45.7 s (assembly per rank 12.9 s -> 4.9 s) | 1.2x | code of 2026-10-12 |
| | CR0 | 10.8 s -> 9.4 s (assembly per rank 2.3 s -> 0.9 s) | 1.15x | code of 2026-10-12 |
| Product assembly, P1 production gCCE2, 30 states, `recursive` | CR60 sigma 4 | 167.9 s -> 162.2 s (assembly per rank 17.3 s -> 10.0 s) | 1.04x | code of 2026-10-12 |
| | CR0 | 31.0 s -> 29.9 s (assembly per rank 3.1 s -> 1.8 s) | 1.04x | code of 2026-10-12 |
| `evolution_update: recursive` (default where it exists) | 16 inputs (NV-13C, h-BN, P1) | | 1.0-1.6x | `exact` |

### Fixed
- Ensemble cCCE: the conjugation of the coherence and the half interval computed in double precision.
- Density-matrix output: a component that divides noise by zero is set to 0 and reported, instead of `nan`.
- `seed` had no effect on sampled states when `nstate` > 0.
- Builds with Intel `icpx`: the adaptive code is compiled without fast floating point, which had folded `isnan`/`isinf`.
- **`hfmedi` (hyperfine-mediated interaction)**: the term is now the complete second-order correction of each branch,
  H_med^(α) = Σ_{β≠α} E_αβ E_αβ† / (ω_α − ω_β) with E_αβ = Σ_i ⟨α|S|β⟩·A_i·I_i: every intermediate eigenstate β of the
  qubit Hamiltonian (it used only the other qubit level, missing for example the |+1⟩ path of the NV |0⟩ branch), every
  ordered pair (i, j) (it summed i < j only, which was not Hermitian for a spin-1 qubit: the CPU used one triangle of it
  and the GPU refused it), and the i = j terms. Hermitian by construction. Results with `hfmedi: true` change.
  An intermediate level degenerate with the branch is skipped when S does not couple them (spin-1 at zero field: +1 and
  -1, which gave `nan`) and refused when it does.
  The h-BN inputs that use it (tutorial 07 `hyperfine_mediated_term`) had an asymmetric `qzfs` typo, now
  `[[-1080000,0,0],[0,-1180000,0],[0,0,2320000]]`; with both, the `check_numerics` case `hbn_hfmedi` moves by up to 0.95.
- **`knight` (Knight shift)**: the term multiplies by the nuclear Zeeman vector −γ_n B instead of the bare field B, and
  the (3|m_s| − 2) factor takes |m_s| of each branch as ⟨S_z²⟩ (it was always 1). Zero for a field along the qubit axis,
  as before; in a tilted field results with `knight: true` change (`natab_knight_shift`: up to 0.98).

### Refused now (inputs that used to be accepted)
- `alphams` at the top of a `qubitfile` (it had been ignored since v1.1.0).
- Arrays shorter than needed (Defect types shorter than `naddspin`, `bathadjust` rows with fewer than 3 numbers).
- An `avaax` axis beyond the defect's `navaax`.
- Combinations a method does not support; `test/code_verification/check_refusals` lists them.

### Removed
- The MKL dependency, `hf_cutoff`, the bundled `zlib/`, the top-level `main.cpp` and `config/` (moved into `src/` and the tools).
