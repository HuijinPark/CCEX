# 04_Diamond_NV_P1Bath

NV⁻ in a P1 electron-spin bath. `bath/14N_<ppm>ppm/` holds the P1 configurations (bath, defect, and the DFT tensor sets), `bath/DiaP1_gyro` the gyromagnetic ratios.

| Calculation | Inputs | What it shows | Recorded run | 2026-09-23 on WSL, 2 ranks, 1 configuration |
|---|---|---|---|---|
| `Population_gCCE2_10ppm/` | ccein.json | gCCE2 population (density matrix, `npulse=0`) of the NV in a 10 ppm P1 bath | 1 conf, 4 ranks | 208 s; the input had a trailing comma (fixed 2026-09-23); the folder has no `job` |
| `SpinEcho_CCE2_1ppm/` | ccein.json | CCE2 spin echo, 1 ppm; the recorded run sampled 50 bath states per configuration (`-N 50`) | 50 conf × 50 states, 128 ranks | 5 s |
| `SpinEcho_CCE2e2n_1ppm/` | ccein.json | CCE2 with the electron-electron-nuclear clusters (CCE2e2n), 1 ppm | 50 conf × 50 states, 128 ranks | 158 s |
| `SpinEcho_CCE4_1ppm/` | ccein.json | CCE4, 1 ppm: the convergence check against CCE2 (heavy) | 50 conf × 50 states, 128 ranks | 25 s |
| `SpinEcho_gCCE2_10ppm_ZFS/` | ccein_with_ZFS.json, ccein_with_ZFS_v2.json, ccein_without_ZFS.json | gCCE2 at 10 ppm and zero field, with and without the NV zero-field splitting | 1 conf × 1 states, 4 ranks | 10 s |
| `SpinEcho_pCCE2_1ppm/` | ccein_sK1.json, ccein_sK2.json, ccein_sK3.json, ccein_sK4.json | pCCE2 with sK = 1..4 spins per centre, 1 ppm | 50 conf × 50 states, 128 ranks | stopped at 635 s |

Each calculation folder holds its `ccein*.json` inputs, `run.sh` and the idm `job`
that submits it (`Population_gCCE2_10ppm/` has no `job`). `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

`Population_gCCE2_10ppm/run.sh` computes the one configuration its input names and ignores `NCONF`.

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
