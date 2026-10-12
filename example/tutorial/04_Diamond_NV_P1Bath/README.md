# 04_Diamond_NV_P1Bath

NV⁻ in a P1 electron-spin bath. `bath/14N_<ppm>ppm/` holds the P1 configurations (bath, defect, and the DFT tensor sets), `bath/DiaP1_gyro` the gyromagnetic ratios. `bath.tar.gz` is the same folder packed.

| Calculation | Inputs | What it shows | Recorded run | 2026-09-23 on WSL, 2 ranks, 1 configuration |
|---|---|---|---|---|
| `Population_gCCE2_10ppm/` | ccein.js | gCCE2 population (density matrix, `npulse=0`) of the NV in a 10 ppm P1 bath | 1 conf, 4 ranks | 208 s; had no job and a trailing comma in the input (fixed 2026-09-23) |
| `SpinEcho_CCE2_1ppm/` | ccein.js | CCE2 spin echo, 1 ppm; the recorded run sampled 50 bath states per configuration (`-N 50`) | 50 conf × 50 states, 128 ranks | 5 s |
| `SpinEcho_CCE2e2n_1ppm/` | ccein.js | CCE2 with the electron-electron-nuclear clusters (CCE2e2n), 1 ppm | 50 conf × 50 states, 128 ranks | 158 s |
| `SpinEcho_CCE4_1ppm/` | ccein.js | CCE4, 1 ppm: the convergence check against CCE2 (heavy) | 50 conf × 50 states, 128 ranks | 25 s |
| `SpinEcho_gCCE2_10ppm_ZFS/` | ccein_with_ZFS.js, ccein_with_ZFS_v2.js, ccein_without_ZFS.js | gCCE2 at 10 ppm and zero field, with and without the NV zero-field splitting (`job` is the field-sweep job the two came from) | 1 conf × 1 states, 4 ranks | 10 s |
| `SpinEcho_pCCE2_1ppm/` | ccein_sK1.js, ccein_sK2.js, ccein_sK3.js, ccein_sK4.js | pCCE2 with sK = 1..4 spins per centre, 1 ppm | 50 conf × 50 states, 128 ranks | stopped at 635 s |

Each calculation folder holds its `ccein*.js` inputs, `run.sh` and the idm `job`
it was recorded with. `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
