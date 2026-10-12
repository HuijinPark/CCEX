# 06_Pulse_Manual

Pulse sequences on the natab and P1 systems. explicit pulse sequences (`sequence` in the input): Hahn echo, WAHUHA, CPMG. The natab calculations read `bath/` here (one configuration; the NV defect file is 02's, same cell), the P1 ones read `../04_Diamond_NV_P1Bath/bath/`.

| Calculation | Inputs | What it shows | Recorded run | 2026-09-23 on WSL, 2 ranks, 1 configuration |
|---|---|---|---|---|
| `CPMG_P1_1ppm/` | ccein_origin.json, ccein_test.json, ccein_wahuha.json | CPMG with eight pulses on the 1 ppm P1 bath: the original input, a test variant and the WAHUHA-style one | 5 conf × 1 states, 128 ranks | 34 s |
| `HahnEcho_natab_ensemble_origin/` | ccein.json | CCE2 Hahn echo written as an explicit pulse sequence, ensemble over configurations | 1 conf, 5 ranks | 6 s |
| `HahnEcho_natab_ensemble_wahuha/` | ccein.json | the same sequence given in the WAHUHA style of the input | 1 conf, 5 ranks | 6 s |
| `HahnEcho_natab_single_origin/` | ccein.json | the explicit Hahn echo with one sampled bath state | 1 conf × 1 states, 5 ranks | 5 s |
| `HahnEcho_natab_single_wahuha/` | ccein.json | the same, WAHUHA-style input | 1 conf × 1 states, 5 ranks | 5 s |
| `SpinEcho_P1_1ppm/` | ccein_origin.json, ccein_wahuha.json | the spin echo on the 1 ppm P1 bath as an explicit sequence, original and WAHUHA-style | 5 conf × 1 states, 128 ranks | 12 s |
| `WAHUHA_P1_1ppm/` | ccein_manual.json, ccein_test.json | WAHUHA on the 1 ppm P1 bath, the manual sequence and a test input | 2 conf × 1 states, 128 ranks | 16 s |
| `WAHUHA_natab_single_ccce/` | ccein.json | the four-pulse WAHUHA sequence, CCE2 | 1 conf × 1 states, 5 ranks | failed (1 s); refused by the current code: method CCE treats every pulse as a pi flip and WAHUHA has 90-degree pulses; kept as recorded, the gcce folder is the one to run |
| `WAHUHA_natab_single_gcce/` | ccein.json | the four-pulse WAHUHA sequence, gCCE2 | 1 conf × 1 states, 5 ranks | 7 s |

Each calculation folder holds its `ccein*.json` inputs, `run.sh` and the idm `job`
that submits it. `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
