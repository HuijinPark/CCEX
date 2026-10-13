# 03_Diamond_multiNV_natab

Several NV centres with a ¹³C pair. two and six NV centres; each calculation folder has its own small `bath` and `state1` files.

| Calculation | Inputs | What it shows | Recorded run | 2026-09-23 on WSL, 2 ranks, 1 configuration |
|---|---|---|---|---|
| `NV2_13C2/` | ccein_NV6_13C2_hahnecho.json, ccein_NV6_13C2_ramsey.json | two NV centres with a 13C pair: gCCE density matrix, Ramsey (ccein_NV6_13C2_ramsey.json) and Hahn echo | 1 conf × 1 states, 1 ranks | 1 s |
| `NV6_13C2/` | ccein_NV6_13C2_hahnecho.json, ccein_NV6_13C2_ramsey.json | six NV centres with a 13C pair: the same two | 1 conf × 1 states, 1 ranks | stopped at 631 s; heavy: the Hahn echo alone took 8 min on 2 ranks and wrote 11 GB of density matrices; the Ramsey run was cut by the 10-minute limit |

Each calculation folder holds its `ccein*.json` inputs, `run.sh` and the idm `job`
that submits it. `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

Here both `run.sh` compute the one fixed `bath` of their folder and ignore `NCONF`.

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
