# 05_Diamond_NV_VpBath

NV⁻ in a Vp defect bath. `bath/<ppm>ppm/` holds the configurations for 10, 12 and 50 ppm.

| Calculation | Inputs | What it shows | Recorded run | 2026-09-23 on WSL, 2 ranks, 1 configuration |
|---|---|---|---|---|
| `SpinEcho_CCE2/` | ccein_10ppm.json, ccein_12ppm.json, ccein_50ppm.json | CCE2 spin echo of the NV in a Vp defect bath at 10, 12 and 50 ppm | 1 conf × 1 states, 128 ranks | 20 s |

Each calculation folder holds its `ccein*.json` inputs, `run.sh` and the idm `job`
it was recorded with. `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
