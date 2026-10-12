# 02_Diamond_NV_natab

NV⁻ in diamond, natural-abundance ¹³C. the basic case; `bath/` holds one configuration of a 150 Å cell (the recorded ensembles used 20, made with `bath/genbath/`), the NV defect file and the ¹³C gyromagnetic ratio.

| Calculation | Inputs | What it shows | Recorded run | 2026-09-23 on WSL, 2 ranks, 1 configuration |
|---|---|---|---|---|
| `ensemble_CCE/` | ccein.js | CCE2 Hahn echo, one run per bath configuration, then the ensemble average | 20 conf, 5 ranks | 5 s |
| `ensemble_Dsj/` | ccein.js, ccein_cce.js | the same with `method=dsj` (ccein.js) beside plain CCE2 (ccein_cce.js) | 1 conf, 5 ranks | 8 s |
| `ensemble_gCCE/` | ccein_multi2.js, ccein_single.js | gCCE2 on one configuration, single qubit (ccein_single.js) and two qubits (ccein_multi2.js), plotted together | 1 conf, 20 ranks | 30 s; ccein_multi2.js had a missing comma (fixed 2026-09-23) |
| `knight_shift/` | ccein_wi_knight.js, ccein_wo_knight.js | CCE2 with the Knight shift on (ccein_wi_knight.js) and off (ccein_wo_knight.js) in a tilted field, plotted together | 1 conf, 20 ranks | 9 s |
| `single-sample_CCE/` | ccein.js | CCE2 with `nstate=20` sampled bath states per configuration, then the average over configurations and states | 20 conf × 20 states, 5 ranks | 58 s |

Each calculation folder holds its `ccein*.js` inputs, `run.sh` and the idm `job`
it was recorded with. `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
