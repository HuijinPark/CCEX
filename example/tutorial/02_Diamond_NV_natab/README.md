# 02_Diamond_NV_natab

NV⁻ in diamond, natural-abundance ¹³C. the basic case; `bath/` holds one configuration of a 150 Å cell (the recorded ensembles used 20, made with `bath/genbath/`), the NV defect file and the ¹³C gyromagnetic ratio.

| Calculation | Inputs | What it shows | Recorded run | 2026-09-23 on WSL, 2 ranks, 1 configuration |
|---|---|---|---|---|
| `ensemble_CCE/` | ccein.json | CCE2 Hahn echo, one run per bath configuration, then the ensemble average | 20 conf, 5 ranks | 5 s |
| `ensemble_Dsj/` | ccein.json, ccein_cce.json | the same with `method=dsj` (ccein.json) beside plain CCE2 (ccein_cce.json) | 1 conf, 5 ranks | 8 s |
| `ensemble_gCCE/` | ccein_multi2.json, ccein_single.json | gCCE2 on one configuration, single qubit (ccein_single.json) and two qubits (ccein_multi2.json), plotted together | 1 conf, 20 ranks | 30 s; ccein_multi2.json had a missing comma (fixed 2026-09-23) |
| `knight_shift/` | ccein_wi_knight.json, ccein_wo_knight.json | CCE2 with the Knight shift on (ccein_wi_knight.json) and off (ccein_wo_knight.json) in a tilted field, plotted together | 1 conf, 20 ranks | 9 s |
| `single-sample_CCE/` | ccein.json | CCE2 with `nstate=20` sampled bath states per configuration, then the average over configurations and states | 20 conf × 20 states, 5 ranks | 58 s |
| `builtin_isotopes/` | ccein_gyrofile.json, ccein_builtin.json | `ensemble_CCE`'s input with the gyro file and without `gyrofile` (¹³C from the built-in isotope table, 6.72828544 against the file's 6.728284 rad/ms/G): the log line saying where ¹³C came from, the largest difference between the curves (1.1e-4, from the 2e-7 relative difference of γ; with the table's value written in a gyro file the outputs are byte-identical), and a plot | added 2026-10-12 (no idm record) | 4 s on 4 ranks |

Each calculation folder holds its `ccein*.json` inputs, `run.sh` and the idm `job`
that submits it. `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
