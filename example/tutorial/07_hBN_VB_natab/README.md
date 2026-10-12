# 07_hBN_VB_natab

V_B⁻ in hexagonal BN. `bath/flat/` holds the bath, the defect, the gyromagnetic ratios and the DFT hyperfine and quadrupole tensors (made with `tools/SpinTensorReader`); the `do_*` scripts beside them are how the tensors were prepared.

| Calculation | Inputs | What it shows | Recorded run | 2026-09-23 on WSL, 2 ranks, 1 configuration |
|---|---|---|---|---|
| `ensemble_CCE/` | ccein_hf0_qd2.js, ccein_hf3_qd0.js, ccein_hf3_qd2.js | CCE2 with the DFT hyperfine and quadrupole tensors: hyperfine order 3 with quadrupole order 0 and 2, and order 0 with 2 | 1 conf, 10 ranks | 17 s |
| `hyperfine_mediated_term/` | ccein_wi_hfmedi.js, ccein_wo_hfmedi.js | CCE2 with and without the hyperfine-mediated term | 1 conf, 20 ranks | 4 s |
| `single-sample_gCCE/` | ccein.js | gCCE2 density matrix of one sampled bath state | 1 conf × 1 states, 10 ranks | 0 s |
| `single-sample_gCCE_DPS/` | ccein_wi_OH.js, ccein_wo_OH.js | gCCE2 with and without the `OH` term (DPS) | 1 conf × 1 states, 10 ranks | 14 s; the zero-field-splitting tensor in both inputs had its second row as [-1180,0,0] (not symmetric, refused as non-Hermitian); fixed to [0,-1180,0] 2026-09-23 |

Each calculation folder holds its `ccein*.js` inputs, `run.sh` and the idm `job`
it was recorded with. `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
