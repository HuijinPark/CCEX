# 07_hBN_VB_natab

V_B⁻ in hexagonal BN. `bath/flat/` holds the bath, the defect and the DFT hyperfine and quadrupole tensors (made with `tools/SpinTensorReader`); the inputs read the hyperfine file `hBN_VB_flat_ATensor_HSE_CellPBE_v3` (tensor format v3, 2026-10-13: the hyperfine of the S = 1 defect, nearest ¹⁴N A_zz 48.6 MHz; until then `..._tag`, read 2x); the gyromagnetic ratios are one level up, `bath/hBN_natab_gyro`, beside the `do_*` scripts that show how the tensors were prepared.

| Calculation | Inputs | What it shows | Recorded run | 2026-10-13 on WSL, 4 ranks, 1 configuration (run.sh with the analyzer) |
|---|---|---|---|---|
| `ensemble_CCE/` | ccein_hf0_qd2.json, ccein_hf3_qd0.json, ccein_hf3_qd2.json | CCE2 with the DFT hyperfine and quadrupole tensors: hyperfine order 3 with quadrupole order 0 and 2, and order 0 with 2 | 1 conf, 10 ranks | 5.1 s |
| `hyperfine_mediated_term/` | ccein_wi_hfmedi.json, ccein_wo_hfmedi.json | CCE2 with and without the hyperfine-mediated term | 1 conf, 20 ranks | 2.2 s |
| `single-sample_gCCE/` | ccein.json | gCCE2 density matrix of one sampled bath state | 1 conf × 1 states, 10 ranks | 0.4 s |
| `single-sample_gCCE_DPS/` | ccein_wi_OH.json, ccein_wo_OH.json | gCCE2 with and without the `OH` term (DPS) | 1 conf × 1 states, 10 ranks | 2.8 s; the zero-field-splitting tensor in both inputs had its second row as [-1180,0,0] (not symmetric, refused as non-Hermitian); fixed to [0,-1180,0] 2026-09-23 |

Each calculation folder holds its `ccein*.json` inputs, `run.sh` and the idm `job`
that submits it. `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
