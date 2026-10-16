# 10_coordinate_rotation

The coordinate-frame rotation (`frame.rotation`). CCEX applies `bfield` in the frame of the bath file, so a bath
written in a frame where the defect axis is not z has to be rotated: by hand, or by this block, which turns every
position about the reference qubit so that `qubit_axis` lands on z (`bath_axis` fixes the azimuth; both default to
+z; the block is on when it is written). `test/code_verification/check_rotation/README.md` describes what moves and
what does not, including the DFT tensor files and the `coordinate_frame` of a `Defect` entry. Every key not written in
an input has its default (`main.out -e`). Times: i7-14700K, WSL, 4 ranks.

`bath/bath_Diamond_NV_natab_cubic_1` is 00_quickstart's bath rotated so that the NV axis points along [111], as a bath
written in the conventional cubic cell would; `bath/make_cubic_frame.py` made it from `00_quickstart/bath/` and says how.

## cCCE2_Hahn_rotation

- **What it shows:** the Hahn echo of 00_quickstart's NV three ways: the bath as generated (NV axis along z,
  `ccein_aligned.json`), the cubic-frame bath with the rotation on (`qubit_axis` [1, 1, 1], `ccein_cubic_rotated.json`),
  and the cubic-frame bath without it (`ccein_cubic_norotation.json`: the field is 54.7° off the NV axis).
- **Keys that matter:** `frame.rotation.qubit_axis` [1, 1, 1]; `configuration.bathfile` (the cubic-frame file).
- **How to run:** `cd cCCE2_Hahn_rotation && NP=4 bash run.sh`: 1.8 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; the rotated run reproduces the aligned one (largest difference 1e-10, T2 1.068
  ms); without the rotation T2 is 1.117 ms and the curve differs by up to 0.91.
- **What to change next:** `bfield` as a vector to tilt the field on purpose; a bath with DFT tensors
  (`frame.rotation.hf_tensor_frame`, `qd_tensor_frame`).
