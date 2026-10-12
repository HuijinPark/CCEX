# 10_Coordinate_rotation

The coordinate-frame rotation (`coordinate_frame_rotation` in the input). CCEX applies
`bfield` along z, so a bath written in a frame where the defect axis is not z has to be
rotated first: either by hand, or by this block, which rotates every position about the
reference qubit so that `qubit_axis` lands on z (`bath_axis` says where z of the file
points; both default to +z). `test/code_verification/check_rotation/README.md` describes what
moves and what does not, including the DFT tensor files and the `coordinate_frame` of a
`Defect` entry.

| Calculation | Inputs | What it shows | Recorded run |
|---|---|---|---|
| `HahnEcho_CCE2/` | ccein_aligned.json, ccein_cubic_rotated.json, ccein_cubic_norotation.json | CCE2 Hahn echo of the NV of 02: the bath as generated (NV axis along z), the same bath written in the cubic crystal frame with the rotation on (`qubit_axis` [1,1,1]), and that cubic-frame bath with the rotation off; the differences from the first curve and a plot | added 2026-09-23; WSL, 2 ranks |

`bath/` holds a copy of 02's natural-abundance bath, its NV defect file and the ¹³C
gyromagnetic ratio, and `bath_Diamond_NV_natab_cubic_1`: that bath rotated so that the
NV axis points along [111]; `bath/make_cubic_frame.py` made it and says how. The
rotated bath with the rotation on reproduces the original to the printed precision
(difference 1e-10); with the rotation off the field sits 54.7° from the NV axis and the
echo is a different curve.

    cd HahnEcho_CCE2 && NP=4 bash run.sh

`example/tutorial/README.md` explains the layout.
