# test/code_verification/

Does the code still do what it did, and what it says. Seven checks and the data
they read; `run.sh` runs the seven in order and exits with the first failure.

| folder | question | run |
|---|---|---|
| `check_analytic/` | does CCEX give the closed-form coherence of one 13C beside the NV (FID and Hahn echo, order 1, ensemble)? the one check not comparing CCEX with CCEX | `bash check_analytic/run.sh` |
| `check_numerics/` | are the numbers the ones earlier builds gave? 25 calculations against this branch's recorded run (exact) and the tags v1.1.0 and v1.2.2 (with the known differences named) | `bash check_numerics/run.sh` |
| `check_rotation/` | does the coordinate-frame rotation do what it says? a unit test and 149 run checks against a hand-rotated input | `bash check_rotation/run.sh` |
| `check_refusals/` | are the 93 inputs that must be refused refused, most before anything is read, and do three that must run run? | `bash check_refusals/run.sh` |
| `check_input_forms/` | do the new input forms say the same physics as the legacy forms? 32 comparisons under one binary, the analyzer reading the outputs, then the gyro, comments, nodiv, evolution_update, gpu_options and autoconv stages | `bash check_input_forms/run.sh` |
| `check_mean_field/` | does `vector_mean_field` give the default z mean field's curve byte for byte where the answer is known (one P1 with its 14N and no other bath spin)? | `bash check_mean_field/run.sh` |
| `check_rdip_ladder/` | does each `rdip_ladder` rung equal a separate run at that rdip (≤ 1e-9), and are the inputs it does not cover refused? about 7 min | `bash check_rdip_ladder/run.sh` |
| `bath_data/` | not a check: the frozen baths, tensors and state files the checks read | — |

Every check needs `bin/main.out` (`bash do_compile.sh wsl-gnu`) and `mpirun`;
`check_numerics` and the analyzer stage need a python3 with numpy. What a
check writes goes to its `_work/`, which is ignored and wiped at each run.
Each folder's README says what it checks and how it decides; `check_mean_field/` and
`check_rdip_ladder/` have none and are described in `../README.md`.
