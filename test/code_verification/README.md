# test/code_verification/

Does the code still do what it did, and what it says. Four checks and the data
they read; `run.sh` runs the four in order and exits with the first failure.

| folder | question | run |
|---|---|---|
| `check_numerics/` | are the numbers the ones earlier builds gave? 25 calculations against this branch's recorded run (exact) and the tags v1.1.0 and v1.2.2 (with the known differences named) | `bash check_numerics/run.sh` |
| `check_rotation/` | does the coordinate-frame rotation do what it says? a unit test and 149 run checks against a hand-rotated input | `bash check_rotation/run.sh` |
| `check_refusals/` | are the 16 inputs that must be refused refused, before anything is read? | `bash check_refusals/run.sh` |
| `check_input_forms/` | do the new input forms say the same physics as the legacy forms? 19 pairs under one binary, then the analyzer reading the outputs | `bash check_input_forms/run.sh` |
| `bath_data/` | not a check: the frozen baths, tensors and state files the checks read | — |

Every check needs `bin/main.out` (`bash do_compile.sh wsl-gnu`) and `mpirun`;
`check_numerics` and the analyzer stage need a python3 with numpy. What a
check writes goes to its `_work/`, which is ignored and wiped at each run.
Each folder's README says what it checks and how it decides.
