# check_refusals/ — the inputs that must be refused are refused

    bash test/code_verification/check_refusals/run.sh     # needs bin/main.out

Sixteen input combinations that make no physical sense together, each run once
to see that CCEX stops with the expected message before it reads any tensor
file. They are edits of the legacy inputs `../check_input_forms/legacy.sh`
writes (so it is run first), on the real systems: the h-BN bath with DFT
tensor files (S4), two NV centres with explicit `intmap` tensors (S5), the P1
defect bath (S3).

| refused | because |
|---|---|
| `hf_readmode 3` / `qd_readmode` with the rotation on and no `hf_tensor_frame` / `qd_tensor_frame` | a tensor file has a frame, and the run must not guess it |
| a tilted `bfield` with the rotation on | the field is applied along z of the computational frame |
| `qubit_position_frame` other than `bath` | only one meaning is defined |
| two qubits without `reference_qubit`, or with a name that is not a qubit | the rotation is about one qubit |
| an `intmap` tensor without its frame | same as the tensor files |
| `magnetic_field_axis`, `defect_axis_reference` without the rotation | keys that belong to the rotation block alone |
| `zfs` given as D and E without `axis`; `axis` without `coordinate_frame`; object-form `zfs` without `coordinate_frame`; E non-zero with one axis | a zero-field splitting needs its axes and their frame |
| a Jahn-Teller axis index out of `1..navaax` | it indexes the per-axis tensors |
| an unknown `unit` | the conversion would be silent |
| `coordinate_frame_rotation` and `bath_coordinate_rotation` both present | the same option under two names |

Since 2026-10-01 it also checks failures that used to end with exit status 0 (or, for the second file in info mode, crash):

| refused | because |
|---|---|
| an h-BN tensor file without its `g-factor___` lines, or without its `etc____` lines | the reader said the condition could not be read and then exited 0 (`READ_Tensor_const`, `READ_Tensor_etc`). These come after a tensor reader has started, so the "before any file is read" condition is not applied to them |
| a result file that cannot be written in full, in `normal` and in `info` mode | `Output_save` and `Output_save_info` checked only that the files opened. With its result files linked to `/dev/full` (every write fails with ENOSPC, as on a full disk), `check_numerics/cases/natab_cce2_hahn` ended with exit 0. Until 2026-10-02 the case used a 4 KB file-size limit (`ulimit -f 4`), which also broke Open MPI's own start-up |
| a `_noDiv` file that cannot be opened, in both modes | a directory in its place; `Output_save_info` checked only the `_wiDiv` file and would have written to a null file |

Since 2026-10-11 it also checks the requests of the `factors` and `except` work modes of method D's regroup
((internal notes) 4절, `-stage3b.md` 1절), written by hand on the P1 input of the regroup
regression (`test/cluster_verification/check_adaptive/regression/inputs/regroup-p1-0.json`). Each has to be refused
before any state is calculated (no response file of a state is written):

| refused | because |
|---|---|
| a factors state that is not a state of the run (0 with `nstate` 5, or 6), a step outside `0..nstep-1`, a step that is not an integer | the factors are one state's point of the run |
| factors or except with `method` gCCE, with `order` 1, with `probes`, with `cluster_contributions`, with `-P` | the modes are defined for cCCE bundle order 2 over an explicit partition, with nothing else in the product |
| `state`/`step` outside factors, `except` outside except | a key a mode does not read is not ignored |
| an empty except list, a term listed twice, a site set that is not a term of the partition's expansion | the driver works X out of the same exponent table; one that differs is a different table |

Since 2026-10-11 it also checks `evolution_update` (`docs/user/InputFile_Backend.md`), on inputs of
`../check_numerics/cases/`, refused before the calculation starts:

| refused | because |
|---|---|
| a value other than `exact` or `recursive` (a string, or a number) | only the two are defined |
| `recursive` given with gCCE `propagator` = `expm` (`p1_manual_hahn`, `evolution` = `matrix`), with `backend` = `cuda` (`p1_cce2_single`) or `auto` (`p1_cce2_ensemble`) | the recurrence exists on the CPU only (on the GPU the run would compute as `exact`), and `expm` forms no phases to carry over. Without the key these paths run `exact` (A29, 2026-10-12). The ensemble cCCE and gCCE `matrix` were refused until A18 (2026-10-11) and run now (`../check_input_forms/`) |

Since 2026-10-12 it also checks the `gpu_options` block (`docs/user/InputFile_Backend.md#tag_gpu_options`, A26), on
`../check_numerics/cases/`: each of the five flat keys it replaced (`gpu_min_dim` ...) is refused with the new form
named, and so are a block that is not an object, an unknown key in it, `packed_max_dim` 3, 8 or a string, `min_dim` 0,
`packed_max_dim` in `bin/main.out` (it needs the GPU binary) and with gCCE (it is cCCE only). The accepted forms are
`../check_input_forms/gpu_options.sh`.

`../check_rotation/run.sh` refuses the same combinations on toy systems; this
folder shows the guards fire on real inputs and before any file is interpreted.
What it writes goes to `../check_input_forms/_work/`, which is ignored.
