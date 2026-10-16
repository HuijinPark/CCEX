# Automatic convergence tags

An `auto_convergence` block in an ordinary input makes `main.out -f input.json` choose the time
window, rdip, rbath and the number of configurations (and states) of a coherence calculation by
itself, then run the production. Every calculation it makes is an ordinary CCEX run in its own
folder under `workdir`. Supported: cCCE with `cluster.mode` standard or pcce (order ≥ 2), gCCE (order ≥ 1), fixed bundles
(`cluster.mode` bundle without a `cluster.adaptive` block), CPU or GPU backends. How each step decides,
the outputs and the exit codes are in `example/tools/AutoConvergence/README.md`; `example/tutorial/09_auto_convergence/`
is a worked input. The runner reads the keys (of the block and of the rest of the
input) without regard to case, as CCEX does, and refuses two keys of one block that differ only in case.

```
mpirun -n 8 main.out -f input.json   # every CCEX run with mpirun -n 8 (rank 0 alone runs the runner)
main.out -f input.json               # every CCEX run on one rank
```

(The block has no `launcher` key, and one is refused; the runs use the `mpirun -n N` main.out was started
with.)

Where a run is (nothing computed) is the tool's own command, run from `example/tools/AutoConvergence` or with that
folder on `PYTHONPATH` (`main.out` itself has no status option):

```
PYTHONPATH=example/tools/AutoConvergence python3 -m autoconv status input.json
```

The automatic run takes only `-f FILE`; other command-line options are refused (set them in the
input). `configuration.nstate` of the input chooses the approach: `0` is the ensemble (only the number of
configurations is chosen), `> 0` samples states (`sampling.test_states` and `max_states`
set how many). `rbath`, `rdip`, `nstep` and `deltat` of the input are not needed; `nstep` is the
number of time points of every run (150 when absent).

The block is grouped by step. A section left out uses its defaults; a key given wins over what
would be found; `null` or an empty list is refused (leave the key out to have it found).
Unknown keys and the runner's earlier flat names are refused with a message that says where each
went.

<a id="tag_ac_configurations"></a>
## The bath configurations

The configurations are the CCEX [configuration](InputFile_General.md#tag_configuration) block's:
`bathfile` (one name with `{k}`), and when given `bath_state_file`, `defect_axis_file` and `defect_state_file` (one name
each, `{k}` too), all in `configuration.path` when it is given. With `configuration.k_range` [first, last] the
configurations are those k (every one must have its files); without it, every k whose bath file exists (`config_7`, not
`config_07`). `bath_state_file` left out with nstate > 0: the states are made per run and removed; `defect_axis_file`
left out with a Defect species that has axis states: the axes are made per run. Each run's input holds one
configuration: the names with k in place of `{k}`, absolute, without `k_range` and `path`.

The block has no section of its own for these file names: a former `auto_convergence.configurations` section is
refused ([the old names](InputFile.md#OldNames) say where each key went), and
`example/tools/InputKeys/convert_input_keys.py` rewrites such an input.

<a id="tag_ac_fixed"></a>
## Fixed values instead of convergence

Every section converges by default. `"converge": false` uses the value given instead, and that
value is then required:

| Section with `"converge": false` | value used | effect |
|---|---|---|
| `time_window` | `value` (+ `unit`) | the window is used as it is; T_d is only measured |
| `rdip` | `value` (Å) | the rdip step and its check at the final rbath are skipped |
| `rbath` | `value` (Å) | the rbath step is skipped; a plateau stops instead of raising rbath |
| `sampling` | `number` (and `states` when nstate > 0) | the sampling step is skipped; the production uses exactly that many |

All four `false` is refused (run CCEX itself).

<a id="tag_ac_keys"></a>
## Keys

| Section | Key | Meaning | Default |
|---|---|---|---|
| (top) | `enabled` | false: `main.out -f` runs the ordinary calculation and the block is not read | true |
| | `resume` | `continue`: reuse what the workdir holds; `fresh`: move it to `<workdir>.old<n>` and start over | `continue` |
| | `workdir` | folder for the runs (`runs/`), the decisions (`decisions.jsonl`) and the results (`result/`) | required |
| `time_window` | `converge`, `value`, `unit` | see above; units `s`, `ms`, `us`, `ns`; `value` may also be `{"value": ..., "unit": ...}` without `unit` beside it ([Units](InputFile.md#Units)) | true, —, `ms` |
| | `threshold` | T_d is where the T2 fit reaches this level | 0.01 |
| | `revival_check_factor` | the test runs go to this × T_d to see a revival | 2.5 |
| | `production_factor` | production window = this × T_d | 1.5 |
| `rdip` | `converge`, `value` | see above | true |
| | `candidates` | rdip candidates (Å); the largest is the reference | found |
| | `num_candidates` | how many to start with when found | 7 |
| | `threshold` | relative difference of the fitted area A = T2 Γ(1+1/n) at which two rdip agree | 0.02 |
| `rbath` | `converge`, `value` | see above | true |
| | `candidates` | rbath candidates (Å) | found |
| | `num_candidates` | how many to start with when found | 4 |
| | `threshold` | relative difference of A at which two rbath agree | 0.02 |
| `t2_fit` | `points` | points of every T2 fit: `peaks` (t = 0 and the local maxima; all points when fewer than 3) or `all` | `peaks` |
| | `floor` | only the points before \|mean\| falls below this are fitted | 0.05 |
| | `n_range` | allowed stretching exponent n of exp(−(t/T2)^n); outside it the fit fails | [0.3, 4.0] |
| `sampling` | `test`, `production` | configuration IDs for the convergence tests and for the production (never the same; within `configuration.k_range` when it is given) | found |
| | `converge`, `number`, `states` | see above | true, —, — |
| | `threshold` | SE(A)/A target of the sampling | 0.05 |
| | `test_number`, `test_states` | configurations and states of the test runs (states for nstate > 0) | 10, 5 |
| | `states_threshold` | the states' share of the sampling error at which the state test stops | 0.05 |
| | `max_states` | the most states the state test computes (with state files, at most the files there) | 30 |
| | `min_configurations` | the fewest production configurations | 10 |
| | `max_configurations` | the most production configurations used | every production file |

<a id="tag_ac_refused"></a>
## Refused together with the block

Other cluster modes (`dsj`, `itb`, `dsjitb`, the adaptive refinement) or quantities than coherence, `rdip_ladder` in the
input, a `cluster.adaptive` block, `cluster.pcce.cache` (the runner gives its own), `cluster.pcce.partition` legacy,
`output.cluster_contributions`, `output.cluster_only`, a nonzero `cluster.selection.max_per_order`, a nonzero
`configuration.bath_shift`, more than one name in `configuration.bathfile` (one bath file per configuration) or a name
without `{k}`, keys in their former places (CCEX refuses them too, see [the old names](InputFile.md#OldNames)), the block's `launcher` and
`configurations`, more than one qubit, a list together with its count, a fixed value with `converge` true, and options
out of range. The input's `output.path` is not used (every run writes into its own folder under `workdir`).

<a id="tag_ac_rdip_ladder"></a>
## rdip_ladder (an internal key)

`rdip_ladder` is a CCEX input key the runner writes into its runs' inputs: a list of 1 to 64 increasing rdip values
computed from one order-2 run (cCCE, gCCE, cluster.mode pcce), each written as `<output.file>_rdip<r>_...` beside the
run's own rdip. It serves the rdip step, which gets every smaller candidate from the run at the largest one; it is not
part of the documented input (left out of the key lists and of `main.out -e`), and an input given to the automatic
convergence may not hold it.

<a id="tag_ac_example"></a>
## Example

```
"configuration": {"gyrofile": "data/gyro_13C", "nstate": 0, "bathfile": ["data/config_{k}"], "k_range": [1, 60]},
"auto_convergence": {
  "workdir": "auto",
  "time_window": {"converge": false, "value": 15, "unit": "us"},
  "rbath": {"converge": false, "value": 288},
  "sampling": {"threshold": 0.03}
}
```
