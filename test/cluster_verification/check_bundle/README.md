# check_bundle/ — the bundle clustering against complete calculations

**These two scripts are historical (2026-09-12..16) and do not run with the current build.** They
set the option names of that time (`adaptive_cluster_max_dim`, `adaptive_group_max_dim`,
`adaptive_connections`, ...; since D24 an unread key near a known one stops the run), run the
binaries `bin/main.cpu.out` / `bin/main.cuda.out` the build no longer writes, and
`check_auto_groups.py` calls the removed oracle `check_adaptive/oracle/adaptive_refine.py`. They are kept
because they are the scripts behind the records listed below and cited in
(internal notes). The current checks of the bundle and adaptive code are
`../check_adaptive/verify.sh` (`unit`, `regress`, `inprocess`) and the unit tests of `../check_pcce/`.

| script (historical) | checked then | needed then |
|---|---|---|
| `check_auto_groups.py` | the automatic connection choice (`bundle.connections=auto`, omitted `rdip`) against the full graph and independent subsets | a baseline run; MPI and CUDA variants |
| `check_baths.py --fixtures DIR --output DIR` | bundles on the full NV and PL6 baths across MPI ranks and backends | `../../code_verification/bath_data/ensemble_data` |

Removed on 2026-10-17 (A96; last in commit 8d656004d: `git show 8d656004d:test/cluster_verification/check_bundle/<name>`):
`check_adaptive.py` (adaptive CCE/gCCE against a baseline of full-bath and subset runs),
`check_auto_connections.py` (`bundle.connections=auto` against the full graph) and `check_refinement.py`
(end-to-end refinement against complete small-bath calculations). Nothing ran them; they set keys refused since A47
(`npulse`, ...), and `../check_adaptive/verify.sh` checks the same code now.

Inputs: `simple_ensemble.json` (distance connections), `checked_ensemble.json`
(`rdip` omitted), `ensemble.json`, `auto_ensemble.json`, `gcce_single.json`, and
the six-spin `bath_6` / `bath_hub_6`. `checked_ensemble.json` with `bath_6` is also the input of the
`checked_merged` stage of `../../code_verification/check_refusals/run.sh` (A96).

Kept from the 2026-09 experiments: the records `baths-20260912.json`, `integration-20260912.json` and `validation-20260912.json`.
