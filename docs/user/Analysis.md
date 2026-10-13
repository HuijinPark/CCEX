# Analysis

What runs on the output of CCEX: two Python tools in `example/analysis/`. Neither
runs CCEX. Both need numpy; the plots and fits also scipy and matplotlib (on WSL:
`PATH=~/anaconda3/bin:$PATH`, and `MPLBACKEND=Agg` without a display).

| Tool | Answers | Entry point |
|---|---|---|
| `Analyzer/` | make the raw `_noDiv` / `_wiDiv` files usable: correct the division errors, average over bath configurations and states, fit T2, plot | `Analyzer/CoherenceAnalyzer.v1/coherence.py` (library) and `run.py` (command line) |
| `ConvTest/` | from which value of a variable (order, `rbath`, number of configurations, ...) the coherence stops changing | `ConvTest/ConvTest/do_convTest.py` |

## The files

One line per time step: the time in ms, a tab, the complex coherence as
`+re+imj` (two real columns are read too). CCEX writes `_noDiv`, the product of
the cluster coherences before the division by the sub-clusters, and `_wiDiv`,
after it; the analyzer writes `_ErrCorr` and `<name>_EnsAvg_conv<N>`. A density
matrix comes as blocks of rows, the first row of a block starting with the time.

## The analyzer

Three steps, and every tutorial's `run.sh` does the first two:

1. **Error correction** (`-err`): CCE divides a cluster's coherence by its
   sub-clusters' and the division blows up near zeros, so the divided curve
   (`_wiDiv`) is repaired against the undivided one (`_noDiv`) with a tolerance
   (`-err_tol`, default 2.0). Only a point where |wiDiv − noDiv| reaches the
   tolerance is changed. The rule is the group's convention: a bad last point is
   first moved to noDiv + tolerance − 1e-6 (so it is no longer bad and stays
   there), a bad first point is set to 0, and any other bad point becomes the
   mean of the previous corrected value and the next good one. Every recorded
   result was made with it; call it, do not rewrite it.
2. **Ensemble average** (`-ea`) over the values of a variable (`-vn`, `-v`): one
   `EnsAvg_conv<N>` file per remaining variable.
3. **Plots and fits**: `-pl` one figure per file, `-pl_m` all curves in one, `-fit`
   a stretched-exponential fit with T2 and the power (`-fit_initp`, `-fit_space`),
   `-pl_i` an intensity plot of coherence against time and a variable.

The file names carry the variables: `-fi` names the pattern with the variable
names in it and `-v` the values, so `-fi CCE2_Diamond_NV_natab_CONF_wiDiv -vn CONF -v 1 2 3`
is three files. `-uc` sets the time unit of the plots (`ms`, `us`, `ns`, `s`),
`-ntn` turns a `nan` into 0, `-fo` names the output. A complete call, from
`example/tutorial/02_Diamond_NV_natab/ensemble_CCE/run.sh`:

```
ANALYZER=example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_CONF_wiDiv -err_d results/ -err_fi CCE2_Diamond_NV_natab_CONF_noDiv -vn "CONF" -v 1 2 3 4 -uc "ms" -err -ntn -fo results/CCE2_Diamond_NV_natab_CONF_ErrCorr
# ensemble average, then a plot with a fit
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_CONF_ErrCorr -vn "CONF" -v 1 2 3 4 -uc "ms" -ea -ntn -fo results/CCE2_Diamond_NV_natab_EnsAvg -pl -fit -fit_initp 2 -fit_space 5
```

`python3 run.py -h` lists every option. For a new figure, reading the curves
with `coherence.py` and drawing them yourself is simpler than the plot classes:

```python
import sys; sys.path.insert(0, 'example/analysis/Analyzer/CoherenceAnalyzer.v1')
from coherence import read_coherence, error_correct, ensemble_average, write_coherence
t, L = read_coherence('results/CCE2_Diamond_NV_natab_EnsAvg_conv1', time_scale=1e3)  # t in us
```

`docs/images/tutorials/make_figures.py` is an example of that.

## The convergence test

`do_convTest.py` takes the curve at the last value of the variable as the
converged reference, computes for every other value the area difference in
percent of the reference area, and calls the first value under the tolerance the
converged point; then it draws the curves and the differences. It takes six
positional arguments (file list, variable values, variable name, unit,
tolerance, output figure) and has no `--help`; the argument table is at the top
of the file. No tutorial calls it now (their `do_convtest` scripts were removed on
2026-09-23). Written from scratch it is twenty lines with `read_coherence`.
