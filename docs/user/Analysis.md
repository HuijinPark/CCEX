# Analysis

What runs on the output of CCEX: the analyzer in `example/analysis/` (it does not
run CCEX). It needs numpy; the plots and fits also scipy and matplotlib (on WSL:
`PATH=~/anaconda3/bin:$PATH`, and `MPLBACKEND=Agg` without a display).

| Tool | Answers | Entry point |
|---|---|---|
| `Analyzer/` | make the raw `_noDiv` / `_wiDiv` files usable: correct the division errors, average over bath configurations and states, fit T2, plot | `Analyzer/CoherenceAnalyzer.v1/coherence.py` (library) and `run.py` (command line) |

## The files

One line per time step: the time in ms, a tab, the complex coherence as
`+re+imj` (two real columns are read too). CCEX writes `_noDiv`, the product of
the cluster coherences before the division by the sub-clusters, and `_wiDiv`,
after it; the analyzer writes `_ErrCorr` and `<name>_EnsAvg_conv<N>`. A density
matrix comes as blocks of rows, the first row of a block starting with the time.

## The analyzer

Three steps (the tutorials' `run.sh` plot the raw `_wiDiv` with `example/tutorial/plot_coherence.py`
and do not call it):

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
names in it and `-v` the values, so `-fi ccein_kCONF_wiDiv -vn CONF -v 1 2 3`
is three files. `-uc` sets the time unit of the plots (`ms`, `us`, `ns`, `s`),
`-ntn` turns a `nan` into 0, `-fo` names the output. A complete call, for the four bath configurations
of `example/tutorial/01_ensemble_and_states/cCCE2_Hahn_ensemble` (`configuration.k_range` [1, 4], written as
`results/ccein_k<k>_wiDiv` and `_noDiv`):

```
ANALYZER=example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi ccein_kCONF_wiDiv -err_d results/ -err_fi ccein_kCONF_noDiv -vn "CONF" -v 1 2 3 4 -uc "ms" -err -ntn -fo results/ccein_kCONF_ErrCorr
# ensemble average, then a plot with a fit
python3 $ANALYZER -d results/ -fi ccein_kCONF_ErrCorr -vn "CONF" -v 1 2 3 4 -uc "ms" -ea -ntn -fo results/ccein_EnsAvg -pl -fit -fit_initp 2 -fit_space 5
```

`python3 run.py -h` lists every option. For a new figure, reading the curves
with `coherence.py` and drawing them yourself is simpler than the plot classes:

```python
import sys; sys.path.insert(0, 'example/analysis/Analyzer/CoherenceAnalyzer.v1')
from coherence import read_coherence, error_correct, ensemble_average, write_coherence
t, L = read_coherence('results/ccein_EnsAvg_conv4', time_scale=1e3)  # t in us
```

`example/tutorial/plot_coherence.py` (the tutorials' figures) is an example of that.

## The convergence test

There is no separate convergence-test script. The automatic convergence (`auto_convergence` block,
`example/tutorial/09_auto_convergence/`) chooses the time window, rdip, rbath and the number of configurations; for one
variable by hand, run the values and compare the curves with `read_coherence`.
