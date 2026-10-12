# example/analysis/

Two Python tools that post-process the coherence files CCEX writes. Neither runs
CCEX. This page is written so that a person or an AI can decide quickly what to
call and what to write themselves.

| Tool | Answers | Entry point | Needs |
|---|---|---|---|
| `Analyzer/` | make the raw `_noDiv`/`_wiDiv` files usable: correct the division errors, average over bath configurations and states, fit T2, plot | `Analyzer/CoherenceAnalyzer.v1/coherence.py` (library) and `run.py` (CLI) | numpy; the CLI and plots also scipy and matplotlib |
| `ConvTest/` | from which value of a variable (order, rbath, number of configurations, ...) the coherence stops changing | `ConvTest/ConvTest/do_convTest.py` (CLI, positional arguments) | numpy, matplotlib, psutil |

On WSL both need `PATH=~/anaconda3/bin:$PATH` (the system python3 has no numpy)
and, without a display, `MPLBACKEND=Agg`.

## Use as-is, do not rewrite

**The error-correction rule** (`coherence.error_correct`). CCE divides a cluster's
coherence by its sub-clusters' and the division blows up near zeros; the
analyzer repairs the divided curve against the undivided one with a tolerance,
by a rule that is the group's convention (`Analyzer/README.md` states it). Every
regression reference and every published curve was made with it, so a
re-implementation that differs even slightly changes what "the same result"
means. Call the function, or the CLI with `-err`, exactly as
`test/code_verification/check_numerics/*/*/job` does.

**Reading and writing the files** (`coherence.read_coherence`,
`coherence.write_coherence`). The format is simple (time in ms, tab, `+re+imj`),
but the reader also handles the two-column real/imag form and a `nan` in the
coherence column, and the writer's 10-decimal format is what every recorded
result was written with.

**The ensemble average** (`coherence.ensemble_average`), so that the `EnsAvg`
files it writes stay comparable with the recorded ones.

## Fine to do yourself

**Plots and fits.** `DataSpace` and `CoherencePlot_v2` draw and fit (stretched
exponential, T2 and power) the way the 2023 figures were made. For a new figure,
reading the curves with `coherence.py` and drawing them directly is simpler than
learning the class.

**The convergence test.** `do_convTest.py` takes the curve at the last value of
the variable as the converged reference, computes for each other value the
area difference in percent of the reference area, and calls the first value
under the tolerance the converged point; it then draws the curves and the
differences. That is twenty lines with `coherence.read_coherence`. Use the
script when a tutorial's `do_convtest` already drives it (they pass six
positional arguments: file list, variable values, variable name, unit,
tolerance, output figure); otherwise re-implement. Its option parser is a custom
one (`pyform.py`), so there is no `--help`; the argument table is at the top of
`do_convTest.py`.

## How the tutorials use them

A `job` in `example/tutorial/` runs CCEX, then `run.py -err` on each
configuration, then `run.py -ea` over all of them; a `do_convtest` next to it
runs `do_convTest.py` over the `ErrCorr` files of several orders or radii. The
regression suite under `test/code_verification/check_numerics/` does the first two
steps and is the reference for the exact calls.
