# example/analysis/

Two Python tools that post-process the coherence files CCEX writes (neither
runs CCEX). This page is written so that a person or an AI can decide quickly what to
call and what to write themselves.

| Tool | Answers | Entry point | Needs |
|---|---|---|---|
| `Analyzer/` | make the raw `_noDiv`/`_wiDiv` files usable: correct the division errors, average over bath configurations and states, fit T2, plot | `Analyzer/CoherenceAnalyzer.v1/coherence.py` (library) and `run.py` (CLI) | numpy; the CLI and plots also scipy and matplotlib |
| `ConvTest/` | from which value of a variable (order, rbath, number of configurations, ...) the coherence stops changing | `ConvTest/ConvTest/do_convTest.py` (CLI, positional arguments) | numpy, matplotlib, psutil |

On WSL the first two need `PATH=~/anaconda3/bin:$PATH` (the system python3 has no numpy)
and, without a display, `MPLBACKEND=Agg`.

## Use as-is, do not rewrite

**The error-correction rule** (`coherence.error_correct`). CCE divides a cluster's
coherence by its sub-clusters' and the division blows up near zeros; the
analyzer repairs the divided curve against the undivided one with a tolerance,
by a rule that is the group's convention (`Analyzer/README.md` states it). Every
tutorial result and every published curve was made with it, so a
re-implementation that differs even slightly changes what "the same result"
means. Call the function, or the CLI with `-err`, exactly as
the tutorials' `run.sh` do.

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
differences. That is twenty lines with `coherence.read_coherence`. The
tutorials' `do_convtest` scripts that drove it were removed on 2026-09-23; it takes six
positional arguments (file list, variable values, variable name, unit,
tolerance, output figure), or re-implement it. `python3 ConvTest/ConvTest/do_convTest.py --help`
lists the arguments and the options `--add` and `--findconv` (the parser in `pyform.py` wraps `argparse`).

## How the tutorials use them

A tutorial's `run.sh` runs CCEX, then `run.py -err` on each
configuration, then `run.py -ea` over all of them, and is the reference for the
exact calls. The `do_convtest` scripts that ran `do_convTest.py` were removed on
2026-09-23; `test/code_verification/check_numerics/` compares raw outputs and does
not call the analyzer.
