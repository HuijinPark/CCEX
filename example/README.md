# example/

Three kinds of things live here: `tutorial/` is what you run, `analysis/` is what
you run on the output, `tools/` is what makes the input.

## tutorial/ — worked examples

Twelve numbered folders, one physical system (or one feature) each; `tutorial/README.md`
has the table. Inside, a calculation is a folder with `ccein*.json` (the input: JSON
with `#` comments), `run.sh` (runs every input with `mpirun` and then the analyzer;
`NP`, `NCONF`, `NSTATE` set the size) and, where it was recorded on idm, a `job` that submits `run.sh` there; the
bath files are in the system's `bath/`. A rerun writes `results/`, which git
ignores; the recorded results are in the archive.

## analysis/ — what runs on the output

| Folder | Does |
|---|---|
| `Analyzer/` | `CoherenceAnalyzer.v1/run.py` reads the `_noDiv` and `_wiDiv` coherences CCEX writes, corrects the division errors (`_ErrCorr`), averages over configurations (`EnsAvg`), fits and plots. The `do_*.sh` beside it are one-line wrappers; `data/` is the example set its `README.md` describes. `test/code_verification/check_input_forms` (`errcorr.sh`) and the tutorials' `run.sh` call `run.py` |
| `ConvTest/` | `do_convTest.py` and `sort.py`: convergence over order and bath size; the tutorials' `do_convtest` scripts that called it were dropped on 2026-09-23 (git history has them) |

## tools/ — what makes the input

| Folder | Does |
|---|---|
| `BathGenerator/` | `BathGenerator.py` with `do_BathGenerator`: a bath file from a POSCAR unit cell (`example/` holds h-BN cells and a hexagonal shape) |
| `IsotopeTable/` | `make_isotope_table.py`: writes `src/input/isotope_table.cpp` (spins and gyromagnetic ratios of the nuclear isotopes CCEX knows without a gyrofile) from the EasySpin database in `SpinTensorReader/` |
| `SpinTensorReader/` | `main.py` with the `spHamilReader/` package: hyperfine and quadrupole tensors from VASP or QE output. `examples/` has five cases, `runningScripts_Exspin/` example runs with their OUTCARs |
| `Runner/` | `MkCCEJobs/mkjobs.sngapp` templates driven by `domkjobs` and `subdomkjobs`: one job folder per variable value. The electron-bath tutorials (04, 05, the P1 runs of 06) were generated with it; their `opt/` and `domkjobs` were dropped on 2026-09-23 and only the inputs and jobs remain. `shortcut/` holds two idm helpers (`QosChanger`, `CdMkDir`) |

## Running a tutorial

`run.sh` finds the binary and the analyzer through `git rev-parse --show-toplevel`,
and an input names its bath files relative to its own directory, so a tutorial runs
from any checkout:

    bash do_compile.sh wsl-gnu
    cd example/tutorial/02_Diamond_NV_natab/ensemble_CCE && NP=4 bash run.sh
