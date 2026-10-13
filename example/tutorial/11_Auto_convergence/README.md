# 11_Auto_convergence

The automatic convergence (the `auto_convergence` block in the input). With the block,
`mpirun -n N main.out -f ccein.json` does not run one calculation: rank 0 starts the runner
(`example/tools/AutoConvergence/`, Python with numpy, scipy and matplotlib), which runs CCEX
itself, each run in its own folder with `mpirun -n N`, and chooses the time window, rdip,
rbath and the number of bath configurations, then computes the production. The rules are
in `example/tools/AutoConvergence/README.md` and (internal notes).

| Calculation | Inputs | What it shows | Recorded run |
|---|---|---|---|
| `ensemble_CCE2/` | ccein.json | CCE2 Hahn echo of the NV in natural-abundance diamond (02's system), ensemble (`nstate` 0), rdip candidates 4–12 Å and rbath candidates 30–45 Å given, everything else found: window 2.58 ms (1.5 × T_d 1.72 ms), rdip 8 Å (against the reference 15 Å, the largest candidate raised by 1.25), rbath 35 Å, 25 production configurations (SE(A)/A 4.0 %, target 5 %), status `converged` | added 2026-10-12; WSL, 4 ranks: 3 min 24 s for the 40 baths, 42 s for the convergence |

`run.sh` first makes the bath configurations, since the convergence needs tens of them and
the repository keeps one: 40 natural-abundance ¹³C baths in a 110 Å cell (they cover rbath up
to 53 Å), made with `tools/BathGenerator` from 02's `bath/genbath/` cell files, configuration
k with the random seed k, so a rerun makes the same files (`results/bath/`; `NBATH` sets how
many). Then it runs the convergence and prints `results/auto/result/summary.txt`. Everything
goes to `results/auto/`:

- `result/summary.txt` (and `result/json/summary.json`): what was found and why, step by step
  (0 the configurations and the split into 10 test and 30 production configurations, 1 the
  window, 2 rdip, 3 rbath, 5 the configurations, 6 the converged curve) and the warnings.
- `result/1_window.png`, `1_window_fit.png`, `2_rdip.png`, `3_rbath.png`, `5_configurations.png`,
  `6_converged_coherence.png`; `result/coherence_mean.txt`, the mean production curve.
- `runs/` one folder per CCEX run, `decisions.jsonl` the decisions, `result/progress.txt` while it runs.

Calling `run.sh` again continues: the baths are kept, finished runs and decisions are reused.
The automatic convergence (`main.out`) exits with 0 when `converged` and 3 when `unresolved`
(the summary says what is missing, e.g. "make larger baths"). `run.sh` only prints that status
("automatic convergence ended with status N") and itself ends with the status of its last
command, so read the status line or `summary.txt`, not the exit status of `run.sh`.

What the recorded result does not say: the combination check (rbath 35, rdip 8 against the
largest computed, rbath 45, rdip 15, on three production configurations) differs by 9 % in the
area, more than the 2 % thresholds of each step (a warning in the summary: each choice is
within its threshold of the largest computed, the combination is not). The 101 time points
over 2.58 ms are 26 µs apart, coarser than the ¹³C precession at 500 G (about 2 µs), so the
fast ESEEM ripple is aliased in the curves (the runner keeps `nstep` of the input and does not
check whether it is enough; raise `nstep` for a resolved curve, at proportionally more cost).

    cd ensemble_CCE2 && NP=4 bash run.sh

`example/tutorial/README.md` explains the layout.
