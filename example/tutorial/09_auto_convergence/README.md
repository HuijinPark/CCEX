# 09_auto_convergence

The automatic convergence (the `auto_convergence` block). With the block, `mpirun -n N main.out -f ccein.json` does not
run one calculation: rank 0 starts the runner (`example/tools/AutoConvergence/`, Python with numpy, scipy and
matplotlib), which runs CCEX itself, each run in its own folder with `mpirun -n N`, chooses the time window, rdip,
rbath and the number of bath configurations, and then computes the production. The rules are in
`example/tools/AutoConvergence/README.md` and `docs/user/InputFile_AutoConvergence.md`. Every key not written in the
input has its default (`main.out -e`).

## cCCE2_Hahn_ensemble

- **What it shows:** 00_quickstart's NV, cCCE2 Hahn echo of the ensemble, with the rdip candidates 4-12 Å and the
  rbath candidates 30-45 Å given and everything else found: the window (1.5 × the decay time), rdip, rbath and how
  many configurations the mean needs (SE(A)/A ≤ 5 %).
- **Keys that matter:** `auto_convergence` {`workdir`, `rdip.candidates`, `rbath.candidates`, `sampling`},
  `configuration.bathfile` ["results/bath/config_{k}"] (every existing k is used), `nstep` 101 (the points of every
  run), `deltat` (only the first window).
- **How to run:** `NP=4 bash run.sh` makes 40 configurations first (`../../make_natab_baths.sh`, 110 Å cells: up to
  rbath 53 Å; `NBATH` sets how many) and runs the convergence: 4 min in all (i7-14700K, WSL, 4 ranks), under 0.2 GB per
  process. Calling it again continues (the baths are kept, finished runs and decisions reused). The automatic
  convergence exits 0 when `converged` and 3 when `unresolved` (the summary says what is missing); `run.sh` prints
  that status and `results/auto/result/summary.txt`.
- **Expected result:** `expected.png` (the runner's `6_converged_coherence.png`); status converged: window 2.58 ms,
  rdip 8 Å, rbath 35 Å, 25 production configurations (σA/A 4.0 %), T2 0.82 ms (stretched-exponential fit, n 1.37).
  `results/auto/result/` also has `1_window.png`, `2_rdip.png`, `3_rbath.png`, `5_configurations.png` and
  `coherence_mean.txt`.
- **What to change next:** leave out `rdip` and `rbath` (the runner proposes candidates from the spin density; more
  runs); raise `nstep`: the 101 points over 2.58 ms are 26 µs apart, coarser than the ¹³C precession at 500 G (about
  2 µs), so the ESEEM ripple is aliased (the runner keeps `nstep` and does not check it).

What the result does not say: the combination check (rbath 35, rdip 8 against the largest computed, rbath 45, rdip 15,
on three production configurations) differs by 9 % in the area, more than the 2 % threshold of each step (a warning in
the summary).
