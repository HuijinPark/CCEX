# 08_Bundle_clustering

The clustering this branch added, on the NV of 02 (`bath/` is a copy of 02's bath, defect and gyro files). `docs/methods/clustering-terminology.md`
fixes the words: a **bundle** is a set of strongly coupled spins that always enters a
cluster whole and is computed exactly inside; **bundle clustering** (`cluster_mode`
`adaptive` in the input, `src/clustering/bundle/`) builds the bundles once, from the
coupling strengths, under the dimension cap `bundle_max_dim`, and then runs CCE with
`order` counting bundles per cluster; **adaptive refinement** (the `Adaptive` block,
`src/clustering/algo_adaptive/`, `docs/methods/adaptive-input-block.md`) repeats that calculation, measures
the error at the bundle boundaries and re-bundles until the curve settles.

| Calculation | Inputs | What it shows | Recorded run |
|---|---|---|---|
| `CCE2_bundle/` | ccein_standard.js, ccein_bundle.js | CCE2 Hahn echo (rbath 30) with ordinary clusters and with bundles of up to three ¹³C (`bundle_max_dim` 8): the `[adaptive]` line that says how the 199 spins were bundled (84 bundles), the largest difference between the two curves, and a plot | added 2026-09-23; WSL, 2 ranks: 1 s and 21 s |
| `CCE2_adaptive/` | ccein.js | the bundled calculation (rbath 20, 2 ms) with an `Adaptive` block, merge policy, `atol` 1e-3: three rounds, 28 → 26 bundles, then "no feasible improving merge" (exit status 2); `results/adaptive/` holds one folder per step and `summary.json` | added 2026-09-23; WSL, 2 ranks: 7 min |

Bundles of one spin each are ordinary CCE, so the two curves of `CCE2_bundle` differ
only by what the bundles change (here 1.7e-2 at most). `Adaptive` needs an explicit
`seed`: every round recomputes the same system.

    cd CCE2_bundle && NP=4 bash run.sh

`example/tutorial/README.md` explains the layout.
