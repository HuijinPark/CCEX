# 12_rdip_ladder

Two keys that make a CCE2 run cheaper without changing what it computes, on the NV of 02
(the inputs read 02's `bath/`).

- `rdip_ladder` (`docs/user/InputFile_General.md#tag_rdip_ladder`): rungs r<sub>1</sub> < r<sub>2</sub> < ...
  below `rdip`. One order-2 run also writes, for every rung, the curve that a separate run at
  `rdip` = r<sub>k</sub> gives (`<outfile>_rdip<r>_wiDiv`/`_noDiv`); the main output stays the run at
  `rdip`. Order 2 of CCE, gCCE and pCCE, `quantity` coherence, standard clusters, `savemode`
  normal; anything else is refused with `Error: rdip_ladder ...`
  (`test/code_verification/check_rdip_ladder/run.sh` lists the refusals). The automatic
  convergence (11) uses it for its rdip step.
- `evolution_update` (`docs/user/InputFile_Backend.md#tag_evolution_update`): `exact` forms the
  phases exp(−iEτ) with sin and cos at every time point, `recursive` at every 16th point and by
  one step's factor in between. `recursive` is the default on the CPU since 2026-10-12; the
  curves agree to rounding, not bit for bit.

| Calculation | Inputs | What it shows | Recorded run |
|---|---|---|---|
| `ensemble_CCE2/` | ccein_ladder.json, ccein_rdip6.json, ccein_exact.json, ccein_recursive.json | CCE2 Hahn echo (rbath 50, rdip 12, 10000 points of 0.2 µs): one run with `rdip_ladder` [4, 6, 8, 10] gives five curves; its rdip-6 rung equals the separate rdip-6 run (difference 0) and its main output the run without the ladder (0). The same input with `evolution_update` `exact` and `recursive`: largest difference 1.0e-10 in `wiDiv` and `noDiv` (one unit of the 10th decimal written), 5.5 s against 4.6 s. A plot of the rungs and their distance from the main output | added 2026-10-12; WSL, 4 ranks, 1 configuration: ladder 6.8 s, rdip 6 1.0 s, exact 5.5 s, recursive 4.6 s (CCEX wall times) |

The ladder run took 6.8 s for five curves; the run at rdip 12 alone took 4.6 s and the separate
rdip-6 run 1.0 s (the other rungs were not run separately here).
`results/CCE2_ladder_rungs.png` shows the five curves and |L<sub>rung</sub> − L<sub>rdip 12</sub>|; its largest
value over the 2 ms is 7.0e-2 for rdip 4 Å, 2.5e-2 for 6 Å, 1.0e-2 for 8 Å and 3.6e-3 for 10 Å.

    cd ensemble_CCE2 && NP=4 NCONF=1 bash run.sh

`example/tutorial/README.md` explains the layout.
