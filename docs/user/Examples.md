# Examples

Eleven worked examples live in `example/tutorial/`, in the order to learn them. Each topic folder holds one folder per
calculation (named `<method><order>_<pulse>_<approach>`) with its inputs (only the keys that matter), a `run.sh` and
the `expected.png` a rerun is compared with; the topic's `README.md` says, per calculation, what it shows, the keys
that matter, how to run it (time, ranks, memory), the expected result and what to change next.

| Topic | README |
|---|---|
| quickstart: NV in natural-abundance ¹³C, cCCE2 Hahn echo, seconds | [00_quickstart](../../example/tutorial/00_quickstart/README.md) |
| bath configurations (`k_range`) and sampled bath states (`nstate`) | [01_ensemble_and_states](../../example/tutorial/01_ensemble_and_states/README.md) |
| Hahn echo, CPMG, a sequence written out | [02_pulses](../../example/tutorial/02_pulses/README.md) |
| an electron-spin bath: P1, Vp | [03_P1](../../example/tutorial/03_P1/README.md) |
| V_B⁻ in h-BN with DFT tensors | [04_hBN](../../example/tutorial/04_hBN/README.md) |
| several NV centres | [05_multiNV](../../example/tutorial/05_multiNV/README.md) |
| cCCE against gCCE, and the cCCE options | [06_cCCE_vs_gCCE](../../example/tutorial/06_cCCE_vs_gCCE/README.md) |
| pCCE, bundles, adaptive refinement | [07_pCCE_and_bundle](../../example/tutorial/07_pCCE_and_bundle/README.md) |
| the GPU backend | [08_GPU](../../example/tutorial/08_GPU/README.md) |
| the automatic convergence | [09_auto_convergence](../../example/tutorial/09_auto_convergence/README.md) |
| the coordinate-frame rotation | [10_coordinate_rotation](../../example/tutorial/10_coordinate_rotation/README.md) |

```
bash do_compile.sh wsl-gnu                        # once; makes bin/main.out
cd example/tutorial/00_quickstart
NP=4 bash run.sh                                  # results/coherence.png; compare with expected.png
```

`example/README.md` has the table "what should I use for my case" (situation → method → example), and
`example/tutorial/README.md` the map from the tutorial folders of earlier versions to these. Each topic's `README.md`
and `expected.png` describe it next to its inputs.
