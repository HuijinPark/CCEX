# Examples

Ten worked examples live in `example/tutorial/`, numbered in the order to read
them: from the smallest hand-written input to the features this branch added.
Each folder is one physical system (or one feature on a system above), and each
calculation inside it is a folder with its inputs (`ccein*.json`, JSON with `#`
comments), one `run.sh` that runs them and the post-processing, and the cluster
`job` the results were first recorded with. The system folder's `README.md` has
the table of its calculations.

## Running one

```
bash do_compile.sh wsl-gnu                        # once; makes bin/main.out
cd example/tutorial/02_Diamond_NV_natab/ensemble_CCE
NP=4 NCONF=1 bash run.sh                      # writes results/
```

`run.sh` runs every input of the folder with `mpirun -n $NP`, then the error
correction, the ensemble average and the plot the recorded run did
(`example/analysis/README.md`). `NCONF` is the number of bath configurations,
and the electron-bath runs also take `NSTATE`, the sampled bath states per
configuration. Into `results/` go `<name>_<conf>[_state<s>]_noDiv` and
`_wiDiv` from CCEX (the coherence without and with the free-evolution
normalisation, or the density matrix), `_ErrCorr` and `EnsAvg_conv<N>` from
the analyzer, and `.png` plots.

The figures below were made on 2026-09-24 on a laptop (WSL, 4 ranks) from the
results of each `run.sh` with one bath configuration, by
`docs/images/tutorials/make_figures.py`; the recorded runs on the cluster used
up to 50 configurations, so these curves are the shape, not the converged
number. A coherence figure shows `|L(t)|`; a density-matrix figure shows the
diagonal of `rho(t)` and its largest off-diagonal element. The time is the one
each system's post-processing uses.

## What run.sh does

Every calculation folder has the same shape of script; this is the one of
`02_Diamond_NV_natab/ensemble_CCE/`:

```
#!/bin/sh
# 02_Diamond_NV_natab/ensemble_CCE: CCE2 Hahn echo, one run per bath configuration, then the ensemble average
#
# Recorded run: 20 bath configuration(s), 5 MPI ranks on idm.
# Here NCONF configurations on NP ranks. Set them to repeat the recorded run.
set -e
ROOT=$(git rev-parse --show-toplevel)
CCEX=$ROOT/bin/main.out
ANALYZER=$ROOT/example/analysis/Analyzer/CoherenceAnalyzer.v1/run.py
NP=${NP:-4}
NCONF=${NCONF:-1}
mkdir -p results

for i in $(seq 1 $NCONF); do
  mpirun -n $NP $CCEX -f ccein.json -I ../bath/bath_Diamond_NV_natab_${i} -o results/CCE2_Diamond_NV_natab_$i > results/process_CCE2_Diamond_NV_natab_$i
done

# --- post-processing (needs numpy, matplotlib; MPLBACKEND=Agg without a display)
# error correction: wiDiv against noDiv
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_CONF_wiDiv -err_d results/ -err_fi CCE2_Diamond_NV_natab_CONF_noDiv -vn "CONF" -v $(seq -s ' ' 1 $NCONF) -uc "ms" -err -ntn -fo results/CCE2_Diamond_NV_natab_CONF_ErrCorr
# ensemble average
python3 $ANALYZER -d results/ -fi CCE2_Diamond_NV_natab_CONF_ErrCorr -vn "CONF" -v $(seq -s ' ' 1 $NCONF) -uc "ms" -ea -ntn -fo results/CCE2_Diamond_NV_natab_EnsAvg -pl -fit -fit_initp 2 -fit_space 5
```

## The ten

| Page | What it shows |
|---|---|
| [01 Arbitrary system](Examples_01.md) | the three input file formats on a hand-written system; gCCE1, density matrix |
| [02 Diamond NV, natural-abundance ¹³C](Examples_02.md) | the basic NV case: ensemble CCE2, sampled bath states, gCCE, dsj, the Knight shift |
| [03 Diamond, several NV centres](Examples_03.md) | several NV centres: the gCCE density matrix, Ramsey and Hahn echo |
| [04 Diamond NV in a P1 bath](Examples_04.md) | an electron-spin bath: CCE2, CCE4 and CCE2e2n, pCCE, gCCE2 with the zero-field splitting, the population |
| [05 Diamond NV in a Vp bath](Examples_05.md) | the same spin echo over Vp concentrations |
| [06 Pulse sequences written out](Examples_06.md) | pulse sequences written out: Hahn echo, CPMG, WAHUHA |
| [07 h-BN V_B⁻](Examples_07.md) | DFT hyperfine and quadrupole tensors, the hyperfine-mediated term, the OH term |
| [08 Bundle clustering](Examples_08.md) | bundle clustering and its adaptive refinement |
| [09 GPU backend](Examples_09.md) | the CUDA backend: cpu, auto and cuda compared |
| [10 Coordinate frame rotation](Examples_10.md) | the coordinate frame rotation: a bath in the crystal frame |
