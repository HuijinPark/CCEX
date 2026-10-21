# 06_cCCE_vs_gCCE

The two main methods and the options of cCCE. cCCE (the default `method`) treats the qubit as a two-level system that
only sets the field each bath spin sees (pure dephasing); gCCE puts the qubit into every cluster's Hamiltonian, so
qubit-bath flip-flops are kept and the qubit can be any spin with its zero-field splitting, at any field. Each folder:
`cd <folder> && NP=4 bash run.sh`, then compare `results/coherence.png` with `expected.png`. Every key not written in
an input has its default (`main.out -e`). Times: i7-14700K, WSL, 4 ranks.

## gCCE2_Hahn_vsCCE

- **What it shows:** 00_quickstart's NV by cCCE2 (`ccein_cCCE.jsonc`) and by gCCE2 with the NV zero-field splitting
  (`ccein_gCCE.jsonc`). At 500 G along the NV axis the decay is the same; the ESEEM ripple differs.
- **Keys that matter:** `method` gCCE, the ZFS as a `Qubit.intmap` self entry (kHz).
- **How to run:** `NP=4 bash run.sh`: 6 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; T2 1.068 ms in both; largest |L| difference 0.49 (the ripple).
- **What to change next:** a tilted or a small field, where the two separate.

![gCCE2_Hahn_vsCCE: expected result](gCCE2_Hahn_vsCCE/expected.png)

Inputs: [`ccein_cCCE.jsonc`](gCCE2_Hahn_vsCCE/ccein_cCCE.jsonc), [`ccein_gCCE.jsonc`](gCCE2_Hahn_vsCCE/ccein_gCCE.jsonc)

## cCCE2_Hahn_dsj

- **What it shows:** `cluster.mode` dsj (the clusters are the disjoint groups of the coupling graph within `rdip`)
  against the standard pairs (`ccein_standard.jsonc`).
- **Keys that matter:** `cluster.mode` dsj.
- **How to run:** `NP=4 bash run.sh`: 1.1 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; T2 1.092 ms (dsj) against 1.068 ms; largest difference 0.064.
- **What to change next:** `cluster.mode` itb or dsjitb (cCCE only).

![cCCE2_Hahn_dsj: expected result](cCCE2_Hahn_dsj/expected.png)

Inputs: [`ccein.jsonc`](cCCE2_Hahn_dsj/ccein.jsonc), [`ccein_standard.jsonc`](cCCE2_Hahn_dsj/ccein_standard.jsonc)

## cCCE2_Hahn_knightShift

- **What it shows:** in a field tilted from the NV axis ([50, 0, 100] G) the qubit's transverse field shifts the ¹³C
  precession; `hamiltonian.cCCE.knight` adds that shift to cCCE (`ccein_wi_knight.jsonc`) against without
  (`ccein_wo_knight.jsonc`).
- **Keys that matter:** `hamiltonian.cCCE.knight` true, `bfield` [50, 0, 100], the ZFS (`Qubit.intmap` self entry).
- **How to run:** `NP=4 bash run.sh`: 1.2 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; T2 1.064 ms with the Knight shift; without it the divided curve (`_wiDiv`)
  diverges (|L| up to 5e5, cut at 1.2 in the plot) - the analyzer's error correction
  (`example/analysis/Analyzer`, `-err`) repairs such curves against `_noDiv`.
- **What to change next:** gCCE on the same input (the qubit is in its Hamiltonian; `hamiltonian.gCCE` has no Knight key).

![cCCE2_Hahn_knightShift: expected result](cCCE2_Hahn_knightShift/expected.png)

Inputs: [`ccein_wi_knight.jsonc`](cCCE2_Hahn_knightShift/ccein_wi_knight.jsonc), [`ccein_wo_knight.jsonc`](cCCE2_Hahn_knightShift/ccein_wo_knight.jsonc)

## gCCE1_Hahn_handwritten

- **What it shows:** the input files written by hand: `bath` (a count line, then x y z name: one ¹³C), `gyro_13C`
  (name, spin, γ in rad/ms/G) and a spin-1 qubit at the origin in the `Qubit.qubit` list (alpha ms 0, beta ms +1);
  gCCE order 1 of the qubit density matrix for one sampled bath state.
- **Keys that matter:** `quantity` dm, `Qubit.qubit` with `Qubit.alpha_state` / `beta_state` `{"ms": [...]}`, `s_range` [1, 1], `seed`.
- **How to run:** `NP=4 bash run.sh`: 0.6 s, 0.05 GB per rank.
- **Expected result:** `expected.png`; |L| stays 1.000 (one ¹³C 37 Å away, order 1).
- **What to change next:** a second ¹³C next to the first and `order` 2; `output.savemode` `cluster_density_matrix`
  for the qubit + bath density matrix of every cluster.

![gCCE1_Hahn_handwritten: expected result](gCCE1_Hahn_handwritten/expected.png)

Input: [`ccein.jsonc`](gCCE1_Hahn_handwritten/ccein.jsonc)

## meCCE2_Hahn_bathRelaxation

- **What it shows:** `method` meCCE, the cCCE with a Lindblad master equation in each cluster: jump operators S₋ and
  S₊ on every ¹³C at γ = 0.01 rad/ms each (`rate` 1.5915e-6 MHz × 2π·10³) relax the bath spins (`ccein_meCCE.jsonc`),
  against cCCE (`ccein_cCCE.jsonc`), on 00_quickstart's system with 500 points of 4 µs. The rates of all the ¹³C add up,
  so even this small rate halves T2. Without jump operators meCCE gives the cCCE curve (checked by
  `test/code_verification/check_master_equation`).
- **Keys that matter:** `method` meCCE, `hamiltonian.meCCE.jump_operators` [{`bath_name` "13C", `operator` "-" or
  "+", `rate` (MHz)}].
- **How to run:** `NP=4 bash run.sh`: 18-20 s (nearly all meCCE), 0.06 GB per rank.
- **Expected result:** `expected.png`; T2 1.068 ms (cCCE) and 0.53 ms (meCCE).
- **What to change next:** `operator` "z" (pure dephasing of the ¹³C); `method` megCCE with a qubit jump
  (`bath_name` "qubit", in `hamiltonian.megCCE`); see [master-equation CCE](../../../docs/methods/master-equation.md).

![meCCE2_Hahn_bathRelaxation: expected result](meCCE2_Hahn_bathRelaxation/expected.png)

Inputs: [`ccein_cCCE.jsonc`](meCCE2_Hahn_bathRelaxation/ccein_cCCE.jsonc), [`ccein_meCCE.jsonc`](meCCE2_Hahn_bathRelaxation/ccein_meCCE.jsonc)
