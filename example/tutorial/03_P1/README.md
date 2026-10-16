# 03_P1

An NV in a bath of electron spins: P1 centres (substitutional N, an electron spin 1/2 whose ¹⁴N hyperfine shifts it by
one of 12 detunings, one per Jahn-Teller axis and ¹⁴N state, the `Defect` block) at 1 and 10 ppm, and Vp defects
(spin 1/2, no `Defect` block). An electron bath couples more strongly than ¹³C, so the times are microseconds and the
bath is sampled in pure states (`nstate`, drawn from `seed`). Each folder: `cd <folder> && NP=4 bash run.sh`, then
compare `results/coherence.png` with `expected.png`. Every key not written in an input has its default
(`main.out -e`). Times: i7-14700K, WSL, 4 ranks.

`bath/` holds 5 configurations at 1 ppm (`14N_1ppm/`, cells for `rbath` 850 Å) and 5 at 10 ppm (`14N_10ppm/`), the
NV position of each (`_Defect`), the P1 gyromagnetic ratio (`DiaP1_gyro`) and the Vp baths (`Vp/`, one configuration
at 10, 12 and 50 ppm); larger sets (100 configurations at 1 to 100 ppm, 99 MB) are in git history.
The recorded production runs used 50 configurations × 50 states on idm (128 ranks).

## cCCE2_Hahn_states

- **What it shows:** the Hahn echo of the NV (ms -1 / 0) in 1 ppm P1, cCCE2, averaged over 2 configurations
  (`k_range`) × 10 sampled states; one state alone gives a curve far from the mean.
- **Keys that matter:** `Qubit.alphams` -1, `bathfile` with `{k}` and `k_range` [1, 2], `nstate` 10, `seed`, `rbath`
  850 Å, `rdip` 620 Å, `deltat` 10 µs × `nstep` 300, the `Defect` block (12 detunings).
- **How to run:** `NP=4 bash run.sh`: 2.2 s, 0.06 GB per rank.
- **Expected result:** `expected.png`; T2 of the mean of 20 curves 295 µs.
- **What to change next:** `k_range` [1, 5] and `nstate` 50 (the production numbers per configuration); `pulse.count` 8.

## cCCE4_Hahn_states

- **What it shows:** the same to order 4, the convergence check of the expansion; `cluster.selection.max_per_order`
  keeps the 5000 strongest three-spin and 10000 strongest four-spin clusters.
- **Keys that matter:** `order` 4, `cluster.selection.max_per_order` [0, 0, 0, 5000, 10000].
- **How to run:** `NP=2 bash run.sh` to stay within 2 GB: 0.84 GB per rank (the cluster list); 34 s on 2 ranks.
- **Expected result:** `expected.png`; T2 of the mean 295 µs, the same as order 2 (converged at order 2 here).
- **What to change next:** `max_per_order` 0 (all clusters; much more memory).

## cCCE2_Hahn_onsite14N

- **What it shows:** each P1 written with its ¹⁴N as an on-site nuclear spin (`n_onsite_spins` 1, four Jahn-Teller
  axes with the ¹⁴N hyperfine and quadrupole tensors), so a cluster of one P1 is its electron and its ¹⁴N.
- **Keys that matter:** the `Defect` species: `n_onsite_spins`, `types`, `spins`, `gyros`, `quadrupole_moments`,
  `n_axis_states` 4, `hyperfine`, `efg`.
- **How to run:** `NP=4 bash run.sh`: 65 s, 0.07 GB per rank.
- **Expected result:** `expected.png`; T2 of the mean 322 µs (against 295 µs with the detunings of `cCCE2_Hahn_states`).
- **What to change next:** the on-site states from a file (`configuration.defect_state_file`).

## gCCE2_Hahn_zfs

- **What it shows:** at zero field the NV must be a spin 1 with its zero-field splitting (D = 2.88 GHz, the `intmap`
  tensor in kHz) in gCCE; without it the ms ±1 levels are degenerate and the divided curve (`_wiDiv`) diverges.
- **Keys that matter:** `method` gCCE, `bfield` 0, the `Qubit.qubit` entry and its `intmap` self tensor (present in
  `ccein_with_ZFS.json`, absent in `ccein_without_ZFS.json`), `nstate` 1, `seed` 8810, 10 ppm, `deltat` 50 ns.
- **How to run:** `NP=4 bash run.sh`: 1.5 s, 0.06 GB per rank.
- **Expected result:** `expected.png`; T2 with the ZFS 2.70 µs; without it |L| reaches 4.6e37 (cut at 1.2 in the
  plot).
- **What to change next:** a small field (`bfield` 10) without the ZFS; more states.

## gCCE2_Ramsey_population

- **What it shows:** the NV prepared in ms = 0 (`alphams` = `betams` = 0), no pulse, the qubit density matrix
  (`quantity` dm) of gCCE2 over 50 ms in 10 ppm P1 at 300 G: the ms = 0 population, which flip-flops with the bath
  would lower.
- **Keys that matter:** `quantity` dm, `alphams` 0, `betams` 0, `pulse.count` 0, `bfield` 300, `deltat` 10 µs ×
  `nstep` 5000.
- **How to run:** `NP=4 bash run.sh`: 76 s, 0.06 GB per rank.
- **Expected result:** `expected.png` (the y range from the data); the population stays within 0.997 to 1.003.
- **What to change next:** a scan of `bfield`: where an NV transition comes into resonance with the bath, the
  population drops.

## cCCE2_Hahn_VpBath

- **What it shows:** the Hahn echo in a bath of Vp defects (spin 1/2 without on-site nuclei: no `Defect` block) at
  10, 12 and 50 ppm, one configuration and one sampled state each.
- **Keys that matter:** `configuration.gyrofile` (`Vp/DiaVp_gyro`), `bathfile`, `nstate` 1, `seed`, `rbath`, `rdip`
  and `deltat` per concentration.
- **How to run:** `NP=4 bash run.sh`: 1.7 s, 0.06 GB per rank.
- **Expected result:** `expected.png`; T2 4.8 µs (10 ppm), 15.6 µs (12 ppm, a longer window), 1.4 µs (50 ppm) for this
  one state.
- **What to change next:** more states (`nstate` 20) before comparing concentrations.
