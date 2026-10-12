# check_numerics/ — does this build still give the numbers an earlier one gave

Twenty-five calculations, cut from the tutorials so that each runs in seconds
(`cases/<name>/ccein.json`, the bath, gyro, tensor, state and axis files they
read in `inputs/`), run with `bin/main.out` and compared against the outputs
recorded from earlier builds on the same machine (`reference/<version>/`).

    bash test/code_verification/check_numerics/run.sh [NP]     # needs bin/main.out (bash do_compile.sh wsl-gnu), mpirun, numpy

One table per reference. The verdict of a case is the largest |difference| over
its `_noDiv` and density-matrix files; the `_wiDiv` coherence is printed but not
judged (it divides by the free evolution and amplifies rounding wherever that is
small). The exit status is 0 only when the newest `reference/head-*` shows no
difference above 1e-9: that reference is this line of development's own, so a
change that moves any number is seen. Against the tagged versions the tolerance
is 1e-4, and `expected.json` names the cases that differ for a known reason,
which the table reports as "expected" instead of FAIL.

## The cases

| name | system | method | approach | pulse | quantity |
|---|---|---|---|---|---|
| natab_cce2_hahn / _ramsey / _cpmg8 | NV + natural-abundance ¹³C, rbath 30 | CCE2 | ensemble (nstate 0) | Hahn / Ramsey (npulse 0) / CPMG-8 | coherence |
| natab_cce2_single | same | CCE2 | single (nstate 2, state files) | Hahn | coherence |
| natab_gcce2_hahn, natab_gcce2_two_qubits | same; the second's `Qubit` block names two NV, but its `qubitfile` wins and the block is never read (one NV, as in every version; found 2026-10-02) | gCCE2 | ensemble | Hahn | coherence |
| natab_dsj2 | same | dsj | ensemble | Hahn | coherence |
| natab_knight_shift | same, tilted field | CCE2 | ensemble | Hahn | coherence |
| natab_manual_hahn | same, the Hahn echo written as a `sequence` | CCE2 | ensemble | manual | coherence |
| p1_cce2_ensemble / _single | NV + P1 1 ppm, rbath 500 (103 spins), JT axes from a file | CCE2 | nstate 0 / nstate 2 | Hahn | coherence |
| p1_cce2e2n, p1_cce4, p1_pcce2_sk2 | same (rbath 400 for CCE2e2n) | CCE2e2n / CCE4 / pCCE2 sK 2 | single | Hahn | coherence |
| p1_gcce2_zfs | NV + P1 10 ppm, B = 0, zero-field splitting | gCCE2 | single | Hahn | coherence |
| p1_gcce2_population | NV + P1 10 ppm | gCCE2 | single | Ramsey | dm |
| p1_cpmg8_manual, p1_manual_hahn | NV + P1 1 ppm, `sequence` | gCCE2 | single | CPMG-8 / Hahn | coherence |
| vp_cce2 | NV + Vp 10 ppm | CCE2 | single | Hahn | coherence |
| hbn_cce2_hf3_qd2 / _hf3_qd0, hbn_hfmedi | V_B⁻ + h-BN, DFT hyperfine and quadrupole tensors | CCE2 | ensemble | Hahn | coherence |
| hbn_gcce2_dm, hbn_gcce2_oh | same | gCCE2 | single | Hahn | dm / coherence (OH term) |
| arbitrary_gcce1_dm | the hand-written three-file system of tutorial 01 | gCCE1 | single | Hahn | dm |

Everything a case draws at random in a tutorial is a file here: the bath states
(`state_1`, `state_2`), the P1 Jahn-Teller axes (`avaax_nav12`, `avaax_nav4`),
the ¹⁴N sublevels (`exstate_*`). Two binaries that read the same files can be
compared; with `seed` alone they cannot, because the random sequence changed on
this branch. `make_cases.py` wrote the cases and these files from the tutorials
and can write them again.

Not in the pool: the two-NV density matrices of tutorial 03 (`NV2_13C2`).
Every version up to 64f0a1a41, and the tutorial's own recorded result, held
`nan` in them (elements the free evolution keeps at zero, divided by zero);
since 64f0a1a41 they are finite, and the case can join the pool once the tags'
`nan` is named in `expected.json`. Also out: the WAHUHA sequences of tutorial
06, whose gCCE coherence grows past 1e4 under every version. Both are
calculations to look at, not references.

## The references

| folder | what | how it was made |
|---|---|---|
| `head-<date>-<commit>/` | this branch at that commit, on WSL, 2 ranks | `bash make_reference.sh current` |
| `v1.1.0/`, `v1.2.2/` | the tagged releases, built from the tag on the same WSL | `bash make_reference.sh v1.2.2` (a worktree of the tag under /tmp, `do_compile.sh wsl-gnu`, removed afterwards) |

v1.0.0 does not build with the current compiler (`src/bath.cpp: taking address
of temporary array`) and has no reference. The two tags give the same numbers
as each other. Against them this branch differs, on purpose, in seven cases
(`expected.json`): ensemble cCCE returned the complex conjugate for an even
pulse count until 8185aa31e (Ramsey and CPMG change the sign of the imaginary
part), deltat was a float until 39d03aa5c (the fast P1 phases and the population
at 5 ms move), and the seeded pCCE partition search draws a different partition.
The density-matrix files of the arbitrary case are indexed `out_0_1` now and
`out_0_0` then; `expected.json` maps the names (`renames`) and `compare.py` compares the values
(1.6e-8 against both tags on 2026-10-02). A case with no output, or with a reference file the run
did not write, fails whatever `expected.json` says.

When a change is meant to move a number, make the new head reference after it
(`bash make_reference.sh current`) and say why in the commit. `run.sh` prints a table for every
folder under `reference/`, but only the newest `head-*` by name decides the exit status; the earlier
`head-*` folders are kept. `run.sh` writes to `_work/`, which is ignored.

**2026-10-02: eight cases re-recorded in the three folders of that date** (`p1_cce2_ensemble`, `p1_cce2_single`,
`p1_cce2e2n`, `p1_cce4`, `p1_cpmg8_manual`, `p1_manual_hahn`, `p1_pcce2_sk2`, `vp_cce2`). Their inputs
wrote `"alphams": -1` at the top of the qubitfile form, which reads `qalphams`: every version ran
them with the NV branch +1. The inputs now say `"qalphams": -1` and the references were run again
with the same command (NP 2; the tags built from the tag), so they hold the branch the inputs
meant. The change against the old references was up to 1.26 (`p1_cce2_single`), 8e-6 for the
ensemble case. Since then CCEX refuses such a key (D24).

