# v1.1.0 equivalence suite

    bash test/equivalence/run.sh        # needs bin/main.out (bash do_compile.sh wsl)

Exit status 0 only if every stage passed. Everything written goes to `_work/`,
which is gitignored.

## What this answers

`test/rotation/` shows that the rotation does what it says. It cannot show that
the branch left v1.1.0 alone, because every comparison it makes is between two
runs of the *same* binary. That is what this suite is for.

## Why it takes three axes

The obvious check -- run the new input under v1.1.0 and compare -- is not
available: v1.1.0's parser rejects every key this branch adds. So the answer is
assembled from comparisons that can each actually be made.

| axis | compares | script | verdict |
|---|---|---|---|
| 1 | the same **legacy** input, v1.1.0 binary vs this one | `axis1.sh` | byte-identical |
| 2 | **legacy form vs new form** of the same physics, one binary | `axis2.sh` | byte-identical, or a reported deviation where a unit conversion sits between them |
| 3 | a **hand-transformed** input with the rotation off vs the original with it on | `../rotation/run.sh` `[A/B]` | numeric agreement |

Axes 1 and 2 compose: if legacy(v1.1.0) == legacy(dev) and legacy(dev) ==
new-form(dev), then new-form(dev) == legacy(v1.1.0), which is the claim. Axis 1
is therefore the load-bearing one, and `run.sh` stops if it fails.

## The systems

Existing fixtures, reduced. Equivalence needs reproducibility, not convergence,
so `rbath` and `nstep` are cut until a case runs in seconds on one core. What
each was cut to, and from what, is in the generator next to the numbers.

| | system | reaches | wall (1 core) |
|---|---|---|---|
| S1a/b | NV + natural-abundance 13C, CCE-2, ensemble and sampled bath | the ordinary CCE path, point-dipole hyperfine | 0.3 s |
| S2a | the same, gCCE-2 | gCCE | 0.35 s |
| S2b/c | NV + P1 10 ppm, CCE-2 and gCCE-2 | `Defect` with `naddspin = 0`, `navaax = 12`, `avaaxfile` | 0.3 / 1.0 s |
| S3 | P1 bath, `Defect` `naddspin = 1` (14N), `navaax = 4` | `gyros` `eqs` `hypf` `efg` `zfs` `detuning`, state and extra-state files | 1.4 s |
| S4 | hBN V_B, `hf_readmode 3` + `qd_readmode 2` | both DFT tensor readers, the vertex boundary, `hf/qd_tensor_frame` | 0.44 s |
| S4b | the same over `Afile` instead of `Afile_vertex` | `CheckBD_Range`, the axis-aligned MinDif/MaxDif box | 0.43 s |
| S5 | two NV, gCCE, explicit `intmap` self/ZFS tensors | `intmap` `tensor_frame`, `reference_qubit`, multi-qubit rotation | 0.29 s |

S3's fixture is written for 24 ranks at `rbath` 850 (469 spins, 80 s). At
`rbath` 360 it keeps 45 spins and 518 pair clusters, which still reaches every
Defect code path. S4's is cut from `rbath` 20 to 6: 3897 spins and 121 s become
101 spins and 0.44 s, well inside the region the A and Q files cover, so both
readers stay on their file path rather than falling back.

S4 and S4b differ only in which boundary test the tensor file selects, and a file
picks exactly one: `Afile_vertex` carries `v1..v8` vertices and reaches
`CheckBD_vertex`, `Afile` carries `MinDif[A]` / `MaxDif[A]` and reaches
`CheckBD_Range`. Without S4b the axis-aligned box -- the one thing that makes a
rotation of the tensor file's own coordinates impossible -- would never run.

S4b needs `hf_ignore_oor = 1` and no quadrupole: both boxes are wider than the set
of atoms their files list, so a spin can sit inside the box with no row of its own.
The HF reader has a flag for that and the QD reader does not, which leaves the QD
range boundary uncovered. So are `qd_readmode` 3 and 4, which upstream also
records as untested rather than verified.

## What is compared

The **raw** `_noDiv` / `_wiDiv`. Not the error-corrected file and not the
ensemble average: both smooth, and a real difference could hide under them.

`errcorr.sh` is the other half of that decision. It runs the analyzer on each
baseline and reports `max|L|` before and after correction, so the suite says out
loud whether the systems it compared are numerically sound. They mostly are;
`S2a` is not -- CCE-2 on that 13C bath produces `max|L| = 4.2` at two points,
which the correction takes back to 1. That is the ordinary CCE division
divergence, it is present in v1.1.0 exactly as it is here, and it does not
disqualify the system as an equivalence fixture. It does mean the raw file is
not physics.

## Unit conversions

`gyros`, `eqs`, `hypf`, `efg`, `zfs` and `detuning` are each given twice: once in
the unit whose scale factor is exactly 1 -- where byte identity with the legacy
array is required -- and once in a second unit, where it is not. Three of the
four second units come back byte-equal anyway. `efg` in `V/angstrom^2` lands
1e-10 away, which is one unit in the last printed place of a file written to ten
decimals.

## Refusals

`reject.sh` re-checks the refused combinations on the physical fixtures rather
than the toy ones, and requires two things of each: the run stops with the
expected message, and neither DFT tensor reader announced itself first. The
second half is what would regress quietly if a guard ever drifted later into the
run -- the refusal would still happen, but only after a file had been read in a
frame it was not written in.
