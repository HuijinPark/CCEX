# coordinate_frame_rotation tests

    bash test/rotation/run.sh        # needs bin/main.out (bash do_compile.sh wsl)

    # legacy baseline for the DFT tensor readers, which nothing else covers
    bash test/rotation/dft_legacy_check.sh <reference-binary> [<test-binary>]

Exit status 0 only if every check passed. Everything the suite writes goes to
`_work/`, which is wiped at the start of each run and is gitignored.

## Where the frame changes

    read options -> CLI flags -> B FIELD ALIGNMENT CHECK
      -> read bath (source frame) -> bathadjust -> rbath/rbathcut selection
      -> HF tensor lookup (source frame) -> QD tensor lookup (source frame)
      -> ROTATE positions and bath-frame tensors/Defect templates together
      -> Defect construction -> clusterization

`Config_validateBfieldAlignment` runs where `Config_resolveEvolution` does,
immediately after the getopt loop, because `-B` lands there and changes the
answer. It never modifies the field -- "bfield" is already in the computational
frame -- it only refuses a transverse component while the rotation is on.

Selection happens before the rotation on purpose: it is a distance to the
reference qubit, and the rotation is taken about that point, so it comes out the
same either way. The tensor lookups have to come first, because those files are
keyed by position in the original frame.

## Which tensors move

Per (bath spin, qubit), not per run. `readHftensorfile` records whether each
tensor it stored was computed here from the source geometry or read from the
file; `applyCoordinateFrameRotation` transforms the first kind always, and the
second kind only when `hf_tensor_frame` says `bath`. The two occur together:

| hf_readmode | what is stored | transformed? |
|---|---|---|
| 0 | point dipole | always |
| 1 | point dipole + isotropic FC (the file's tensor is discarded) | always; `hf_tensor_frame` does not apply |
| 2 / 3, matched | DFT tensor from the file | only when `hf_tensor_frame = bath` |
| 2 / 3, fallback | point dipole (spin not covered by the file) | always |

`qd_readmode` keeps a single flag: `readQdtensorfile` has no computed-here
fallback to distinguish. Modes 3 and 4 share this post-read path but have no
fixture here, so they are untested, not verified.

## Why it is split in two

`test_rotation_unit.cpp` decides everything that does not need a calculation --
the matrix, the fixed point, what the transform preserves. It links against
`obj/*.o` from the main build, so it always tests the same code the binary runs.
(`main.cpp` is compiled at link time by the top-level Makefile rather than into
`obj/`, which is why the test can define `rank` / `verbosity` / `nprocess`
itself without a duplicate-symbol clash.)

`run.sh` covers what only a real run can show: that a disabled rotation changes
nothing, that the bath order and labels and sampled states survive, that all
bath files get the same transform, that the input files are not touched, and
that the refused combinations really exit non-zero with the message they claim.

`check_geometry.py` compares an unrotated log against a rotated one. It rebuilds
R from the convention (ez = q, ex = q x b, ey = ez x ex, as ROWS) instead of
holding a table of expected numbers, so a sign flip in the C++ fails here rather
than being copied into the expected values. Positions come out of the BathArray
report at 3 decimals -- that, and nothing else, sets its tolerances.

## The test systems

**Nuclear bath.** One qubit at `[10, 20, 30]` (deliberately not the origin, so
that "r0 stays fixed" is a real check), 6 `13C` in `bath_a` and 3 more in
`bath_b`, all inside `rbath = 10`. One spin sits at `r0 + (1,1,1)`, exactly on
the [111] axis, so after the rotation it must land at `r0 + (0, 0, sqrt(3))`.

**P1 defect bath.** The same 6 positions as `P1` electron spins, with a `Defect`
section: `naddspin = 0`, `navaax = 12`, a scalar per-axis `detuning`, and
`avaax_p1` pinning the assignments so no RNG is involved. It is checked twice:

- an **identity** rotation must leave the run byte-identical, which isolates the
  feature itself from any geometry change;
- a **[111]** rotation must move the bath (the coherence has to differ) while
  `paxes[i]` and the `detuning[iax]` table they index into stay exactly the same.

Per-spin detuning is `detuning[paxes[i]]`, so those two together pin it down;
CCEX does not print the per-spin value anywhere reachable at `naddspin = 0`.

`naddspin = 0` is deliberate. `rxyzs`, `hypf` and `efg` only exist for added
subspins. This legacy P1 entry has no `axis` or `coordinate_frame`, so it keeps
the backward-compatible `qubit` default; its scalar detuning is frame-independent.

**NV defect bath.** Six spin-1 `NV_0` bath spins use three deterministically
selected frozen configurations. One physical `axis` belongs to `NV_0`; the three
`navaax` entries are configurations, not three orientations. The suite describes
the same [111] NV five ways: shared `{D,E,unit}`, shared tensor, the old three-entry
tensor list, D in GHz, and an already transformed qubit-frame axis. It also expresses
the same indexed detuning as a legacy MHz list, a unit-aware MHz object, and a GHz
object. Their numerical coherence outputs must agree. The suite checks that a shared
tensor is broadcast and that a bath-frame Defect is reported as transformed.

**Unit-aware on-site Defect.** One spin-1/2 `DUT` main spin has one on-site
spin-1 `14N`, so all four newly wrapped quantities reach a real Hamiltonian:
`gyros`, `eqs`, `hypf`, and `efg`. A legacy input is compared with canonical
unit objects, then one unit is changed at a time. The suite requires numerical
equivalence for MHz/GHz hyperfine, radkHz/G versus MHz/T gyro, 1e-30 m^2 versus
mbarn quadrupole moment, and Hartree/Bohr^2 versus V/angstrom^2 EFG. Invalid
units, value counts, and missing tensor frames are also required to fail.

`bath_a`, `bath_b` and `bath_p1` are toy configurations written for this suite.
They are not physical diamond lattices and are not meant to be.

## bfield and defect_axis_reference

With the rotation on, the computational z axis is the qubit axis and "bfield" is
read in that frame: `[0,0,500]` is 500 G along +qubit_axis, `[0,0,-500]` the same
magnitude the other way, and both are accepted. A transverse component is
refused. `check_srcfield.py` reads the field back out of the log in both frames
and requires `R^T [0,0,Bz] = Bz * normalize(qubit_axis)` -- the expected value is
built from qubit_axis alone, so it is not just re-deriving what the code did.

`defect_axis_reference: "qubit_axis"` is validated and logged and nothing else.
The P1 test runs the same defect system with and without it and requires
byte-identical coherence plus identical paxes and detuning, so "changes nothing"
is checked rather than asserted.

## The A/B equivalence test

The strongest check here, and the only one that can catch a frame-dependent
quantity nothing transforms. `make_manual_frame.py` writes the same physical
system out by hand in the qubit-aligned frame -- bath positions about r0, tensor
file positions and boundary vertices about the origin (they are relative to the
qubit), tensor blocks as R T R^T -- and that input is run with no rotation option
at all. It has to agree with the original input run with the rotation on.

Numerical agreement is the criterion, not byte identity: the two reach the same
numbers by different floating-point routes, one parsing what the other computes.
Byte identity is reported when it happens, and at the time of writing it does,
for the point-dipole system and for the DFT-tensor one alike.

## Multi-qubit

The rotation moves every qubit and every bath spin with one R about one r0 -- the
position of `reference_qubit`, which therefore does not move. That is what makes
the source-frame rbath / rbathcut selection valid for several qubits too: nothing
changes any relative distance.

MQ/C holds `bath_axis` and `qubit_axis` fixed -- so `R` is the same in both runs --
and moves only the rotation ORIGIN, from NV0 to NV1. The two coordinate sets then
differ by the rigid translation `(I-R)(r0-r1)`, so distances, tensors and
coherence are unchanged while the printed coordinates are not. The test checks
both halves: the first would break if the key were ignored, the second if the
index were hard-coded to 0.

This is not a general claim that changing `reference_qubit` never changes the
answer. If the new reference qubit has a different physical axis, `qubit_axis`
has to change with it; `R` is then different and the translation relation does
not apply.

`Qubit.intmap` tensors carry their own provenance, and the rotation reads it
rather than guessing from what kind of tensor it is:

| provenance | source | transformed? |
|---|---|---|
| `INTMAP_DEFAULT_ZERO` | never written | no (zero either way) |
| `INTMAP_AUTO_SOURCE_GEOMETRY` | `setIntmap_dipAuto`, from the source qubit positions | always |
| `INTMAP_EXPLICIT_BATH` | explicit, `tensor_frame: "bath"` | always |
| `INTMAP_EXPLICIT_QUBIT` | explicit, `tensor_frame: "qubit"` | never |
| `INTMAP_EXPLICIT_UNSPECIFIED` | explicit, no `tensor_frame` | fatal at nqubit>1; kept at nqubit==1 (legacy `qzfs`) |

Each explicit `Qubit.intmap` entry may also carry `unit: "Hz"`, `"kHz"`,
`"MHz"`, or `"GHz"`. An absent unit retains the legacy kHz interpretation.
The suite feeds the same non-zero self tensor through all five spellings and
requires identical internal intmap tensors. Unknown and non-string units must
fail before the simulation begins.

`validateCoordinateFrameRotationInputs` (reader.cpp) runs in main right after
`Config_validateBfieldAlignment` -- every option source in, no reader started --
and refuses what cannot be made consistent. The validation tests check the
message AND that neither `Read the Hyperfine interaction from DFT inputfile` nor
`Read the Quadrupole interaction from DFT inputfile` appears in the log, which is
the part that would silently regress if the check drifted later in the run.

External HF/QD tensor files stay single-qubit. That restriction is **not new**:
`readHftensorfile` and `readQdtensorfile` have always refused `nqubit > 1`, and
both guards are still there. The rotation just reports it before a file is
opened.

## What is still refused

Multi-qubit runs with an external HF or QD tensor file (a pre-existing
restriction, see above), `qubit_position_frame: "qubit"`, an explicit intmap
tensor with no `tensor_frame` at `nqubit > 1`, and a tilted magnetic field
alongside the rotation. cCCE remains single-qubit -- `simulator_cce.cpp` is
untouched -- so multi-qubit coherence is tested with gCCE.

## What is read as written, and what is not

Four different things get lumped together as "the ZFS"; they are not the same:

| input | frame | transformed? |
|---|---|---|
| `Defect` entry with `coordinate_frame: "bath"` | source bath/crystal frame | axis/rxyzs and hypf/efg/zfs transformed |
| `Defect` entry with `coordinate_frame: "qubit"` (also the legacy default) | computational frame | no |
| `Defect.detuning` | scalar frequency | never |
| `bfield` | computational-frame input | never |
| legacy top-level `qzfs` | read as written, for backward compatibility | never |
| `Qubit.intmap` self/ZFS tensor | whatever its `tensor_frame` says | `"bath"` yes, `"qubit"` no |

`qzfs` belongs to the single-qubit `qubitfile` path and has no `tensor_frame`
sub-key to declare, which is why it is kept rather than transformed. A ZFS given
through `Qubit.intmap` does have one, and is transformed when it says `"bath"` --
so "the ZFS is not rotated" is only true of the legacy spelling.

A `Defect` section written in the crystal frame must declare
`coordinate_frame: "bath"`. `applyCoordinateFrameRotation` then transforms its
normalized axis and relative `rxyzs` with R, and its `hypf`, `efg`, and `zfs`
with R T R^T before Defect construction. Scalar `detuning` is unchanged. A
legacy entry without either a new `axis` or object-form `zfs` defaults to
`"qubit"`, preserving its old result.
