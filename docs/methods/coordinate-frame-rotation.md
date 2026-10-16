# The coordinate frame rotation

What the `frame.rotation` block of the input file does (the rotation is on when the block is given, `enabled`
defaulting to true), in the order the run does it: what is transformed, what is left alone, what is refused, and what
each of its keys means. The short form is the [frame.rotation](../user/InputFile_General.md#tag_coordinate_frame_rotation)
entry of the input-file reference.

Where it happens in the run. The magnetic field is settled first, before any file is opened. The bath is
then read and the rbath selection (with its inner radius) is applied, and the HF and QD tensor files are
looked up -- by POSITION, against un-moved coordinates, because the position columns of
those files are written in the source bath frame. Tensor COMPONENTS are a separate
matter: they sit in whichever basis hf_tensor_frame / qd_tensor_frame declares. <br>
Only then does the frame change, in one step: every qubit position and every bath-spin
position is rotated, and the tensors that are still in the bath frame -- HF, QD, and the
Qubit.intmap entries whose provenance says so -- are transformed with them, while tensors
already in the computational frame are left as they are. Every Defect entry whose
<tt>coordinate_frame</tt> is <tt>"bath"</tt> is normalized in this same step, before its
on-site spins are constructed; clusterization follows. <br>

The transformation is <b>passive</b>: it builds the rotation R with <br>
<p align="center">R · <i>qubit_axis</i>/|<i>qubit_axis</i>| = [0, 0, 1]</p>
and applies r<sub>new</sub> = r<sub>0</sub> + R (r<sub>old</sub> − r<sub>0</sub>) to every qubit AND every bath spin
alike, with r<sub>0</sub> the position of the REFERENCE QUBIT -- the one named by
"reference_qubit", or the only qubit when that key is omitted. The reference qubit has
r<sub>old</sub> = r<sub>0</sub> and so does not move; everything else moves rigidly around it. It is
NOT the active rotation that would carry bath_axis onto qubit_axis. <br>
The frame is
<p align="center">e<sub>z</sub> = q̂, &nbsp;&nbsp; e<sub>x</sub> = (q̂ × b̂)/|q̂ × b̂|, &nbsp;&nbsp; e<sub>y</sub> = e<sub>z</sub> × e<sub>x</sub></p>
with e<sub>x</sub>, e<sub>y</sub>, e<sub>z</sub> as the ROWS of R (q̂ and b̂ the unit qubit_axis and bath_axis). bath_axis therefore only fixes the
azimuth of x and y about the new z; it does not change where z points. <br>

Both vectors are given in the bath file's own Cartesian frame. <b>Transformed:</b> every
Qubit.xyz, every bath-spin position, the bath-frame HF / QD tensors stored next to them
(every point-dipole tensor CCEX computed itself, plus the tensor-file ones where
hf_tensor_frame / qd_tensor_frame say "bath"), and the [Qubit.intmap](../user/InputFile_Qubit.md#tag_intmap)
tensors whose provenance calls for it -- the automatic qubit-pair point-dipole tensor,
and every explicit entry with [tensor_frame](../user/InputFile_Qubit.md#tag_tensor_frame) = "bath", self/ZFS
entries included. <br>
For each Defect entry with [coordinate_frame](../user/InputFile_Defect.md#tag_df_coordinate_frame) = "bath",
its [axis](../user/InputFile_Defect.md#tag_df_axis), nuclear_positions displacement vectors, and hyperfine / efg / zfs tensors
are transformed by the same R. The displacement is a relative vector, so it receives
R r without the rotation origin. Defect detuning is scalar and is not transformed.
<br>

<b>Not transformed:</b> explicit Qubit.intmap entries with tensor_frame = "qubit"; the
"zfs" of the single-qubit [Qubit.file](../user/InputFile_Qubit.md#tag_Qubit_file) form, which has no tensor_frame to declare and is read as a computational-frame tensor; each Defect entry whose coordinate_frame is "qubit"; scalar
Defect detuning; and [bfield](../user/InputFile_General.md#tag_bfield). These are read as written in the
computational frame (or are frame-independent). Legacy Defect entries without a
coordinate_frame also default to "qubit" to preserve their previous numerical meaning.
<br>

For "bfield" that has a concrete meaning: with the rotation on, the computational z axis
IS the qubit axis, so
``` "bfield" : [0.0, 0.0, 500.0] ```
is 500 G along +qubit_axis -- along the physical [111] for an NV -- and [0,0,-500] is the
same magnitude the other way. Both signs are accepted; what is required is the LINE, not
the direction. A scalar "bfield" means [0,0,bz] and works the same way. <br>
A field with any transverse component is refused: a tilted field is not supported
alongside the rotation in this version. The check runs once the input file has been read, against the final field, and compares
|B<sub>⊥</sub>| = √(B<sub>x</sub>² + B<sub>y</sub>²) against 10<sup>−12</sup> max(1, |B|) -- a tight
bound, because nothing rotates these numbers and anything larger is a component the user
meant. A zero field is always allowed. <br>
The run log prints the field both ways round, so the physical direction is visible:
R<sup>T</sup> [0, 0, B<sub>z</sub>] = B<sub>z</sub> · <i>qubit_axis</i>/|<i>qubit_axis</i>|. <br>

Refused: collinear bath_axis and qubit_axis, because the cross product then leaves the
azimuth undefined -- the one exception is bath_axis = qubit_axis = [0,0,1], which is the
identity. Also refused, with more than one qubit: an external HF or QD tensor file (a
restriction that predates this option -- see hf_tensor_frame below), a missing
reference_qubit or qubit_position_frame, and an explicit intmap entry with no
tensor_frame. <br>
Multi-qubit coherence is validated with method = gCCE. cCCE remains single-qubit and
src/calculation/simulator_cce.cpp is untouched by this option. <br>

"reference_qubit" names the qubit the rotation is taken about, by its exact
[Qubit](../user/InputFile_Qubit.md#tag_qubit) name. Its position is r<sub>0</sub>, and it is the one thing that
does not move. With a single qubit it is optional -- that qubit is used, and one supplied
through [Qubit.file](../user/InputFile_Qubit.md#tag_Qubit_file) is named "q0". <b>With more than one qubit it is
required</b>, there being no sensible default. It must match an existing qubit exactly,
and duplicate qubit names are refused. Any index may be named, not just the first. <br>
When bath_axis and qubit_axis -- and therefore R -- are held fixed, changing only
reference_qubit changes the rotation ORIGIN. The two coordinate sets then differ by the
rigid translation (I − R)(r<sub>0</sub> − r<sub>1</sub>), so distances, interaction tensors and coherence
are unchanged; only the printed coordinates move. <br>
If reference_qubit is changed to a qubit with a DIFFERENT physical axis, qubit_axis must
be changed to that qubit's actual source-frame axis as well. R is then different
too, and the relation above no longer holds. <br>

"qubit_position_frame" states the frame every [Qubit.xyz](../user/InputFile_Qubit.md#tag_xyz) is written in.
One value for the whole QubitArray -- per-qubit frames are not accepted. The only
supported value is <tt>"bath"</tt>: the qubit positions use the same Cartesian basis AND
the same origin as the bath files, so a single R and r<sub>0</sub> move every qubit
and every bath spin together. <br>
<tt>"qubit"</tt> is refused. It would put the qubit positions in the computational frame
while the bath files are still read in the source frame, and readBathfiles would then
compare the two in its rbath selection, its minimum-distance bookkeeping and
its point-dipole tensors. <br>
<b>Required with more than one qubit when the rotation is on.</b> With a single qubit it may be
omitted and is treated as <tt>"bath"</tt>, so existing inputs are unaffected. <br>

The old top-level names of this block are refused with the new place ([old names](../user/InputFile.md#OldNames)). <br>

"hf_tensor_frame" and "qd_tensor_frame" say which basis the COMPONENTS of an external
tensor file are written in. "bath" means the original bath Cartesian basis, so each
tensor is transformed as T → R T R<sup>T</sup>; "qubit" means they are
already in the computational basis and are left alone. They say nothing about the
POSITION columns of those files, which are always read in the original, unrotated frame
-- that is what the lookup matches against, and the files themselves are never touched.
<br>
The hyperfine {"dipole": "point", "contact": true} is a special case: it keeps only the isotropic Fermi-contact term from
the file and rebuilds the anisotropic part as a point-dipole tensor from the source
geometry, so what is stored is bath-frame whatever the key says. hf_tensor_frame is
therefore neither required nor used there; giving it is allowed and reported as
not applicable. <br>
With "dipole": "dft" the two kinds can occur together: a spin the tensor file does not
cover falls back to a point-dipole tensor built from the source geometry. Those fallbacks
are always transformed, even under hf_tensor_frame = "qubit", while the matched
file tensors are left alone. <br>
The quadrupole structures "strained" and "bubble" go through the same post-read transformation as "plain". The test
suite exercises their two-file reader with distinct wiDef and woDef boundaries and
checks that each boundary comes from the requested file. "plain" additionally has the
full hand-transformed versus automatic-rotation tensor/coherence comparison. <br>
With more than one qubit, an external HF or QD tensor file is refused outright, before
the file is opened. <b>This is not a new restriction introduced by the rotation:</b>
readHftensorfile and readQdtensorfile have always refused more than one qubit, because a tensor
file is written around a single defect and which qubit it belongs to is undefined. The
rotation only reports it earlier and in its own terms; both original guards remain in
place. A per-qubit tensor-file schema is left for later work. <br>

There is no default. With the rotation on and a [hyperfine](../user/InputFile_General.md#tag_hyperfine) "dipole" of "dft"
(or a [quadrupole](../user/InputFile_General.md#tag_quadrupole) that is not "off"), the corresponding key MUST be
given:
transforming a tensor twice, or not at all, produces something that still looks like a
tensor, so it is the one mistake nothing downstream can catch. With the point dipole alone the key
has no numerical effect -- the point-dipole tensors are built from the source geometry
and are always transformed. <br>
