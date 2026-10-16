# Qubit tags
The `Qubit` block describes the qubits (old key names are in the table of [old names](InputFile.md#OldNames)):
[qubit](#tag_qubit), one entry per qubit (the number of qubits is the length of the list), each with its position
written out ([xyz](#tag_xyz)) or read from a file ([file](#tag_Qubit_file)). Beside the list: [intmap](#tag_intmap),
[alpha_state](#tag_alpha_state), [beta_state](#tag_beta_state) and [psi0](#tag_psi0). <br>

<b>Example:</b> ``` "Qubit" : {"qubit": [{"file": "bath/bath_Diamond_NV_natab_defect"}], "alpha_state": {"ms": [1]}, "beta_state": {"ms": [0]}} ```

Before 2.0.0 one qubit from a file was a form of its own, `"Qubit": {"file": ..., "spin": ..., "name": ...}`; it is the
entry `"qubit": [{"file": ..., "spin": ..., "name": ...}]` now, which also takes `gyro` and `detuning`, and
`Qubit.file`, `Qubit.spin` and `Qubit.name` are refused with that form. The
Overhauser field of the bath on the qubits is
[hamiltonian.common.qubit_bath_mean_field](InputFile_General.md#tag_qubit_bath_mean_field).

<details>
<summary><b>Example:</b> Click here to show example code</summary>
```
{
    "Qubit" : {
                "qubit" : [
                            {
                                "name" : "Q1",
                                "spin" : 0.5,
                                "gyro" : -17608.597050,
                                "xyz"  : [1.0, 2.0, 3.0],
                                "detuning" : 0.0
                            },
                            {
                                "name" : "Q2",
                                "spin" : 1.0,
                                "gyro" : -17608.597050,
                                "xyz"  : [4.0, 5.0, 6.0],
                                "detuning" : 0.0
                            }
                ],
                "intmap" : [
                            {
                                # zero-field splitting of Q2
                                "between" : ["Q2", "Q2"],
                                "tensor"  : [[1.0, 2.0, 3.0],
                                            [4.0, 5.0, 6.0], 
                                            [7.0, 8.0, 9.0]]
                            },
                            {
                                "between" : ["Q1", "Q2"],
                                "tensor"  : [[1.0, 2.0, 3.0],
                                            [4.0, 5.0, 6.0], 
                                            [7.0, 8.0, 9.0]]
                            }
                ],
                "alpha_state" : {"ms": [0.5, 0]},
                "beta_state"  : {"vector": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]},
                "psi0" : [0.0, 0.70710678, 0.0, 0.70710678, 0.0, 0.0]
              }
}
```
The first `intmap` entry is the zero-field splitting of Q2 (`between` names the same qubit twice); the
second is the Q1–Q2 coupling. `alpha_state` gives one m<sub>s</sub> per qubit in the order of `qubit`:
|m<sub>1</sub> = +1/2, m<sub>2</sub> = 0⟩. The `vector` of `beta_state`, and `psi0`, have the length of the whole qubit
space, here 2 × 3 = 6, ordered as kron(Q1, Q2) with the last qubit varying fastest: `beta_state` is |−1/2, +1⟩, the same
state as `{"ms": [-0.5, 1]}`.
</details>

- [Qubit](#tag_Qubit) 
    - [qubit](#tag_qubit): [name](#tag_name), [spin](#tag_spin), [gyro](#tag_gyro), [xyz](#tag_xyz) or
      [file](#tag_Qubit_file), [detuning](#tag_detuning)
    - [intmap](#tag_intmap): [between](#tag_between), [tensor](#tag_tensor), [unit](#tag_intmap_unit), [frame](#tag_intmap_frame)
    - [alpha_state](#tag_alpha_state), [beta_state](#tag_beta_state)
    - [psi0](#tag_psi0)

<a id="tag_Qubit"></a>
## Qubit
<b>Description:</b><br>
Qubit information. Without a `Qubit` block (or without a non-empty `qubit` list) the run stops. <br>

<b>Data Type:</b> Dictionary

<b>Default:</b> -

<b>Example:</b> 
    ``` "Qubit" : { "qubit" : [ { "name" : "Q1", "spin" : 0.5, "gyro" : -17608.597050, "xyz" : [1.0, 2.0, 3.0], "detuning" : 0.0 }, 
                                 { "name" : "Q2", "spin" : 1.0, "gyro" : -17608.597050, "xyz" : [4.0, 5.0, 6.0], "detuning" : 0.0 } ],
                    "alpha_state" : {"ms": [0.5, 0]}, "beta_state" : {"ms": [-0.5, 1]} } ```

<a id="tag_Qubit_file"></a>
## file
<b>Description:</b><br>
In place of [xyz](#tag_xyz): the file the position of this qubit is read from, the first line holding three numbers
x y z (Angstrom), in the frame of the bath file (e.g. a defect file of `example/tools/BathGenerator`). Relative to the
working directory ([configuration.path](InputFile_General.md#tag_configuration_path) does not apply to it). The other
keys of the entry are as for any qubit (name `q0` by default for the first, spin 1, the electron's gyro, detuning 0).
A file that does not open stops the run; an entry with both `file` and `xyz` is refused. Before 2.0.0 this was the
`Qubit.file` form (with `Qubit.spin`, `Qubit.name`). <br>

<pre>
    File format is the following :
    [x] [y] [z]
</pre>
<pre>
    example)
    10.0 12.0 13.0    - x,y,z position of qubit in cartesian, angstrom
</pre>

<b>Data Type:</b> String (File path)

<b>Default:</b> - (the entry gives `xyz` or `file`)

<b>Example:</b> ``` "Qubit" : {"qubit": [{"file": "defect", "name": "NV", "spin": 1.0}], "alpha_state": {"ms": [1]}} ```

<b>Parent-tags:</b> [qubit](#tag_qubit)

<a id="tag_qubit"></a>
## qubit
<b>Description:</b><br>
The qubits, one entry each: name, spin number, gyromagnetic ratio, position (`xyz`, or `file`) and detuning. The
number of qubits is the length of the list (at least one); [alpha_state](#tag_alpha_state) and
[beta_state](#tag_beta_state) give their levels in this order. <br>

<b>Data Type:</b> Dictionary List

<b>Default:</b> -

<b>Example:</b> 
    ``` "qubit" : [ { "name" : "Q1", "spin" : 0.5, "gyro" : -17608.597050, "xyz" : [1.0, 2.0, 3.0], "detuning" : 0.0 }, 
                    { "name" : "Q2", "spin" : 1.0, "gyro" : -17608.597050, "xyz" : [4.0, 5.0, 6.0], "detuning" : 0.0 } ] ```

<b>Sub-tags:</b>

- [name](#tag_name) : Qubit name
    - [spin](#tag_spin) : Spin number
    - [gyro](#tag_gyro) : Gyromagnetic ratio \(rad/ms/G\)
    - [xyz](#tag_xyz)   : x,y,z position of qubit in cartesian, angstrom
    - [file](#tag_Qubit_file) : in place of xyz, a file holding the position
    - [detuning](#tag_detuning) : Detuning frequency \(kHz\)

<a id="tag_name"></a>
## name
<b>Description:</b><br>
Name of the qubit. <br>
<b>Data Type:</b> String

<b>Default:</b> "q%d" with the 0-based position in the list, e.g. "q0", "q1", "q2", ...

<b>Example:</b> ``` "name" : "Q1" ```

<b>Parent-tags:</b> [qubit](#tag_qubit)

<a id="tag_spin"></a>
## spin
<b>Description:</b><br>
Spin number of the qubit. <br>
<b>Data Type:</b> Float

<b>Default:</b> 1.0

<b>Example:</b> ``` "spin" : 0.5 ```

<b>Parent-tags:</b> [qubit](#tag_qubit)

<a id="tag_gyro"></a>
## gyro
<b>Description:</b><br>
Gyromagnetic ratio of the qubit. <br>
<b>Data Type:</b> Double, or `{"value": ..., "unit": ...}` with the units of
[Defect gyros](InputFile_Defect.md#tag_df_gyros): `radkHz/G` (= `rad/ms/G`), `kHz/G`, `MHz/G`, `MHz/T`, `Hz/T`,
`radMHz/T` (= `rad/us/T`), `rad/s/T`

<b>Unit:</b> rad/ms/G

<b>Default:</b> -17608.597050

<b>Example:</b> ``` "gyro" : -17608.597050 ```

<b>Parent-tags:</b> [qubit](#tag_qubit)

<a id="tag_xyz"></a>
## xyz
<b>Description:</b><br>
x,y,z position of the qubit in cartesian coordinates. <br>
<b>Data Type:</b> Double List, or `{"value": ..., "unit": ...}` with the unit `A` (angstrom) or `nm` ([Units](InputFile.md#Units))

<b>Unit:</b> angstrom

<b>Default:</b> -

<b>Example:</b> ``` "xyz" : [1.0, 2.0, 3.0] ```

<b>Parent-tags:</b> [qubit](#tag_qubit)

<a id="tag_detuning"></a>
## detuning
<b>Description:</b><br>
Detuning frequency of the qubit: adds detuning × S<sub>z</sub> to its Hamiltonian, shifting its levels. The phase of the
coherence turns faster while |L| stays the same, and in gCCE the resonances of the qubit with the bath spins move. For
a static shift of the qubit levels that the other keys do not give (00_quickstart's NV, gCCE1 Ramsey: |L| with 0 and
1000 kHz differ by 3.5e-6). <br>
<b>A plain number is kHz here, and MHz in [Defect.species[].detuning](InputFile_Defect.md#tag_df_detuning)</b> (the same
name, a factor 1000 apart): a plain number other than 0 is noted in the log, `Note: Qubit.qubit[i].detuning x is read in
kHz (= y MHz ...)`; `{"value": ..., "unit": "MHz"}` names the unit. `[x]` is refused: write `x`. <br>
<b>Data Type:</b> Double, or `{"value": ..., "unit": ...}` with the unit `Hz`, `kHz`, `MHz` or `GHz` ([Units](InputFile.md#Units))

<b>Unit:</b> kHz

<b>Default:</b> 0.0

<b>Example:</b> ``` "detuning" : 0.0 ```

<b>Parent-tags:</b> [qubit](#tag_qubit)

<a id="tag_intmap"></a>
## intmap
<b>Description:</b><br>
Interaction map between qubits. <br>
Defines the interaction map between qubits, including the self interaction of a qubit \(e.g. zero-field splitting\). <br>
It is optional: a self entry that is not given is the zero tensor, and a pair of different qubits that is not given gets the
point-dipole tensor CCEX computes from their `xyz` and `gyro`. A self entry names the qubit twice (`["q0", "q0"]` for
an unnamed first qubit). <br>
<b>Data Type:</b> Dictionary List

<b>Unit:</b> Selected independently for each entry by its optional
[unit](#tag_intmap_unit); kHz when omitted (a Defect zfs is MHz: a nonzero tensor without unit is noted in the log).

<b>Default:</b> zero self terms; point-dipole tensors between different qubits

<b>Example:</b> 
    ``` intmap : [ { "between" : ["q1","q1"],
                       "unit"    : "MHz",
                       "tensor"  : [[0.0, 0.0, 0.0],
                                    [0.0, 0.0, 0.0], 
                                    [0.0, 0.0, 0.0]] }, 
                     { ... } ] 
    ```

<b>Sub-tags:</b>

- [between](#tag_between) : Interaction qubit pair
    - [tensor](#tag_tensor)   : Interaction tensor
    - [unit](#tag_intmap_unit) : Frequency unit of this tensor
    - [frame](#tag_intmap_frame) : Frame the tensor components are written in

<a id="tag_between"></a>
## between
<b>Description:</b><br>
Names of an interacting qubit pair. <br>
This tag is required to define the interaction map for a specific qubit pair. <br>
For zero-field splitting interactions, the qubit names must be identical. \(e.g. `["q1","q1"]` \) <br>
<b>Data Type:</b> String List

<b>Default:</b> -

<b>Example:</b> ``` "between" : ["q1","q2"] ```

<b>Related tags:</b> [tensor](#tag_tensor)

<b>Parent-tags:</b> [intmap](#tag_intmap)

<a id="tag_tensor"></a>
## tensor
<b>Description:</b><br>
Interaction tensor. <br>
This tag is required to define the interaction tensor for a specific qubit pair. <br>
This tensor is used to calculate the interaction between qubits. <br>
<b>Data Type:</b> Double 2D List \(3x3\)

<b>Unit:</b> Given by the sibling [unit](#tag_intmap_unit) key; kHz when omitted.

<b>Default:</b> `[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]`

<b>Example:</b> ``` "tensor" : [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]] ```

<b>Related tags:</b> [between](#tag_between)

<b>Parent-tags:</b> [intmap](#tag_intmap)

<a id="tag_alpha_state"></a>
## alpha_state
<b>Description:</b><br>
The alpha state |α⟩ of the qubits, the one place it is given (per-qubit `alphams` and the Qubit.file form's
`Qubit.alphams` until then, now refused with this form). An object with exactly one of: <br>

- `"ms"`: [m<sub>1</sub>, m<sub>2</sub>, ...], one m<sub>s</sub> per qubit in the order of [qubit](#tag_qubit) (a list
  also for one qubit, `[m]`): |α⟩ is the product of the |m<sub>s</sub>⟩ states. A list of another length is refused.
- `"index"`: n, the n-th eigenstate of the qubit Hamiltonian (zero-field splitting, Zeeman, intmap) sorted by
  **descending** energy, starting at 0. That order changes with the field: an NV at 500 G
  along its axis orders m<sub>s</sub> = +1, −1, 0, so `alpha_state {"index": 2}`, `beta_state {"index": 1}` is the
  0 ↔ −1 transition, not +1 ↔ 0. To name a sublevel, use `"ms"` instead. The run log prints which m<sub>s</sub> state
  |α⟩ and |β⟩ mostly are (the lines `psia is mostly`, `psib is mostly`). The eigenvector's phase is fixed: its
  component of largest magnitude (the first of equal ones, within 1e-12 relative) is made real and positive, so that
  psi0 = |α⟩ + |β⟩ does not depend on the phase the eigensolver happened to pick (gCCE with a non-axial ZFS, E ≠ 0,
  depends on that relative phase). An index outside 0 .. (dimension − 1) is refused.
- `"vector"`: the state itself over the whole qubit space (kron of the qubits, the last one varying fastest), each
  component a number (real) or `[re, im]`, normalized when read.

More than one of them in one state, or none, is refused. Each side is decided on its own: `beta_state` may be
given another way or not at all. <br>
The m<sub>s</sub> of each qubit come from `"ms"`. [pulse.per_qubit](InputFile_Pulse.md#tag_pulse_per_qubit), which turns each
qubit between its own two levels, needs both sides as `"ms"` and refuses `"index"` and `"vector"` (before 2.0.0 a default level was taken there). <br>
<b>Data Type:</b> Object `{"ms": [...]}`, `{"index": n}` or `{"vector": [...]}`

<b>Default:</b> `{"ms": [1, 1, ...]}` (m<sub>s</sub> +1 on every qubit)

<b>Example:</b> ``` "alpha_state" : {"ms": [1]} ```, ``` "alpha_state" : {"index": 0} ```, ``` "alpha_state" : {"vector": [1.0, 0.0, 0.0]} ```,
``` "alpha_state" : {"vector": [[0.7071, 0], [0, 0.7071], 0]} ```

<b>Related tags:</b> [beta_state](#tag_beta_state) [psi0](#tag_psi0)

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_beta_state"></a>
## beta_state
<b>Description:</b><br>
The beta state |β⟩ of the qubits, as [alpha_state](#tag_alpha_state): `{"ms": [...]}`, `{"index": n}`
or `{"vector": [...]}`, exactly one (per-qubit `betams` and `Qubit.betams` before 2.0.0). <br>
<b>Data Type:</b> Object `{"ms": [...]}`, `{"index": n}` or `{"vector": [...]}`

<b>Default:</b> `{"ms": [0, 0, ...]}` (m<sub>s</sub> 0 on every qubit)

<b>Example:</b> ``` "beta_state" : {"ms": [0]} ```, ``` "beta_state" : {"vector": [0.0, 1.0, 0.0]} ```

<b>Related tags:</b> [alpha_state](#tag_alpha_state)

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_psi0"></a>
## psi0
<b>Description:</b><br>
Initial state of the qubit in matrix representation, symbol:  |psi0⟩ <br>
This tag is used to define a specific initial state of the qubit. <br>
If this tag is not provided, the initial state defaults to the normalized sum |psi0⟩ = (|alpha⟩ + |beta⟩) / ‖|alpha⟩ + |beta⟩‖. <br>
Read by gCCE and megCCE (their initial qubit state); cCCE and meCCE work from |alpha⟩ and |beta⟩ alone, and report a
written psi0 as not used. A given vector is used as written, not normalized: give it normalized. Use it to start the
qubits in another superposition than the default. Each component is a number (real) or `[re, im]`. <br>
<b>Data Type:</b> List of numbers or `[re, im]` pairs

<b>Default:</b> -

<b>Example:</b> ``` "psi0" : [1.0, 0.0, 0.0] ```

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_intmap_unit"></a>
## intmap unit
<b>Description:</b><br>
Frequency unit of the [tensor](#tag_tensor) in this Qubit.intmap entry. Each entry
may choose its own unit. Supported values are <tt>"Hz"</tt>, <tt>"kHz"</tt>,
<tt>"MHz"</tt>, and <tt>"GHz"</tt>, written with their case (`"mhz"` is refused with the table's spelling). The tensor is converted
to the internal rad*kHz representation while parsing. <br>
Without this key the tensor numbers are kHz. The code, not the input author, supplies the factor of 2π. <br>
This key belongs to Qubit.intmap entries only; a qubit's zero-field splitting is an intmap self entry (`Qubit.zfs`
of the old file form, always kHz, until then). A nonzero tensor without `unit` is noted in the log as read in kHz. <br>

<b>Data Type:</b> String

<b>Default:</b> <tt>"kHz"</tt>

<b>Example:</b> ``` "unit" : "MHz" ```

<b>Parent-tags:</b> [intmap](#tag_intmap)

<a id="tag_intmap_frame"></a>
## intmap frame
<b>Description:</b><br>
The frame the components of this [intmap](#tag_intmap) tensor are written in (`tensor_frame` before 2.0.0, with `"bath"`
for `"file"`). Read with [frame.rotation](InputFile_General.md#tag_frame_rotation) on, where every written
intmap entry needs it, with one qubit or more (before 2.0.0 a single qubit's entry without it was read as written); with the
rotation off it is not used (reported). <br>

- <tt>"file"</tt>  : the bath file's frame. The tensor is turned with the positions, T → R T R<sup>T</sup>.
- <tt>"qubit"</tt> : the computational frame, whose z is `qubit_axis`. The tensor is used as written.

<b><tt>"qubit"</tt> does not mean an individual qubit's local frame.</b> There is one computational frame for the whole
run. A qubit whose crystallographic orientation differs from the first one is described by ITS OWN ZFS tensor expressed
in that common frame; with <tt>"file"</tt> the same R turns it, and the off-diagonal elements that appear are correct
and must not be diagonalized away. <br>
The rule is the same for a self/ZFS entry (<tt>between: ["NV0","NV0"]</tt>) and for a pair entry
(<tt>between: ["NV0","NV1"]</tt>): the frame is taken from this key, never guessed from which kind of tensor it is. An
entry that is not given at all is either the zero tensor or, off the diagonal, the point-dipole tensor CCEX computes from
the qubit positions; the latter is always turned. <br>

<b>Data Type:</b> String

<b>Default:</b> none (required with the rotation on)

<b>Example:</b> ``` "frame" : "file" ```

<b>Parent-tags:</b> [intmap](#tag_intmap)
