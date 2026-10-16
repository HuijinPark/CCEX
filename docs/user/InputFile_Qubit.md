# Qubit tags
The `Qubit` block describes the qubits, in one of two forms (old key names are in the table of
[old names](InputFile.md#OldNames)):

<b>Form 1.</b> [file](#tag_Qubit_file): one qubit, its position read from a file, with
[spin](#tag_Qubit_spin), [zfs](#tag_Qubit_zfs), [alphams](#tag_Qubit_alphams) and [betams](#tag_Qubit_betams) beside it. <br>

<b>Example:</b> ``` "Qubit" : {"file": "bath/bath_Diamond_NV_natab_defect", "alphams": 1, "betams": 0} ```

<b>Form 2.</b> [qubit](#tag_qubit): one or more qubits written out; the number of qubits is the length of the list. Beside the list: [intmap](#tag_intmap), [alpha_state](#tag_alpha_state),
[beta_state](#tag_beta_state) and [psi0](#tag_psi0). <br>

A block with `file` and any of `qubit`, `intmap`, `alpha_state`, `beta_state`, `psi0` is refused: give one form. The
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
                                "detuning" : 0.0,
                                "alphams" : 0.5,
                                "betams"  : -0.5
                            },
                            {
                                "name" : "Q2",
                                "spin" : 1.0,
                                "gyro" : -17608.597050,
                                "xyz"  : [4.0, 5.0, 6.0],
                                "detuning" : 0.0,
                                "alphams" : 0,
                                "betams"  : 1
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
                "alpha_state" : {"vector": [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]},
                "beta_state"  : {"vector": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]},
                "psi0" : [0.0, 0.70710678, 0.0, 0.70710678, 0.0, 0.0]
              }
}
```
The first `intmap` entry is the zero-field splitting of Q2 (`between` names the same qubit twice); the
second is the Q1–Q2 coupling. The `vector` of `alpha_state` and `beta_state`, and `psi0`, sit inside the `Qubit` block
and have the length of the whole qubit space, here 2 × 3 = 6, ordered as kron(Q1, Q2) with the last qubit varying
fastest: `alpha_state` is |m<sub>1</sub> = +1/2, m<sub>2</sub> = 0⟩ and `beta_state` is |−1/2, +1⟩, the same states that
`alphams`/`betams` name.
</details>

- [Qubit](#tag_Qubit) 
- Form 1: [file](#tag_Qubit_file), [spin](#tag_Qubit_spin), [zfs](#tag_Qubit_zfs), [alphams](#tag_Qubit_alphams), [betams](#tag_Qubit_betams)
- Form 2:
    - [qubit](#tag_qubit): [name](#tag_name), [spin](#tag_spin), [gyro](#tag_gyro), [xyz](#tag_xyz), [detuning](#tag_detuning),
      [alphams](#tag_alphams), [betams](#tag_betams)
    - [intmap](#tag_intmap): [between](#tag_between), [tensor](#tag_tensor), [unit](#tag_intmap_unit), [tensor_frame](#tag_tensor_frame)
    - [alpha_state](#tag_alpha_state), [beta_state](#tag_beta_state)
    - [psi0](#tag_psi0)

<a id="tag_Qubit"></a>
## Qubit
<b>Description:</b><br>
Qubit information, in one of the two forms above. Without a `Qubit` block (or with neither `file` nor a non-empty
`qubit` list) the run stops. <br>

<b>Data Type:</b> Dictionary

<b>Default:</b> -

<b>Example:</b> 
    ``` "Qubit" : { "qubit" : [ { "name" : "Q1", "spin" : 0.5, "gyro" : -17608.597050, "xyz" : [1.0, 2.0, 3.0], "detuning" : 0.0, "alphams" : 0.5, "betams" : -0.5 }, 
                                 { "name" : "Q2", "spin" : 1.0, "gyro" : -17608.597050, "xyz" : [4.0, 5.0, 6.0], "detuning" : 0.0, "alphams" : 0, "betams" : 1 } ] } ```

<a id="tag_Qubit_file"></a>
## Qubit.file
<b>Description:</b><br>
File path of the position of a single qubit, relative to the working directory
([configuration.path](InputFile_General.md#tag_configuration_path) does not apply to it). The qubit is named `q0`, its gyromagnetic
ratio is the electron's, its detuning 0; [spin](#tag_Qubit_spin), [zfs](#tag_Qubit_zfs), [alphams](#tag_Qubit_alphams)
and [betams](#tag_Qubit_betams) sit beside `file` in the block. An unreadable file stops the run. <br>

<pre>
    File format is the following :
    [x] [y] [z]
</pre>
<pre>
    example)
    10.0 12.0 13.0    - x,y,z position of qubit in cartesian, angstrom
</pre>

<b>Data Type:</b> String (File path)

<b>Default:</b> -

<b>Example:</b> ``` "Qubit" : {"file": "filepath", "alphams": 1, "betams": 0} ```

<a id="tag_Qubit_zfs"></a>
## Qubit.zfs
<b>Description:</b><br>
The zero-field-splitting tensor of the qubit of [Qubit.file](#tag_Qubit_file), 3 × 3. It has no
`unit` or `tensor_frame` sub-key: the values are kHz and are read as written, never rotated by
[frame.rotation](InputFile_General.md#tag_coordinate_frame_rotation). (The `Qubit.qubit` form gives it as an
[intmap](#tag_intmap) entry.) <br>

<b>Data Type:</b> Double List (3 × 3)

<b>Unit:</b> kHz

<b>Default:</b> zero tensor

<b>Example:</b> ``` "Qubit" : {"file": "defect", "zfs": [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]]} ```

<a id="tag_Qubit_spin"></a>
## Qubit.spin
<b>Description:</b><br>
The spin quantum number of the qubit of [Qubit.file](#tag_Qubit_file). It is also the S a
hyperfine tensor file of format v3 has to name. <br>

<b>Data Type:</b> Double

<b>Default:</b> `1.0`

<b>Example:</b> ``` "Qubit" : {"file": "defect", "spin": 1.0} ```

<a id="tag_Qubit_alphams"></a>
## Qubit.alphams
<b>Description:</b><br>
The m_s value of the qubit level used as |0⟩ (alpha), with [Qubit.file](#tag_Qubit_file). <br>

<b>Data Type:</b> Double

<b>Default:</b> `1.0`

<b>Example:</b> ``` "Qubit" : {"file": "defect", "alphams": -1} ```

<a id="tag_Qubit_betams"></a>
## Qubit.betams
<b>Description:</b><br>
The m_s value of the qubit level used as |1⟩ (beta), with [Qubit.file](#tag_Qubit_file). <br>

<b>Data Type:</b> Double

<b>Default:</b> `0.0`

<b>Example:</b> ``` "Qubit" : {"file": "defect", "betams": 0} ```

<a id="tag_qubit"></a>
## qubit
<b>Description:</b><br>
The qubits, one entry each: name, spin number, gyromagnetic ratio, position, detuning and the m_s of the two levels. The
number of qubits is the length of the list (at least one). <br>

<b>Data Type:</b> Dictionary List

<b>Default:</b> -

<b>Example:</b> 
    ``` "qubit" : [ { "name" : "Q1", "spin" : 0.5, "gyro" : -17608.597050, "xyz" : [1.0, 2.0, 3.0], "detuning" : 0.0, "alphams" : 0.5, "betams" : -0.5 }, 
                    { "name" : "Q2", "spin" : 1.0, "gyro" : -17608.597050, "xyz" : [4.0, 5.0, 6.0], "detuning" : 0.0, "alphams" : 0, "betams" : 1 } ] ```

<b>Sub-tags:</b>

- [name](#tag_name) : Qubit name
    - [spin](#tag_spin) : Spin number
    - [gyro](#tag_gyro) : Gyromagnetic ratio \(rad/ms/G\)
    - [xyz](#tag_xyz)   : x,y,z position of qubit in cartesian, angstrom
    - [detuning](#tag_detuning) : Detuning frequency \(kHz\)
    - [alphams](#tag_alphams)   : qubit sublevel state
    - [betams](#tag_betams)     : qubit sublevel state

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
<b>Data Type:</b> Double

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
Detuning frequency of the qubit. <br>
<b>Data Type:</b> Double (or `[x]`), or `{"value": ..., "unit": ...}` with the unit `Hz`, `kHz`, `MHz` or `GHz` ([Units](InputFile.md#Units))

<b>Unit:</b> kHz

<b>Default:</b> 0.0

<b>Example:</b> ``` "detuning" : 0.0 ```

<b>Parent-tags:</b> [qubit](#tag_qubit)

<a id="tag_alphams"></a>
## alphams
<b>Description:</b><br>
Alpha state of the qubit. <br>
Represents one of the effective spin states of the qubit. <br>
The alpha state |α⟩ is the product of the qubits' alphams unless [alpha_state](#tag_alpha_state) gives it. When
alpha_state is an index that picks a state disagreeing with an alphams written in the input, the run stops. <br>
<b>Data Type:</b> Float

<b>Default:</b> 1.0

<b>Example:</b> ``` "alphams" : 1.0 ```

<b>Related-tags:</b> [betams](#tag_betams) [alpha_state](#tag_alpha_state)

<b>Parent-tags:</b> [qubit](#tag_qubit)

<a id="tag_betams"></a>
## betams
<b>Description:</b><br>
Beta state of the qubit. <br>
One of the effective spin states of the qubit. <br>
The beta state |β⟩ is the product of the qubits' betams unless [beta_state](#tag_beta_state) gives it. When
beta_state is an index that picks a state disagreeing with a betams written in the input, the run stops. <br>
<b>Data Type:</b> Float

<b>Default:</b> 0.0

<b>Example:</b> ``` "betams" : 0.0 ```

<b>Related-tags:</b> [alphams](#tag_alphams) [beta_state](#tag_beta_state)

<b>Parent-tags:</b> [qubit](#tag_qubit)

<a id="tag_intmap"></a>
## intmap
<b>Description:</b><br>
Interaction map between qubits. <br>
Defines the interaction map between qubits, including the self interaction of a qubit \(e.g. zero-field splitting\). <br>
It is optional: a self entry that is not given is the zero tensor, and a pair of different qubits that is not given gets the
point-dipole tensor CCEX computes from their `xyz` and `gyro`. <br>
<b>Data Type:</b> Dictionary List

<b>Unit:</b> Selected independently for each entry by its optional
[unit](#tag_intmap_unit); kHz when omitted.

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
    - [tensor_frame](#tag_tensor_frame) : Basis the tensor components are written in

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
The alpha state |α⟩ of the qubits, given explicitly instead of from the qubits' [alphams](#tag_alphams). An object with exactly one of: <br>

- `"index"`: n, the n-th eigenstate of the qubit Hamiltonian (zero-field splitting, Zeeman, intmap) sorted by
  **descending** energy, starting at 0. That order changes with the field: an NV at 500 G
  along its axis orders m<sub>s</sub> = +1, −1, 0, so `alpha_state {"index": 2}`, `beta_state {"index": 1}` is the
  0 ↔ −1 transition, not +1 ↔ 0. To name a sublevel, write [alphams](#tag_alphams)/[betams](#tag_betams) instead. The run
  log prints which m<sub>s</sub> state |α⟩ and |β⟩ mostly are (the lines `psia is mostly`, `psib is mostly`), and the run stops when an index and an explicitly written
  alphams/betams name different states. An index outside 0 .. (dimension − 1) is refused.
- `"vector"`: the state itself, real numbers over the whole qubit space (kron of the qubits,
  the last one varying fastest), normalized when read.

Both `index` and `vector` in one state, or neither, is refused. Each side is decided on its own: `beta_state` may be
given another way or not at all, and a side that is not given is the product of the qubits' `alphams` (`betams`). <br>
<b>Data Type:</b> Object `{"index": n}` or `{"vector": [...]}`

<b>Default:</b> - (the product of the qubits' alphams)

<b>Example:</b> ``` "alpha_state" : {"index": 0} ```, ``` "alpha_state" : {"vector": [1.0, 0.0, 0.0]} ```

<b>Related tags:</b> [beta_state](#tag_beta_state) [alphams](#tag_alphams) [psi0](#tag_psi0)

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_beta_state"></a>
## beta_state
<b>Description:</b><br>
The beta state |β⟩ of the qubits, as [alpha_state](#tag_alpha_state): `{"index": n}`
or `{"vector": [...]}`, not both; not given, the product of the qubits' [betams](#tag_betams). <br>
<b>Data Type:</b> Object `{"index": n}` or `{"vector": [...]}`

<b>Default:</b> - (the product of the qubits' betams)

<b>Example:</b> ``` "beta_state" : {"vector": [0.0, 1.0, 0.0]} ```

<b>Related tags:</b> [alpha_state](#tag_alpha_state) [betams](#tag_betams)

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_psi0"></a>
## psi0
<b>Description:</b><br>
Initial state of the qubit in matrix representation, symbol:  |psi0⟩ <br>
This tag is used to define a specific initial state of the qubit. <br>
If this tag is not provided, the initial state defaults to the normalized sum |psi0⟩ = (|alpha⟩ + |beta⟩) / ‖|alpha⟩ + |beta⟩‖. <br>
<b>Data Type:</b> Double List

<b>Default:</b> -

<b>Example:</b> ``` "psi0" : [1.0, 0.0, 0.0] ```

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_intmap_unit"></a>
## intmap unit
<b>Description:</b><br>
Frequency unit of the [tensor](#tag_tensor) in this Qubit.intmap entry. Each entry
may choose its own unit. Supported values are <tt>"Hz"</tt>, <tt>"kHz"</tt>,
<tt>"MHz"</tt>, and <tt>"GHz"</tt>, matched case-insensitively. The tensor is converted
to the internal rad*kHz representation while parsing. <br>
Without this key the tensor numbers are kHz. The code, not the input author, supplies the factor of 2π. <br>
This key belongs to Qubit.intmap entries only. [Qubit.zfs](#tag_Qubit_zfs), read with the single-qubit
[Qubit.file](#tag_Qubit_file) form, has no unit sub-key and is always kHz. <br>

<b>Data Type:</b> String

<b>Default:</b> <tt>"kHz"</tt>

<b>Example:</b> ``` "unit" : "MHz" ```

<b>Parent-tags:</b> [intmap](#tag_intmap)

<a id="tag_tensor_frame"></a>
## tensor_frame
<b>Description:</b><br>
Which Cartesian basis the components of this [intmap](#tag_intmap) tensor are written
in. Only meaningful together with
[frame.rotation](InputFile_General.md#tag_coordinate_frame_rotation). <br>

- <tt>"bath"</tt>  : the original bath Cartesian basis. The tensor is transformed as
  T → R T R<sup>T</sup> along with the positions.
- <tt>"qubit"</tt> : the common computational basis, the one in which the reference
  qubit's axis is +z. The tensor is left exactly as written.

<b><tt>"qubit"</tt> does not mean an individual qubit's local frame.</b> There is one
computational frame for the whole run, fixed by `frame.rotation.reference_qubit`. A qubit whose
crystallographic orientation differs from the reference one is described by ITS OWN ZFS
tensor expressed in that common frame; with <tt>"bath"</tt> the same global R
transforms it, and the off-diagonal elements that appear are correct and must not be
diagonalized away. <br>

The rule is the same for a self/ZFS entry (<tt>between: ["NV0","NV0"]</tt>) and for a
pair entry (<tt>between: ["NV0","NV1"]</tt>): the frame is taken from this key, never
guessed from which kind of tensor it is. <br>
With the rotation on and more than one qubit, every explicit intmap entry MUST carry it;
omitting it is an error, not a default. With a single qubit the key is optional, and without
it the tensor is left as written. <br>
An entry that is not given at all is either the zero tensor or, off the diagonal, the
point-dipole tensor CCEX computes from the source qubit geometry; the latter is always
transformed. <br>

<b>Data Type:</b> String

<b>Default:</b> absent. Fatal with the rotation on and more than one qubit; otherwise the
tensor is left as written.

<b>Example:</b> ``` tensor_frame : "bath" ```

<b>Parent-tags:</b> [intmap](#tag_intmap)
