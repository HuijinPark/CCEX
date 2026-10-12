# Qubit-related tags
This tags are used for adjusting Qubit Hamiltonian and qubit inital state.
The format is the following:

<b>Case1.</b> [qubitfile](InputFile_General.md#tag_qubitfile) tag <br>
  : Only single-qubit systems are supported. <br>
    Note that this tag takes priority over the [Qubit](#tag_Qubit) tag. <br>

<b>Example:</b> ``` "qubitfile"  : "file-path" ```

<b>Case2.</b> [Qubit](#tag_Qubit) tag <br>
  : Both multiple qubit and single system are supported. <br>

<details>
<summary><b>Example:</b> Click here to show example code</summary>
```
{
    "Qubit" : {
                "nqubit" : 2,
                "qubit" : [
                            {
                                "name" : "Q1",
                                "spin" : 0.5,
                                "gyro" : -17608.597050,
                                "xyz"  : [1.0, 2.0, 3.0],
                                "detuning" : 0.0,
                                "alphams" : 0,
                                "betams"  : 1,
                            },
                            {
                                "name" : "Q2",
                                "spin" : 1.0,
                                "gyro" : -17608.597050,
                                "xyz"  : [4.0, 5.0, 6.0],
                                "detuning" : 0.0,
                                "alphams" : 0,
                                "betams"  : 1,
                            }
                ],
                "intmap" : [
                            {
                                # Zero field splitting
                                "between" : ["Q1", "Q1"], 
                                "tensor"  : [[1.0, 2.0, 3.0],
                                            [4.0, 5.0, 6.0], 
                                            [7.0, 8.0, 9.0]],
                            },
                            {
                                "between" : ["Q1", "Q2"],
                                "tensor"  : [[1.0, 2.0, 3.0],
                                            [4.0, 5.0, 6.0], 
                                            [7.0, 8.0, 9.0]],
                            }
                ]
              },
    "psia" : [1.0, 2.0, 3.0],
    "psib" : [4.0, 5.0, 6.0],
    "psi0" : [7.0, 8.0, 9.0],
    "overhaus" : true,
}
```
</details>

- [qubitfile](InputFile_General.md#tag_qubitfile)
- [qzfs](#tag_qzfs)
- [qspin](#tag_qspin)
- [qalphams](#tag_qalphams) 
- [qbetams](#tag_qbetams)
- [Qubit](#tag_Qubit) 
- [nqubit](#tag_nqubit)
- [qubit](#tag_qubit)
- [name](#tag_name)         
- [spin](#tag_spin)
- [gyro](#tag_gyro)         
- [xyz](#tag_xyz)          
- [detuning](#tag_detuning)     
- [alphams](#tag_alphams) 
- [betams](#tag_betams)
- [intmap](#tag_intmap)
- [between](#tag_between)      
- [tensor](#tag_tensor)       
- [psia](#tag_psia)  
- [psib](#tag_psib)  
- [psi0](#tag_psi0)
- [overhaus](#tag_overhaus)     
- [alphaidx](#tag_alphaidx) 
- [betaidx](#tag_betaidx)

<a id="tag_Qubit"></a>
## Qubit
<b>Description:</b><br>
Qubit information. <br>
This tag is required to give the details of qubit information. <br>

<b>Data Type:</b> Dictionary

<b>Default:</b> -

<b>Example:</b> 
    ``` Qubit : { "nqubit" : 2, 
                    "qubit" : [ { "name" : "Q1", "spin" : 0.5, "gyro" : -17608.597050, "xyz" : [1.0, 2.0, 3.0], "detuning" : 0.0, "alphams" : 0, "betams" : 1 }, 
                                 { "name" : "Q2", "spin" : 1.0, "gyro" : -17608.597050, "xyz" : [4.0, 5.0, 6.0], "detuning" : 0.0, "alphams" : 0, "betams" : 1 } ] } ```

<Related variable:</b> [QubitArray struct](#QubitArray), [Qubit struct](#Qubit)

<b>Sub-tags:</b>

- [nqubit](#tag_nqubit) 
    - [qubit](#tag_qubit) 
    - [intmap](#tag_intmap) 
    - [psia](#tag_psia)
    - [psib](#tag_psib)
    - [psi0](#tag_psi0)
    - [overhaus](#tag_overhaus)
    - [alphaidx](#tag_alphaidx)
    - [betaidx](#tag_betaidx)

<a id="tag_qzfs"></a>
## qzfs
<b>Description:</b><br>
The zero-field-splitting tensor of the qubit, 3 × 3. Only on the single-qubit path that [qubitfile](InputFile_General.md#tag_qubitfile) opens; the `Qubit` block has its own keys for the same quantities. It has no `unit` or
`tensor_frame` sub-key: the values are kHz and are read as written, never rotated by
[coordinate_frame_rotation](InputFile_General.md#tag_coordinate_frame_rotation). <br>

<b>Data Type:</b> Double List (3 × 3)

<b>Unit:</b> kHz

<b>Default:</b> zero tensor

<b>Example:</b> ``` "qzfs" : [[-960000, 0, 0], [0, -960000, 0], [0, 0, 1920000]] ```

<b>Related variable:</b> [QubitArray](#QubitArray)::intmap

<a id="tag_qspin"></a>
## qspin
<b>Description:</b><br>
The spin quantum number of the qubit. Only on the single-qubit path that [qubitfile](InputFile_General.md#tag_qubitfile) opens; the `Qubit` block has its own keys for the same quantities. <br>

<b>Data Type:</b> Double

<b>Default:</b> `1.0`

<b>Example:</b> ``` "qspin" : 1.0 ```

<b>Related variable:</b> [Qubit](#Qubit)::spin

<a id="tag_qalphams"></a>
## qalphams
<b>Description:</b><br>
The m_s value of the qubit level used as |0⟩ (alpha). Only on the single-qubit path that [qubitfile](InputFile_General.md#tag_qubitfile) opens; the `Qubit` block has its own keys for the same quantities. <br>

<b>Data Type:</b> Double

<b>Default:</b> `1.0`

<b>Example:</b> ``` "qalphams" : -1 ```

<b>Related variable:</b> [Qubit](#Qubit)::alpha

<a id="tag_qbetams"></a>
## qbetams
<b>Description:</b><br>
The m_s value of the qubit level used as |1⟩ (beta). Only on the single-qubit path that [qubitfile](InputFile_General.md#tag_qubitfile) opens; the `Qubit` block has its own keys for the same quantities. <br>

<b>Data Type:</b> Double

<b>Default:</b> `0.0`

<b>Example:</b> ``` "qbetams" : 0 ```

<b>Related variable:</b> [Qubit](#Qubit)::beta

<a id="tag_nqubit"></a>
## nqubit
<b>Description:</b><br>
Number of qubits in the system. <br>
The number of qubit information in tag [Qubit](#tag_Qubit) should be the same as this value. <br>
<b>Data Type:</b> Integer

<b>Default:</b> -

<b>Example:</b> ``` nqubit : 2 ```

<b>Related variable:</b> [QubitArray](#QubitArray)::nqubit

<a id="tag_qubit"></a>
## qubit
<b>Description:</b><br>
Qubit information details, e.g. qubit name, spin number, gyro magnetic ratio, etc.. <br>
This tag is required to give the details of qubit information. <br>

<b>Data Type:</b> Dictionary List

<b>Default:</b> -

<b>Example:</b> 
    ``` qubit : [ { "name" : "Q1", "spin" : 0.5, "gyro" : -17608.597050, "xyz" : [1.0, 2.0, 3.0], "detuning" : 0.0, "alphams" : 0, "betams" : 1 }, 
                    { "name" : "Q2", "spin" : 1.0, "gyro" : -17608.597050, "xyz" : [4.0, 5.0, 6.0], "detuning" : 0.0, "alphams" : 0, "betams" : 1 } ] ```

<b>Related variable:</b> [Qubit struct](#Qubit)

<b>Sub-tags:</b>

- [name](#tag_name) : Bath spin name
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

<b>Default:</b> "q%d" e.g. "q1", "q2", "q3", ...

<b>Example:</b> ``` name : "Q1" ```

<b>Related variable:</b> [Qubit](#Qubit)::name

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_spin"></a>
## spin
<b>Description:</b><br>
Spin number of the qubit. <br>
<b>Data Type:</b> Float

<b>Default:</b> 1.0

<b>Example:</b> ``` spin : 0.5 ```

<b>Related variable:</b> [Qubit](#Qubit)::spin

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_gyro"></a>
## gyro
<b>Description:</b><br>
Gyromagnetic ratio of the qubit. <br>
<b>Data Type:</b> Double

<b>Unit:</b> rad/ms/G

<b>Default:</b> -17608.597050

<b>Example:</b> ``` gyro : -17608.597050 ```

<b>Related variable:</b> [Qubit](#Qubit)::gyro

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_xyz"></a>
## xyz
<b>Description:</b><br>
x,y,z position of the qubit in cartesian coordinates. <br>
<b>Data Type:</b> Double List

<b>Unit:</b> angstrom

<b>Default:</b> -

<b>Example:</b> ``` xyz : [1.0, 2.0, 3.0] ```

<b>Related variable:</b> [Qubit](#Qubit)::xyz

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_detuning"></a>
## detuning
<b>Description:</b><br>
Detuning frequency of the qubit. <br>
<b>Data Type:</b> Double

<b>Unit:</b> kHz

<b>Default:</b> 0.0

<b>Example:</b> ``` detuning : 0.0 ```

<b>Related variable:</b> [Qubit](#Qubit)::detuning

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_alphams"></a>
## alphams
<b>Description:</b><br>
Alpha state of the qubit. <br>
Represents one of the effective spin states of the qubit. <br>
Priority: [psia](#tag_psia) > [alphaidx](#tag_alphaidx) > [alphams](#tag_alphams). <br>
This tag is ultimately applied to set [QubitArray](#QubitArray)::psia \(i.e. [psia](#tag_psia)\), unless it is overridden by a higher-priority tag <br>
<b>Data Type:</b> Float

<b>Default:</b> 1.0

<b>Example:</b> ``` alphams : 1.0 ```

<b>Related-tags:</b> [betams](#tag_betams) [psia](#tag_psia) [alphaidx](#tag_alphaidx)

<b>Related variable:</b> [Qubit](#Qubit)::alphams [QubitArray](#QubitArray)::psia

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_betams"></a>
## betams
<b>Description:</b><br>
Beta state of the qubit. <br>
One of the effective spin states of the qubit. <br>
Priority: [psib](#tag_psib) > [betaidx](#tag_betaidx) > [betams](#tag_betams). <br>
This tag is ultimately applied to set [QubitArray](#QubitArray)::psib \(i.e. [psib](#tag_psib)\), unless it is overridden by a higher-priority tag <br>
<b>Data Type:</b> Float

<b>Default:</b> 0.0

<b>Example:</b> ``` betams : 0.0 ```

<b>Related-tags:</b> [alphams](#tag_alphams) [psib](#tag_psib) [betaidx](#tag_betaidx)

<b>Related variable:</b> [Qubit](#Qubit)::betams [QubitArray](#QubitArray)::psib

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_intmap"></a>
## intmap
<b>Description:</b><br>
Interaction map between qubits. <br>
This tag is required to define the interaction map between qubits, including its self interaction \(e.g. Zero-field splitting\). <br>
<b>Data Type:</b> Dictionary List

<b>Unit:</b> Selected independently for each entry by its optional
[unit](#tag_intmap_unit); kHz when omitted.

<b>Default:</b> -

<b>Example:</b> 
    ``` intmap : [ { "between" : ["q1","q1"],
                       "unit"    : "MHz",
                       "tensor"  : [[0.0, 0.0, 0.0],
                                    [0.0, 0.0, 0.0], 
                                    [0.0, 0.0, 0.0]] }, 
                     { ... } ] 
    ```

<b>Related variable:</b> [QubitArray](#QubitArray)::intmap

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

<b>Example:</b> ``` between : ["q1","q2"] ```

<b>Related tags:</b> [tensor](#tag_tensor)

<b>Related variable:</b> [QubitArray](#QubitArray)::intmap

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

<b>Example:</b> ``` tensor : [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]] ```

<b>Related tags:</b> [between](#tag_between)

<b>Related variable:</b> [QubitArray](#QubitArray)::intmap

<b>Parent-tags:</b> [intmap](#tag_intmap)

<a id="tag_psia"></a>
## psia
<b>Description:</b><br>
Alpha state of the qubit in matrix representation.<br>
One of the effective spin states of the qubit in matrix representation, simbol : |alpha⟩ . <br>
Priority: [psia](#tag_psia) > [alphaidx](#tag_alphaidx) > [alphams](#tag_alphams). <br>
<b>Data Type:</b> Double List

<b>Default:</b> - \(If not entered, alphams state would be the \)

<b>Example:</b> ``` psia : [1.0, 0.0, 0.0] ```

<b>Related tags:</b> [psib](#tag_psib) [alphams](#tag_alphams) [alphaidx](#tag_alphaidx)

<b>Related variable:</b> [QubitArray](#QubitArray)::psia

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_psib"></a>
## psib
<b>Description:</b><br>
Beta state of the qubit in matrix representation. <br>
One of the effective spin states of the qubit in matrix representation, simbol : |beta⟩ . <br>
Priority: [psib](#tag_psib) > [betaidx](#tag_betaidx) > [betams](#tag_betams). <br>
<b>Data Type:</b> Double List

<b>Default:</b> -

<b>Example:</b> ``` psib : [0.0, 1.0, 0.0] ```

<b>Related tags:</b> [psia](#tag_psia) [betams](#tag_betams) [betaidx](#tag_betaidx)

<b>Related variable:</b> [QubitArray](#QubitArray)::psib

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_psi0"></a>
## psi0
<b>Description:</b><br>
Initial state of the qubit in matrix representation, symbol:  |psi0⟩ <br>
This tag is used to define a specific initial state of the qubit. <br>
If this tag is not provided, the initial state defaults to  |psi0⟩ =  |alpha⟩ +  |beta⟩ . <br>
<b>Data Type:</b> Double List

<b>Default:</b> -

<b>Example:</b> ``` psi0 : [1.0, 0.0, 0.0] ```

<b>Related variable:</b> [QubitArray](#QubitArray)::psi0

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_overhaus"></a>
## overhaus
<b>Description:</b><br>
Overhauser field of the qubit. <br>
This tag is required to compute spin dynamics including the Overhauser field of the qubit. <br>
<b>Data Type:</b> Bool

<b>Default:</b> false

<b>Example:</b> ``` overhaus : true ```

<b>Related variable:</b> [QubitArray](#QubitArray)::overhaus

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_alphaidx"></a>
## alphaidx
<b>Description:</b><br>
Alpha state of the qubit. (Adiabatic state of the qubit system) <br>
Use this tag to assign the eigenvector of the diagonalized qubit Hamiltonian as the initial qubit state. <br>
Priority: [psia](#tag_psia) > [alphaidx](#tag_alphaidx) > [alphams](#tag_alphams). <br>
The index counts the eigenstates of the qubit Hamiltonian (zero-field splitting, Zeeman, intmap) sorted by **descending** energy, starting at 0. That order changes with the field: an NV at 500 G along its axis orders m<sub>s</sub> = +1, −1, 0, so `alphaidx: 2, betaidx: 1` is the 0 ↔ −1 transition, not +1 ↔ 0. To name a sublevel, write [alphams](#tag_alphams)/[betams](#tag_betams) instead. The run log prints which m<sub>s</sub> state psia and psib mostly are, and the run stops when alphaidx/betaidx and an explicitly written alphams/betams name different states (since 2026-10-03; until then the indices silently won). <br>
This tag is ultimately applied to set [QubitArray](#QubitArray)::psia \(i.e. [psia](#tag_psia)\), unless it is overridden by a higher-priority tag <br>
<b>Data Type:</b> Int

<b>Default:</b> -1 \(= None\)

<b>Example:</b> ``` alphaidx : 0 ```

<b>Related tags:</b> [betaidx](#tag_betaidx) [psia](#tag_psia) [alphams](#tag_alphams)

<b>Related variable:</b> [QubitArray](#QubitArray)::_alphaidx [QubitArray](#QubitArray)::psia

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_betaidx"></a>
## betaidx
<b>Description:</b><br>
Beta state of the qubit. (Adiabatic state of the qubit system) <br>
Use this tag to assign the eigenvector of the diagonalized qubit Hamiltonian as the initial qubit state. <br>
Priority: [psib](#tag_psib) > [betaidx](#tag_betaidx) > [betams](#tag_betams). <br>
The index counts the eigenstates of the qubit Hamiltonian (zero-field splitting, Zeeman, intmap) sorted by **descending** energy, starting at 0. That order changes with the field: an NV at 500 G along its axis orders m<sub>s</sub> = +1, −1, 0, so `alphaidx: 2, betaidx: 1` is the 0 ↔ −1 transition, not +1 ↔ 0. To name a sublevel, write [alphams](#tag_alphams)/[betams](#tag_betams) instead. The run log prints which m<sub>s</sub> state psia and psib mostly are, and the run stops when alphaidx/betaidx and an explicitly written alphams/betams name different states (since 2026-10-03; until then the indices silently won). <br>
This tag is ultimately applied to set [QubitArray](#QubitArray)::psib \(i.e. [psib](#tag_psib)\), unless it is overridden by a higher-priority tag <br>
<b>Data Type:</b> Int

<b>Default:</b> -1 \(= None\)

<b>Example:</b> ``` betaidx : 0 ```

<b>Related tags:</b> [alphams](#tag_alphams) [psib](#tag_psib) [betams](#tag_betams)

<b>Related variable:</b> [QubitArray](#QubitArray)::_betaidx [QubitArray](#QubitArray)::psib

<b>Parent-tags:</b> [Qubit](#tag_Qubit)

<a id="tag_intmap_unit"></a>
## intmap unit
<b>Description:</b><br>
Frequency unit of the [tensor](#tag_tensor) in this Qubit.intmap entry. Each entry
may choose its own unit. Supported values are <tt>"Hz"</tt>, <tt>"kHz"</tt>,
<tt>"MHz"</tt>, and <tt>"GHz"</tt>, matched case-insensitively. The tensor is converted
to the internal rad*kHz representation while parsing. <br>
Omitting this key preserves the legacy behavior: the tensor numbers are interpreted as
kHz. The code, not the input author, supplies the factor of $2\pi$. <br>
This key belongs to Qubit.intmap entries only. The legacy top-level [qzfs](#tag_qzfs), which is read by the single-qubit [qubitfile](InputFile_General.md#tag_qubitfile) path, has no
unit sub-key and is always kHz. <br>

<b>Data Type:</b> String

<b>Default:</b> <tt>"kHz"</tt>

<b>Example:</b> ``` "unit" : "MHz" ```

<b>Parent-tags:</b> [intmap](#tag_intmap)

<a id="tag_tensor_frame"></a>
## tensor_frame
<b>Description:</b><br>
Which Cartesian basis the components of this [intmap](#tag_intmap) tensor are written
in. Only meaningful together with
[coordinate_frame_rotation](InputFile_General.md#tag_coordinate_frame_rotation). <br>

- <tt>"bath"</tt>  : the original bath Cartesian basis. The tensor is transformed as
      $ T \rightarrow R\,T\,R^{T} $ along with the positions.
    - <tt>"qubit"</tt> : the common computational basis, the one in which the reference
      qubit's axis is +z. The tensor is left exactly as written.

<b><tt>"qubit"</tt> does not mean an individual qubit's local frame.</b> There is one
computational frame for the whole run, fixed by `reference_qubit.` A qubit whose
crystallographic orientation differs from the reference one is described by ITS OWN ZFS
tensor expressed in that common frame; with <tt>"bath"</tt> the same global $R$
transforms it, and the off-diagonal elements that appear are correct and must not be
diagonalized away. <br>

The rule is the same for a self/ZFS entry (<tt>between: ["NV0","NV0"]</tt>) and for a
pair entry (<tt>between: ["NV0","NV1"]</tt>): the frame is taken from this key, never
guessed from which kind of tensor it is. <br>
With the rotation on and more than one qubit, every explicit intmap entry MUST carry it;
omitting it is an error, not a default. With a single qubit the key stays optional and an
input without it behaves exactly as before -- the tensor is left as written. <br>
An entry that is not given at all is either the zero tensor or, off the diagonal, the
point-dipole tensor CCEX computes from the source qubit geometry; the latter is always
transformed. <br>

<b>Data Type:</b> String

<b>Default:</b> absent. Fatal with the rotation on and nqubit > 1; otherwise the
tensor is left as written.

<b>Example:</b> ``` tensor_frame : "bath" ```

<b>Parent-tags:</b> [intmap](#tag_intmap)
