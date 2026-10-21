# Defect tags
The `Defect` block describes the paramagnetic defects of the bath: the list of defect species, the model (old key names are in the table of
[old names](InputFile.md#OldNames)). The files of their axes and states are bath-configuration keys:
[configuration.defect_axis_file](InputFile_General.md#tag_defect_axis_file) and
[configuration.defect_state_file](InputFile_General.md#tag_defect_state_file).

<a id="DefectModel"></a>
**The bath model.** A bath spin whose name in the bath file is a species' `name` is a paramagnetic defect: an electron
spin (its spin and γ from the gyro file) with a set of indexed *axis states*, `1 .. n_axis_states`. An axis state stands
for what the electron's environment can be: a principal (Jahn-Teller) axis of a P1 centre, or a frozen state of its
nucleus. Each defect spin is given one axis state, from
[configuration.defect_axis_file](InputFile_General.md#tag_defect_axis_file) or at random from `configuration.seed`, and
its tensors are those of that state. A nucleus on the defect site enters in one of two ways:
- through the electron's [detuning](#tag_df_detuning) only: the hyperfine coupling taken as a static shift
  A·m<sub>I</sub> of the electron (one axis state per axis and m<sub>I</sub>, as the twelve detunings of a P1 with
  <sup>14</sup>N), with `n_onsite_spins` 0;
- as an on-site spin of its own (`n_onsite_spins` > 0, with [types](#tag_df_types), [spins](#tag_df_spins),
  [gyros](#tag_df_gyros) and per axis state its [position](#tag_df_nuclear_positions), [hyperfine](#tag_df_hyperfine)
  and [efg](#tag_df_efg)), which joins its electron in every cluster the defect is in.

With [configuration.s_range](InputFile_General.md#tag_s_range) the bath is sampled in product states: the
electrons' states from `bath_state_file` or `seed`, the on-site spins' from `defect_state_file` or `seed`, and the static
field of the spins outside a cluster enters as a mean field
([bath_bath_mean_field](InputFile_General.md#tag_bath_bath_mean_field)). The Hamiltonian and the approximations are in
[CCE methods](../methods/cce.md); `example/tutorial/03_P1/` runs the P1 bath both ways.

<a id="tag_df_forms"></a>
**Two forms of the values**: a key written as a list `[[axis state, name, value], ...]` holds a value per axis
state; the object form holds one value for every axis state (`zfs`: `{"D", "axis", "unit"}` or `{"tensor", "unit"}`) or
wraps the same list with a unit (`hyperfine`, `efg`, `detuning`, `gyros`, `quadrupole_moments`: `{"values": [...],
"unit": ...}`, where `"value"` is taken too, the word the other keys use; both in one object are refused). A unit is
written with its case, as in the tables below (`"mhz"` is refused with the table's `"MHz"`; it was compared without
regard to case). A plain frequency here is MHz, where a qubit's ([Qubit.qubit[].detuning](InputFile_Qubit.md#tag_detuning),
[Qubit.intmap](InputFile_Qubit.md#tag_intmap)) is kHz.

<details><summary><b>Example:</b> Click here to show example code</summary>
```
{
"configuration"  : {
    "bathfile"          : ["./bath_DiaP1_1ppm_1"],
    "defect_axis_file"  : ["./bathJT_DiaP1_1ppm_1"],
    "defect_state_file" : ["./stateEx_DiaP1_1ppm_"]
},
"Defect"         : {
    "species"    : [
		{
			"name"	: "NV0",
			"n_onsite_spins"  	: 1,
			"types" 	: ["15N"],
			"spins" 	: [0.5],   
			"gyros" 	: [-2.71261804],
			"n_axis_states"    	: 4,
			"nuclear_positions"      : [ [1, "15N", [  1.36404 ,   -0.78843 ,   -0.56598 ]],
		                     [2, "15N", [ -1.36801 ,   -0.79025 ,   -0.56701 ]],
		                     [3, "15N", [  0.00000 ,    1.57603 ,   -0.56511 ]],
		                     [4, "15N", [  0.00000 ,   -0.00733 ,    1.67433 ]]],
			"hyperfine" 		: [  [1, "15N", [ -10.33219 ,    1.01563 ,    1.50162 ,    1.01563 ,   -5.18317 ,   -0.16410 ,    1.50162 ,   -0.16410 ,   -6.50270]],  
			  				 [2, "15N", [  -7.34999 ,   -2.73741 ,   -0.89293 ,   -2.73741 ,   -8.16537 ,   -1.21839 ,   -0.89293 ,   -1.21839 ,   -6.50270]], 
			  				 [3, "15N", [  -5.59086 ,    1.72177 ,   -0.60870 ,    1.72177 ,   -9.92450 ,    1.38249 ,   -0.60870 ,    1.38249 ,   -6.50270]],
			  				 [4, "15N", [  -6.01395 ,    0.00000 ,    0.00000 ,    0.00000 ,   -5.59086 ,   -1.82620 ,    0.00000 ,   -1.82620 ,  -10.41325]]]
		}
                       ]
}
}
```
</details>

- [Defect](#tag_Defect)       
- [species](#tag_df_species): [name](#tag_df_name), [n_onsite_spins](#tag_df_n_onsite_spins), [n_axis_states](#tag_df_n_axis_states),
  [frame](#tag_df_frame), [types](#tag_df_types), [spins](#tag_df_spins),
  [gyros](#tag_df_gyros), [quadrupole_moments](#tag_df_quadrupole_moments), [nuclear_positions](#tag_df_nuclear_positions),
  [hyperfine](#tag_df_hyperfine), [efg](#tag_df_efg), [zfs](#tag_df_zfs), [detuning](#tag_df_detuning)

<a id="tag_Defect"></a>
## Defect
<b>Description:</b><br>
Detailed information of the paramagnetic \(PM\) defect bath: `species`, the list of defect species. <br>
<b>Data Type:</b> Dictionary `{"species": [...]}`

<b>Default:</b> -

<b>Example:</b>
    ```
    "Defect" : {"species" : [{
        "name" : "NV_0", "n_onsite_spins" : 0, "n_axis_states" : 3,
        "frame" : "file",
        "zfs" : {"D" : 2870.0, "axis" : [1.0, 1.0, 1.0], "unit" : "MHz"},
        "detuning" : {
            "values" : [[1,"e",-2.14], [2,"e",0.0], [3,"e",2.14]],
            "unit" : "MHz"
        }
    }]}
    ```

<a id="tag_df_species"></a>
## Defect.species
<b>Description:</b><br>
The defect species, one entry each. A bath spin whose name is a species' `name` is that defect. <br>
<b>Data Type:</b> Dictionary List

<b>Default:</b> no `Defect` block: no bath spin is a defect (every bath spin is a plain spin with its gyro-file entry).

<b>Sub-tags:</b>

- [name](#tag_df_name) : Defect name
    - [n_onsite_spins](#tag_df_n_onsite_spins) : Number of on-site nuclear spins of the PM defect
    - [n_axis_states](#tag_df_n_axis_states) : Number of indexed axis states of the defect
    - [frame](#tag_df_frame) : Frame of the Defect vectors and tensors
    - [types](#tag_df_types) : On-site nuclear spin types
    - [spins](#tag_df_spins) : Spin numbers of the on-site nuclear spins
    - [gyros](#tag_df_gyros) : Gyromagnetic ratio of the on-site nuclear spins
    - [quadrupole_moments](#tag_df_quadrupole_moments) : eQ values of the on-site nuclear spins
    - [nuclear_positions](#tag_df_nuclear_positions) : Relative x,y,z position of on-site nuclear spins
    - [hyperfine](#tag_df_hyperfine) : Hyperfine interaction tensor
    - [efg](#tag_df_efg) : Electric field gradient tensor
    - [zfs](#tag_df_zfs) : Zero-field splitting tensor
    - [detuning](#tag_df_detuning) : Detuning frequency

<a id="tag_df_name"></a>
## name
<b>Description:</b><br>
Name of the paramagnetic defect. <br>
<b>Data Type:</b> String

<b>Default:</b> -

<b>Example:</b> ``` "name" : "P1" ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_n_onsite_spins"></a>
## n_onsite_spins
<b>Description:</b><br>
Number of on-site (nuclear) spins of the paramagnetic defect. <br>
Each on-site spin (with its [types](#tag_df_types), [spins](#tag_df_spins), [gyros](#tag_df_gyros) and the tensors per
axis state) joins its electron in every cluster the defect is in, so a cluster of one P1 is its electron and its
<sup>14</sup>N. With 0 the nucleus enters only through the electron's [detuning](#tag_df_detuning); write the on-site
spins when its own dynamics matters (`example/tutorial/03_P1/cCCE2_Hahn_onsite14N/` against `cCCE2_Hahn_states/`). <br>
<b>Data Type:</b> Int

<b>Default:</b> 0

<b>Example:</b> ``` "n_onsite_spins" : 0 ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_n_axis_states"></a>
## n_axis_states
<b>Description:</b><br>
Number of indexed axis states available for the paramagnetic defect: an index can be a principal axis, or a
frozen nuclear state (for example three m<sub>I</sub> detunings). A shared object-form zfs is broadcast to indices 1
through n_axis_states. <br>

The indexed entries of [nuclear_positions](#tag_df_nuclear_positions), [hyperfine](#tag_df_hyperfine),
[efg](#tag_df_efg), zfs and [detuning](#tag_df_detuning) must address an
axis state in <b>1 .. n_axis_states</b>; index 0 is reserved for the main spin and is
written by CCEX itself. An index above n_axis_states is refused. <br>

<b>Data Type:</b> Int

<b>Default:</b> 0

<b>Example:</b> ``` "n_axis_states" : 4 ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_axis"></a>
## axis (zfs.axis)
<b>Description:</b><br>
The symmetry axis of the D form of [zfs](#tag_df_zfs), written inside it: `"zfs": {"D": ..., "axis": [x, y, z], "unit":
...}`, in the species' [frame](#tag_df_frame), normalized when read. Before 2.0.0 it was a key of the species,
`Defect.species[].axis`, which only the D form used; that key is refused with the new place (the converter moves it into
a D-form zfs, and drops it elsewhere). <br>

<b>Data Type:</b> Double List(3)

<b>Default:</b> - (required with `D`)

<b>Example:</b> ``` "zfs" : {"D": 2870.0, "axis": [1.0, 1.0, 1.0], "unit": "MHz"} ```

<b>Parent-tags:</b> [zfs](#tag_df_zfs)

<a id="tag_df_frame"></a>
## Defect frame
<b>Description:</b><br>
The frame this species' nuclear_positions displacement vectors, hyperfine / efg / zfs tensor components and zfs axis
are written in (`coordinate_frame` before 2.0.0, with `"bath"` for `"file"`). Read with
[frame.rotation](InputFile_General.md#tag_frame_rotation) on, where a species that writes `nuclear_positions`,
`hyperfine`, `efg` or `zfs` needs it (before 2.0.0 every species did, also a P1 of detunings alone, which has no
direction; before 2.0.0 an absent one meant `"qubit"`); otherwise, and with the rotation off, it is not used (reported).
Detuning is a scalar frequency and is unaffected. <br>

| Option | Description |
|---|---|
| `"file"` | the bath file's frame (e.g. the crystal frame): turned with the positions |
| `"qubit"` | the computational frame (z = `qubit_axis`): used as written |

<b>Data Type:</b> String

<b>Default:</b> none (required with the rotation on)

<b>Example:</b> ``` "frame" : "file" ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_types"></a>
## types
<b>Description:</b><br>
Name of on-site nuclear spins of the paramagnetic defect. <br>
The number of types should be the same as [n_onsite_spins](#tag_df_n_onsite_spins). <br>
<b>Data Type:</b> String List

<b>Default:</b> -

<b>Example:</b> ``` "types" : ["1H", "13C"] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_spins"></a>
## spins
<b>Description:</b><br>
Spin numbers of the on-site nuclear spins of the paramagnetic defect. <br>
The number of spins should be the same as [n_onsite_spins](#tag_df_n_onsite_spins). <br>
<b>Data Type:</b> Float List

<b>Default:</b> -

<b>Example:</b> ``` "spins" : [0.5, 1.0] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_gyros"></a>
## gyros
<b>Description:</b><br>
Gyromagnetic ratio of the on-site nuclear spins of the paramagnetic defect. <br>
The number of gyros should be the same as [n_onsite_spins](#tag_df_n_onsite_spins). <br>
The unit-aware form is
``` "gyros" : {"values" : [3.0777051853], "unit" : "MHz/T"} ```
Plain frequency units such as MHz/T and kHz/G mean γ/2π; units beginning
with <tt>rad</tt> are angular gyromagnetic ratios. Values are normalized internally to
rad kHz/G, numerically the same as rad/ms/G. The original bare list remains supported
in that internal unit. <br>
<b>Data Type:</b> Object {values (or value), unit}, or legacy Double List

<b>Units:</b> radkHz/G (default), rad/ms/G, kHz/G, MHz/G, MHz/T, Hz/T,
radMHz/T, rad/us/T, or rad/s/T

<b>Default:</b> -

<b>Example:</b> ``` "gyros" : {"values" : [3.0777051853], "unit" : "MHz/T"} ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_quadrupole_moments"></a>
## quadrupole_moments
<b>Description:</b><br>
Quadrupole moment \(eQ\) values of the on-site nuclear spins of the paramagnetic defect. <br>
The number of eQ values should be the same as [n_onsite_spins](#tag_df_n_onsite_spins). <br>
The stored quantity follows the existing area convention: one internal unit is
10<sup>−30</sup> m² = 1 fm² = 10 mbarn. The unit-aware form is
``` "quadrupole_moments" : {"values" : [20.44], "unit" : "mbarn"} ```
and is normalized to that internal unit. The original bare list remains supported in
10<sup>−30</sup> m². <br>
<b>Data Type:</b> Object {values (or value), unit}, or legacy Double List

<b>Units:</b> 1e-30 m^2 (default), fm^2, m^2, cm^2, barn, or mbarn

<b>Default:</b> [0.0 , ..]

<b>Example:</b> ``` "quadrupole_moments" : {"values" : [20.44], "unit" : "mbarn"} ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_nuclear_positions"></a>
## nuclear_positions
<b>Description:</b><br>
Relative x,y,z position of on-site nuclear spins of the paramagnetic defect over the electron spin position. <br>
The number of `nuclear_positions` should be the same as [n_onsite_spins](#tag_df_n_onsite_spins) × [n_axis_states](#tag_df_n_axis_states) <br>
<b>Data Type:</b> 2D List

<b>Unit:</b> angstrom

<b>Default:</b> [0.0, 0.0, 0.0]

<b>Example:</b> ``` "nuclear_positions" : [[1, "15N", [0.0, 0.0, 0.0]],
                                  [2, "1H", [0.0, 0.0, 0.0]]] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_hyperfine"></a>
## hyperfine
<b>Description:</b><br>
The hyperfine interaction between on-site nuclear spins and the electron spin associated with a paramagnetic defect. <br>
The number of `hyperfine` should be the same as [n_onsite_spins](#tag_df_n_onsite_spins) × [n_axis_states](#tag_df_n_axis_states) <br>
The unit-aware form wraps the same indexed entries in <tt>values</tt>. Each tensor may
be a flat row-major 9-number array or a nested 3x3 array:
```
"hyperfine" : {
    "values" : [[1, "14N", [[81.312,0,0], [0,81.312,0], [0,0,114.0264]]]],
    "unit" : "MHz"
}
```
With the rotation on, the species' [frame](#tag_df_frame) says its
Cartesian tensor basis. Object values are normalized to MHz. The
original bare indexed list remains supported and is interpreted as MHz. <br>
<b>Data Type:</b> Object {values (or value), unit}, or legacy 2D List

<b>Units:</b> Hz, kHz, MHz (default), or GHz

<b>Default:</b> -

<b>Example:</b> ``` "hyperfine" : [[1, "15N", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]],
                                 [2, "1H", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]]] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_efg"></a>
## efg
<b>Description:</b><br>
The electric field gradient tensors of the paramagnetic defect. <br>
The number of `efg` should be the same as [n_onsite_spins](#tag_df_n_onsite_spins) × [n_axis_states](#tag_df_n_axis_states) <br>
The unit-aware form follows the same indexed tensor syntax as hyperfine:
```
"efg" : {
    "values" : [[1, "14N", [[0.521910,0,0], [0,0.522771,-0.000086], [0,-0.000086,-1.044681]]]],
    "unit" : "Hartree/Bohr^2"
}
```
With the rotation on, the species' [frame](#tag_df_frame) says its basis. Object
values are normalized to Hartree/Bohr^2 before the existing eQ–EFG quadrupole conversion.
The original bare indexed list remains supported in Hartree/Bohr^2. <br>
<b>Data Type:</b> Object {values (or value), unit}, or legacy 2D List

<b>Units:</b> Hartree/Bohr^2 (default; aliases Ha/Bohr^2, Ha/a0^2, a.u.),
V/m^2, or V/angstrom^2

<b>Default:</b> -

<b>Example:</b> ``` "efg" : [[1, "14N", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]],
                               [2, "1H", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]]] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_zfs"></a>
## zfs
<b>Description:</b><br>
The zero-field splitting of the electron spin associated with a paramagnetic defect.
Three input forms are accepted. <br>

The compact axial form constructs one Cartesian tensor from its normalized [axis](#tag_df_axis) n̂:
<p align="center">T = D (n̂ n̂<sup>T</sup> − I/3).</p>
It is broadcast to every axis state. It has no E: an axis fixes z but not the transverse x/y basis, so an E other than
0 needs the explicit tensor (`E`, which took 0 only, is refused, and so is a species `axis` beside the zfs).
```
"frame" : "file",
"zfs" : { "D" : 2870.0, "axis" : [1.0, 1.0, 1.0], "unit" : "MHz" }
```

The shared tensor form accepts either a nested 3x3 matrix or a flat row-major 9-number
array, and also broadcasts it to every axis state:
```
"frame" : "file",
"zfs" : {
    "tensor" : [[0.0, 956.6667, 956.6667],
                [956.6667, 0.0, 956.6667],
                [956.6667, 956.6667, 0.0]],
    "unit" : "MHz"
}
```

The per-index list form is in MHz. Here a separate 3x3 tensor may be supplied for every axis state:
```
"zfs" : [[1, "e", [0.0, 956.6667, 956.6667, 956.6667, 0.0, 956.6667, 956.6667, 956.6667, 0.0]],
         [2, "e", [0.0, 956.6667, 956.6667, 956.6667, 0.0, 956.6667, 956.6667, 956.6667, 0.0]]]
```

<b>Data Type:</b> Object {D, axis, unit}, Object {tensor, unit}, or legacy 2D List

<b>Example:</b> one per form (with the rotation on every form needs the species' `frame`):
    ```
    "zfs" : {"D": 2870.0, "axis": [1.0, 1.0, 1.0], "unit": "MHz"}
    "zfs" : {"tensor": [[0.0, 956.6667, 956.6667], [956.6667, 0.0, 956.6667], [956.6667, 956.6667, 0.0]], "unit": "MHz"}
    "zfs" : [[1, "e", [0.0, 956.6667, 956.6667, 956.6667, 0.0, 956.6667, 956.6667, 956.6667, 0.0]]]
    ```

<b>Unit:</b> Object forms accept Hz, kHz, MHz (default), or GHz and are converted to
MHz while parsing. The legacy indexed form is MHz. The Hamiltonian later converts the
stored MHz values to angular kHz (rad/ms).

<b>Default:</b> zero tensor

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_detuning"></a>
## detuning
<b>Description:</b><br>
Detuning frequencies that will be applied to the electron spin associated with the paramagnetic defect. <br>
The number of `detuning` should be the same as [n_axis_states](#tag_df_n_axis_states) <br>
This tag is used to approximate the spin interaction tensors related to the paramagnetic defect, <br>
e.g., <b>S</b> · <b>A</b> · <b>I</b> ≈ Δ m<sub>I</sub> S<sub>z</sub>. <br>
The unit-aware form wraps the same indexed values in an object:
```
"detuning" : {
    "values" : [[1, "e", -2.14], [2, "e", 0.0], [3, "e", 2.14]],
    "unit" : "MHz"
}
```
The original bare indexed list remains supported and is interpreted as MHz (a qubit's
[detuning](InputFile_Qubit.md#tag_detuning) is kHz). Missing axis states retain the zero default. <br>
<b>Data Type:</b> Object {values (or value), unit}, or legacy 2D List

<b>Unit:</b> The object form accepts Hz, kHz, MHz (default), or GHz and converts to
MHz while parsing. The legacy indexed form is MHz. The Hamiltonian later converts the
stored MHz values to angular kHz (rad/ms).

<b>Default:</b> zero for every unspecified axis state

<b>Example:</b> ``` "detuning" : [[1, "e",  1.0],
                                     [2, "e",  2.0],
                                     [3, "e",  3.0],
                                     [4, "e",  4.0]] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)
