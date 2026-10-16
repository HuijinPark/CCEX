# Defect tags
The `Defect` block describes the paramagnetic defects of the bath: the list of defect species, the model (old key names are in the table of
[old names](InputFile.md#OldNames)). The files of their axes and states are bath-configuration keys:
[configuration.defect_axis_file](InputFile_General.md#tag_defect_axis_file) and
[configuration.defect_state_file](InputFile_General.md#tag_defect_state_file).

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
- [species](#tag_df_species): [name](#tag_df_name), [n_onsite_spins](#tag_df_n_nuclear_spins), [n_axis_states](#tag_df_n_configurations),
  [axis](#tag_df_axis), [coordinate_frame](#tag_df_coordinate_frame), [types](#tag_df_types), [spins](#tag_df_spins),
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
        "axis" : [1.0, 1.0, 1.0], "coordinate_frame" : "bath",
        "zfs" : {"D" : 2870.0, "E" : 0.0, "unit" : "MHz"},
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

<b>Sub-tags:</b>

- [name](#tag_df_name) : Defect name
    - [n_onsite_spins](#tag_df_n_nuclear_spins) : Number of on-site nuclear spins of the PM defect
    - [n_axis_states](#tag_df_n_configurations) : Number of indexed axis states of the defect
    - [axis](#tag_df_axis) : Physical symmetry-axis direction
    - [coordinate_frame](#tag_df_coordinate_frame) : Basis of Defect vectors and tensors
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

<a id="tag_df_n_nuclear_spins"></a>
## n_onsite_spins
<b>Description:</b><br>
Number of on-site (nuclear) spins of the paramagnetic defect. <br>
<b>Data Type:</b> Int

<b>Default:</b> 0

<b>Example:</b> ``` "n_onsite_spins" : 0 ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_n_configurations"></a>
## n_axis_states
<b>Description:</b><br>
Number of indexed axis states available for the paramagnetic defect: an index can be a principal axis, or a
frozen nuclear state (for example three m<sub>I</sub> detunings). A shared object-form zfs is broadcast to indices 1
through n_axis_states. <br>

The indexed entries of [nuclear_positions](#tag_df_nuclear_positions), [hyperfine](#tag_df_hyperfine),
[efg](#tag_df_efg), zfs and [detuning](#tag_df_detuning) must address a
configuration in <b>1 .. n_axis_states</b>; index 0 is reserved for the main spin and is
written by CCEX itself. An index above n_axis_states is refused. <br>

<b>Data Type:</b> Int

<b>Default:</b> 0

<b>Example:</b> ``` "n_axis_states" : 4 ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_axis"></a>
## axis
<b>Description:</b><br>
Physical symmetry-axis direction of this defect, expressed in
[coordinate_frame](#tag_df_coordinate_frame). The vector is normalized when read.
It is required by the D/E form of [zfs](#tag_df_zfs). It is optional metadata for
the explicit tensor form and is normally omitted for a scalar-detuning-only P1 model.
Because an axis alone does not define transverse x and y directions, a non-zero E must
be supplied with the explicit tensor form. <br>

<b>Data Type:</b> Double List(3)

<b>Default:</b> absent

<b>Example:</b> ``` "axis" : [1.0, 1.0, 1.0] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_coordinate_frame"></a>
## Defect coordinate_frame
<b>Description:</b><br>
Basis in which this Defect's axis, nuclear_positions displacement vectors, and hyperfine / efg / zfs
tensor components are written. <tt>"bath"</tt> means the source bath/crystal basis;
when frame.rotation is given they are re-expressed in the computational
qubit frame. <tt>"qubit"</tt> means they are already in that frame and remain unchanged.
Detuning is a scalar frequency and is unaffected. <br>
This key is required when <tt>axis</tt> or an object-form <tt>zfs</tt>, <tt>hyperfine</tt>, or
<tt>efg</tt> is used. Without it the components are taken as written in the computational
frame (<tt>"qubit"</tt>). <br>

<b>Data Type:</b> String: "bath" or "qubit"

<b>Default:</b> "qubit"

<b>Example:</b> ``` "coordinate_frame" : "bath" ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_types"></a>
## types
<b>Description:</b><br>
Name of on-site nuclear spins of the paramagnetic defect. <br>
The number of types should be the same as [n_onsite_spins](#tag_df_n_nuclear_spins). <br>
<b>Data Type:</b> String List

<b>Default:</b> -

<b>Example:</b> ``` "types" : ["1H", "13C"] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_spins"></a>
## spins
<b>Description:</b><br>
Spin numbers of the on-site nuclear spins of the paramagnetic defect. <br>
The number of spins should be the same as [n_onsite_spins](#tag_df_n_nuclear_spins). <br>
<b>Data Type:</b> Float List

<b>Default:</b> -

<b>Example:</b> ``` "spins" : [0.5, 1.0] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_gyros"></a>
## gyros
<b>Description:</b><br>
Gyromagnetic ratio of the on-site nuclear spins of the paramagnetic defect. <br>
The number of gyros should be the same as [n_onsite_spins](#tag_df_n_nuclear_spins). <br>
The unit-aware form is
``` "gyros" : {"values" : [3.0777051853], "unit" : "MHz/T"} ```
Plain frequency units such as MHz/T and kHz/G mean γ/2π; units beginning
with <tt>rad</tt> are angular gyromagnetic ratios. Values are normalized internally to
rad kHz/G, numerically the same as rad/ms/G. The original bare list remains supported
in that internal unit. <br>
<b>Data Type:</b> Object {values, unit}, or legacy Double List

<b>Units:</b> radkHz/G (default), rad/ms/G, kHz/G, MHz/G, MHz/T, Hz/T,
radMHz/T, rad/us/T, or rad/s/T

<b>Default:</b> -

<b>Example:</b> ``` "gyros" : {"values" : [3.0777051853], "unit" : "MHz/T"} ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_quadrupole_moments"></a>
## quadrupole_moments
<b>Description:</b><br>
Quadrupole moment \(eQ\) values of the on-site nuclear spins of the paramagnetic defect. <br>
The number of eQ values should be the same as [n_onsite_spins](#tag_df_n_nuclear_spins). <br>
The stored quantity follows the existing area convention: one internal unit is
10<sup>−30</sup> m² = 1 fm² = 10 mbarn. The unit-aware form is
``` "quadrupole_moments" : {"values" : [20.44], "unit" : "mbarn"} ```
and is normalized to that internal unit. The original bare list remains supported in
10<sup>−30</sup> m². <br>
<b>Data Type:</b> Object {values, unit}, or legacy Double List

<b>Units:</b> 1e-30 m^2 (default), fm^2, m^2, cm^2, barn, or mbarn

<b>Default:</b> [0.0 , ..]

<b>Example:</b> ``` "quadrupole_moments" : {"values" : [20.44], "unit" : "mbarn"} ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_nuclear_positions"></a>
## nuclear_positions
<b>Description:</b><br>
Relative x,y,z position of on-site nuclear spins of the paramagnetic defect over the electron spin position. <br>
The number of `nuclear_positions` should be the same as [n_onsite_spins](#tag_df_n_nuclear_spins) × [n_axis_states](#tag_df_n_configurations) <br>
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
The number of `hyperfine` should be the same as [n_onsite_spins](#tag_df_n_nuclear_spins) × [n_axis_states](#tag_df_n_configurations) <br>
The unit-aware form wraps the same indexed entries in <tt>values</tt>. Each tensor may
be a flat row-major 9-number array or a nested 3x3 array:
```
"hyperfine" : {
    "values" : [[1, "14N", [[81.312,0,0], [0,81.312,0], [0,0,114.0264]]]],
    "unit" : "MHz"
}
```
It requires the Defect's [coordinate_frame](#tag_df_coordinate_frame), because its
Cartesian tensor basis must be explicit. Object values are normalized to MHz. The
original bare indexed list remains supported and is interpreted as MHz. <br>
<b>Data Type:</b> Object {values, unit}, or legacy 2D List

<b>Units:</b> Hz, kHz, MHz (default), or GHz

<b>Default:</b> -

<b>Example:</b> ``` "hyperfine" : [[1, "15N", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]],
                                 [2, "1H", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]]] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_efg"></a>
## efg
<b>Description:</b><br>
The electric field gradient tensors of the paramagnetic defect. <br>
The number of `efg` should be the same as [n_onsite_spins](#tag_df_n_nuclear_spins) × [n_axis_states](#tag_df_n_configurations) <br>
The unit-aware form follows the same indexed tensor syntax as hyperfine:
```
"efg" : {
    "values" : [[1, "14N", [[0.521910,0,0], [0,0.522771,-0.000086], [0,-0.000086,-1.044681]]]],
    "unit" : "Hartree/Bohr^2"
}
```
It requires the Defect's [coordinate_frame](#tag_df_coordinate_frame). Object
values are normalized to Hartree/Bohr^2 before the existing eQ–EFG quadrupole conversion.
The original bare indexed list remains supported in Hartree/Bohr^2. <br>
<b>Data Type:</b> Object {values, unit}, or legacy 2D List

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

The compact axial form constructs one Cartesian tensor from the normalized Defect axis:
<p align="center">T = D (n̂ n̂<sup>T</sup> − I/3).</p>
It is broadcast to every indexed configuration. E may be omitted or set to zero. A
non-zero E is refused because one axis fixes z but not the transverse x/y basis; use an
explicit tensor in that case.
```
"axis" : [1.0, 1.0, 1.0],
"coordinate_frame" : "bath",
"zfs" : { "D" : 2870.0, "E" : 0.0, "unit" : "MHz" }
```

The shared tensor form accepts either a nested 3x3 matrix or a flat row-major 9-number
array, and also broadcasts it to every indexed configuration:
```
"coordinate_frame" : "bath",
"zfs" : {
    "tensor" : [[0.0, 956.6667, 956.6667],
                [956.6667, 0.0, 956.6667],
                [956.6667, 956.6667, 0.0]],
    "unit" : "MHz"
}
```

The per-index list form is in MHz. Here a separate 3x3 tensor may be supplied for every configuration:
```
"zfs" : [[1, "e", [0.0, 956.6667, 956.6667, 956.6667, 0.0, 956.6667, 956.6667, 956.6667, 0.0]],
         [2, "e", [0.0, 956.6667, 956.6667, 956.6667, 0.0, 956.6667, 956.6667, 956.6667, 0.0]]]
```

<b>Data Type:</b> Object {D, E, unit}, Object {tensor, unit}, or legacy 2D List

<b>Unit:</b> Object forms accept Hz, kHz, MHz (default), or GHz and are converted to
MHz while parsing. The legacy indexed form is MHz. The Hamiltonian later converts the
stored MHz values to angular kHz (rad/ms).

<b>Default:</b> zero tensor

<b>Parent-tags:</b> [Defect.species](#tag_df_species)

<a id="tag_df_detuning"></a>
## detuning
<b>Description:</b><br>
Detuning frequencies that will be applied to the electron spin associated with the paramagnetic defect. <br>
The number of `detuning` should be the same as [n_axis_states](#tag_df_n_configurations) <br>
This tag is used to approximate the spin interaction tensors related to the paramagnetic defect, <br>
e.g., <b>S</b> · <b>A</b> · <b>I</b> ≈ Δ m<sub>I</sub> S<sub>z</sub>. <br>
The unit-aware form wraps the same indexed values in an object:
```
"detuning" : {
    "values" : [[1, "e", -2.14], [2, "e", 0.0], [3, "e", 2.14]],
    "unit" : "MHz"
}
```
The original bare indexed list remains supported and is interpreted as MHz. Missing
configuration indices retain the zero default. <br>
<b>Data Type:</b> Object {values, unit}, or legacy 2D List

<b>Unit:</b> The object form accepts Hz, kHz, MHz (default), or GHz and converts to
MHz while parsing. The legacy indexed form is MHz. The Hamiltonian later converts the
stored MHz values to angular kHz (rad/ms).

<b>Default:</b> zero for every unspecified configuration

<b>Example:</b> ``` "detuning" : [[1, "e",  1.0],
                                     [2, "e",  2.0],
                                     [3, "e",  3.0],
                                     [4, "e",  4.0]] ```

<b>Parent-tags:</b> [Defect.species](#tag_df_species)
