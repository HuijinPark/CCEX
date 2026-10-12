# DefectArray-related tags
This tags are used for adjusting paramagnetic defect bath information.
The format is the following:

<details><summary><b>Example:</b> Click here to show example code</summary>
```
{
"Defect"         : [
		{
			"dfname"	: "NV0",
			"naddspin" 	: 1,
			"types" 	: ["15N"],
			"spins" 	: [0.5],   
			"gyros" 	: [-2.71261804],
			"navaax" 	: 4,
			"rxyzs"      : [ [1, "15N", [  1.36404 ,   -0.78843 ,   -0.56598 ]],
		                     [2, "15N", [ -1.36801 ,   -0.79025 ,   -0.56701 ]],
		                     [3, "15N", [  0.00000 ,    1.57603 ,   -0.56511 ]],
		                     [4, "15N", [  0.00000 ,   -0.00733 ,    1.67433 ]]],
			"hypf" 		: [  [1, "15N", [ -10.33219 ,    1.01563 ,    1.50162 ,    1.01563 ,   -5.18317 ,   -0.16410 ,    1.50162 ,   -0.16410 ,   -6.50270]],  
			  				 [2, "15N", [  -7.34999 ,   -2.73741 ,   -0.89293 ,   -2.73741 ,   -8.16537 ,   -1.21839 ,   -0.89293 ,   -1.21839 ,   -6.50270]], 
			  				 [3, "15N", [  -5.59086 ,    1.72177 ,   -0.60870 ,    1.72177 ,   -9.92450 ,    1.38249 ,   -0.60870 ,    1.38249 ,   -6.50270]],
			  				 [4, "15N", [  -6.01395 ,    0.00000 ,    0.00000 ,    0.00000 ,   -5.59086 ,   -1.82620 ,    0.00000 ,   -1.82620 ,  -10.41325]]]
		},
        {
            ...
        }
                       ]
}
```
</details>

- [Defect](#tag_Defect)       
- [dfname](#tag_df_dfname)       
- [apprx](#tag_df_apprx)        
- [naddspin](#tag_df_naddspin)     
- [navaax](#tag_df_navaax)       
- [axis](#tag_df_axis)
- [coordinate_frame](#tag_df_coordinate_frame)
- [types](#tag_df_types)        
- [spins](#tag_df_spins)        
- [gyros](#tag_df_gyros)        
- [eqs](#tag_df_eqs)          
- [rxyzs](#tag_df_rxyzs)        
- [hypf](#tag_df_hypf)         
- [efg](#tag_df_efg)          
- [zfs](#tag_df_zfs)          
- [detuning](#tag_df_detuning)

<a id="tag_Defect"></a>
## Defect
<b>Description:</b><br>
Detailed information of paramagnetic \(PM\) defect bath can be provided by this tag.
<b>Data Type:</b> Dictionary List

<b>Default:</b> -

<b>Example:</b>
    ```
    "Defect" : [{
        "dfname" : "NV_0", "naddspin" : 0, "navaax" : 3,
        "axis" : [1.0, 1.0, 1.0], "coordinate_frame" : "bath",
        "zfs" : {"D" : 2870.0, "E" : 0.0, "unit" : "MHz"},
        "detuning" : {
            "values" : [[1,"e",-2.14], [2,"e",0.0], [3,"e",2.14]],
            "unit" : "MHz"
        }
    }]
    ```

<b>Related variable:</b> [DefectArray struct](#DefectArray) [Defect struct](#Defect)

<b>Sub-tags:</b>

- [dfname](#tag_df_dfname) : Defect name
    - [apprx](#tag_df_apprx)   : Approximation type
    - [naddspin](#tag_df_naddspin) : Number of on-site nuclear spins of the PM defect
    - [navaax](#tag_df_navaax) : Number of indexed defect configurations
    - [axis](#tag_df_axis) : Physical symmetry-axis direction
    - [coordinate_frame](#tag_df_coordinate_frame) : Basis of Defect vectors and tensors
    - [types](#tag_df_types) : On-site nuclear spin types
    - [spins](#tag_df_spins) : Spin numbers of the on-site nuclear spins
    - [gyros](#tag_df_gyros) : Gyromagnetic ratio of the on-site nuclear spins
    - [eqs](#tag_df_eqs) : eQ values of the on-site nuclear spins
    - [rxyzs](#tag_df_rxyzs) : Relative x,y,z position of on-site nuclear spins
    - [hypf](#tag_df_hypf) : Hyperfine interaction tensor
    - [efg](#tag_df_efg) : Electric field gradient tensor
    - [zfs](#tag_df_zfs) : Zero-field splitting tensor
    - [detuning](#tag_df_detuning) : Detuning frequency

<a id="tag_df_dfname"></a>
## dfname
<b>Description:</b><br>
Name of the paramagnetic defect. <br>
<b>Data Type:</b> String

<b>Default:</b> -

<b>Example:</b> ``` dfname : "P1" ```

<b>Related variable:</b> [Defect](#Defect)::dfname

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_apprx"></a>
## apprx
<b>Description:</b><br>
If `apprx` is set to `false,` the electron and nuclear spins of the paramagnetic defect are clusterized separately. <br>
This corresponds to a conventional CCE simulation for all spins. <br>
If `apprx` is set to `true,` the electron and nuclear spins are treated as one set, <br>
resulting in a simulation based on the CCE-NeNn method. <br>
<b>Data Type:</b> Bool

<b>Default:</b> true

<b>Example:</b> ``` apprx : true ```

<b>Related variable:</b> [Defect](#Defect)::apprx

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_naddspin"></a>
## naddspin
<b>Description:</b><br>
Number of on-site nuclear spins of the paramagnetic defect. <br>
<b>Data Type:</b> Int

<b>Default:</b> 0

<b>Example:</b> ``` naddspin : 0 ```

<b>Related variable:</b> [DefectArray](#DefectArray)::naddspins [Defect](#Defect)::naddspin

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_navaax"></a>
## navaax
<b>Description:</b><br>
Number of indexed configurations available for the paramagnetic defect. The historical
name is retained for compatibility, but an index can represent a frozen nuclear state
(for example three $m_I$ detunings) rather than a different physical axis. A shared
object-form zfs is broadcast to indices 1 through navaax. <br>

The indexed entries of [rxyzs](#tag_df_rxyzs), [hypf](#tag_df_hypf),
[efg](#tag_df_efg), zfs and [detuning](#tag_df_detuning) must address a
configuration in <b>1 .. navaax</b>; index 0 is reserved for the main spin and is
written by CCEX itself. <br>
<b>Changed since v1.0.0:</b> an index of navaax+1 is now refused. It used to be
accepted, and wrote one element past the end of every one of those arrays, so an
input that relied on it was already producing undefined results. This is the one
change here that can stop an input which previously ran to completion. <br>

<b>Data Type:</b> Int

<b>Default:</b> 0

<b>Example:</b> ``` navaax : 4 ```

<b>Related variable:</b> [Defect](#Defect)::navaax [DefectArray](#DefectArray)::paxes

<b>Parent-tags:</b> [Defect](#tag_Defect)

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

<b>Related variable:</b> [Defect](#Defect)::axis

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_coordinate_frame"></a>
## Defect coordinate_frame
<b>Description:</b><br>
Basis in which this Defect's axis, rxyzs displacement vectors, and hypf / efg / zfs
tensor components are written. <tt>"bath"</tt> means the source bath/crystal basis;
when coordinate_frame_rotation is enabled they are re-expressed in the computational
qubit frame. <tt>"qubit"</tt> means they are already in that frame and remain unchanged.
Detuning is a scalar frequency and is unaffected. <br>
This key is required when <tt>axis</tt> or an object-form <tt>zfs</tt>, <tt>hypf</tt>, or
<tt>efg</tt> is used. For a legacy Defect without those features it defaults to
<tt>"qubit"</tt>, preserving the old behavior. <br>

<b>Data Type:</b> String: "bath" or "qubit"

<b>Default:</b> "qubit" for legacy entries

<b>Example:</b> ``` "coordinate_frame" : "bath" ```

<b>Related variable:</b> [Defect](#Defect)::coordinate_frame

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_types"></a>
## types
<b>Description:</b><br>
Name of on-site nuclear spins of the paramagnetic defect. <br>
The number of types should be the same as [naddspin](#tag_df_naddspin). <br>
<b>Data Type:</b> String List

<b>Default:</b> -

<b>Example:</b> ``` types : ["1H", "13C"] ```

<b>Related variable:</b> [Defect](#Defect)::types

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_spins"></a>
## spins
<b>Description:</b><br>
Spin numbers of the on-site nuclear spins of the paramagnetic defect. <br>
The number of spins should be the same as [naddspin](#tag_df_naddspin). <br>
<b>Data Type:</b> Float List

<b>Default:</b> -

<b>Example:</b> ``` spins : [0.5, 1.0] ```

<b>Related variable:</b> [Defect](#Defect)::spins

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_gyros"></a>
## gyros
<b>Description:</b><br>
Gyromagnetic ratio of the on-site nuclear spins of the paramagnetic defect. <br>
The number of gyros should be the same as [naddspin](#tag_df_naddspin). <br>
The unit-aware form is
``` "gyros" : {"values" : [3.0777051853], "unit" : "MHz/T"} ```
Plain frequency units such as MHz/T and kHz/G mean $\gamma/2\pi$; units beginning
with <tt>rad</tt> are angular gyromagnetic ratios. Values are normalized internally to
rad kHz/G, numerically the same as rad/ms/G. The original bare list remains supported
in that internal unit. <br>
<b>Data Type:</b> Object {values, unit}, or legacy Double List

<b>Units:</b> radkHz/G (default), rad/ms/G, kHz/G, MHz/G, MHz/T, Hz/T,
radMHz/T, rad/us/T, or rad/s/T

<b>Default:</b> -

<b>Example:</b> ``` "gyros" : {"values" : [3.0777051853], "unit" : "MHz/T"} ```

<b>Related variable:</b> [Defect](#Defect)::gyros

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_eqs"></a>
## eqs
<b>Description:</b><br>
Quadrupole moment \(eQ\) values of the on-site nuclear spins of the paramagnetic defect. <br>
The number of eQ values should be the same as [naddspin](#tag_df_naddspin). <br>
The stored quantity follows the existing area convention: one internal unit is
$10^{-30}\,\mathrm{m^2}=1\,\mathrm{fm^2}=10\,\mathrm{mbarn}$. The unit-aware form is
``` "eqs" : {"values" : [20.44], "unit" : "mbarn"} ```
and is normalized to that internal unit. The original bare list remains supported in
$10^{-30}\,\mathrm{m^2}$. <br>
<b>Data Type:</b> Object {values, unit}, or legacy Double List

<b>Units:</b> 1e-30 m^2 (default), fm^2, m^2, cm^2, barn, or mbarn

<b>Default:</b> [0.0 , ..]

<b>Example:</b> ``` "eqs" : {"values" : [20.44], "unit" : "mbarn"} ```

<b>Related variable:</b> [Defect](#Defect)::eqs

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_rxyzs"></a>
## rxyzs
<b>Description:</b><br>
Relative x,y,z position of on-site nuclear spins of the paramagnetic defect over the electron spin position. <br>
The number of `rxyzs` should be the same as \" [naddspin](#tag_df_naddspin) x [navaax](#tag_df_navaax) \" <br>
<b>Data Type:</b> 2D List

<b>Unit:</b> angstrom

<b>Default:</b> [0.0, 0.0, 0.0]

<b>Example:</b> ``` rxyzs : [[1, "15N", [0.0, 0.0, 0.0]],
                                  [2, "1H", [0.0, 0.0, 0.0]]] ```

<b>Related variable:</b> [Defect](#Defect)::rxyzs

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_hypf"></a>
## hypf
<b>Description:</b><br>
The hyperfine interaction between on-site nuclear spins and the electron spin associated with a paramagnetic defect. <br>
The number of `hypf` should be the same as \" [naddspin](#tag_df_naddspin) x [navaax](#tag_df_navaax) \" <br>
The unit-aware form wraps the same indexed entries in <tt>values</tt>. Each tensor may
be a flat row-major 9-number array or a nested 3x3 array:
```
"hypf" : {
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

<b>Example:</b> ``` hypf : [[1, "15N", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]],
                                 [2, "1H", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]]] ```

<b>Related variable:</b> [Defect](#Defect)::hypf

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_efg"></a>
## efg
<b>Description:</b><br>
The electric field gradient tensors of the paramagnetic defect. <br>
The number of `efg` should be the same as \" [naddspin](#tag_df_naddspin) x [navaax](#tag_df_navaax) \" <br>
The unit-aware form follows the same indexed tensor syntax as hypf:
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

<b>Example:</b> ``` efg : [[1, "14N", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]],
                               [2, "1H", [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]]] ```

<b>Related variable:</b> [Defect](#Defect)::efg

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_zfs"></a>
## zfs
<b>Description:</b><br>
The zero-field splitting of the electron spin associated with a paramagnetic defect.
Three input forms are accepted. <br>

The compact axial form constructs one Cartesian tensor from the normalized Defect axis:
$$T = D\left(\hat n\hat n^T-I/3\right).$$
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

The original per-index syntax remains supported and is interpreted in MHz exactly as
before. Here a separate 3x3 tensor may be supplied for every configuration:
```
"zfs" : [[1, "e", [0.0, 956.6667, 956.6667, 956.6667, 0.0, 956.6667, 956.6667, 956.6667, 0.0]],
         [2, "e", [0.0, 956.6667, 956.6667, 956.6667, 0.0, 956.6667, 956.6667, 956.6667, 0.0]]]
```

<b>Data Type:</b> Object {D, E, unit}, Object {tensor, unit}, or legacy 2D List

<b>Unit:</b> Object forms accept Hz, kHz, MHz (default), or GHz and are converted to
MHz while parsing. The legacy indexed form is MHz. The Hamiltonian later converts the
stored MHz values to angular kHz (rad/ms).

<b>Default:</b> zero tensor

<b>Related variable:</b> [Defect](#Defect)::zfs

<b>Parent-tags:</b> [Defect](#tag_Defect)

<a id="tag_df_detuning"></a>
## detuning
<b>Description:</b><br>
Detuning frequencies that will be applied to the electron spin associated with the paramagnetic defect. <br>
The number of `detuning` should be the same as \" [navaax](#tag_df_navaax) \" <br>
This tag is used to approximate the spin interaction tensors related to the paramagnetic defect, <br>
e.g., $$\vec{S} \cdot \mathbf{A} \cdot \vec{I} \approx \Delta m_I S_z$$. <br>
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

<b>Example:</b> ``` detuning : [[1, "e",  1.0],
                                     [2, "e",  2.0],
                                     [3, "e",  3.0],
                                     [4, "e",  4.0]] ```

<b>Related variable:</b> [Defect](#Defect)::detuning

<b>Parent-tags:</b> [Defect](#tag_Defect)
