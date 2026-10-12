# General tags
Tags in this section is about simulatation parameters such as file path, numerical parameters.

- [method](#tag_method)       
- [quantity](#tag_quantity)     
- [order](#tag_order)        
- [bfield](#tag_bfield)       
- [rbath](#tag_rbath)        
- [rdip](#tag_rdip)         
- [deltat](#tag_deltat)       
- [nstep](#tag_nstep)        
- [rbathcut](#tag_rbathcut)     
- [rdipcut](#tag_rdipcut)       
- [rdip_ladder](#tag_rdip_ladder)
- [nstate](#tag_nstate)        
- [vector_mean_field](#tag_vector_mean_field)
- [seed](#tag_seed)       
- [qubitfile](#tag_qubitfile)   
- [gyrofile](#tag_gyrofile)      
- [bathfile](#tag_bathfile)      
- [bathadjust](#tag_bathadjust)
- [coordinate_frame_rotation](#tag_coordinate_frame_rotation)
- [defect_axis_reference](#tag_defect_axis_reference)
- [avaaxfile](#tag_avaaxfile)
- [statefile](#tag_statefile)                     
- [exstatefile](#tag_exstatefile)                  
- [DefectTotSpin](#tag_DefectTotSpin)                
- [CorrTotSpin](#tag_CorrTotSpin)                  
- [hf_readmode](#tag_hf_readmode)                       
- [hf_tensorfile](#tag_hf_tensorfile)                
-                    
- [hf_ignore_oor](#tag_hf_ignore_oor)                
- [qd_readmode](#tag_qd_readmode)                  
- [qd_tensorfile](#tag_qd_tensorfile)                
- [qd_cellpara](#tag_qd_cellpara)                     
- [qd_tensorfile_woqubit](#tag_qd_tensorfile_woqubit)
- [hfmedi](#tag_hfmedi)
- [knight](#tag_knight)


<a id="tag_method"></a>
## method
<b>Description:</b><br>
CCE Methods to calculate the dynamics of the system.

You can specify this value via the "-m" command-line option. <br>
<b>Data Type:</b> String

<b>Options:</b><br>
| Option  | Description |
|---------|-------------------------------------|
| `CCE`   | Conventional CCE method             |
| `gCCE`  | Generalized CCE method              |
| `pCCE`  | Partition CCE method                |
| `Dsj`   | Disjoint method                     |
| `Itb`   | Inter bath method                   |
| `Dsjitb`| Disjoint + Inter bath method        |

<b>Default:</b> `CCE`

<b>Example:</b> ``` method : CCE ```

<b>Related variable:</b> [Config](#Config)::order

<a id="tag_quantity"></a>
## quantity
<b>Description:</b><br>
Quantity to calculate the dynamics of the system.

You can specify this value via the "-q" command-line option. <br>
<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|--------|-------------------------------------|
| `coherence` | Coherence function             |
| `dm`        | Density matrix                 |

With `dm` (gCCE), the reduced density matrix of the qubits is written as blocks of rows, one
block per time step, and the cluster expansion is applied element by element. An element the
free evolution keeps at zero (a level the initial state does not touch) is left at zero; when
clusters would put a value into it, the run prints a `[dm]` line with the largest such value,
since the elementwise product cannot carry it. Until 2026-09-24 those elements came out nan. <br>
<b>Unit:</b> arbitrary unit (a.u.)

<b>Default:</b> `coherence`

<b>Example:</b> ``` quantity : coherence ```

<b>Related variable:</b> [Config](#Config)::quantity

<a id="tag_order"></a>
## order
<b>Description:</b><br>
Order of the CCE method.

<b>Data Type:</b> Integer

<b>Default:</b> -

<b>Example:</b> ``` order : 2 ```
<b>Related variable:</b> [Config](#Config)::order

<a id="tag_bfield"></a>
## bfield
<b>Description:</b><br>
Magnetic field.

You can specify a magnetic field along the z-axis via the "-B" command-line option. <br>
<b>Data Type:</b> Float or Float List

<b>Unit:</b> gauss (G)

<b>Default:</b> 0.0

<b>Example:</b>

``` bfield : 0.0 ```

``` bfield : [0.0, 0.0, 0.0] ```

<b>Related variable:</b> [Config](#Config)::bfield

<a id="tag_rbath"></a>
## rbath
<b>Description:</b><br>
Bath radius

<b>Data Type:</b> Float

<b>Unit:</b> angstrom

<b>Default:</b> -

<b>Example:</b> ``` rbath : 0.0 ```

<b>Related variable:</b> [Config](#Config)::rbath

<a id="tag_rdip"></a>
## rdip
<b>Description:</b><br>
Dipolar radius

<b>Data Type:</b> Float

<b>Unit:</b> angstrom

<b>Default:</b> -

<b>Example:</b> ``` rdip : 0.0 ```

<b>Related variable:</b> [Config](#Config)::rdip

<a id="tag_deltat"></a>
## deltat
<b>Description:</b><br>
Time interval.

<b>Data Type:</b> Float

<b>Unit:</b> ms

<b>Default:</b> -

<b>Example:</b> ``` deltat : 0.0 ```

<b>Related variable:</b> [Config](#Config)::deltat

<a id="tag_nstep"></a>
## nstep
<b>Description:</b><br>
Number of time steps.

<b>Data Type:</b> Integer

<b>Default:</b> -

<b>Example:</b> ``` nstep : 10 ```

<b>Related variable:</b> [Config](#Config)::nstep

<a id="tag_rbathcut"></a>
## rbathcut
<b>Description:</b><br>
Bath inner-cutoff radius.

<b>Data Type:</b> Float

<b>Unit:</b> angstrom

<b>Default:</b> 0.0

<b>Example:</b> ``` rbathcut : 0.0 ```

<b>Related variable:</b> [Config](#Config)::rbathcut

<a id="tag_rdipcut"></a>
## rdipcut
<b>Description:</b><br>
Dipolar inner-cutoff radius.

<b>Data Type:</b> Float

<b>Unit:</b> angstrom

<b>Default:</b> 0.0

<b>Example:</b> ``` rdipcut : 0.0 ```

<b>Related variable:</b> [Config](#Config)::rdipcut

<a id="tag_rdip_ladder"></a>
## rdip_ladder
<b>Description:</b><br>
Extra rdip values (rungs) computed from the same run. With rungs r<sub>1</sub> < r<sub>2</sub> < ... below `rdip`,
the run also writes, for every rung, the coherence that a separate run at `rdip` = r<sub>k</sub> gives
(`<outfile>_rdip<r_k>_state<s>_wiDiv`/`_noDiv`, or `<outfile>_rdip<r_k>_wiDiv`/`_noDiv` for `nstate` 0; r<sub>k</sub> is
printed with `%g`). The main output is the run at `rdip` as before. The automatic convergence uses it for its rdip step.

Order 2 only, of `method` CCE, gCCE (every path: ensemble, `evolution` vector and matrix) or pCCE (both
[pcce_connections](InputFile_Cluster.md#tag_pcce_connections) rules; the partition does not depend on `rdip`, only which
partition pairs are connected), with `quantity` coherence, `cluster_mode` standard, `savemode` normal, every pair within
`rdip` as a cluster (no `nk` cut) and without `cluster_contributions`. Anything else is refused with
`Error: rdip_ladder ...`. Each rung is its separate run up to the rounding of the order of the products
(`test/code_verification/check_rdip_ladder`: equal to within 1e-9; on a curve that diverges far above 1, a difference
in the last bit of the double, as between two rank counts of one run).

<b>Data Type:</b> Array of floats, positive, strictly increasing, each between `rdipcut` and `rdip` (1 to 64 values)

<b>Unit:</b> angstrom

<b>Default:</b> - (off)

<b>Example:</b> ``` "rdip_ladder" : [200, 300, 450] ```

<b>Related variable:</b> [Config](#Config)::rdip_ladder

<a id="tag_nstate"></a>
## nstate
<b>Description:</b><br>
Number of state templates to be generated.

You can specify this value via the "-N" command-line option. <br>
<b>Data Type:</b> Integer

<b>Options:</b><br>
| Option | Description |
|--------|-------------------------------------|
| `0`    | No state template will be generated (Ensemble approach) |
| `>0`   | State template will be generated (Single-sample approach) <br> Mean-field effect will be automatically added |

<b>Default:</b> 0

<b>Example:</b> ``` nstate : 10 ```

<b>Related Tags:</b> [seed](#tag_seed)

<b>Related variable:</b> [Config](#Config)::nstate

<a id="tag_vector_mean_field"></a>
## vector_mean_field
<b>Description:</b><br>
With [nstate](#tag_nstate) &gt; 0 each cluster feels the static field of every bath spin outside it, from that spin's sampled
state m<sub>j</sub>. By default only the z part is kept: J<sub>zz</sub> m<sub>j</sub> I<sub>z</sub>. With `true` the whole
vector J·(0, 0, m<sub>j</sub>) = J[:, z] m<sub>j</sub> is used, so the x and y parts enter too; a cluster removes its own
members' field (all three parts) and keeps their full coupling tensor instead. The x and y parts matter at first order for a
bath spin whose own z field nearly cancels (a close neighbour's J<sub>zz</sub> m<sub>j</sub> against its Zeeman term, as
in clustered P1 baths) or whose axis the qubit tilts; there the z-only field can make the cluster expansion diverge
((internal notes), 강건성 6–13). The transverse field is added and removed only between main bath
spins: every mean-field term involving a defect sub-spin (e.g. a P1's <sup>14</sup>N) stays z only, while the couplings
inside a cluster keep their full tensors. Added 2026-10-03; off by default, so earlier results are unchanged. Against
exact 10-spin subsystems it removed the divergence of a broken state but was not more accurate in every normal region,
so it is not the default (same document, 12–13). <br>
<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "vector_mean_field" : true ```

<b>Related Tags:</b> [nstate](#tag_nstate)

<b>Related variable:</b> [Config](#Config)::vector_mean_field

<a id="tag_seed"></a>
## seed
<b>Description:</b><br>
Seed for random number generator. <br>
If you want to generate the same state templates, you should set the same seed value.

<b>Data Type:</b> Integer

<b>Default:</b> Time(NULL)

<b>Example:</b> ``` seed : 1241798 ```

<b>Related Tags:</b> [nstate](#tag_nstate)

<b>Related variable:</b> [Config](#Config)::seed

<a id="tag_qubitfile"></a>
## qubitfile
<b>Description:</b><br>
File path of the qubit spatial configuration. <br>
This tag is required when you use the qubit configuration file. <br>
Also, this only supports for a single qubit system. <br>

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

<b>Example:</b> ``` qubitfile : "filepath" ```

<b>Related variable:</b> [Config](#Config)::qubitfile

<b>Related Tags:</b> [QubitArrayTags](#QubitArrayTags)

<a id="tag_gyrofile"></a>
## gyrofile
<b>Description:</b><br>
File name including bath spin information, e.g. gyromagnetic ratio, spin number <br>

<pre>
    File format is the following :
    [Name] [SpinNumber] [Gyromagnetic ratio (rad/ms/G)]
</pre>
<pre>
    example)
    13C 0.5 6.728284
    P1  0.5 -17608.597050
</pre>

Optional. A bath species the gyrofile does not name (or every species, without a gyrofile) is taken
from the built-in isotope table: every nuclear isotope (mass number > 0) with a positive spin of the EasySpin isotope
database (magnetic moments from N. Stone, IAEA INDC(NDS)-0658, 2014), named mass number + symbol
(`1H`, `11B`, `13C`, `14N`, `29Si`, ...), with γ = g<sub>n</sub> μ<sub>N</sub>/ħ (CODATA 2018). The table is
`src/input/isotope_table.cpp`, written by `example/tools/IsotopeTable/make_isotope_table.py`. A gyrofile
entry wins over the table (a value that differs from it by more than 0.1 % is noted in the log);
species that are not nuclear isotopes (`P1`, `Vp`, an electron spin) must be in the gyrofile. A name
found in neither stops the run. The log says where each species came from. <br>

<b>Data Type:</b> String (File path)

<b>Default:</b> - (the built-in isotope table)

<b>Example:</b> ``` gyrofile : "filepath" ```

<b>Related variable:</b> [Config](#Config)::gyrofile

<a id="tag_bathfile"></a>
## bathfile
<b>Description:</b><br>
File name of the bath spatial configuration. <br>
You can specify this value via the "-I" command-line option. <br>

<pre>
    File format is the following :
    [x] [y] [z] [Name]
</pre>
<pre>
    example)
    1.0 0.0 0.0 0.0   - Total line number 
    1.0 2.0 3.0 13C   - x,y,z in cartesian, angstrom
    4.0 5.0 6.0  P1   - Spin name should be the same to the name in [gyrofile](#tag_gyrofile)
</pre>

<b>Data Type:</b> String (File path) or String List

<b>Default:</b> -

<b>Example:</b> 
    ``` bathfile : "filepath" ```

``` bathfile : ["filepath1", "filepath2", ... ] ```

<b>Related variable:</b> [Config](#Config)::bathfile

<a id="tag_bathadjust"></a>
## bathadjust
<b>Description:</b><br>
Adjust the bath configuration. <br>
If you want to adjust the bath configuration, you can set the bathadjust value. <br>
Usually, this tag is required for multiple qubit calculation that wouldb have different bath configuration <br>

<b>Data Type:</b> Double List

<b>Default:</b> -

<b>Example:</b> 
    ``` bathadjust : [1.0, 2.0, 3.0] ```

``` bathadjust : [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], ... [7.0, 8.0, 9.0]] ```

<b>Related variable:</b> [Config](#Config)::bathadjust

<a id="tag_coordinate_frame_rotation"></a>
## coordinate_frame_rotation
<b>Description:</b><br>
Re-express the qubits and the bath in the frame whose z axis is the reference qubit's
symmetry axis. <br>
Use it when the bath file is written in the crystal frame (e.g. diamond [001]) but the
qubit is quantized along another direction (e.g. an NV along [111]). <br>
More than one qubit is supported: every Qubit.xyz and every bath-spin position receives
the same $R$ and the same $r_0$. <br>
The keys of the block are "enabled", "bath_axis", "qubit_axis", "reference_qubit",
"qubit_position_frame", "hf_tensor_frame" and "qd_tensor_frame". What is transformed,
what is refused and what each key means are in `docs/methods/coordinate-frame-rotation.md`. <br>

<b>Data Type:</b> Object { Boolean, Double List(3), Double List(3), String x4 }

<b>Default:</b> ``` {"enabled": false, "bath_axis": [0,0,1], "qubit_axis": [0,0,1]} ```
An absent section, or "enabled" : false, leaves the run numerically identical to a build
without this option.

<b>Example:</b>
    ```
    "coordinate_frame_rotation" : {
        "enabled"         : true,
        "bath_axis"       : [0.0, 0.0, 1.0],
        "qubit_axis"      : [1.0, 1.0, 1.0],
        "reference_qubit" : "q0",
        "qubit_position_frame" : "bath",
        "hf_tensor_frame" : "bath",
        "qd_tensor_frame" : "bath"
    }
    ```

<b>Related variable:</b> [Config](#Config)::rot_enabled, [Config](#Config)::rot_bath_axis, [Config](#Config)::rot_qubit_axis, [Config](#Config)::rot_reference_qubit

<a id="tag_defect_axis_reference"></a>
## defect_axis_reference
<b>Description:</b><br>
Name the axis a future Jahn-Teller classification should measure defect axes against. <br>
<b>This key changes no number that this version computes.</b> It is read, validated and
written to the run log, and nothing else: it neither performs JT assignment nor changes
the Defect Hamiltonian. Setting it does not move the coherence, the avaax assignments or
the detuning table. This is separate from each Defect's own
[coordinate_frame](InputFile_Defect.md#tag_df_coordinate_frame): bath-frame Defect vectors and tensors
are still converted by coordinate_frame_rotation regardless of this reference key. <br>

The only accepted value is the string "qubit_axis", meaning the central qubit's symmetry
axis. That axis is defined only by
[coordinate_frame_rotation](#tag_coordinate_frame_rotation), so the key requires that
section to be present and enabled; otherwise it is an error rather than a no-op. <br>
Anything else is refused: another string, an array, a boolean. <br>

With the rotation on, the run log reports what was recorded and says plainly that it is
not applied: <br>
<pre>
      Defect axis reference           : qubit_axis
      Reference axis, source frame    : [ normalize(qubit_axis) ]
      Reference axis, calc frame      : [   0, 0, 1 ]
      Status                          : recorded for future JT-axis classification;
                                        not applied to the Defect Hamiltonian
</pre>

<b>Data Type:</b> String

<b>Default:</b> absent -- no reference axis is recorded.

<b>Example:</b>
    ``` "defect_axis_reference" : "qubit_axis" ```

<b>Related variable:</b> [Config](#Config)::defect_axis_reference

<a id="tag_avaaxfile"></a>
## avaaxfile
<b>Description:</b><br>
File path of the principal axis of the bath spin. <br>
This tag might require during paramagnetic defect bath simulation. <br>
During PM defect bath simulation, if you don\'t give this file, the principal axes are automatically assigned. <br>

You can specify this value via the "-a" command-line option. <br>
<pre>
    File format is the following :
    [Principal axis index]
</pre>
<pre>
    example)
    5.0 - The number of total line
    1.0 - principal axis index
    2.0 -  :
    3.0
    3.0
</pre>

<b>Data Type:</b> String (File path)

<b>Default:</b> -

<b>Example:</b> ``` avaaxfile : "filepath" ```

<b>Related tags:</b> [Defect](InputFile_Defect.md#tag_Defect)

<b>Related variable:</b> [Config](#Config)::avaaxfile

<a id="tag_statefile"></a>
## statefile
<b>Description:</b><br>
File name of the bath state template. <br>
This tag might require when you want to obtain the result for a certain set of initial bath state. <br>
If the file path of the statefile doesn\'t exist, the program will automatically generate the state template. <br>
You can specify this value via the "-s" command-line option. <br>

<pre>
    File format is the following :
    [SpinSubLevel]
</pre>
<pre>
    example)
     5.0   - The number of total line
    -0.5   - Spin sublevel when S=1/2
     0.5   -    :
    -0.5
    -0.5
</pre>

<b>Data Type:</b> String (File path)

<b>Default:</b> -

<b>Example:</b> ``` statefile : "filepath" ```

<b>Related tags:</b> [nstate](#tag_nstate)

<b>Related variable:</b> [Config](#Config)::statefile

<a id="tag_exstatefile"></a>
## exstatefile
<b>Description:</b><br>
File name of the additional bath spin state template. <br>
This tag might require when you want to calculate paramagnetic defect bath. <br>
Each bath state represent the additional bath spin state. <br>
If the file path of the exstatefile doesn\'t exist, the program will automatically generate the state template. <br>
You can specify this value via the "-S" command-line option. <br>

<pre>
    File format is the following :
    [SpinSubLevel]
</pre>
<pre>
    example)
     5.0   - The number of total line
    -1.0   - Spin sublevel when S=1
     0.0   -    :
     1.0
    -1.0
</pre>

<b>Data Type:</b> String (File path)

<b>Default:</b> -

<b>Example:</b> ``` exstatefile : "filepath" ```

<b>Related tags:</b> [Defect](InputFile_Defect.md#tag_Defect)

<b>Related variable:</b> [Config](#Config)::exstatefile

<a id="tag_DefectTotSpin"></a>
## DefectTotSpin
<b>Description:</b><br>
Total spin number of the qubit. <br>
(Required when you use Hyperfine tensor file obtained by QUANTUM ESPRESSO) <br>

<b>Data Type:</b> Double

<b>Default:</b> 1.0

<b>Example:</b> ``` DefectTotSpin : 1.0 ```

<b>Related Tags:</b> [hf_tensorfile](#tag_hf_tensorfile) [CorrTotSpin](#tag_CorrTotSpin)

<b>Related variable:</b> [Config](#Config)::DefectTotSpin

<a id="tag_CorrTotSpin"></a>
## CorrTotSpin
<b>Description:</b><br>
Spin number correction <br>
(Required when you use Hyperfine tensor file obtained by QUANTUM ESPRESSO) <br>

<b>Data Type:</b> Double

<b>Default:</b> 0.0

<b>Example:</b> ``` CorrTotSpin : 0.0 ```

<b>Related Tags:</b> [hf_tensorfile](#tag_hf_tensorfile) [DefectTotSpin](#tag_DefectTotSpin)

<b>Related variable:</b> [Config](#Config)::CorrTotSpin

<a id="tag_hf_readmode"></a>
## hf_readmode
<b>Description:</b><br>
Interaction tensor between qubit and bath spins can be determined by this option. <br>
If you want to use the hyperfine tensor obtained from DFT calculation, you should set this option. <br>

<b>Options:</b><br>
| Option | Description |
|--------|-------------------------------------|
| `0`    | No DFT hyperfine tensor <br> (Point-dipole approximation) |
| `1`    | Fermi-contact term + point-dipole tensor |
| `2`    | DFT dipolar tensor                       |
| `3`    | Fermi-contact term + DFT dipole tensor   |

<b>Changed since v1.0.0:</b> a value outside 0..3 is now refused. It used to be
accepted -- the range check could never fire -- and [hf_tensorfile](#tag_hf_tensorfile) was then never read and no point-dipole tensor was built either, so
the run completed with every hyperfine tensor left unset. Only an input that was
already meaningless is affected, but it now stops instead of returning a number. <br>

<b>Default:</b> 0

<b>Example:</b> ``` hf_readmode : 0 ```

<b>Related Tags:</b> [hf_tensorfile](#tag_hf_tensorfile) [hf_ignore_oor](#tag_hf_ignore_oor)

<b>Related variable:</b> [Config](#Config)::hf_readmode

<a id="tag_hf_tensorfile"></a>
## hf_tensorfile
<b>Description:</b><br>
The hyperfine tensor file, read when [hf_readmode](#tag_hf_readmode) is not 0. It is
written by `example/tools/SpinTensorReader` from a DFT run and holds, after a header, one
line per species with its g-factor and its default tensor (`etc____`, used for atoms outside
the range), the range of relative positions the file covers (`MinDif[A]`, `MaxDif[A]`), and one
row per atom: species, position relative to the defect (Å), the isotropic part (MHz) and the
nine tensor components (MHz). <br>
A bath spin is matched to a row by its position relative to the qubit. (`hf_cutoff`, read
beside these keys until 2026-09-24, never had an effect and is ignored if present.) A spin outside the
range gets a point-dipole tensor; a spin inside the range with no row is decided by
[hf_ignore_oor](#tag_hf_ignore_oor). The positions are in the bath frame; the frame of the
components is `hf_tensor_frame` of [coordinate_frame_rotation](#tag_coordinate_frame_rotation). <br>
<pre>
AtomName   r_1i_X  r_1i_Y  r_1i_Z  iso[MHz]  Axx Axy Axz Ayx Ayy Ayz Azx Azy Azz
10B  g-factor___:  0.600215 [No-unit]       one line per species
10B  etc____:      ---  ---  ---  0.00000    the tensor used outside the range
MinDif[A]:  -11.31  -6.53  0.00              the range of relative positions the file covers
MaxDif[A]:   10.06   6.53  0.00
B   -10.05626  -4.36351  0.00000  0.00000  0.06098 ... -0.04   one row per atom
</pre>

<b>Data Type:</b> String (file path)

<b>Default:</b> - (required when `hf_readmode` is not 0)

<b>Example:</b> ``` "hf_tensorfile" : "../bath/flat/hBN_VB_flat_ATensor_HSE_CellPBE_tag" ```

<b>Related Tags:</b> [hf_readmode](#tag_hf_readmode) [hf_ignore_oor](#tag_hf_ignore_oor)

<b>Related variable:</b> [Config](#Config)::hf_tensorfile

<a id="tag_hf_ignore_oor"></a>
## hf_ignore_oor
<b>Description:</b><br>
What to do with a bath spin whose position relative to the qubit lies inside the range the
hyperfine tensor file covers but has no row in the file. Spins outside the range always get
a point-dipole tensor. <br>

<b>Data Type:</b> Integer

<b>Options:</b><br>
| Option | Description |
|---|---|
| `0` | stop the run and print the position that was not found |
| `1` | give the spin a point-dipole tensor and go on |

<b>Default:</b> `0`

<b>Example:</b> ``` "hf_ignore_oor" : 0 ```

<b>Related Tags:</b> [hf_tensorfile](#tag_hf_tensorfile) [hf_readmode](#tag_hf_readmode)

<b>Related variable:</b> [Config](#Config)::hf_ignore_oor

<a id="tag_qd_readmode"></a>
## qd_readmode
<b>Description:</b><br>
Where the quadrupole (electric field gradient) tensors of the bath spins come from. With a
file, a bath spin inside the range the file covers takes the row at its position relative to
the defect (a spin inside the range with no row stops the run); a spin outside the range takes
the species default of the file (`etc`) in mode 2, or, in modes 3 and 4, the tensor of the
structure without the defect from [qd_tensorfile_woqubit](#tag_qd_tensorfile_woqubit), with
the relative position folded by [qd_cellpara](#tag_qd_cellpara) along every axis whose cell
parameter is not zero. <br>

<b>Data Type:</b> Integer

<b>Options:</b><br>
| Option | Description |
|---|---|
| `0` | no quadrupole interaction |
| `1` | removed (was: experimental constants); refused |
| `2` | one DFT file computed with the defect, [qd_tensorfile](#tag_qd_tensorfile) |
| `3` | strained structure: [qd_tensorfile](#tag_qd_tensorfile) with the defect and [qd_tensorfile_woqubit](#tag_qd_tensorfile_woqubit) without it; [qd_cellpara](#tag_qd_cellpara) required |
| `4` | bubble structure: the same two files and cell parameters |

<b>Default:</b> `0`

<b>Example:</b> ``` "qd_readmode" : 2 ```

<b>Related Tags:</b> [qd_tensorfile](#tag_qd_tensorfile) [qd_tensorfile_woqubit](#tag_qd_tensorfile_woqubit) [qd_cellpara](#tag_qd_cellpara)

<b>Related variable:</b> [Config](#Config)::qd_readmode

<a id="tag_qd_tensorfile"></a>
## qd_tensorfile
<b>Description:</b><br>
The quadrupole tensor file computed with the defect, read when [qd_readmode](#tag_qd_readmode)
is 2, 3 or 4. The same layout as [hf_tensorfile](#tag_hf_tensorfile), with the electric field
gradient components (Hartree/Bohr²) in place of the hyperfine ones and one `eQ___` line per
species (its quadrupole moment) in place of the g-factor. Written by
`example/tools/SpinTensorReader`. <br>

<b>Data Type:</b> String (file path)

<b>Default:</b> - (required when `qd_readmode` is 2, 3 or 4)

<b>Example:</b> ``` "qd_tensorfile" : "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE" ```

<b>Related Tags:</b> [qd_readmode](#tag_qd_readmode) [qd_tensorfile_woqubit](#tag_qd_tensorfile_woqubit) [qd_cellpara](#tag_qd_cellpara)

<b>Related variable:</b> [Config](#Config)::qd_tensorfile

<a id="tag_qd_cellpara"></a>
## qd_cellpara
<b>Description:</b><br>
Lattice parameter of the strained or bubble structure. <br>
This tag is required when you use the quadrupole tensor obtained from DFT calculation. <br>

<b>Data Type:</b> Double List

<b>Default:</b> -

<b>Unit:</b> angstrom

<b>Example:</b> ``` "qd_cellpara" : [1.0, 2.0, 3.0] ```

<b>Related Tags:</b> [qd_readmode](#tag_qd_readmode) [qd_tensorfile](#tag_qd_tensorfile) [qd_tensorfile_woqubit](#tag_qd_tensorfile_woqubit)

<b>Related variable:</b> [Config](#Config)::qd_cellpara

<a id="tag_qd_tensorfile_woqubit"></a>
## qd_tensorfile_woqubit
<b>Description:</b><br>
The quadrupole tensor file of the same structure computed without the defect, in the layout
of [qd_tensorfile](#tag_qd_tensorfile). Read when [qd_readmode](#tag_qd_readmode) is 3 or 4;
it supplies the tensors of the spins outside the range of the file with the defect. <br>

<b>Data Type:</b> String (file path)

<b>Default:</b> - (required when `qd_readmode` is 3 or 4)

<b>Example:</b> ``` "qd_tensorfile_woqubit" : "../bath/flat/hBN_VB_flat_QTensor_HSE_CellPBE" ```

<b>Related Tags:</b> [qd_readmode](#tag_qd_readmode) [qd_tensorfile](#tag_qd_tensorfile) [qd_cellpara](#tag_qd_cellpara)

<b>Related variable:</b> [Config](#Config)::qd_tensorfile_woqubit

<a id="tag_hfmedi"></a>
## hfmedi
<b>Description:</b><br>
CCE only. Add the hyperfine-mediated interaction between bath spins, the second-order
coupling through the qubit (N. Zhao et al., Nature Nanotechnology 6, 242 (2011), Supplementary, Eq. 9).
CCEX adds the full second-order correction, for each qubit branch α on its own:
H_med^(α) = Σ_{β≠α} E_αβ E_αβ† / (ω_α − ω_β), E_αβ = Σ_i ⟨α|S|β⟩·A_i·I_i, where α is the eigenstate of the
qubit Hamiltonian closest to the branch state and β runs over every other eigenstate (every intermediate state).
Expanded, it holds both orders of every pair (i, j) and the i = j terms of each bath spin; it is Hermitian.
(Before 2.0.0 it held one pair order and the other qubit level only.)
`example/tutorial/07_hBN_VB_natab/hyperfine_mediated_term/` runs with and without it. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "hfmedi" : true ```

<b>Related variable:</b> [Config](#Config)::hfmedi

<a id="tag_knight"></a>
## knight
<b>Description:</b><br>
CCE only. Add the Knight-field shift of the bath spins that the qubit induces under a
transverse magnetic field (same reference as [hfmedi](#tag_hfmedi), SI Eq. 14):
H_K^(α) = −γ_n B·(g − 1)·I = −(3|m_s| − 2) (γ_e/Δ) B·Ã·I for each bath spin, with Δ = 3/2 of the zz entry of `qzfs`,
Ã the hyperfine tensor with its z row set to zero, and |m_s| = ⟨S_z²⟩ of the branch's eigenstate.
It is zero when the field is along z. (Before 2.0.0 the term used the bare field B instead of −γ_n B and |m_s| = 1
for both branches.)
`example/tutorial/02_Diamond_NV_natab/knight_shift/` runs with and without it in a tilted field. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "knight" : true ```

<b>Related variable:</b> [Config](#Config)::knight
