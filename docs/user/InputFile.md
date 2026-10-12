# Input File
<a id="InputFile"></a>

Every key of the input file, in seven groups. Each group is a page of its own, and
every key in the lists below links to its entry there.

<a id="HowToWrite"></a>
## How to write the input file

The input file consists of key-value pairs where each key represents a specific setting.
Each parameter should be written in the following format:

```
{
    parameter_name_1 : value_1,
    parameter_name_2 : value_2,
                :           
    parameter_name_N : value_N
}
```

A `#` or `!` outside a JSON string starts a comment that runs to the end of the line, so a whole line or the rest of a
line after a value can be a comment. Inside a string (`"#note"`, `"outfile": "run#1"`) they are ordinary characters
(until 2026-10-08 the line was cut there too and such an input was refused).

```
# a comment line
{
    "method" : "CCE",       # a comment after a value
    ! a comment line of the other kind
    "outfile" : "run#1"
}
```

A key that is left out takes its default. A key that is present has to have the type its
entry gives: `"knight": "true"`, a quoted number, `"nstep": 10.5` or a string inside a
number array stops the run with the key named. Until 2026-10-02 such a value was taken as
missing and the default was used without a word. Two older forms are still read: a Bool
may be written `0` or `1`, and a qubit's `detuning` may be written `[x]`.

A key that no part of the run reads is reported (`Warning: input key ... is not read`), and
when it is one or two letters away from a key that was looked up and is absent from the same
block, the run stops: `"nstates"` would leave `nstate` at its default. In the `qubitfile` form
the branch keys are `qalphams`/`qbetams`; a top-level `alphams` there was ignored until
2026-10-02 and is now refused. When both `qubitfile` and a `Qubit` block are given, the
`qubitfile` is used and the block is not read (reported as unread).

<a id="GeneralTag"></a>
## General tags
The keys of this group are described in [InputFile_General.md](InputFile_General.md).

Tags in this section is about simulatation parameters such as file path, numerical parameters.

- [method](InputFile_General.md#tag_method)       
- [quantity](InputFile_General.md#tag_quantity)     
- [order](InputFile_General.md#tag_order)        
- [bfield](InputFile_General.md#tag_bfield)       
- [rbath](InputFile_General.md#tag_rbath)        
- [rdip](InputFile_General.md#tag_rdip)         
- [deltat](InputFile_General.md#tag_deltat)       
- [nstep](InputFile_General.md#tag_nstep)        
- [rbathcut](InputFile_General.md#tag_rbathcut)     
- [rdipcut](InputFile_General.md#tag_rdipcut)       
- [rdip_ladder](InputFile_General.md#tag_rdip_ladder)
- [nstate](InputFile_General.md#tag_nstate)        
- [vector_mean_field](InputFile_General.md#tag_vector_mean_field)
- [seed](InputFile_General.md#tag_seed)       
- [qubitfile](InputFile_General.md#tag_qubitfile)   
- [gyrofile](InputFile_General.md#tag_gyrofile)      
- [bathfile](InputFile_General.md#tag_bathfile)      
- [bathadjust](InputFile_General.md#tag_bathadjust)
- [coordinate_frame_rotation](InputFile_General.md#tag_coordinate_frame_rotation)
- [defect_axis_reference](InputFile_General.md#tag_defect_axis_reference)
- [avaaxfile](InputFile_General.md#tag_avaaxfile)
- [statefile](InputFile_General.md#tag_statefile)                     
- [exstatefile](InputFile_General.md#tag_exstatefile)                  
- [DefectTotSpin](InputFile_General.md#tag_DefectTotSpin)                
- [CorrTotSpin](InputFile_General.md#tag_CorrTotSpin)                  
- [hf_readmode](InputFile_General.md#tag_hf_readmode)                       
- [hf_tensorfile](InputFile_General.md#tag_hf_tensorfile)                
- [hf_ignore_oor](InputFile_General.md#tag_hf_ignore_oor)                
- [qd_readmode](InputFile_General.md#tag_qd_readmode)                  
- [qd_tensorfile](InputFile_General.md#tag_qd_tensorfile)                
- [qd_cellpara](InputFile_General.md#tag_qd_cellpara)                     
- [qd_tensorfile_woqubit](InputFile_General.md#tag_qd_tensorfile_woqubit)
- [hfmedi](InputFile_General.md#tag_hfmedi)
- [knight](InputFile_General.md#tag_knight)

<a id="QubitArrayTags"></a>
## Qubit-related tags
The keys of this group are described in [InputFile_Qubit.md](InputFile_Qubit.md).

This tags are used for adjusting Qubit Hamiltonian and qubit inital state.
The format is the following:

<b>Case1.</b> [qubitfile](InputFile_General.md#tag_qubitfile) tag <br>
  : Only single-qubit systems are supported. <br>
    Note that this tag takes priority over the [Qubit](InputFile_Qubit.md#tag_Qubit) tag. <br>

<b>Example:</b> ``` "qubitfile"  : "file-path" ```

<b>Case2.</b> [Qubit](InputFile_Qubit.md#tag_Qubit) tag <br>
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
- [qzfs](InputFile_Qubit.md#tag_qzfs)
- [qspin](InputFile_Qubit.md#tag_qspin)
- [qalphams](InputFile_Qubit.md#tag_qalphams) 
- [qbetams](InputFile_Qubit.md#tag_qbetams)
- [Qubit](InputFile_Qubit.md#tag_Qubit) 
- [nqubit](InputFile_Qubit.md#tag_nqubit)
- [qubit](InputFile_Qubit.md#tag_qubit)
- [name](InputFile_Qubit.md#tag_name)         
- [spin](InputFile_Qubit.md#tag_spin)
- [gyro](InputFile_Qubit.md#tag_gyro)         
- [xyz](InputFile_Qubit.md#tag_xyz)          
- [detuning](InputFile_Qubit.md#tag_detuning)     
- [alphams](InputFile_Qubit.md#tag_alphams) 
- [betams](InputFile_Qubit.md#tag_betams)
- [intmap](InputFile_Qubit.md#tag_intmap)
- [between](InputFile_Qubit.md#tag_between)      
- [tensor](InputFile_Qubit.md#tag_tensor)       
- [psia](InputFile_Qubit.md#tag_psia)  
- [psib](InputFile_Qubit.md#tag_psib)  
- [psi0](InputFile_Qubit.md#tag_psi0)
- [overhaus](InputFile_Qubit.md#tag_overhaus)     
- [alphaidx](InputFile_Qubit.md#tag_alphaidx) 
- [betaidx](InputFile_Qubit.md#tag_betaidx)

<a id="ClusterTags"></a>
## Cluster-related tags
The keys of this group are described in [InputFile_Cluster.md](InputFile_Cluster.md).

- [order](InputFile_General.md#tag_order)        
- [method](InputFile_General.md#tag_method)       
- [addsubclus](InputFile_Cluster.md#tag_addsubclus)   
- [nk](InputFile_Cluster.md#tag_nk)           
- [sK](InputFile_Cluster.md#tag_sK)           
- [max_trial](InputFile_Cluster.md#tag_max_trial)    
- [max_iter](InputFile_Cluster.md#tag_max_iter)     
- [kmeans_pp](InputFile_Cluster.md#tag_kmeans_pp)    
- [iter_detail](InputFile_Cluster.md#tag_iter_detail)
- [cluster_mode](InputFile_Cluster.md#tag_cluster_mode)
- [bundle_max_dim](InputFile_Cluster.md#tag_bundle_max_dim)
- [bundle_merge_cutoff](InputFile_Cluster.md#tag_bundle_merge_cutoff)
- [bundle_connections](InputFile_Cluster.md#tag_bundle_connections)
- [bundle_connection_budget](InputFile_Cluster.md#tag_bundle_connection_budget)
- [expansion_max_clusters](InputFile_Cluster.md#tag_expansion_max_clusters)
- [expansion_max_work](InputFile_Cluster.md#tag_expansion_max_work)
- [adaptive_work_file](InputFile_Cluster.md#tag_adaptive_work_file)
- [Adaptive](InputFile_Cluster.md#tag_Adaptive)
- [pcce_partition](InputFile_Cluster.md#tag_pcce_partition)
- [pcce_connections](InputFile_Cluster.md#tag_pcce_connections)
- [pcce_cache_mode](InputFile_Cluster.md#tag_pcce_cache_mode)
- [pcce_cache_file](InputFile_Cluster.md#tag_pcce_cache_file)
- [pcce_detail_limit](InputFile_Cluster.md#tag_pcce_detail_limit)

<a id="BackendTags"></a>
## Backend tags
The keys of this group are described in [InputFile_Backend.md](InputFile_Backend.md).

Where and how the cluster dynamics is computed, and what the run reports about it.

- [backend](InputFile_Backend.md#tag_backend)
- [gpu_min_dim](InputFile_Backend.md#tag_gpu_min_dim)
- [gpu_batch_steps](InputFile_Backend.md#tag_gpu_batch_steps)
- [gpu_batch_clusters](InputFile_Backend.md#tag_gpu_batch_clusters)
- [gpu_cluster_batch_mb](InputFile_Backend.md#tag_gpu_cluster_batch_mb)
- [probe_memory_gb](InputFile_Backend.md#tag_probe_memory_gb)
- [gpu_cache_mb](InputFile_Backend.md#tag_gpu_cache_mb)
- [propagator](InputFile_Backend.md#tag_propagator)
- [evolution](InputFile_Backend.md#tag_evolution)
- [evolution_update](InputFile_Backend.md#tag_evolution_update)
- [gcce_refinement](InputFile_Backend.md#tag_gcce_refinement)
- [gcce_refinement_max_sweeps](InputFile_Backend.md#tag_gcce_refinement_max_sweeps)
- [profile](InputFile_Backend.md#tag_profile)
- [cluster_contributions](InputFile_Backend.md#tag_cluster_contributions)
- [cluster_contribution_max_mb](InputFile_Backend.md#tag_cluster_contribution_max_mb)

<a id="PulseTags"></a>
## Pulse-related tags
The keys of this group are described in [InputFile_Pulse.md](InputFile_Pulse.md).

- [npulse](InputFile_Pulse.md#tag_npulse)       
- [pulsename](InputFile_Pulse.md#tag_pulsename)    
- [sequence](InputFile_Pulse.md#tag_sequence)
- [pulseiter](InputFile_Pulse.md#tag_pulseiter)
- [pulse_time](InputFile_Pulse.md#tag_pulse_time)
- [detuning_factor](InputFile_Pulse.md#tag_detuning_factor)

<a id="DefectArrayTags"></a>
## DefectArray-related tags
The keys of this group are described in [InputFile_Defect.md](InputFile_Defect.md).

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

- [Defect](InputFile_Defect.md#tag_Defect)       
- [dfname](InputFile_Defect.md#tag_df_dfname)       
- [apprx](InputFile_Defect.md#tag_df_apprx)        
- [naddspin](InputFile_Defect.md#tag_df_naddspin)     
- [navaax](InputFile_Defect.md#tag_df_navaax)       
- [axis](InputFile_Defect.md#tag_df_axis)
- [coordinate_frame](InputFile_Defect.md#tag_df_coordinate_frame)
- [types](InputFile_Defect.md#tag_df_types)        
- [spins](InputFile_Defect.md#tag_df_spins)        
- [gyros](InputFile_Defect.md#tag_df_gyros)        
- [eqs](InputFile_Defect.md#tag_df_eqs)          
- [rxyzs](InputFile_Defect.md#tag_df_rxyzs)        
- [hypf](InputFile_Defect.md#tag_df_hypf)         
- [efg](InputFile_Defect.md#tag_df_efg)          
- [zfs](InputFile_Defect.md#tag_df_zfs)          
- [detuning](InputFile_Defect.md#tag_df_detuning)

<a id="OutputTags"></a>
## Output tags
The keys of this group are described in [InputFile_Output.md](InputFile_Output.md).

- [savemode](InputFile_Output.md#tag_savemode)
- [outfile](InputFile_Output.md#tag_outfile)
- [save_nodiv](InputFile_Output.md#tag_save_nodiv)
