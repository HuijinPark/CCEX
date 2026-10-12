# Backend tags
Where and how the cluster dynamics is computed, and what the run reports
about it: the CPU or CUDA backend, the gCCE propagation, the timing and the
per-cluster diagnostics.

- [backend](#tag_backend)
- [gpu_options](#tag_gpu_options) (min_dim, packed_max_dim, evolution_batch_size, evolution_batch_steps, cluster_batch_mb, cache_mb)
- [probe_memory_gb](#tag_probe_memory_gb)
- [propagator](#tag_propagator)
- [evolution](#tag_evolution)
- [evolution_update](#tag_evolution_update)
- [gcce_refinement](#tag_gcce_refinement)
- [gcce_refinement_max_sweeps](#tag_gcce_refinement_max_sweeps)
- [profile](#tag_profile)
- [cluster_contributions](#tag_cluster_contributions)
- [cluster_contribution_max_mb](#tag_cluster_contribution_max_mb)

<a id="tag_backend"></a>
## backend
<b>Description:</b><br>
Where the dynamics of each cluster is computed. <br>
`cuda` and `auto` need the CUDA build, `bin/main_gpu.out` from `bash do_compile.sh gpu`
(docs/methods/gpu.md); in a CPU-only binary `auto` runs everything on the CPU and `cuda`
is refused. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `cpu` | every cluster on the CPU |
| `cuda` | CUDA for every CCE, pCCE and gCCE cluster; the run stops if the device or the build is missing |
| `auto` | CUDA for the clusters whose Hilbert dimension is at least [gpu_options](#tag_gpu_options).`min_dim`, the CPU for the rest |

<b>Default:</b> `cpu`

<b>Example:</b> ``` "backend" : "auto" ```

<b>Related variable:</b> [Config](#Config)::backend

<a id="tag_gpu_options"></a>
## gpu_options
<b>Description:</b><br>
The GPU settings, in one block (since 2026-10-12, A26; the flat keys before it — `gpu_min_dim`, `gpu_batch_steps`,
`gpu_batch_clusters`, `gpu_cluster_batch_mb`, `gpu_cache_mb` — are refused with a message naming the new form, e.g.
`gpu_min_dim is now gpu_options.min_dim`). Every key is optional; an unknown key is refused. Without the block every
default below holds; the outputs are those of the earlier flat defaults except that the GPU build with `backend` `cuda`
or `auto` now runs the small cCCE clusters in the packed kernel (`packed_max_dim`, below). <br>

<b>Data Type:</b> Object

<b>Keys:</b><br>
| Key | Was | Type, default | Meaning |
|---|---|---|---|
| `min_dim` | `gpu_min_dim` | positive integer, `256` | with `backend` = `auto`, the smallest cluster Hilbert dimension that goes to the GPU; smaller clusters stay on the CPU. Ignored by `cpu` and `cuda` |
| `packed_max_dim` | (new) | `0`, `2` or `4`; `4` in `bin/main_gpu.out` with `backend` `cuda` or `auto` and `method` CCE, else `0` | the packed kernel (below): cCCE clusters of bath dimension up to this value on the GPU, many at once. `0` is off |
| `evolution_batch_size` | `gpu_batch_clusters` | positive integer, `8` | the largest number of consecutive clusters of one dimension evolved together on the GPU (`backend` `cuda`/`auto`); `1` evolves one cluster at a time |
| `evolution_batch_steps` | `gpu_batch_steps` | positive integer, `64` | the largest number of time points a GPU evolution processes together |
| `cluster_batch_mb` | `gpu_cluster_batch_mb` | positive integer, `64` (MiB per rank) | the evolution-buffer budget that limits how many clusters one GPU batch may hold |
| `cache_mb` | `gpu_cache_mb` | nonnegative integer, `64` (MiB per rank) | GPU memory that freed allocations are kept in for reuse; `0` disables the cache |

<b>The packed kernel (`packed_max_dim`).</b> The `cuda`/`auto` backends send a cluster to the GPU on its own, which
for the dimension-2 and 4 clusters of a P1 bath costs far more than the CPU (A12, A26 step 0: about 400 times slower).
With `packed_max_dim` = 4 (or 2), every cCCE cluster of bath dimension 2 and 4 (or 2) is instead computed by one GPU
thread — both eigensystems by Jacobi rotations and the whole time chain — thousands of clusters per kernel. Each rank
still builds its clusters' Hamiltonians and multiplies their coherences into the product as before; one process per node
(node rank 0, in a thread of its own) hands the packed clusters of all the node's ranks to the GPU through node shared
memory, so the node holds one CUDA context. Larger clusters follow `backend` as before. Single state and ensemble, every
pulse sequence with at most 4 distinct interval lengths (Hahn, CPMG, Ramsey, most `Manual` ones; others stay on the
CPU), `evolution_update` `exact` and `recursive` (the kernel follows the run's value). Needs `bin/main_gpu.out` and
`method` = `CCE` (refused otherwise). The results agree with the CPU within the GPU tolerance 1e-8 + 2e-6 |CPU|, not bit
for bit (the eigensolver differs); the standing check is `test_small_cluster_cuda` in `make gpu-test`.
<b>Default (user, 2026-10-12):</b> on (`4`) in `bin/main_gpu.out` with `backend` = `cuda` or `auto` and `method` CCE
(quantity coherence); off (`0`) with `backend` = `cpu`, with any other method, and always in the CPU build `bin/main.out`,
where nothing goes to the GPU and only an explicit nonzero value is refused. An explicit value works with any `backend`
(`cpu` included: then the packed clusters are the only GPU work). The run log prints the value and why
(`chosen by : default (backend=auto, method=CCE)`, `default (CPU build)`, ... or `given`). <br>
<b>When it pays.</b> Starting it costs a fixed part of a second (the CUDA context and page-locking the node shared memory,
done while the bath is read and the clusters are made; until it is up the ranks compute on the CPU), so a run of a few
seconds gains nothing. Measured on idm4-gpu (32 ranks + one RTX PRO 6000, job 691, the production P1 cCCE2 inputs, whole
run against the CPU binary on the same node, median of AB BA AB): CR60 σ4 30 states 70.7 → 44.2 s (0.62), its ensemble
6.5 → 3.9 s (0.59), CR0 30 states 13.6 → 11.4 s (0.84), but the CR0 ensemble, a 1.8 s run, 1.80 → 1.71 s (0.95);
(internal notes). In that job the 12 small fixture runs (under 1 s) finished before the GPU was
up and ran entirely on the CPU. <br>

<b>Example:</b> ``` "gpu_options" : {"min_dim": 256, "packed_max_dim": 4} ```

<b>Related variables:</b> [Config](#Config)::gpu_min_dim, gpu_packed_max_dim, gpu_batch_clusters, gpu_batch_steps,
gpu_cluster_batch_mb, gpu_cache_mb

<a id="tag_propagator"></a>
## propagator
<b>Description:</b><br>
gCCE only. How the free evolution of a cluster is propagated. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `eigen` | one eigendecomposition of the cluster Hamiltonian; every step is exp(-i H τ) = M diag(exp(-i λ τ)) M^† |
| `expm` | a matrix exponential per time step, the original CCEX path |

<b>Default:</b> `eigen`

<b>Example:</b> ``` "propagator" : "eigen" ```

<b>Related variable:</b> [Config](#Config)::propagator

<a id="tag_evolution"></a>
## evolution
<b>Description:</b><br>
gCCE only. What is propagated. `vector` needs a pure initial state (`nstate` > 0),
`propagator` = `eigen` and `method` = `gCCE`. When the key is absent the run picks `vector`
where those hold and `matrix` otherwise; an explicit `vector` that cannot run stops with an
error, and an explicit value with any other method is refused. The two agree to rounding but
are not bit-identical, so a converged series should not switch between them. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `vector` | propagate the state vector; matrix-vector products, O(n²) per step |
| `matrix` | propagate the density matrix ρ(t) = U ρ₀ U^†, O(n³) per step; works for every configuration |

<b>Default:</b> `vector` where it can run, else `matrix`

<b>Example:</b> ``` "evolution" : "matrix" ```

<b>Related variable:</b> [Config](#Config)::evolution

<a id="tag_evolution_update"></a>
## evolution_update
<b>Description:</b><br>
How the phases exp(-i E τ) of the eigenvalues are formed at each time point of the uniform grid (`deltat`). `exact`
computes them with sin and cos at every point. `recursive` does so at points 0, 16, 32, … and in between multiplies the
previous point's phases by one step's factor exp(-i E f `deltat`) (f the interval's fraction of the sequence), so the
rounding of at most 15 products accumulates before the next exact point (A16 (b), 2026-10-11). <br>
`recursive` exists on the CPU time evolution of cCCE (`method` CCE, pCCE, dsj, itb or dsjitb; one
sampled state, `nstate` > 0, and the ensemble, `nstate` = 0, since 2026-10-11 (A18)), in the packed kernel of
`bin/main_gpu.out` (`gpu_options.packed_max_dim`, cCCE clusters of bath dimension 2 and 4) and on the CPU time evolution
of gCCE with `propagator` = `eigen` (`evolution` = `vector`, and `matrix` since A18). Every pulse sequence these paths
run is covered: the ensemble cCCE carries over the phases of its τ = t/2 (Hahn), t (Ramsey), t/(2 `npulse`) and 2τ
(CPMG) or of each `Manual` interval, the others those of each interval of the sequence. The GPU large-cluster path
(`bin/main_gpu.out`: every other cluster with `backend` = `cuda`, the clusters of `gpu_options.min_dim` and above with
`auto`) computes every phase with sin and cos, as `exact`. <br>
Without the key (since 2026-10-12, A29; per path in the GPU build since A36) each path takes `recursive` where it has
it and `exact` elsewhere, without a message: with `backend` = `cpu`, and in the CPU build `bin/main.out` with `cpu` or
`auto` (all on the CPU there), the value is `recursive` (gCCE `expm`: `exact`); in `bin/main_gpu.out` with `auto` the
clusters on the CPU (below `min_dim`) and the packed kernel run `recursive` and the clusters on the GPU `exact`; with
`cuda` the packed kernel runs `recursive` and every other cluster `exact` on the GPU (with the packed kernel off, as for
gCCE and pCCE, the value is `exact`). The run log prints the value and why, per path in the GPU build
(`chosen by : default (supported path)`, `default (clusters on the CPU and packed kernel: recursive; clusters on the GPU
(backend=auto): exact)`, `default (backend=cuda: every cluster on the GPU, exact there)`, … or `given`). Until A36
(2026-10-12) `backend` `cuda` or `auto` made the whole run `exact`; `"evolution_update": "exact"` gives those runs, and
the runs before A29, bit for bit. <br>
A given `recursive` is refused rather than run as `exact` where a cluster may reach a path without it: gCCE with
`propagator` = `expm` (no phases to carry over), and in `bin/main_gpu.out` `backend` = `cuda` or `auto` with a method
the GPU takes (CCE, pCCE, gCCE `eigen`) — leave the key out there to get `recursive` on the CPU clusters and the packed
kernel. `backend` = `cuda` in the CPU build is refused in any case. <br>
Measured on the production P1 inputs (CR0 and CR60 σ4, cCCE2 and gCCE2, 30 states; idm job 17217, 64 ranks): the run
took 0.63–0.76 of the time of the build before A16 (geometric mean 0.70; this includes the byte-identical A16 (a)
changes, and against those alone the same job gives about 0.73–0.76 for cCCE2 and 0.82–0.84 for gCCE2). The `wiDiv`
curves moved by at most 2.2e-10 for cCCE2 (the outputs have 10 decimals) and 1.3e-8 for gCCE2 at points with |L| ≤ 1
((internal notes)). The results are not bit-identical to `exact`, so a converged
series should not switch between the two. Must be a JSON string. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `exact` | sin and cos at every time point |
| `recursive` | sin and cos every 16th time point, the previous point's phases times one step's factor in between |

<b>Default:</b> `recursive` where it exists (CPU cCCE, the packed kernel, CPU gCCE with `propagator` = `eigen`), `exact` elsewhere (the GPU large-cluster path, gCCE `expm`); per path in `bin/main_gpu.out` since A36 (the whole run `exact` with `backend` `cuda`/`auto` until then; `exact` everywhere until A29, 2026-10-12)

<b>Example:</b> ``` "evolution_update" : "exact" ```

<b>Related variable:</b> [Config](#Config)::evolution_update

<a id="tag_probe_memory_gb"></a>
## probe_memory_gb
<b>Description:</b><br>
The probe step of the adaptive clustering (`cluster_mode` = `adaptive`): the memory, in GB per node, that the probes
running at once may hold. Since 2026-10-02 the probes of a bath state go to the ranks in turn and the ranks do not wait
for each other between states, so the large probes of several states can run together; a probe waits while the probes
already running on its node would pass this budget, and a probe larger than the budget runs alone. A probe is estimated
on the path it takes: on the CPU 6 dense dim × dim complex matrices (16 bytes each), one more per pulse for gCCE; on the
GPU what the device allocates for the time batch and evolution; with `savemode` = `allfull` also every time's full
density matrix. Probes under 64 MiB do not count. 0 takes 80 % of the memory the job may use on the node (its RAM, or
less if a cgroup limits the job) and, when the run may put clusters on the GPU, of the smallest GPU. The run prints the budget (`[probe-budget]`) and how long
each rank waited. <br>

<b>Data Type:</b> Double (nonnegative)

<b>Default:</b> `0` (automatic)

<b>Example:</b> ``` "probe_memory_gb" : 60 ```

<b>Related variable:</b> [Config](#Config)::probe_memory_gb

<a id="tag_gcce_refinement"></a>
## gcce_refinement
<b>Description:</b><br>
gCCE with `propagator` = `eigen`: refine the eigenvalues and eigenvectors the solver returns by
compensated correction sweeps before they are used. Must be a JSON boolean. <br>
Off by default since 2026-10-02: on the NV–¹³C gCCE4 case the refinement was made for, and on P1 gCCE3, it
changed the error-corrected final coherence by at most 6.7e-7 (idm job 16869,
(internal notes); the raw `wiDiv` differs by more at points where it diverges with either
setting — up to 12.7 on NV–¹³C and 2.75e11 on P1 gCCE3 — so the predeclared 1e-6 bound on every raw point was not met), while it took most of the setup time of a gCCE
cluster (dim 768: about 4 s of 4.8 s on the CPU). Turn it on where the bath energies are very small beside the
qubit's and the curve is read to better than 1e-6. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false` (`true` until 2026-10-02)

<b>Example:</b> ``` "gcce_refinement" : true ```

<b>Related variable:</b> [Config](#Config)::gcce_refinement

<a id="tag_gcce_refinement_max_sweeps"></a>
## gcce_refinement_max_sweeps
<b>Description:</b><br>
The largest number of refinement sweeps per cluster. A cluster that has not converged
within it stops the run with an error naming its dimension. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `32`

<b>Example:</b> ``` "gcce_refinement_max_sweeps" : 32 ```

<b>Related variable:</b> [Config](#Config)::gcce_refinement_max_sweeps

<a id="tag_profile"></a>
## profile
<b>Description:</b><br>
Print timings at the end of the run. Every rank writes one `[profile]` line per method,
backend and cluster dimension with the number of clusters and the seconds spent building,
setting up and evolving them; with CUDA also `[gpu-batch]` and `[gpu-memory]` lines. In a plain CCE/gCCE run
rank 0 also prints one `[phase]` line per rank: the wall seconds of the state loop split into state preparation,
cluster creation, cluster calculation, product assembly, freeing, waiting at the barriers, the reduction and saving.
Must be a JSON boolean. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "profile" : true ```

<b>Related variable:</b> [Config](#Config)::profile

<a id="tag_cluster_contributions"></a>
## cluster_contributions
<b>Description:</b><br>
Record the contribution of every cluster to the coherence over all time points
(cCCE and gCCE, `quantity` = `coherence`). For each bath state the run writes
`<outfile>.contributions.state-<s>.tsv` (one row per cluster: sites, bundles, production
power, dimension, the largest |log L̃|, undefined points, worst step) and once
`<outfile>.contributions.sites.tsv` (site, name, spin, gyromagnetic ratio, position).
A diagnostic; the raw table per rank is bounded by
[cluster_contribution_max_mb](#tag_cluster_contribution_max_mb). Must be a JSON boolean. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "cluster_contributions" : true ```

<b>Related variable:</b> [Config](#Config)::cluster_contributions

<a id="tag_cluster_contribution_max_mb"></a>
## cluster_contribution_max_mb
<b>Description:</b><br>
The memory per rank the contribution table of `cluster_contributions` may take. The run
stops before the propagation if the table would exceed it. <br>

<b>Data Type:</b> Integer (positive)

<b>Unit:</b> MiB per rank

<b>Default:</b> `256`

<b>Example:</b> ``` "cluster_contribution_max_mb" : 512 ```

<b>Related variable:</b> [Config](#Config)::cluster_contribution_max_mb
