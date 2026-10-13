# Cluster-related tags
- [order](InputFile_General.md#tag_order)        
- [method](InputFile_General.md#tag_method)       
- [addsubclus](#tag_addsubclus)   
- [nk](#tag_nk)           
- [sK](#tag_sK)           
- [max_trial](#tag_max_trial)    
- [max_iter](#tag_max_iter)     
- [kmeans_pp](#tag_kmeans_pp)    
- [iter_detail](#tag_iter_detail)
- [cluster_mode](#tag_cluster_mode)
- [bundle_max_dim](#tag_bundle_max_dim)
- [bundle_merge_cutoff](#tag_bundle_merge_cutoff)
- [bundle_connections](#tag_bundle_connections)
- [bundle_connection_budget](#tag_bundle_connection_budget)
- [expansion_max_clusters](#tag_expansion_max_clusters)
- [expansion_max_work](#tag_expansion_max_work)
- [adaptive_work_file](#tag_adaptive_work_file)
- [Adaptive](#tag_Adaptive)
- [pcce_partition](#tag_pcce_partition)
- [pcce_connections](#tag_pcce_connections)
- [pcce_cache_mode](#tag_pcce_cache_mode)
- [pcce_cache_file](#tag_pcce_cache_file)
- [pcce_detail_limit](#tag_pcce_detail_limit)

<a id="tag_addsubclus"></a>
## addsubclus
<b>Description:</b><br>
For CCE simulations with order ≥ 3, enabling this option allows all sub-clusters of the highest order to be included. <br>
This leads to more stable and robust results, and may also reduce computation time. <br>
<b>Data Type:</b> Bool

<b>Default:</b> true

<b>Example:</b> ``` addsubclus : true ```

<b>Related variable:</b> [Cluster](#Cluster)::addsubclus

<a id="tag_nk"></a>
## nk
<b>Description:</b><br>
Number of clusters for each CCE order. <br>
This is useful for CCE simulations with order ≥ 3, leading tons of clusters <br>
The list has exactly `order + 1` entries and the index is the cluster order: `nk[n]` limits the clusters of order n
(entry 0 is the order-0 slot and is kept, usually 0). A 0 means no limit for that order, so all zeros include every
cluster. A list of another length is refused. <br>
<b>Data Type:</b> Int List (length `order + 1`)

<b>Default:</b> all entries 0 (no limit)

<b>Example:</b> ``` nk : [0, 0, 0, 30000] ``` (order 3: at most 30000 clusters of three spins)

<b>Related variable:</b> [Cluster](#Cluster)::nk

<a id="tag_sK"></a>
## sK
<b>Description:</b><br>
pCCE: the number of bath spins in one partition. The bath is divided into partitions of
`sK` spins by a k-means search, and a pCCE cluster of order n is n partitions, n × `sK`
spins. Required with `method` = `pCCE`; refused with `cluster_mode` = `adaptive`.
`example/tutorial/04_Diamond_NV_P1Bath/SpinEcho_pCCE2_1ppm/` runs `sK` 1 to 4. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> - (required)

<b>Example:</b> ``` "sK" : 3 ```

<b>Related variable:</b> [Cluster](#Cluster)::sK

<a id="tag_max_trial"></a>
## max_trial
<b>Description:</b><br>
pCCE: the number of partition trials, each a k-means run from its own starting centres;
the assignment with the best score is kept. With [pcce_partition](#tag_pcce_partition) =
`deterministic` the trials are split over the ranks in a fixed way, so the partition follows
`seed`. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `30000`

<b>Example:</b> ``` "max_trial" : 1000 ```

<b>Related variable:</b> [Cluster](#Cluster)::max_trial

<a id="tag_max_iter"></a>
## max_iter
<b>Description:</b><br>
pCCE: the largest number of k-means iterations in one partition trial. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `30000`

<b>Example:</b> ``` "max_iter" : 1000 ```

<b>Related variable:</b> [Cluster](#Cluster)::max_iter

<a id="tag_kmeans_pp"></a>
## kmeans_pp
<b>Description:</b><br>
pCCE: start every k-means trial from k-means++ centres instead of random ones. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `true`

<b>Example:</b> ``` "kmeans_pp" : true ```

<b>Related variable:</b> [Cluster](#Cluster)::kmeans_pp

<a id="tag_iter_detail"></a>
## iter_detail
<b>Description:</b><br>
pCCE: rank 0 prints the details of the partition trials, at most
[pcce_detail_limit](#tag_pcce_detail_limit) of them. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "iter_detail" : true ```

<b>Related variable:</b> [Cluster](#Cluster)::iter_detail

<a id="tag_cluster_mode"></a>
## cluster_mode
<b>Description:</b><br>
How the bath is divided into clusters. With `adaptive`, strongly coupled spins are tied into
bundles first; [order](InputFile_General.md#tag_order) then counts bundles per cluster and
[bundle_max_dim](#tag_bundle_max_dim) (required) caps the Hilbert dimension of one bundle. With
[rdip](InputFile_General.md#tag_rdip) given, bundles are connected by distance; with `rdip` omitted,
connections are chosen by a coherence check and the weak ones are dropped. The fixed-size
search keys `sK`, `max_trial`, `max_iter`, `kmeans_pp`, `iter_detail`, `pcce_partition` and
`pcce_cache_file` are refused in this mode. Adding the [Adaptive](#tag_Adaptive) block turns
the adaptive refinement on top of it. The words are defined in
[Clustering terminology](../methods/clustering-terminology.md), the method in [Adaptive cluster mode](../methods/adaptive-cluster-mode.md). <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `standard` | clusters of up to `order` spins within `rdip`, as before |
| `adaptive` | bundle clustering: clusters of up to `order` bundles |

<b>Default:</b> `standard`

<b>Example:</b> ``` "cluster_mode" : "adaptive" ```

<b>Related variable:</b> [Config](#Config)::cluster_mode

<a id="tag_bundle_max_dim"></a>
## bundle_max_dim
<b>Description:</b><br>
The largest Hilbert dimension one bundle may have. Required with `cluster_mode` = `adaptive`.
The dimension of a cluster follows as the product of its bundle dimensions (gCCE: times the
qubit dimension); a cluster that cannot be handled is refused during the run. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> - (required)

<b>Example:</b> ``` "bundle_max_dim" : 8 ```

<b>Related variable:</b> [Config](#Config)::bundle_max_dim

<a id="tag_bundle_merge_cutoff"></a>
## bundle_merge_cutoff
<b>Description:</b><br>
The smallest coupling strength for which two groups are merged into one bundle. <br>

<b>Data Type:</b> Double (nonnegative)

<b>Unit:</b> rad/ms

<b>Default:</b> `0`

<b>Example:</b> ``` "bundle_merge_cutoff" : 1.0 ```

<b>Related variable:</b> [Config](#Config)::bundle_merge_cutoff

<a id="tag_bundle_connections"></a>
## bundle_connections
<b>Description:</b><br>
How bundles are connected into clusters when `cluster_mode` = `adaptive`. Omit it together
with `rdip` to get the coherence-checked connections (the run log then says `checked`; that
word is not a value to write). <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `distance` | pairs within `rdip`; `rdip` is required |
| `auto` | the legacy budget mode: weak pairs are left out up to [bundle_connection_budget](#tag_bundle_connection_budget), which is required |

<b>Default:</b> `distance`

<b>Example:</b> ``` "bundle_connections" : "auto" ```

<b>Related variable:</b> [Config](#Config)::bundle_connections

<a id="tag_bundle_connection_budget"></a>
## bundle_connection_budget
<b>Description:</b><br>
With `bundle_connections` = `auto`, the allowed score of the left-out pairs: the largest
time times the sum of their coupling weights must stay below it. <br>

<b>Data Type:</b> Double (nonnegative)

<b>Default:</b> -

<b>Example:</b> ``` "bundle_connection_budget" : 0.01 ```

<b>Related variable:</b> [Config](#Config)::bundle_connection_budget

<a id="tag_expansion_max_clusters"></a>
## expansion_max_clusters
<b>Description:</b><br>
The largest number of clusters the expansion of a bundle partition may list.
`order` may not exceed it. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `200000`

<b>Example:</b> ``` "expansion_max_clusters" : 200000 ```

<b>Related variable:</b> [Config](#Config)::expansion_max_clusters

<a id="tag_expansion_max_work"></a>
## expansion_max_work
<b>Description:</b><br>
The work limit of the expansion and of the dimension search of a bundle partition. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `20000000`

<b>Example:</b> ``` "expansion_max_work" : 20000000 ```

<b>Related variable:</b> [Config](#Config)::expansion_max_work

<a id="tag_adaptive_work_file"></a>
## adaptive_work_file
<b>Description:</b><br>
The plan/probe protocol file the adaptive refinement driver writes for each step it runs.
The driver sets it itself; a hand-written input does not need it. <br>

<b>Data Type:</b> String (path)

<b>Default:</b> -

<b>Example:</b> ``` "adaptive_work_file" : "./work.json" ```

<b>Related variable:</b> [Config](#Config)::adaptive_work_file

<a id="tag_Adaptive"></a>
## Adaptive
<b>Description:</b><br>
Turns on the adaptive refinement: instead of one calculation the run repeats calculating,
measuring the error and correcting the bundles until `atol` and `rtol` are met or `max_rounds`
is reached. Needs `cluster_mode` = `adaptive`, `bundle_max_dim` and an explicit
[seed](InputFile_General.md#tag_seed). The keys are the driver's options with `_` for `-`;
`output`, `atol` and `rtol` are required, and every step goes into `output` (one folder per
step and `summary.json`). The full key table with defaults is in [the `Adaptive` block](../methods/adaptive-input-block.md).
The block was called `Refine` until 2026-09-23; that name is refused. <br>

<b>Data Type:</b> Object

<b>Options:</b><br>
| Option | Description |
|---|---|
| `policy` | `merge`, the only policy and the default. `error-gated`, `budgeted`, `auto-partition` and `cost-preserving` were removed on 2026-10-10 and `auto-order` on 2026-09-25; naming one is refused |
| `output` | folder for the steps and the summary (required) |
| `atol`, `rtol` | absolute and relative tolerance of the coherence (required) |
| `max_rounds` | the largest number of rounds, default 8 |

<b>Default:</b> - (absent: one calculation)

<b>Example:</b> ``` "Adaptive" : { "policy" : "merge", "output" : "./results/adaptive", "atol" : 1e-3, "rtol" : 0, "max_rounds" : 4 } ```


<a id="tag_pcce_partition"></a>
## pcce_partition
<b>Description:</b><br>
pCCE: how the partition of the bath is searched. Refused with `cluster_mode` = `adaptive`. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `deterministic` | the partition follows `seed`; trials are distributed deterministically over the ranks |
| `legacy` | the original random-number and trial distribution |

<b>Default:</b> `deterministic`

<b>Example:</b> ``` "pcce_partition" : "deterministic" ```

<b>Related variable:</b> [Config](#Config)::pcce_partition

<a id="tag_pcce_connections"></a>
## pcce_connections
<b>Description:</b><br>
pCCE: which distance decides that two partitions are connected, so that they can be
in one cluster of order 2 or more. `rdip` and `rdipcut` bound that distance under
either rule. Refused with `cluster_mode` = `adaptive`, whose bundles are connected by
`bundle_connections`. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `centers` | the distance between the two k-means centres, the centroids of the partitions. The strength that orders the clusters is the point-dipole T<sub>zz</sub> between the centroids. |
| `spins` | any spin of the one partition within (`rdipcut`, `rdip`) of any spin of the other, the rule of the adaptive bundles. The strength is the sum of \|T<sub>zz</sub>\| over the spin pairs of the two partitions. |

With `sK` = 1 the two rules are the same. Under `centers` two partitions whose spins
are close can stay unconnected when their centroids are far apart, and the more so
the larger `sK`; `spins` connects them. The rank-0 log says which rule ran and how
many partition pairs it connected. The comparison of the two rules on the P1 bath of tutorial 04
is (internal notes): at that tutorial's `rdip` (620 Å) the two rules
differ less than pCCE2 differs from CCE4, and `spins` costs 27-65 % more; at an `rdip` comparable
to the partition size, `spins` is the closer of the two to CCE4.

<b>Default:</b> `centers` (the rule of every earlier version)

<b>Example:</b> ``` "pcce_connections" : "spins" ```

<b>Related variable:</b> [Config](#Config)::pcce_connections

<a id="tag_pcce_cache_mode"></a>
## pcce_cache_mode
<b>Description:</b><br>
pCCE with `pcce_partition` = `deterministic`: reuse a stored partition instead of searching again. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `off` | search every time |
| `readwrite` | reuse a matching stored partition, otherwise search and store it |
| `read` | require a matching stored partition; stop if there is none |

<b>Default:</b> `off`

<b>Example:</b> ``` "pcce_cache_mode" : "readwrite" ```

<b>Related variable:</b> [Config](#Config)::pcce_cache_mode

<a id="tag_pcce_cache_file"></a>
## pcce_cache_file
<b>Description:</b><br>
The partition cache file of `pcce_cache_mode`, relative to the working directory or absolute;
its directory must exist. Refused with `cluster_mode` = `adaptive`. <br>

<b>Data Type:</b> String (path)

<b>Default:</b> -

<b>Example:</b> ``` "pcce_cache_file" : "./results/partition.cache" ```

<b>Related variable:</b> [Config](#Config)::pcce_cache_file

<a id="tag_pcce_detail_limit"></a>
## pcce_detail_limit
<b>Description:</b><br>
pCCE with [iter_detail](#tag_iter_detail): rank 0 prints the details of at most this many
trials; `0` prints none. <br>

<b>Data Type:</b> Integer (nonnegative)

<b>Default:</b> `32`

<b>Example:</b> ``` "pcce_detail_limit" : 8 ```

<b>Related variable:</b> [Config](#Config)::pcce_detail_limit
