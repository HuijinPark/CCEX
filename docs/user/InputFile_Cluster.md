# Cluster tags
How the bath is cut into clusters: the `cluster` block - the [mode](#tag_cluster_mode), the
`selection` every mode reads, and the blocks of the `pcce` and `bundle` modes and of the adaptive refinement. `order`
and `method` stay at the top. Older key names and places are in the table of [InputFile.md](InputFile.md#OldNames).

```
"cluster": {"mode": "pcce", "selection": {"add_subclusters": true, "max_per_order": [0, 0, 0]},
            "pcce": {"partition_size": 2}}
```

- [order](InputFile_General.md#tag_order), [method](InputFile_General.md#tag_method) (top level)
- [cluster.mode](#tag_cluster_mode)
- cluster.selection: [add_subclusters](#tag_add_subclusters), [max_per_order](#tag_max_clusters_per_order)
- [cluster.pcce](#tag_pcce):[partition_size](#tag_pcce_partition_size), [trials](#tag_pcce_trials),
  [max_iterations](#tag_pcce_max_iterations), [kmeans_pp](#tag_pcce_kmeans_pp), [print_trials](#tag_pcce_print_trials),
  [print_limit](#tag_pcce_print_limit), [partition](#tag_pcce_partition), [connections](#tag_pcce_connections),
  [cache](#tag_pcce_cache) (`mode`, `file`)
- [cluster.bundle](#tag_bundle): [max_dim](#tag_bundle_max_dim), [merge_cutoff](#tag_bundle_merge_cutoff),
  [connections](#tag_bundle_connections), [connection_budget](#tag_bundle_connection_budget),
  [max_clusters](#tag_bundle_max_clusters), [max_work](#tag_bundle_max_work)
- [cluster.adaptive](#tag_adaptive) ([probe_memory_gb](InputFile_Backend.md#tag_probe_memory_gb), and its internal
  [work_file](#tag_adaptive_work_file))

<a id="tag_cluster_mode"></a>
## cluster.mode
<b>Description:</b><br>
How the bath is divided into clusters. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description | method |
|---|---|---|
| `standard` | clusters of up to `order` spins within `rdip` | cCCE, gCCE, meCCE, megCCE |
| `pcce` | partition CCE: the bath is cut into partitions of [pcce.partition_size](#tag_pcce_partition_size) spins, a cluster is up to `order` partitions | cCCE, gCCE, meCCE, megCCE |
| `bundle` | bundle clustering: strongly coupled spins are tied into bundles first, a cluster is up to `order` bundles | cCCE, gCCE |
| `dsj` | disjoint clusters | cCCE |
| `itb` | inter-bath clusters, the pair interaction between disjoint clusters | cCCE |
| `dsjitb` | disjoint + inter-bath (order 2 gives the same result as standard CCE2) | cCCE |

`pcce` runs with cCCE, gCCE, meCCE and megCCE (the gCCE and master-equation forms
on the CPU only, and [pcce.cache](#tag_pcce_cache) stays with cCCE). gCCE with `dsj`, `itb` or `dsjitb` is refused (not
implemented yet); meCCE and megCCE with `bundle`, `dsj`, `itb` or `dsjitb`
likewise. Old method names are refused with their current form ([old names](InputFile.md#OldNames)).
The [pcce](#tag_pcce) block belongs to `pcce`, the [bundle](#tag_bundle) block to `bundle`, the [adaptive](#tag_adaptive)
block to `pcce` and `bundle`; each is refused with another mode. `"adaptive"` as a value is refused with the name `bundle`.

With `bundle`, [order](InputFile_General.md#tag_order) counts bundles per cluster and
[bundle.max_dim](#tag_bundle_max_dim) (required) caps the Hilbert dimension of one bundle. With
[rdip](InputFile_General.md#tag_rdip) given, bundles are connected by distance; with `rdip` omitted,
connections are chosen by a coherence check and the weak ones are dropped. Adding the [adaptive](#tag_adaptive) block turns
the adaptive refinement on top of it. The words: a *site* is one bath spin (or one defect that cannot be split); a
*bundle* is a set of sites that always enters a cluster whole and is treated exactly inside; a *cluster* is the union of
up to `order` connected bundles, the unit whose Hamiltonian is built and propagated; a *partition* is the whole assignment
of every site to exactly one bundle. Bundle clustering builds one partition and computes it once; the adaptive refinement
measures the error and rebuilds the partition. With every bundle one site, `bundle` is standard CCE. <br>

<b>Default:</b> `standard`

<b>Example:</b> ``` "cluster" : {"mode": "bundle", "bundle": {"max_dim": 8}} ```

<a id="tag_add_subclusters"></a>
## cluster.selection.add_subclusters
<b>Description:</b><br>
For CCE simulations with order ≥ 3, enabling this option allows all sub-clusters of the highest order to be included
 <br>
This leads to more stable and robust results, and may also reduce computation time. <br>
<b>Data Type:</b> Boolean (`true` or `false`)

<b>Default:</b> `true`

<b>Example:</b> ``` "cluster" : {"selection": {"add_subclusters": true}} ```

<a id="tag_max_clusters_per_order"></a>
## cluster.selection.max_per_order
<b>Description:</b><br>
Number of clusters for each CCE order. <br>
This is useful for CCE simulations with order ≥ 3, leading tons of clusters <br>
The list has exactly `order + 1` entries and the index is the cluster order: entry n limits the clusters of order n
(entry 0 is the order-0 slot and is kept, usually 0). A 0 means no limit for that order, so all zeros include every
cluster. A list of another length is refused. <br>
<b>Data Type:</b> Int List (length `order + 1`)

<b>Default:</b> all entries 0 (no limit)

<b>Example:</b> ``` "cluster" : {"selection": {"max_per_order": [0, 0, 0, 30000]}} ``` (order 3: at most 30000 clusters of three spins)

<a id="tag_pcce"></a>
## cluster.pcce
<b>Description:</b><br>
The settings of `cluster.mode` = `pcce`. Refused with any
other `cluster.mode`. `example/tutorial/07_pCCE_and_bundle/pCCE2_Hahn_states/` runs `partition_size` 1 to 3 (4 on request). <br>

<b>Data Type:</b> Object

<b>Example:</b> ``` "cluster" : {"mode": "pcce", "pcce": {"partition_size": 2, "trials": 1000, "max_iterations": 1000, "connections": "spins"}} ```

<a id="tag_pcce_partition_size"></a>
## cluster.pcce.partition_size
<b>Description:</b><br>
The number of bath spins in one partition. The bath is divided into partitions of
`partition_size` spins by a k-means search, and a pCCE cluster of order n is n partitions, n × `partition_size`
spins. Required with `cluster.mode` = `pcce`. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> - (required)

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 3}} ```

<a id="tag_pcce_trials"></a>
## cluster.pcce.trials
<b>Description:</b><br>
The number of partition trials, each a k-means run from its own starting centres;
the assignment with the best score is kept. With [pcce.partition](#tag_pcce_partition) =
`deterministic` the trials are split over the ranks in a fixed way, so the partition follows
`seed`. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `30000`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "trials": 1000}} ```

<a id="tag_pcce_max_iterations"></a>
## cluster.pcce.max_iterations
<b>Description:</b><br>
The largest number of k-means iterations in one partition trial. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `30000`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "max_iterations": 1000}} ```

<a id="tag_pcce_kmeans_pp"></a>
## cluster.pcce.kmeans_pp
<b>Description:</b><br>
Start every k-means trial from k-means++ centres instead of random ones. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `true`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "kmeans_pp": true}} ```

<a id="tag_pcce_print_trials"></a>
## cluster.pcce.print_trials
<b>Description:</b><br>
Rank 0 prints the details of the partition trials, at most [print_limit](#tag_pcce_print_limit) of them <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `false`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "print_trials": true}} ```

<a id="tag_pcce_print_limit"></a>
## cluster.pcce.print_limit
<b>Description:</b><br>
With [print_trials](#tag_pcce_print_trials): rank 0 prints the details of at most this many
trials; `0` prints none. <br>

<b>Data Type:</b> Integer (nonnegative)

<b>Default:</b> `32`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "print_trials": true, "print_limit": 8}} ```

<a id="tag_pcce_partition"></a>
## cluster.pcce.partition
<b>Description:</b><br>
How the partition of the bath is searched. <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `deterministic` | the partition follows `seed`; trials are distributed deterministically over the ranks |
| `legacy` | the original random-number and trial distribution |

<b>Default:</b> `deterministic`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "partition": "deterministic"}} ```

<a id="tag_pcce_connections"></a>
## cluster.pcce.connections
<b>Description:</b><br>
Which distance decides that two partitions are connected, so that they can be
in one cluster of order 2 or more. The inner and the outer `rdip` bound that
distance under either rule. (`cluster.mode` bundle connects its bundles by [bundle.connections](#tag_bundle_connections).) <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `centers` | the distance between the two k-means centres, the centroids of the partitions. The strength that orders the clusters is the point-dipole T<sub>zz</sub> between the centroids. |
| `spins` | any spin of the one partition between the inner and the outer `rdip` of any spin of the other, the rule of the bundles. The strength is the sum of \|T<sub>zz</sub>\| over the spin pairs of the two partitions. |

With `partition_size` = 1 the two rules are the same. Under `centers` two partitions whose spins
are close can stay unconnected when their centroids are far apart, and the more so
the larger `partition_size`; `spins` connects them. The rank-0 log says which rule ran and how
many partition pairs it connected. The comparison of the two rules on the P1 bath of tutorial `03_P1`
is (internal notes): at that tutorial's `rdip` (620 Å) the two rules
differ less than pCCE2 differs from CCE4, and `spins` costs 27-65 % more; at an `rdip` comparable
to the partition size, `spins` is the closer of the two to CCE4.

<b>Default:</b> `centers` (the rule of every earlier version)

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "connections": "spins"}} ```

<a id="tag_pcce_cache"></a>
## cluster.pcce.cache
<b>Description:</b><br>
Reuse a stored partition instead of searching again. Needs
[pcce.partition](#tag_pcce_partition) = `deterministic` and method `cCCE` (refused with gCCE, meCCE and megCCE). `mode`: <br>

| mode | Description |
|---|---|
| `off` | search every time |
| `readwrite` | reuse a matching stored partition, otherwise search and store it |
| `read` | require a matching stored partition; stop if there is none |

`file`: the partition cache file, relative to the working directory or absolute; its directory must exist. Required
with `mode` other than `off`. With [configuration.k_range](InputFile_General.md#tag_configuration_k_range) `{k}`
in it is the configuration's k, and a `file` without `{k}` is refused (one cache file per configuration, e.g.
`"partition_{k}.cache"`). <br>

<b>Data Type:</b> Object `{"mode": ..., "file": ...}`

<b>Default:</b> `{"mode": "off"}`

<b>Example:</b> ``` "cluster" : {"pcce": {"partition_size": 2, "cache": {"mode": "readwrite", "file": "./results/partition.cache"}}} ```

<a id="tag_bundle"></a>
## cluster.bundle
<b>Description:</b><br>
The settings of `cluster.mode` = `bundle`. Refused
with any other `cluster.mode`. `max_dim` is required. <br>

<b>Data Type:</b> Object

<b>Example:</b> ``` "cluster" : {"mode": "bundle", "bundle": {"max_dim": 8, "merge_cutoff": 1.0}} ```

<a id="tag_bundle_max_dim"></a>
## cluster.bundle.max_dim
<b>Description:</b><br>
The largest Hilbert dimension one bundle may have. Required with `cluster.mode` =
`bundle`. The dimension of a cluster follows as the product of its bundle dimensions (gCCE: times the
qubit dimension); a cluster that cannot be handled is refused during the run. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> - (required)

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8}} ```

<a id="tag_bundle_merge_cutoff"></a>
## cluster.bundle.merge_cutoff
<b>Description:</b><br>
The smallest coupling strength for which two groups are merged into one bundle. <br>

<b>Data Type:</b> Double (nonnegative)

<b>Unit:</b> rad/ms

<b>Default:</b> `0`

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "merge_cutoff": 1.0}} ```

<a id="tag_bundle_connections"></a>
## cluster.bundle.connections
<b>Description:</b><br>
How bundles are connected into clusters. Omit it together
with `rdip` to get the coherence-checked connections (the run log then says `checked`; that
word is not a value to write). <br>

<b>Data Type:</b> String

<b>Options:</b><br>
| Option | Description |
|---|---|
| `distance` | pairs within `rdip`; `rdip` is required |
| `auto` | the legacy budget mode: weak pairs are left out up to [connection_budget](#tag_bundle_connection_budget), which is required |

<b>Default:</b> `distance`

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "connections": "auto", "connection_budget": 0.01}} ```

<a id="tag_bundle_connection_budget"></a>
## cluster.bundle.connection_budget
<b>Description:</b><br>
With `connections` = `auto`, the allowed score of the left-out pairs: the largest
time times the sum of their coupling weights must stay below it. <br>

<b>Data Type:</b> Double (nonnegative)

<b>Default:</b> -

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "connections": "auto", "connection_budget": 0.01}} ```

<a id="tag_bundle_max_clusters"></a>
## cluster.bundle.max_clusters
<b>Description:</b><br>
The largest number of clusters the expansion of a bundle partition may list.
`order` may not exceed it. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `200000`

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "max_clusters": 200000}} ```

<a id="tag_bundle_max_work"></a>
## cluster.bundle.max_work
<b>Description:</b><br>
The work limit of the expansion and of the dimension search of a bundle partition. <br>

<b>Data Type:</b> Integer (positive)

<b>Default:</b> `20000000`

<b>Example:</b> ``` "cluster" : {"bundle": {"max_dim": 8, "max_work": 20000000}} ```

<a id="tag_adaptive"></a>
## cluster.adaptive
<b>Description:</b><br>
Turns on the adaptive refinement: instead of one calculation the run repeats calculating,
measuring the error and correcting the bundles until `atol` and `rtol` are met or `max_rounds`
is reached. Belongs to `cluster.mode` `pcce` or `bundle` (refused with another mode); the refinement driver itself runs
with `cluster.mode` = `bundle`, `bundle.max_dim` and an explicit [seed](InputFile_General.md#tag_seed). The keys are the
driver's options with `_` for `-`; `output`, `atol` and `rtol` are required, and every step goes into `output` (one folder
per step and `summary.json`). All keys with their defaults are in [the keys of cluster.adaptive](#tag_adaptive_keys).
A top-level `adaptive` block, or the block under an older name, is refused (with the place it belongs). Its
[probe_memory_gb](InputFile_Backend.md#tag_probe_memory_gb) is the memory of the probes, read by every step. <br>

<b>Data Type:</b> Object

<b>Options:</b><br>
| Option | Description |
|---|---|
| `policy` | `merge`, the only policy and the default. The removed policies `error-gated`, `budgeted`, `auto-partition`, `cost-preserving` and `auto-order` are refused |
| `output` | folder for the steps and the summary (required) |
| `atol`, `rtol` | absolute and relative tolerance of the coherence (required) |
| `max_rounds` | the largest number of rounds, default 8 |
| `probe_memory_gb` | the memory of the probes running at once on a node ([its entry](InputFile_Backend.md#tag_probe_memory_gb)) |

<b>Default:</b> - (absent: one calculation)

<b>Example:</b> ``` "cluster" : {"adaptive": { "policy" : "merge", "output" : "./results/adaptive", "atol" : 1e-3, "rtol" : 0, "max_rounds" : 4 }} ```

<a id="tag_adaptive_keys"></a>
## The keys of cluster.adaptive
<b>Description:</b><br>
Every key of the block. The names are the refinement driver's options with `_` for `-` (the stand-alone
`ccex-adaptive` takes the same options as `--max-rounds` etc.), read without regard to case. A `true/false` key is on
only when `true`. An unknown key, a value of the wrong type or a missing required key stops the run before anything is
computed (`Error: cluster.adaptive.<key> ...`). `input`, `binary`, `cwd`, `mpiexec` and `ranks` cannot be given: the
input is this file, the binary is the running process and the ranks are those of `mpirun`. The steps go into `output`:
one folder per step (`0001-plan`, `0002-probe`, ...) and `summary.json`; the bath is read once, since the steps run in the
same process. The `input.json` a step writes has no `adaptive` block. `configuration.seed` is required (every round
recomputes the same bath states). <br>

| Key | Type | Default | Meaning |
|---|---|---|---|
| `output` | path | required | folder of the steps and of `summary.json` |
| `atol` | real | required | absolute tolerance of the coherence |
| `rtol` | real | required | relative tolerance of the coherence (`atol` + `rtol` > 0) |
| `policy` | string | `merge` | the only policy: merge the bundles whose merge reduces the local error most |
| `max_rounds` | integer | `8` | the largest number of rounds |
| `max_patches` | integer | `64` | the largest number of patches (bundle neighbourhoods) checked in one round |
| `max_probes` | integer | `4096` | the largest number of probe calculations |
| `max_candidates` | integer | `16` | the largest number of merge candidates per round |
| `neighbors` | integer | `8` | the number of neighbouring bundles a patch takes |
| `cache_mb` | integer | `256` | memory (MB) for the cached probe curves |
| `timeout` | real | `600` | time limit (s) of one probe run started as a separate process by the stand-alone driver (not used inside `main.out`) |
| `raw_cache` | path | - | file in which probe results are kept and reused across runs |
| `probe_batch_clusters` | integer | `0` | clusters per probe batch; `0` derives it from `cache_mb` |
| `merge_priority` | string | `gain-per-cost` | order of the merges: `gain-per-cost` or `worst-error` |
| `max_total_work` | real | `8000000000000.0` | work limit of the whole run (D³ × states × time points) |
| `prefer_physical` | true/false | off | prefer candidates whose curve stays physical (magnitude of L at most 1) |
| `hub_candidates` | true/false | off | add candidates built around a hub bundle (needs `prefer_physical`) |
| `regroup` | true/false | off | instead of merge's selection, regroup the bundles where the magnitude of L exceeds 1 (cCCE order 2 only; not with `prefer_physical` or `hub_candidates`) |
| `regroup_max_regions` | integer | `64` | with `regroup`: regions checked per round |
| `regroup_max_candidates` | integer | `256` | with `regroup`: candidates per region |
| `regroup_factors` | string | `probes` | with `regroup`: `probes` or `factors` (one `factors` step instead of probes for the culprit pairs) |
| `regroup_whole` | string | `probes` | with `regroup`: `probes` or `except` (needs an x87 80-bit long double build) |
| `probe_memory_gb` | real | `0` (automatic) | the memory of the probes running at once on a node ([its entry](InputFile_Backend.md#tag_probe_memory_gb)) |

<b>Example:</b> ``` "cluster" : {"mode": "bundle", "bundle": {"max_dim": 8}, "adaptive": {"output": "./adaptive_out", "atol": 1e-8, "rtol": 2e-6}} ```

<a id="tag_adaptive_work_file"></a>
## cluster.adaptive.work_file
<b>Description:</b><br>
The plan/probe protocol file the adaptive refinement driver writes for each step it runs. The driver sets it itself, in a step's input, in the `cluster.adaptive` block (beside the user's
`probe_memory_gb` only); a hand-written input does not need it, and it is refused beside the refinement options. <br>

<b>Data Type:</b> String (path)

<b>Default:</b> -

<b>Example:</b> ``` "cluster" : {"adaptive": {"work_file": "./work.json"}} ```
