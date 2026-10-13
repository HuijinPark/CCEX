# Output tags
- [savemode](#tag_savemode) (`mode`, `nodiv`)
- [outfile](#tag_outfile)

<a id="tag_savemode"></a>
## savemode
<b>Description:</b><br>
What is written, and whether `_noDiv` is written beside `_wiDiv`. A string gives the mode (with `nodiv` true); the block
`{"mode": ..., "nodiv": true|false}` gives both (a missing `mode` is `"coherence"`, a missing `nodiv` is `true`). <br>
<b>mode:</b>

- `"coherence"` : The coherence of each state: "output_state{state}_wiDiv/noDiv". <br>
- `"coherence_tagged"` : Like `"coherence"` (the same file names), but every time row also carries the zx and zz
                 components of the first bath spin's hyperfine tensor and the first bath file's name, so that results of
                 many bath files can be told apart after collecting them. <br>
- `"cluster_coherence"` : The result and, for each cluster and each state, its coherence:
                 "output_{cluster}_state{state}_wiDiv/noDiv". <br>
- `"cluster_density_matrix"` : Everything `"cluster_coherence"` writes, and with gCCE also, for every cluster,
                 the density matrix of the whole qubit + cluster space at every time step, as `n` × `n` blocks (`n` the
                 full dimension), to "output_{cluster}_full_state{state}_wiDiv/noDiv" (the cluster's bath indices as in
                 `"cluster_coherence"`; "output_full_state{state}_..." for the zeroth cluster, the qubits alone). The
                 `_wiDiv` and `_noDiv` files hold the same matrix. This costs a matrix of the full dimension per time step
                 and cluster in memory and on disk, so it is meant for small systems (tutorial 01). cCCE has no qubit +
                 bath density matrix: there it writes what `"cluster_coherence"` writes. <br>

<b>nodiv:</b> Write `_noDiv` beside `_wiDiv`. `_wiDiv` is the CCE coherence; `_noDiv` is the product of the positive-power
cluster terms only (no division), which never diverges but is **not** a coherence. It is read by the group analyzer's
error correction (`example/analysis/Analyzer`, \|wiDiv − noDiv\| ≥ 2 marks a divergent point), by most tutorials' `run.sh`,
by the standing checks and by the paper scripts, so it stays on by default. `false` writes `_wiDiv` alone (half the
files), for runs whose analysis does not use it. Must be a JSON boolean. <br>

<b>Names before 2026-10-13</b> (refused, with a message naming the new form):

| old | new |
|---|---|
| `"savemode": "normal"` | `"savemode": "coherence"` |
| `"savemode": "info"` | `"savemode": "coherence_tagged"` |
| `"savemode": "all"` | `"savemode": "cluster_coherence"` |
| `"savemode": "allfull"` | `"savemode": "cluster_density_matrix"` |
| `"save_nodiv": false` | `"savemode": {"mode": "coherence", "nodiv": false}` |

<b>Data Type:</b> String or object (`mode`: string, `nodiv`: boolean)

<b>Default:</b> "coherence" (`nodiv` true)

<b>Example:</b> ``` "savemode" : "cluster_coherence" ```, ``` "savemode" : {"mode": "coherence", "nodiv": false} ```

<b>Related variable:</b> [Output](#Output)::savemode, [Output](#Output)::save_nodiv

<a id="tag_outfile"></a>
## outfile
<b>Description:</b><br>
Output file name prefix. <br>
"outfile_..._wiDiv" is the final result of the coherence function. <br>
"outfile_..._noDiv" is the result without dividing the final coherence function into those obtained from its sub-clusters. <br>
You can specify this value via the "-o" command-line option.
<b>Data Type:</b> String

<b>Default:</b> -

<b>Example:</b> ``` outfile : "filepath" ```

<b>Related variable:</b> [Output](#Output)::outfile
