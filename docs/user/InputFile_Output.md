# Output tags
- [savemode](#tag_savemode)
- [outfile](#tag_outfile)
- [save_nodiv](#tag_save_nodiv)

<a id="tag_savemode"></a>
## savemode
<b>Description:</b><br>
The output file can be saved in different ways. <br>
<b>Options:</b>

- `"all"` : Saves all coherence functions for each cluster and each state. <br>
                The file name will be in the format: "output_{cluster}_state{state}_wiDiv/noDiv" <br>
   - `"normal"` : Saves coherence functions for each state only. <br>
                   The file name will be in the format: "output_state{state}_wiDiv/noDiv" <br>

<b>Data Type:</b> String

<b>Default:</b> "normal"

<b>Example:</b> ``` savemode : "all" ```

<b>Related variable:</b> [Output](#Output)::savemode

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

<a id="tag_save_nodiv"></a>
## save_nodiv
<b>Description:</b><br>
Write `_noDiv` beside `_wiDiv`. `_wiDiv` is the CCE coherence; `_noDiv` is the product of the positive-power cluster terms
only (no division), which never diverges but is **not** a coherence. It is read by the group analyzer's error correction
(`example/analysis/Analyzer`, \|wiDiv − noDiv\| ≥ 2 marks a divergent point), by most tutorials' `run.sh`, by the standing
checks and by the paper scripts, so it stays on by default. `false` writes `_wiDiv` alone (half the files), for runs
whose analysis does not use it. Must be a JSON boolean. <br>

<b>Data Type:</b> Boolean

<b>Default:</b> `true`

<b>Example:</b> ``` "save_nodiv" : false ```

<b>Related variable:</b> [Output](#Output)::save_nodiv
