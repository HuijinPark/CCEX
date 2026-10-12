# docs/

What to read, by what you are trying to do. Dated records are in git history (last row);
everything here is current.

| Folder | For | Read |
|---|---|---|
| `user/` | running CCEX | `GettingStarted.md` (build, `mpirun ... -f ccein.js`, the command-line flags), `InputFile.md` (the input file and its seven groups of keys, one page per group: `InputFile_General.md`, `_Qubit`, `_Cluster`, `_Backend`, `_Pulse`, `_Defect`, `_Output`, each key with its type, unit and default), `Overview.md` (what the code computes and how it is organised), `Analysis.md` (the analyzer and the convergence test, from `example/analysis/`), `Tools.md` (the bath generator, the tensor reader and the job runner, from `example/tools/`), `Examples.md` (the ten tutorials, one page each: `Examples_01.md` … `Examples_10.md`, with the figures in `../images/tutorials/`). Converted from the Doxygen pages on 2026-09-24; the tutorials themselves are `example/tutorial/` |
| `methods/` | understanding a method this branch added | `clustering-terminology.md` (site / bundle / cluster / partition; bundle clustering vs adaptive refinement), `adaptive-cluster-mode.md` (bundle clustering: `cluster_mode=adaptive`, `bundle_max_dim`, with or without `rdip`), `adaptive-input-block.md` (the `Adaptive` block that turns on the refinement, and its keys), `gpu.md` (the CUDA backend: `backend`, `gpu_min_dim`, `gpu_batch_steps`, what runs on the GPU), `coordinate-frame-rotation.md` (the `coordinate_frame_rotation` block: what the rotation transforms, leaves alone and refuses, and its keys) |
| `doxygen/` | the API reference | `do_doxy` (or `doxygen Doxyfile`) runs Doxygen on `include/`, `src/` and `user/` (`mdfilter.sh` strips the page part of a `Page.md#anchor` link, which Doxygen cannot resolve but GitHub can); the output (`html/`, `latex/`, `xml/`) is written here and not committed. The look is the doxygen-awesome theme (two `.css` files, MIT) with `ccex.css` on top; `header.html` is Doxygen's own template plus one line that loads `ccex.js`, which opens the tree entry of the page being shown |
| `images/` | figures the pages above use | |
| (records) | what was decided and measured, by date | the 87 dated records of 2026-09-11 to 2026-09-23 (and the pre-branch `TODO`) were taken out of the tree on 2026-09-24; they are in git history up to commit 83fe451d5: `git show 83fe451d5:docs/archive/<name>`, `git ls-tree 83fe451d5 docs/archive/` lists them |

The standing checks are described in `test/README.md`, the source layout in
`src/README.md`, the examples in `example/README.md`.
