# 01_Arbitrary_system

A hand-written system. a bath, a gyromagnetic-ratio file and a qubit file written by hand, to show the file formats.

| Calculation | Inputs | What it shows | Recorded run | 2026-09-23 on WSL, 2 ranks, 1 configuration |
|---|---|---|---|---|
| `./` | ccein.json | gCCE order 1 on a hand-written three-file system; the density matrix (`quantity=dm`) of one bath state | 1 conf × 1 states, 1 ranks | 0 s |

Each calculation folder holds its `ccein*.json` inputs, `run.sh` and the idm `job`
that submits it. `run.sh` runs every input on `NP` ranks for `NCONF` bath
configurations (default 4 and 1) and then the post-processing that the recorded
run did, all into `results/`:

    cd <calculation> && NP=4 NCONF=1 bash run.sh

`example/tutorial/README.md` explains the layout; the recorded results are in
the archive.
