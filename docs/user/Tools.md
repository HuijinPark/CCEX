# Tools

What makes the inputs CCEX reads: three tools in `example/tools/`. A fourth,
`AutoConvergence/`, is the automatic convergence: an input with an `auto_convergence`
block run as `main.out -f input.json` (without mpirun) chooses the time window, rdip,
rbath and the configurations itself; see `example/tools/AutoConvergence/README.md`.

| Tool | Makes | Needs |
|---|---|---|
| `BathGenerator/` | a bath file, site positions and isotopes, from a POSCAR unit cell | numpy, psutil, mpi4py with an MPI |
| `SpinTensorReader/` | the hyperfine (`Afile`) and quadrupole (`Qfile`) tensor files, from a VASP OUTCAR or a Quantum ESPRESSO output | numpy, pandas, psutil |
| `Runner/` | the job folders of the tutorials on the cluster, one per value of a variable | the idm Slurm cluster |

## BathGenerator

Repeats the unit cell to the requested size and draws the isotope of every site
with the probabilities in a `prob` file; options for a 2D material, a defect, an
inserted relaxed cell and a deformation. The isotopes are drawn without a fixed
seed, so two runs give the same sites with different isotopes.

```
mpirun -n 4 python3 BathGenerator.py -o OUT -i POSCAR -p PROB --size 30 --Defect B --opt2D
```

`-h` lists everything; `do_BathGenerator` has three worked cases. The baths of
the tutorials were made with this copy. A newer generator lives at
<https://github.com/CQML-SKKU/BathGeneratorSpinCluster>.

## SpinTensorReader

Reads the hyperfine and electric-field-gradient tensors a DFT code wrote for
every nucleus around the defect and writes them in the layout the reference
describes under [hf_tensorfile](InputFile_General.md#tag_hf_tensorfile) and
[qd_tensorfile](InputFile_General.md#tag_qd_tensorfile): the species lines,
the range the file covers and one row per atom with its position relative to
the defect and the nine components.

```
python3 main.py VASP POSCAR DEFECT --A Afile --Q Qfile [--version v1|v2] [--rangeformat vertex]
```

`spHamilReader/` is the package (`base/`, `vasp/`, `qe/`), `main.py` the
command line (`-h`). `examples/01..05` are five worked cases with their inputs
and the tensor files they produced, each with a `do` that reproduces it; the
tensors of the h-BN tutorial (07) came from it.

## Runner

Makes the Slurm job folders the tutorials were recorded with: one per value of
a variable (`order`, ppm, `rbath`, ...) from the templates in
`MkCCEJobs/mkjobs.sngapp/`, driven by a tutorial's `domkjobs` and `opt/` files
(`./domkjobs -h`). It is written for the group's idm cluster and its module
names; on another machine the tutorials' `run.sh` do the same job by hand.
