# Tools

The tools in `example/tools/`: two that make the inputs CCEX reads (BathGenerator, SpinTensorReader), the generator of
the built-in isotope table (IsotopeTable), the input-key tools (InputKeys), the writer of the state and axis files of an
input (StateFiles) and the runner of the automatic convergence
(AutoConvergence): an input with an `auto_convergence` block run as `main.out -f input.json` (or
`mpirun -n N main.out -f input.json`, every run then on N ranks) chooses the time window, rdip, rbath and the
configurations itself ([Automatic convergence tags](InputFile_AutoConvergence.md)).

| Tool | Makes | Needs |
|---|---|---|
| `BathGenerator/` | a bath file, site positions and isotopes, from a POSCAR unit cell | numpy, psutil, mpi4py with an MPI |
| `SpinTensorReader/` | the hyperfine (`Afile`) and quadrupole (`Qfile`) tensor files, from a VASP OUTCAR or a Quantum ESPRESSO output | numpy, pandas, psutil |
| `IsotopeTable/` | `src/input/isotope_table.cpp`, the isotopes CCEX knows without a gyrofile | python3 |
| `InputKeys/` | an input in the current keys from an old one; a copy of an input with keys set | python3 |
| `StateFiles/` | the `bath_state_file`, `defect_axis_file` and `defect_state_file` files an input names, drawn by a fixed rule | python3 with numpy |
| `AutoConvergence/` | the automatic convergence (started by `main.out`), and the analysis of existing runs | python3 (3.10 or newer) with numpy, scipy, matplotlib |

## BathGenerator

Repeats the unit cell to the requested size and draws the isotope of every site
with the probabilities in a `prob` file; options for a 2D material, a defect, an
inserted relaxed cell and a deformation. The isotopes are drawn without a fixed
seed, so two runs give the same sites with different isotopes.

```
mpirun -n 4 python3 BathGenerator.py -o OUT -i POSCAR -p PROB --size 30 --Defect B --opt2D
```

<a id="tag_prob_file"></a>
**The `prob` file** gives, for every isotope that can sit on a site, three lines: its name (element symbol with the
mass number), its nuclear spin and its abundance in %. The isotopes of one element (the name without its digits, which
has to be an atom name of the POSCAR) are drawn together; when their abundances add up to less than 100 %, the rest of
the sites of that element hold no spin and are left out of the bath (more than 100 % is refused). The natural-abundance
diamond file of the tutorials (`example/tutorial/00_quickstart/bath/genbath/prob_natab`) is

```
13C
0.5
1.1
```

so 1.1 % of the carbon sites hold a ¹³C. For h-BN (`example/tutorial/04_hBN/bath/prob_hBN`) the file lists ¹⁰B (3,
19.9), ¹¹B (3/2, 80.1), ¹⁴N (1, 99.632) and ¹⁵N (1/2, 0.368).

| Option | Meaning |
|---|---|
| `-o OUT` | the output bath file; the defect position goes to `OUT_defect` |
| `-i POSCAR` | the unit cell of the pure crystal (or, with `--convert`, a bath file) |
| `-p PROB` | the `prob` file |
| `--size L` | the edge of the bath cell in Å (a cube; with `--opt2D` the z extent is the unit cell's) |
| `--Defect NAME` | the atom name of the defect site: of the sites with that name, the one closest to `--Defect_center` is the defect, the centre of the bath, written to `OUT_defect` |
| `--Defect_center x y z` | the fractional point (each between 0 and 1) the defect site is chosen near (default 0.5 0.5 0.5) |
| `-I POSCAR` | a relaxed supercell with the defect, inserted around it |
| `-D POSCAR ... --axis AX ...` | deformed cells and their axes (`x`, `y`, `z`: plane shape; `xy`, `xz`, `yz`: line shape) |
| `--opt2D` | a two-dimensional bath (one layer, as h-BN) |
| `--convert` | convert a bath file (`-i`) to a POSCAR (`-o`) to look at it |
| `--test` | also write a POSCAR of the bath, to check it |

`-h` lists everything; `do_BathGenerator` has four worked cases (three make a bath, one converts a bath to POSCAR). The baths of
the tutorials were made with this copy (`example/tutorial/make_natab_baths.sh` calls it for the 13C baths of 01 and 09). A newer generator lives at
<https://github.com/CQML-SKKU/BathGeneratorSpinCluster>.

## SpinTensorReader

Reads the hyperfine and electric-field-gradient tensors a DFT code wrote for
every nucleus around the defect and writes them in the layout the reference
describes under [hyperfine.file](InputFile_General.md#tag_hyperfine_file) and
[quadrupole.file](InputFile_General.md#tag_quadrupole_file): the species lines,
the range the file covers and one row per atom with its position relative to
the defect and the nine components.

```
python3 main.py VASP POSCAR DEFECT --A OUTCAR --Q OUTCAR [--spin S] [--rangeformat vertex]
python3 main.py QE   POSCAR DEFECT --A hyperfine.out --spin S --Q efg.out [--rangeformat vertex]
```

`--rangeformat` chooses how a file states the region it covers (the DFT supercell around the defect): `vertex` (the
default) writes the eight corners of the supercell relative to the defect (lines `v1` to `v8`, Å), `coord` the smallest
and largest relative coordinates along x, y and z (lines `MinDif[A]`, `MaxDif[A]`). A bath spin inside that region must
have a row in the file ([hyperfine.missing_row](InputFile_General.md#tag_hyperfine_missing_row)).

Every hyperfine file (`Afile`) is in tensor format v3: line 1 `VASP` or
`QUANTUMESPRESSO`, line 2 `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = ..., S = ...`,
and the numbers are the hyperfine of the spin-S defect, which CCEX reads as written. VASP's OUTCAR is already per
S (VASP divides by the total spin), so the numbers are the OUTCAR's and S comes from the OUTCAR's total magnetization
(a `--spin` given must agree). Quantum ESPRESSO's GIPAW `hyperfine.out` is the S = 1/2 coupling, so `--spin S` is
required and the file holds it x 1/(2S). `--version v1` (the OUTCAR x 2, an older format) is refused and there is no `--fidentify`. The quadrupole file (`Qfile`) does not depend on the spin and is
unchanged; the spin is read (and `--spin` checked) only when an `Afile` is made, so `--Q` alone takes an OUTCAR with
magnetization 0 or none (an EFG run need not be spin polarized).

`spHamilReader/` is the package (`base/`, `vasp/`, `qe/`), `main.py` the
command line (`-h`). `examples/01..05` are five worked cases with their inputs
and the tensor files they produced, each with a `do` that reproduces it; the
tensors of the h-BN tutorial (04_hBN) came from it.

<a id="tag_tool_isotopetable"></a>
## IsotopeTable

`make_isotope_table.py` writes `src/input/isotope_table.cpp`, the nuclear isotopes (spin and gyromagnetic ratio in
rad/ms/G) CCEX uses for a bath species the gyrofile does not name, or for every species when the input has no gyrofile
([gyrofile](InputFile_General.md#tag_gyrofile)). Its source is the EasySpin isotope database in
`example/tools/SpinTensorReader/spHamilReader/spin_database/isotopedata.txt`; γ = g<sub>n</sub> μ<sub>N</sub> / ħ. To
regenerate the table after the database changes, run from the repository root and rebuild CCEX:

```
python3 example/tools/IsotopeTable/make_isotope_table.py
```

<a id="tag_tool_inputkeys"></a>
## InputKeys

Two scripts for input files:

```
python3 example/tools/InputKeys/convert_input_keys.py FILE...         # rewrite inputs written with old keys, in place
python3 example/tools/InputKeys/convert_input_keys.py --check FILE... # list only; exit 1 while old keys remain
python3 example/tools/InputKeys/set_input_keys.py IN OUT KEY=VALUE ...
```

`convert_input_keys.py` reads `.json` and `.jsonc` inputs, keeps the comments (`//`, `#` or `!`) and lists what it
changed; what it maps is the table of [old names](InputFile.md#OldNames). `set_input_keys.py` writes a copy of `IN`
with keys set: `KEY` is a key or a dotted
path into a block (`configuration.s_range=[1,5]`; the block has to exist in `IN`), matched without regard to case; `VALUE` is
read as JSON when it is JSON (`50`, `true`, `["bath_1"]`) and as a string otherwise; `null` removes the key (`configuration.s_range=null`: the ensemble). `OUT` is plain JSON (the comments
are not kept). The tutorials' `run.sh` use it to run one input with other values.

<a id="tag_tool_statefiles"></a>
## StateFiles

`make_state_files.py` writes the state and axis files a CCEX input names, so that sampled states and Defect axes come from
files of a fixed rule rather than from CCEX's own random numbers (the files then also serve other codes and other runs):

```
python3 example/tools/StateFiles/make_state_files.py INPUT [--force]
```

It reads the input's `configuration` block - `bathfile`, `k_range` (every k, `{k}` in the names), `s_range` (the states,
`{s}` in a state file name or the number appended), `seed`, `path`, `gyrofile` - and its `Defect.species` (`name`,
`n_axis_states`, `n_onsite_spins`, `spins`), and writes for every k, every bath file whose name in the key is not `""` and
every s of `s_range` the files of the three keys present ([file formats](InputFile_General.md#tag_configuration_files)):

| Key | Each line after the first (one per line of the bath file) | Drawn with numpy |
|---|---|---|
| `bath_state_file` | m of the spin, uniform over its 2S+1 levels (S from the gyrofile, else the built-in isotope table) | `default_rng([seed, k, s])` |
| `defect_axis_file` | the axis state, uniform over 1..`n_axis_states` of its Defect species; 0 for another species (CCEX reads an axis file only when every spin is a Defect species with axis states) | `default_rng([seed, k, 0])` |
| `defect_state_file` | the m of each on-site spin of its Defect species (uniform over 2I+1, `spins`), space-separated; 0 for another species | `default_rng([seed, k, s, 1])` |

The first line of each file is the bath file's first line. `seed` is the number when it is one, its entry for k when it is
a list (one per configuration of `k_range`), 1 when absent; without `k_range` k is 0. The names are read from the input
file's folder (and `configuration.path`). A file that exists is not overwritten without `--force`; a species neither in the
gyrofile nor in the built-in table, or a file key without the Defect species it needs, is refused. The bath-state and axis
rules are those of the automatic convergence, which makes the files it needs with these functions (one rule;
`example/tools/AutoConvergence/examples/nv13c/make_states.py` too). Tests: `python3 -m unittest discover -s
example/tools/StateFiles/tests` (one test runs `bin/main.out` on the files and checks that CCEX reads every one).

<a id="tag_tool_autoconvergence"></a>
## AutoConvergence

The runner of the automatic convergence. `main.out` starts it when the input has an `auto_convergence` block
([Automatic convergence tags](InputFile_AutoConvergence.md)); `CCEX_PYTHON` picks the Python interpreter. Its own
commands, run from `example/tools/AutoConvergence` (or with it on `PYTHONPATH`), compute nothing:
`python3 -m autoconv status INPUT.json` says where a run is, and `analyze` and `judge` judge existing runs
([Commands](InputFile_AutoConvergence.md#tag_ac_commands)); their `--config` is the input with the `auto_convergence`
block, read as the runner reads it. The state and axis files it makes for a run are written by
[StateFiles](#tag_tool_statefiles)' functions.

## Many runs

There is no job-folder generator: `configuration.k_range` runs many bath configurations from one input, and the
automatic convergence chooses the values a convergence study needs.
