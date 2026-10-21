# example/tools/

Six tools: four that make what CCEX reads (`BathGenerator/`, `SpinTensorReader/`, `IsotopeTable/` its built-in isotope table, and `StateFiles/` the state and axis files of an input), the runner of the automatic convergence (`AutoConvergence/`), which runs CCEX itself, and the input-key tools (`InputKeys/`: the converter of inputs written with older keys, and a writer of an input copy with keys set).

| Tool | Makes | Run | Needs |
|---|---|---|---|
| `BathGenerator/` | a bath file (site positions and isotopes) from a POSCAR unit cell, drawing isotopes with the probabilities in a `prob` file; options for 2D materials, a defect, an inserted relaxed cell, deformation | `mpirun -n N python3 BathGenerator.py -o OUT -i POSCAR -p PROB --size 30 --Defect B --opt2D` (`-h` lists everything; `do_BathGenerator` has four cases: three that make a bath and one that converts a configuration to POSCAR) | numpy, psutil, mpi4py with an MPI |
| `SpinTensorReader/` | the hyperfine (`Afile`) and quadrupole (`Qfile`) tensor files CCEX reads, from a VASP OUTCAR or a QE output | `python3 main.py VASP\|QE POSCAR DEFECT --A OUT --Q OUT [--spin S] [--rangeformat vertex] ...` (`-h`; `--spin` required for QE `--A`); `examples/01..05` each have a `do` | numpy, pandas, psutil |
| `IsotopeTable/` | `src/input/isotope_table.cpp`, the 154 nuclear isotopes (spin, gyromagnetic ratio) CCEX uses for a species no gyrofile names, from the EasySpin database in `SpinTensorReader/` | `python3 example/tools/IsotopeTable/make_isotope_table.py` from the repository root | python3 |
| `StateFiles/` | the `bath_state_file`, `defect_axis_file` and `defect_state_file` files a CCEX input names (`configuration` block: `bathfile`, `k_range`, `s_range`, `seed`, the `{k}`/`{s}` names; `Defect.species`), m uniform over the 2S+1 levels from numpy `default_rng([seed, k, s])`, axes uniform over 1..`n_axis_states` from `default_rng([seed, k, 0])`, on-site states from `default_rng([seed, k, s, 1])` - the rule the automatic convergence makes its own files with (`docs/user/Tools.md#tag_tool_statefiles`) | `python3 example/tools/StateFiles/make_state_files.py INPUT [--force]`; tests `python3 -m unittest discover -s example/tools/StateFiles/tests` | python3 with numpy |
| `AutoConvergence/` | the automatic convergence: `main.out -f input.json` with an `auto_convergence` block (or `mpirun -n N main.out -f input.json`: rank 0 runs the runner and every CCEX run gets `mpirun -n N`; a `launcher` key in the block is refused) finds the time window, rdip, rbath and the configurations and runs the production on the baths you prepare (it says which baths to add when they are too small or too few); `PYTHONPATH=example/tools/AutoConvergence python3 -m autoconv status input.json` shows where it is. `CCEX_PYTHON` picks the interpreter, `CCEX_AUTOCONV_DIR` the tool folder; it may be a shared installation (owned by you, root or the owner of `main.out`), but neither it nor any folder above it may be writable by others (other users, or a group other than the owner's own private group; `chmod go-w`). Also `python -m autoconv analyze` for existing runs (`AutoConvergence/README.md`) | `main.out -f input.json`; tests `python -m unittest discover -s tests` | numpy, scipy, matplotlib |
| `InputKeys/` | an input in the current keys from one written before (`method` cCCE/gCCE with `cluster_mode`, the `hamiltonian`, `pcce`, `bundle`, `pulse`, `configuration`, `execution`, `cluster`, `output` and `frame` blocks, ...; the table is [Old names](../../docs/user/InputFile.md#OldNames)); `.json` or `.jsonc`, comments (`//`, `#`, `!`) stay, a file already in the new keys is left alone | `python3 example/tools/InputKeys/convert_input_keys.py FILE...` rewrites in place and lists what changed; `--check FILE...` lists only (exit 1 while old keys remain); `--stdout FILE` prints the new text; `python3 example/tools/InputKeys/set_input_keys.py IN OUT KEY=VALUE ...` writes a copy with keys set (dotted for a block: `configuration.s_range=[1,5]`, `configuration.defect_axis_file=["Random1"]`; VALUE as JSON when it is JSON), what the tutorials' `run.sh` do per run, as the command line has no `-I`, `-N`, `-s`, `-S`, `-a` | python3 |

## BathGenerator

A newer version of this generator lives at
<https://github.com/CQML-SKKU/BathGeneratorSpinCluster>. This copy is the one the
tutorials' baths were made with and is kept for that reason.

The isotopes are drawn without a fixed seed, so two runs give the same sites with
different isotopes. On this WSL the anaconda `mpi4py` aborts at interpreter exit
("double free or corruption") even for a bare `from mpi4py import MPI`; the bath
file is complete before that and the tool is not at fault.

## SpinTensorReader

`spHamilReader/` is the package (`base/`, `vasp/`, `qe/`),
`main.py` the command line. `examples/` are five worked cases with their inputs and
the tensor files they produced; `do` in each reproduces them (`04.Diamond_N2p/do`
takes two arguments, the JT index and the case, e.g. `bash do 4 1-2`). The
package's internal imports are package-relative, and every
example was checked to reproduce its tracked `Afile`/`Qfile` byte for byte.

**Tensor format v3.** Every `Afile` starts with the tag (`VASP` / `QUANTUMESPRESSO`) and the provenance
line `# CCEX-tensor v3: values = hyperfine of the spin-S defect (MHz), source = <VASP OUTCAR|QE hyperfine.out>, S = <S>`,
and holds the hyperfine of the spin-S defect, which CCEX reads x 1 (`docs/user/InputFile_General.md#tag_hyperfine_file`).
VASP: the OUTCAR numbers (VASP divides by the total spin), S from the OUTCAR magnetization (`--spin`, if given, must
agree). QE: GIPAW's `hyperfine.f90` uses the bare spin density with the S = 1/2 prefactor, so `--spin S` is required
and the file holds `hyperfine.out` x 1/(2S). `--version v1` is refused: it wrote the OUTCAR x 2, which is the defect's
value only for S = 1 and only when CCEX then divides by 2S; with a `VASP` line CCEX read it 2x (the h-BN tutorial file,
and example 03's MoS2 `Afile`). `--fidentify` is gone (the lines are always written). `Qfile`s (EFG, independent of the
spin) are unchanged. The examples were re-run in format v3: 01 (QE, S = 1) nearest ¹⁴N A_zz 50.4 MHz (PBEsol); 02 (VASP)
48.8 MHz without the core term (`Afile_nocore`, = the removed `Afile_v1` / 2) and 35.4 MHz with it (the default
`--addhfcore True`); 03 and 05 S = 1/2; 04 v2 only (the v1 outputs removed). Every `Qfile`, and every v2 `Afile` apart
from its two new lines, is byte-identical.

**Wide values.** The `Afile`/`Qfile` fields are `%10.5f` between tabs, and CCEX's tensor reader split a
line at spaces only: a value of 10 characters or more (-123.45678, 1234.56789) ran into the field before it, and the
row was dropped (CCEX then stopped with "doesn't exist in A-file") or read shifted. The writers now put a space in front
of every value, and CCEX splits at any blank, so old files with such values are read too (`runningScripts_Exspin/`'s
NVHn `Afile`s have an iso of -141 MHz). `tests/check_wide_fields.py` (run by the `tensor_format` stage of
`test/code_verification/check_input_forms`) writes values up to -12345.6 with both writers, in the new and in the old
format, and checks that CCEX reads them as written. Examples 01-05 were re-run: every `Afile`/`Qfile` has the same
numbers as before (only the added spaces). The unused copies `v2_vasp/`, `v2_qe/` (not imported anywhere) were removed.
Python 3.12 prints SyntaxWarnings for a few regular expressions; they are harmless.

`runningScripts_Exspin/` holds two example runs (N2p, NVHn) with their OUTCARs and
a `do_summary` written for the author's own machine (absolute paths).

## Runner (removed)

`Runner/` (the `MkCCEJobs` job templates driven by `domkjobs`: one job folder per value of a variable, which then
called `analysis/ConvTest`) was removed before 2.0.0; `configuration.k_range` (one input, many
bath configurations) and the automatic convergence (`AutoConvergence/`) took its place. git history has it
(`git show c3c862972:example/tools/Runner/domkjobs`).
