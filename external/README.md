# external/

Third-party code kept in the tree, so that every server (idm, idm4-gpu, nurion,
WSL) builds without looking for packages. Nothing here is patched.

| Folder | What | Version | Used by |
|---|---|---|---|
| `eigen/` | Eigen, header-only linear algebra (`-I./external/eigen`) | 3.4.1 | the core (Hamiltonians, propagation) and `src/calculation/gpu/` |
| `uthash/` | uthash, a hash table in one header (`uthash.h`) | 2.3.0 | `src/clustering/cluster_hash.cpp` |
| `sqlite/` | the SQLite amalgamation (`sqlite3.c`, `sqlite3.h`), compiled once into `obj/sqlite3.o` | 3.53.4 | the probe cache of the adaptive clustering (`src/clustering/algo_adaptive/base_run.cpp`) |

Not in this folder but third-party as well: **cJSON 1.7.17** (MIT license, in the files), the JSON parser of the input
file, kept as `src/input/json.cpp` and `include/input/json.h` (`CJSON_VERSION_*`). Whether it differs from the upstream
release has not been checked (2026-10-02, review 7eed38 R5). The key tracking of D24 wraps it from `src/input/option.cpp`
and does not touch these two files. To update it: replace both, keep the version macros in step, and run
`test/code_verification/run.sh` (the refusal and input-form checks read every input form).

Licenses are in each folder (`eigen/COPYING.*`, `uthash/LICENSE`; SQLite is
public domain). `sqlite/VENDORING.md` records where the amalgamation came from,
its hashes, and why it is carried rather than linked from the system.

To update one: replace the folder with the new release, change the version in
this table (and in `VENDORING.md` for SQLite), rebuild with `do_compile.sh`,
and run `test/cluster_verification/check_adaptive/verify.sh unit` and `regress`.
