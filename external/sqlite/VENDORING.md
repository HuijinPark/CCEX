# SQLite amalgamation, vendored for `bin/ccex-refine`

The probe cache (`src/clustering/algo_adaptive/base_run.cpp`) is a SQLite file. It was first
shared with the Python implementation, whose machines linked different SQLite versions
(3.26 to 3.50), so the port carried one of its own. The Python is gone (2026-09-20);
the copy stays because it keeps the build self-contained on every server - no system
`libsqlite3-dev` to look for - and the file format has been stable since 3.0.

- Source: <https://sqlite.org/2026/sqlite-amalgamation-3530400.zip> (SQLite **3.53.4**),
  downloaded 2026-09-16.
- The archive's published hash is **SHA3-256** `628a44cfe82c66aed1ccbbe85a562d2e33ebe64b3288981ed76285612227934e`,
  which the download matches. (Its SHA-256 is `1e71ddf93849c6a6ecf58b827c0692073d2dd7ee40196158068f7b29f422e87d`
  - worth writing down, because comparing that against the published value looks like a
  mismatch until you notice the algorithm.)
- Kept here, with SHA-256 of each file:

```
sqlite3.c  b1dd5d74ec7f29055a6684fa06fb3c2f6821c87dd38f9a458dfd2e8a1db28189  9515341 bytes
sqlite3.h  919e7f2e8ed1d8f56ac17b412b8971c76aa5d1a879752cc6058f75e7d5910e1d  690838 bytes
```

`shell.c` and `sqlite3ext.h` are not kept: the port links the library, not the shell,
and loads no extensions.

## Rebuilding

Download the amalgamation for the version you want, check its SHA3-256 against the
download page, and copy `sqlite3.c` and `sqlite3.h` here. What must hold afterwards is
that `test/cluster_verification/check_adaptive/verify.sh regress` still passes: it reads the caches recorded under
`test/cluster_verification/check_adaptive/regression/`.
