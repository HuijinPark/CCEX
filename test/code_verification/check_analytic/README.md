# check_analytic/ — CCEX against a closed form

    bash test/code_verification/check_analytic/run.sh     # bin/main.out, or CCEX_BIN; CCEX_MPI as for the others

One 13C beside the NV (`inputs/`: qubit at (10, 10, 10) Å, the 13C at (13, 11.5, 14) Å,
B = 500 G along z), order 1, which is exact for one bath spin, ensemble average over
its state, 400 points to 0.2 ms. `analytic.py` writes the point-dipole hyperfine
out itself and compares the free-induction and Hahn-echo coherence with their closed
forms (its docstring has them); the tolerance, 2e-10, is the ten decimals CCEX
writes. A self-test runs the echo with a gyro 1 % off and requires the check to fail.

Every other suite under `code_verification/` compares CCEX with CCEX (another input
form, another frame, a recorded run), so an error common to both sides passes there:
a wrong time scale in the free evolution, a wrong constant in the dipole formula.
This one does not, and it is small and machine-independent, so the CI workflow runs
it (D22, 2026-10-02). Recorded on WSL (Intel MPI, Make and CMake builds): deviation
5.0e-11 for both pulse counts; the 1 % wrong gyro gives 0.37; reading the time
column as tau instead of the total time gives 0.37 (echo) and 1.96 (free induction).

`analytic.py` exits 1 on a numerical mismatch and 2 when the output is missing, cannot be
read or does not hold the requested grid (400 points, k * 0.0005 ms); the self-test counts
only a clean CCEX run that `analytic.py` rejects with 1. The qubit file form takes the branch
keys as `qalphams`/`qbetams`.

Exit status 0 only if all three pass. Everything written goes to `_work/`.
