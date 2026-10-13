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
only a clean CCEX run that `analytic.py` rejects with 1. The inputs use the `Qubit.file` form with the branch
keys `Qubit.alphams`/`Qubit.betams` (`qalphams`/`qbetams` until 2026-10-13).

More runs, not closed forms, cover `hfmedi` at zero field (review 20261012-194049-af8e42 R1): with an axial
ZFS the levels +1 and -1 are degenerate but not coupled by S, and the run has to finish with finite values (it wrote NaN);
without ZFS the levels 0 and +1 are degenerate and coupled, and the run has to be refused. Review 20261012-233058-452278:
the same degenerate levels with a zero hyperfine (gyro 0) have no coupling, so the run has to finish with finite values
(P1-02: the check looked at S alone and refused it); and ZFS `diag(0,0,1)` and the same plus 1e13 I (an energy shift
only) have to get the same decision, both accepted and finite (P1-03: the degeneracy tolerance followed the largest
level, not their spread, and refused the shifted one).

Two runs check the Hamiltonian switches added on 2026-10-13 (A47), each against what the switch must give, not against
another CCEX run alone. `hamiltonian.cCCE.delta_A` false: two ¹³C on the NV axis 1.54 Å apart at 14.5 Å (`inputs/bath_pair`;
axial tensors, so one spin alone gives L = 1 in a Hahn echo) decay with delta_A true (max |L−1| 0.53) and give L = 1 to
every written decimal with false, and order 1 (clusters of one spin keep their own tensor) is byte for byte the run
without the key (two checks). `hamiltonian.common.bath_bath_mean_field` "off": two electron-like spins off the axis
(`inputs/bath_off`, one fixed state each, order 1) give the product of each spin alone in its bath (1.5e-10, the written
decimals), while the default "z" differs by 0.44 (one check).

Exit status 0 only if all ten pass. Everything written goes to `_work/`.
