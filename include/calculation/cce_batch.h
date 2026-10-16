/**
 * @file cce_batch.h
 * @brief Small cCCE clusters (bath dimension 2 and 4) evolved several at once on the CPU, one cluster per SIMD lane.
 */
#ifndef CCEX_CCE_BATCH_H
#define CCEX_CCE_BATCH_H
// A30: the time evolution of small cCCE clusters (bath dimension 2 and 4) on the CPU, several clusters at once - one
// cluster per SIMD lane (4 with AVX2, 8 with AVX-512) - with, per cluster, the arithmetic of the scalar path of
// calCoherenceCceInto (simulator_cce.cpp, simulator_cce_recursive.cpp) in the same order, so that the coherences are the
// same bits. The setup of each cluster (Hamiltonians, eigensystems, M, start vectors) stays the scalar path's own
// (cceSetup); only what follows it is batched.
//
// The kernel (cce_batch_kernel.cpp, compiled with -ffp-contract=off: what it fuses is what it writes) reproduces how
// Eigen and the compiler compute the scalar path, which depends on the compiler and the instruction set: which Eigen
// packets a size-4 object uses (AVX-512 fuses their products), how a gemv sums, which scalar std::complex products
// the compiler fuses. It knows the builds it was checked on (variant()); any other build has lanes() == 0 and every
// cluster goes through the scalar path. test/code_verification/check_numerics/test_cce_batch.cpp compares the two bit for
// bit on dumped clusters (make cpu-test).
#include <complex>

namespace ccebatch {
constexpr int maxLanes = 8;
constexpr int maxIntervals = 16;

// What the kernel evolves: one state with any pulse plan of at most maxIntervals intervals, or the ensemble Hahn echo.
struct Plan {
    bool ensemble = false, recursive = false;
    int intervals = 0, nstep = 0;      // one state: npulse+1 intervals
    double deltat = 0.0;
    double fraction[maxIntervals];     // Pulse sequence[p][2]
    int same[maxIntervals];            // Pulse sequence_indices[p]
};

// One window of lanes, structure of arrays: [element][lane]. Element (i,k) of a matrix at i*D+k.
struct alignas(64) Pack {
    double Mre[16][maxLanes], Mim[16][maxLanes];           // M = M_up^dagger M_down
    double Are[16][maxLanes], Aim[16][maxLanes];           // M^dagger (one state)
    double r0re[4][maxLanes], r0im[4][maxLanes];           // psi^dagger M_down (one state)
    double c0re[4][maxLanes], c0im[4][maxLanes];           // M_up^dagger psi (one state)
    double Eu[4][maxLanes], Ed[4][maxLanes];               // the eigenvalues
};

// The lanes of this build's kernel, 0 when the build is not one the kernel was checked on.
int lanes();
// The build the kernel takes the arithmetic of ("g++ 13 AVX2", ...), or "none".
const char* variant();
// Whether this build's kernel evolves this plan for bath dimension D (2 or 4).
bool covers(const Plan& plan, int D);
// The coherences of the lanes: out[l][0..nstep) for lane l; a null out[l] (a padded lane) is not written.
void evolve(const Plan& plan, int D, const Pack& pack, std::complex<double>* const* out);
}
#endif
