#ifndef CCEX_WARP_CLUSTER_H
#define CCEX_WARP_CLUSTER_H
// A40: the warp kernel - one cluster per warp of 32 threads, its matrices in shared memory (warp_cluster_gpu.cu) - for
// clusters of total dimension up to 32 that the packed kernel (one thread per cluster, small_cluster_math.h) cannot hold
// in registers: gCCE with a qubit of dimension 3 (6, 9, 12, 18, 21, 24, 27) and cCCE of bath dimension 8 and 16, up to
// gpu_options.warp_max_dim (default 27 in the GPU build with backend auto or cuda; 0 off). This header is the part both
// sides see: the time plan, the sizes of a cluster's input and output, which dimensions are compiled. Plain C++, no CUDA.
//
// What one cluster computes:
//   gCCE   the coherence (or, quantity dm, the reduced qubit density matrix) of calCoherenceGcceInto: H = V diag(E) V^dagger
//          (parallel cyclic Jacobi, the lower triangle read as Eigen does), then for every time point the evolution=vector
//          chain of propagateStateGcce - exp(-iH tau_0) psi0, then pulse p and exp(-iH tau_p) for every pulse - and
//          psia^dagger Psi Psi^dagger psib with Psi(i,b) the state as a qdim x bdim matrix. A single state starts from
//          psi0 = qubit psi0 (x) bath psi0; the ensemble (evolution=matrix, rho0 = qubit rho0 (x) 1/bdim) is the same chain
//          averaged over the bdim starts qubit psi0 (x) e_m - the density-matrix path's result in other arithmetic.
//   cCCE   calCoherenceCceInto's coherence, as small_cluster_math.h computes it (the same chain, the Hahn-ensemble identity).
// Phases: exact always (sin and cos at every time point; the gCCE ones with the tail of refinedPhases), whatever the run's
// evolution_update - the step-0 measurement (job 698) found the warp kernel 3 times slower with recursive phases.
#include "small_cluster_math.h"

namespace warpcluster {

using smallcluster::C;
constexpr int maxDim = 32;          // total dimension (gCCE: qubit x bath; cCCE: bath)
constexpr int maxQubitDim = 4;
constexpr int maxIntervals = smallcluster::maxIntervals;
enum Mode { modeCce = 0, modeGcce = 1 };

struct Plan {
    int mode;                       // modeCce or modeGcce
    int ensemble;                   // gCCE: nstate = 0 (evolution=matrix); cCCE: nstate = 0
    int dm;                         // gCCE quantity dm: qdim x qdim values per time point (column-major), else 1
    int qdim;                       // gCCE qubit dimension (cCCE: 1)
    smallcluster::Plan time;        // intervals, distinct fractions, nstep, deltat (cCCE: as smallClusterPlan makes it;
                                    // gCCE: the fractions of the pulse sequence, as propagateStateGcce reads them)
    C pulse[maxIntervals-1][maxQubitDim*maxQubitDim];   // gCCE: pulse p on the qubit, qdim x qdim column-major
    C psia[maxQubitDim], psib[maxQubitDim], qpsi0[maxQubitDim];   // gCCE: the qubit states (qpsi0: ensemble starts)
};

// Values (complex) of one cluster's input: gCCE H (D x D, column-major) then psi0 (D) for a single state; cCCE Halpha,
// Hbeta (D x D each) then psi (D) for a single state.
inline long inputSize(const Plan& p, int dim){
    return (p.mode == modeGcce ? 1L : 2L)*dim*dim + (p.ensemble ? 0 : dim);
}
// Values per time point of one cluster's output
inline int outputPer(const Plan& p){ return p.mode == modeGcce && p.dm ? p.qdim*p.qdim : 1; }
// Whether the kernel is compiled for this total dimension. gCCE: a qubit of dimension 3 with a bath of dimension 2, 3, 4,
// 6, 7, 8 or 9 - the products of spin multiplicities 2, 3, 4 and 7 up to 30 (5 and 10 would need a spin-2 bath spin;
// they are left out, as is a qubit of another dimension); cCCE: bath dimension 8 and 16.
inline bool compiled(int mode, int qdim, int dim){
    if (mode == modeCce) return dim == 8 || dim == 16;
    return qdim == 3 && (dim == 6 || dim == 9 || dim == 12 || dim == 18 || dim == 21 || dim == 24 || dim == 27);
}

} // namespace warpcluster
#endif
