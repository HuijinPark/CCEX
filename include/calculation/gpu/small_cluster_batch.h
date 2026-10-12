#ifndef CCEX_SMALL_CLUSTER_BATCH_H
#define CCEX_SMALL_CLUSTER_BATCH_H
// A26 step 1: small cCCE clusters (bath dimension 2 and 4) computed on the GPU in large batches, one thread per cluster
// (small_cluster_gpu.cu): the packed kernel, "gpu_options": {"packed_max_dim": 2 or 4} (0, the default: off, and nothing of
// it runs). The rest of the run is unchanged.
//
// calculate() asks get(ic) for each local cluster. The first time a window of clusters is reached, the rank builds their
// Hamiltonians (cceClusterHamiltonians, the CPU path's own function) and packs them into its part of a node shared-memory
// window; before waiting for that window it packs the next one, so that the card computes one window while the rank builds
// the next and assembles the one before. One process per node drives the card: node rank 0 runs a feeder thread that
// collects the packed windows of all ranks of the node, computes them in combined batches and writes the coherences back
// into the shared memory. Every rank keeps its clusters and does the rest of its work (Hamiltonians, the product, the MPI
// reduction) as before; only the eigensystems and the time evolution move to the card. One CUDA context per node instead of
// one per rank: 32 ranks starting CUDA at once took 3.5 s and their teardown about 2 s more on idm4-gpu (A26 step 1, job 683).
// A cluster the batch does not take (another dimension, a Hamiltonian that is not finite and Hermitian, a Jacobi that did
// not converge, and the clusters met while the card is still starting) gets nullptr, and calculate()
// computes it on the CPU as before.
#include "calculation/simulator.h"
#include "calculation/gpu/small_cluster_math.h"
#include <complex>

class SmallClusterBatch {
public:
    // Collective (every rank calls it at the same place). Active with gpu_options.packed_max_dim 2 or 4, method CCE with one
    // qubit, and a pulse plan of at most 4 distinct intervals. Asked for in a build without CUDA or with another method the
    // run stops in smallClusterWarmUp; without a usable card it stops when the first window is waited for.
    SmallClusterBatch(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, bool isGCCE);
    ~SmallClusterBatch();
    SmallClusterBatch(const SmallClusterBatch&) = delete;
    SmallClusterBatch& operator=(const SmallClusterBatch&) = delete;
    bool on() const { return impl != nullptr; }
    // the local clusters of order n (localclusters[n]: [0][0]-1 clusters, each [power, spins...]) of the current state
    void beginOrder(int n, int** clusters, int ncluster);
    // cluster ic's nstep coherences and its dimension when the GPU computed it, else nullptr
    const std::complex<double>* get(int ic, int& dim);
    // collective, after the last state: stops the feeder, releases the shared memory and prints one [gpu-small] line by
    // rank 0 (counts and seconds over the ranks); the batch is off afterwards
    void report();
private:
    struct Impl;
    Impl* impl = nullptr;
};

// Collective, early in the run (computeBackendInit, before the bath is read and the clusters are made): when the packed
// kernel is asked for, every rank allocates its part of the node shared memory and node rank 0 starts the CUDA context
// (and page-locks that memory) in a thread of its own, so that the feeder finds it ready (about 0.2 s on
// idm4-gpu, otherwise in the first state). smallClusterWarmDown (computeBackendFinish) joins it if no batch took it.
void smallClusterWarmUp(const Config* cnf);
void smallClusterWarmDown();

// The time plan of the run as the batch uses it (interval fractions, distinct slots); false if it has more than
// 4 distinct intervals or more than maxIntervals. Both builds.
bool smallClusterPlan(Config* cnf, Pulse* pls, smallcluster::Plan& plan, bool& ensemble, bool& recursive);
// Test data for the GPU check (test_small_cluster_cuda): with CCEX_SMALL_CLUSTER_DUMP=<file>, rank 0 appends its first
// CCEX_SMALL_CLUSTER_DUMP_N (default 400) cCCE clusters of dimension 2 and 4 - plan, Halpha, Hbeta, psi and the CPU
// coherence just computed. Any build; nothing happens without the variable.
void smallClusterDump(QubitArray* qa, BathArray* cluster, Config* cnf, Pulse* pls, const MatrixXcd* result);
#endif
