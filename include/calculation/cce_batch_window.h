#ifndef CCEX_CCE_BATCH_WINDOW_H
#define CCEX_CCE_BATCH_WINDOW_H
// A30: calculate()'s side of the batched small-cluster cCCE on the CPU (cce_batch.h). When cluster ic is asked for and is
// not in the current window, the window moves to ic, ic+1, ... (this rank's local clusters in their order): each cluster's
// bath array is made, and a cluster of bath dimension 2 or 4 gets its Hamiltonians and setup (cceClusterHamiltonians,
// cceSetup - the scalar path's own) and a lane; full groups of lanes are evolved together, the rest of each dimension at
// the window's end. calculate() then takes the clusters one by one in their order, so the assembly sees the same factors in
// the same order. A cluster the batch does not take (another dimension, a setup with a value that is not finite) keeps its
// bath array for the scalar path.
#include "calculation/simulator.h"
#include "calculation/cce_batch.h"
#include <climits>
#include <complex>
#include <vector>

// The kernel's plan for this run; false when the kernel does not evolve it (an ensemble other than the Hahn echo, more
// than ccebatch::maxIntervals intervals).
bool ccebatchPlan(Config* cnf, Pulse* pls, ccebatch::Plan& plan);
// Lane l of p from a cluster's setup (cceSetup's ws); false when a value is not finite (that cluster goes the scalar way).
bool ccebatchPackLane(const CceWorkspace& ws, bool ensemble, int D, ccebatch::Pack& p, int l);
// Evolves lanes 0..count-1 of p into out[0..count) (nstep values each); the other lanes repeat lane count-1.
void ccebatchEvolvePadded(const ccebatch::Plan& plan, int D, ccebatch::Pack& p, int count, std::complex<double>* const* out);

// A51: whether the kernel may be used in this process: the build has one (ccebatch::lanes() > 0) and, checked once at the
// first call, a few made-up clusters (bath dimension 2 and 4; one state - Ramsey, Hahn echo, a 3-interval CPMG - and the
// ensemble Hahn echo; exact and recursive; a full and a padded group of lanes) come out of the batched path with the bits of
// the scalar path (cceSetup/cceSetupSmall + cceEvolve). The kernel reproduces the arithmetic of the builds it was checked
// on (cce_batch.h), but the scalar path's fused products also follow flags it cannot see (-ffp-contract=off, -march's
// tuning; review 20261013-cpu R1, CI run 38366397340); where they differ the batch is off and every cluster goes the scalar
// way, so the results are the scalar path's everywhere. About a millisecond. With MPI initialised every rank calls it at
// the same place the first time (collective: the ranks take the AND of their checks, so that all decide alike); rank 0
// prints one line when it is off for any reason but the build. why (may be null): the reason when false.
bool ccebatchUsable(const char** why = nullptr);
// A48: whether CceBatchWindow takes this run's clusters of bath dimension 2 (dimOk[0]) and 4 (dimOk[1]), whatever the GPU
// batch does: cCCE (not gCCE) with one qubit, quantity coherence, a plan the kernel evolves (ccebatchPlan), not backend
// cuda (every cluster on the card), not propagator sparse or the CUDA cluster path for that dimension, neither
// CCEX_SMALL_CLUSTER_DUMP (the dump writes the scalar path's clusters) nor CCEX_CPU_BATCH=0 (A/B), and ccebatchUsable.
// The checks are the same on every rank, so that ccebatchUsable is reached by all or none. false with why otherwise.
bool ccebatchRunTakes(QubitArray* qa, Config* cnf, Pulse* pls, bool isGCCE, ccebatch::Plan& plan, bool dimOk[2], const char** why);

class CceBatchWindow {
public:
    // On when ccebatchRunTakes (collective the first time, see ccebatchUsable). In the GPU build it runs beside the packed
    // and warp kernels (A48): it is asked only for the clusters the GPU batch leaves to the CPU (get's end). Off when one
    // group of lanes' coherences alone would pass 4 MiB (long time grids). The storage is made by the first window
    // (review 20261013-221715 X5). Off: nothing is allocated and get() returns nullptr.
    CceBatchWindow(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, bool isGCCE);
    ~CceBatchWindow();
    CceBatchWindow(const CceBatchWindow&) = delete;
    CceBatchWindow& operator=(const CceBatchWindow&) = delete;
    bool on() const { return enabled; }
    // the local clusters of order n of the current state (clusters[ic] = [power, spins...], 1 <= ic <= ncluster);
    // use false: this order's clusters all go the scalar way (savemode cluster_coherence/allfull, except, the rdip ladder,
    // cluster contributions)
    void beginOrder(int n, int** clusters, int ncluster, bool use);
    // cluster ic's nstep coherences and its dimension when the batch computed it, else nullptr - and then *made is its
    // bath array when the window made it (the caller computes and frees it), else nullptr. A window that has to move to ic
    // stops before end (A48: SmallClusterBatch::cpuEnd, so that it takes no cluster the GPU batch will compute).
    const std::complex<double>* get(int ic, int& dim, BathArray** made, int end = INT_MAX);
private:
    void fill(int ic, int end);
    void allocate(int n);
    void run(int k);
    void drop();
    QubitArray* qa; BathArray* ba; DefectArray* dfa; Config* cnf;
    bool enabled = false, use = false, dimOk[2] = {false,false};
    int B = 0, W = 0, capacity = 0, nstep = 0, nqubit = 1, order = 0, ncluster = 0, first = 1, last = 0;
    int** clusters = nullptr;
    ccebatch::Plan plan;
    CceWorkspace ws;
    std::vector<std::complex<double>> results;      // W x nstep, by position in the window
    std::vector<int> dims;                           // by position: 2 or 4 when batched, 0 when not
    std::vector<BathArray*> made;                    // by position: the bath array of a cluster left to the scalar path
    // A42: by position, the Hamiltonians and state of a batched cluster, and by 2*position+s (s 0 Halpha, 1 Hbeta) the
    // Jacobi eigensystems (values [4], vectors [16]) and whether they converged
    std::vector<MatrixXcd> Ha, Hb, psiAt;
    std::vector<double> buildAt, eigW;
    std::vector<std::complex<double>> eigV;
    std::vector<char> eigOk;
    double eigEach[2] = {0.0,0.0};                   // profile: the eigensolves of the last group, per cluster
    ccebatch::Pack* pack = nullptr;                  // [2]: dimension 2 and 4
    int filled[2] = {0,0};
    int slot[2][ccebatch::maxLanes];                 // window position of each lane
    double build[2][ccebatch::maxLanes], setup[2][ccebatch::maxLanes];
};
#endif
