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
#include <complex>
#include <vector>

// The kernel's plan for this run; false when the kernel does not evolve it (an ensemble other than the Hahn echo, more
// than ccebatch::maxIntervals intervals).
bool ccebatchPlan(Config* cnf, Pulse* pls, ccebatch::Plan& plan);
// Lane l of p from a cluster's setup (cceSetup's ws); false when a value is not finite (that cluster goes the scalar way).
bool ccebatchPackLane(const CceWorkspace& ws, bool ensemble, int D, ccebatch::Pack& p, int l);
// Evolves lanes 0..count-1 of p into out[0..count) (nstep values each); the other lanes repeat lane count-1.
void ccebatchEvolvePadded(const ccebatch::Plan& plan, int D, ccebatch::Pack& p, int count, std::complex<double>* const* out);

class CceBatchWindow {
public:
    // On when the build has a kernel (ccebatch::lanes() > 0) and the run is one it covers: cCCE with one qubit, quantity
    // coherence, one state with a pulse plan of at most ccebatch::maxIntervals intervals or the ensemble Hahn echo; the
    // packed GPU batch off (gpuBatch false), backend cpu in a CUDA build, and CCEX_SMALL_CLUSTER_DUMP not set (the dump
    // writes the scalar path's clusters), and CCEX_CPU_BATCH=0 not set (A/B). Off: nothing is allocated and get() returns nullptr.
    CceBatchWindow(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, bool isGCCE, bool gpuBatch);
    ~CceBatchWindow();
    CceBatchWindow(const CceBatchWindow&) = delete;
    CceBatchWindow& operator=(const CceBatchWindow&) = delete;
    bool on() const { return enabled; }
    // the local clusters of order n of the current state (clusters[ic] = [power, spins...], 1 <= ic <= ncluster);
    // use false: this order's clusters all go the scalar way (savemode cluster_coherence/allfull, except, the rdip ladder,
    // cluster contributions)
    void beginOrder(int n, int** clusters, int ncluster, bool use);
    // cluster ic's nstep coherences and its dimension when the batch computed it, else nullptr - and then *made is its
    // bath array when the window made it (the caller computes and frees it), else nullptr
    const std::complex<double>* get(int ic, int& dim, BathArray** made);
private:
    void fill(int ic);
    void run(int k);
    void drop();
    QubitArray* qa; BathArray* ba; DefectArray* dfa; Config* cnf;
    bool enabled = false, use = false, dimOk[2] = {false,false};
    int B = 0, W = 0, nstep = 0, nqubit = 1, order = 0, ncluster = 0, first = 1, last = 0;
    int** clusters = nullptr;
    ccebatch::Plan plan;
    CceWorkspace ws;
    std::vector<std::complex<double>> results;      // W x nstep, by position in the window
    std::vector<int> dims;                           // by position: 2 or 4 when batched, 0 when not
    std::vector<BathArray*> made;                    // by position: the bath array of a cluster left to the scalar path
    ccebatch::Pack* pack = nullptr;                  // [2]: dimension 2 and 4
    int filled[2] = {0,0};
    int slot[2][ccebatch::maxLanes];                 // window position of each lane
    double build[2][ccebatch::maxLanes], setup[2][ccebatch::maxLanes];
};
#endif
