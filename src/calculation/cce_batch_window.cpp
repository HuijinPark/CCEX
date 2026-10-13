// A30: the window of small cCCE clusters calculate() hands to the batched kernel (cce_batch_window.h).
#include "../../include/calculation/cce_batch_window.h"
#include "../../include/calculation/small_hamiltonian.h"
#include "../../include/calculation/sparse_propagator.h"
#include "../../include/calculation/gpu/compute_backend.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>

bool ccebatchPlan(Config* cnf, Pulse* pls, ccebatch::Plan& plan){
    plan = ccebatch::Plan{};
    plan.ensemble = Config_getNstate(cnf) == 0;
    plan.recursive = strcasecmp(cnf->evolution_update,"recursive") == 0;   // as calCoherenceCceInto decides
    plan.nstep = Config_getNstep(cnf);
    plan.deltat = (double)Config_getDeltat(cnf);
    const int npulse = Pulse_getNpulse(pls);
    if (plan.ensemble){
        // the Hahn echo only: the CPMG/Ramsey/Manual ensemble branches stay scalar
        plan.intervals = npulse+1;
        return npulse == 1 && strcasecmp(Pulse_getPulsename(pls),"Manual") != 0;
    }
    if (npulse+1 > ccebatch::maxIntervals) return false;
    plan.intervals = npulse+1;
    double** sequence = Pulse_getSequence(pls);
    int* same = Pulse_getSequenceIndices(pls);
    for (int p=0; p<=npulse; ++p){ plan.fraction[p] = sequence[p][2]; plan.same[p] = same[p]; }
    return true;
}

bool ccebatchPackLane(const CceWorkspace& ws, bool ensemble, int D, ccebatch::Pack& p, int l){
    bool finite = true;
    auto put = [&](double (*re)[ccebatch::maxLanes], double (*im)[ccebatch::maxLanes], int e, doublec z){
        re[e][l] = z.real(); im[e][l] = z.imag();
        finite = finite && std::isfinite(z.real()) && std::isfinite(z.imag());
    };
    for (int i=0; i<D; ++i)
        for (int k=0; k<D; ++k){
            put(p.Mre,p.Mim,i*D+k,ws.M(i,k));
            if (!ensemble) put(p.Are,p.Aim,i*D+k,ws.Madj(i,k));
        }
    for (int k=0; k<D; ++k){
        if (!ensemble){ put(p.r0re,p.r0im,k,ws.rowStart(k)); put(p.c0re,p.c0im,k,ws.colStart(k)); }
        p.Eu[k][l] = ws.valUp(k).real(); p.Ed[k][l] = ws.valDown(k).real();
        finite = finite && std::isfinite(p.Eu[k][l]) && std::isfinite(p.Ed[k][l]);
    }
    return finite;
}

CceBatchWindow::CceBatchWindow(QubitArray* qa_, BathArray* ba_, DefectArray* dfa_, Config* cnf_, Pulse* pls, bool isGCCE,
                               bool gpuBatch) : qa(qa_), ba(ba_), dfa(dfa_), cnf(cnf_){
    B = ccebatch::lanes();
    if (B == 0 || isGCCE || gpuBatch || getenv("CCEX_SMALL_CLUSTER_DUMP")) return;
    if (getenv("CCEX_CPU_BATCH") && !strcmp(getenv("CCEX_CPU_BATCH"),"0")) return;   // A/B: the scalar path for every cluster
    if (QubitArray_getNqubit(qa) != 1 || strcasecmp(Config_getQuantity(cnf),"coherence")) return;
#ifdef CCEX_USE_CUDA
    if (strcasecmp(cnf->backend,"cpu")) return;          // the CUDA paths keep their clusters
#endif
    if (!ccebatchPlan(cnf,pls,plan)) return;
    for (int k=0; k<2; ++k){
        const int d = k == 0 ? 2 : 4;
        dimOk[k] = ccebatch::covers(plan,d) && !sparsePropagatorFor(cnf,d) && !computeUseCuda(cnf,d);
    }
    if (!dimOk[0] && !dimOk[1]) return;
    nstep = plan.nstep;
    nqubit = QubitArray_getNqubit(qa);
    // the coherences of a window: 256 clusters, fewer when they would pass 4 MiB, at least 4 groups of lanes
    W = int(std::max<long>(4L*B, std::min<long>(256, (4L<<20)/(16L*std::max(1,nstep)))));
    results.resize(size_t(W)*nstep);
    dims.assign(W,0);
    made.assign(W,nullptr);
    pack = static_cast<ccebatch::Pack*>(std::aligned_alloc(64, 2*sizeof(ccebatch::Pack)));
    std::memset(static_cast<void*>(pack), 0, 2*sizeof(ccebatch::Pack));
    enabled = true;
}

CceBatchWindow::~CceBatchWindow(){
    drop();
    std::free(pack);
}

void CceBatchWindow::drop(){
    for (BathArray*& b : made) if (b){ BathArray_freeAll(b); b = nullptr; }
    first = 1; last = 0;
}

void CceBatchWindow::beginOrder(int n, int** c, int count, bool u){
    if (!enabled) return;
    drop();
    order = n; clusters = c; ncluster = count; use = u;
}

const std::complex<double>* CceBatchWindow::get(int ic, int& dim, BathArray** out){
    *out = nullptr;
    if (!enabled || !use) return nullptr;
    if (ic < first || ic > last) fill(ic);
    const int at = ic-first;
    if (dims[at]){ dim = dims[at]; return &results[size_t(at)*nstep]; }
    *out = made[at]; made[at] = nullptr;
    return nullptr;
}

void CceBatchWindow::fill(int ic){
    drop();
    first = ic; last = std::min(ncluster, ic+W-1);
    filled[0] = filled[1] = 0;
    const bool profile = cnf->profile;
    for (int j=first; j<=last; ++j){
        const int at = j-first;
        dims[at] = 0;
        BathArray* b = createBathArray(clusters[j]+1, order, ba, dfa, nqubit);
        const int d = BathArray_dim(b);
        const int k = d == 2 ? 0 : d == 4 ? 1 : -1;
        if (k < 0 || !dimOk[k]){ made[at] = b; continue; }
        const double t0 = profile ? MPI_Wtime() : 0.0;
        cceClusterHamiltonians(qa,b,cnf,ws.Halpha,ws.Hbeta);
        const double t1 = profile ? MPI_Wtime() : 0.0;
        MatrixXcd psi0;
        if (!plan.ensemble) psi0 = BathArray_Psi0(b);
        cceSetup(ws, plan.ensemble ? nullptr : &psi0);
        const int lane = filled[k];
        if (!ccebatchPackLane(ws,plan.ensemble,d,pack[k],lane)){ made[at] = b; continue; }   // the scalar path decides
        BathArray_freeAll(b);
        slot[k][lane] = at; dims[at] = d;
        if (profile){ build[k][lane] = t1-t0; setup[k][lane] = MPI_Wtime()-t1; }
        if (++filled[k] == B) run(k);
    }
    for (int k=0; k<2; ++k) if (filled[k]) run(k);
}

void ccebatchEvolvePadded(const ccebatch::Plan& plan, int D, ccebatch::Pack& p, int count, std::complex<double>* const* out){
    double (*rows)[ccebatch::maxLanes] = reinterpret_cast<double (*)[ccebatch::maxLanes]>(&p);
    const size_t nrows = sizeof(ccebatch::Pack)/sizeof(rows[0]);
    for (int l=count; l<ccebatch::lanes(); ++l) for (size_t r=0; r<nrows; ++r) rows[r][l] = rows[r][count-1];
    std::complex<double>* lanes[ccebatch::maxLanes] = {};
    for (int l=0; l<count; ++l) lanes[l] = out[l];
    ccebatch::evolve(plan,D,p,lanes);
}

// Evolves the filled lanes of dimension 2 (k 0) or 4 (k 1)
void CceBatchWindow::run(int k){
    const int count = filled[k], d = k == 0 ? 2 : 4;
    std::complex<double>* out[ccebatch::maxLanes] = {};
    for (int l=0; l<count; ++l) out[l] = &results[size_t(slot[k][l])*nstep];
    const double t0 = cnf->profile ? MPI_Wtime() : 0.0;
    ccebatchEvolvePadded(plan,d,pack[k],count,out);
    if (cnf->profile){
        const double each = (MPI_Wtime()-t0)/count;
        for (int l=0; l<count; ++l) computeProfile(cnf,"cce","cpu-batch",d,build[k][l],setup[k][l],each);
    }
    filled[k] = 0;
}
