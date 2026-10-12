// A26 step 1: the host side of the batched small-cluster cCCE on the GPU (small_cluster_batch.h says how it is used).
#include "../../../include/calculation/gpu/small_cluster_batch.h"
#include "../../../include/calculation/hamiltonian_check.h"
#include "../../../include/calculation/small_hamiltonian.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
#include "../../../include/calculation/gpu/small_cluster_gpu.h"
#endif

using smallcluster::C;
using smallcluster::Plan;

bool smallClusterPlan(Config* cnf, Pulse* pls, Plan& plan, bool& ensemble, bool& recursive){
    plan = Plan{};
    const int npulse = Pulse_getNpulse(pls);
    if (npulse+1 > smallcluster::maxIntervals) return false;
    ensemble = Config_getNstate(cnf)==0;
    recursive = !strcasecmp(cnf->evolution_update,"recursive");
    // The ensemble CPMG/Hahn/Ramsey branch of calCoherenceCceInto takes tau = t/(2 npulse) at both ends and twice that in
    // between, Ramsey t; a single state and a Manual sequence take the sequence's fractions.
    const bool manual = !strcasecmp(Pulse_getPulsename(pls),"Manual");
    double** sequence = Pulse_getSequence(pls);
    plan.intervals = npulse+1;
    plan.nstep = Config_getNstep(cnf);
    plan.deltat = (double)Config_getDeltat(cnf);
    for (int p=0; p<=npulse; ++p){
        double f;
        if (ensemble && !manual) f = npulse==0 ? 1.0 : (p==0 || p==npulse) ? 1.0/(2.0*npulse) : 1.0/npulse;
        else f = sequence[p][2];
        int slot = -1;
        for (int s=0; s<plan.distinct; ++s) if (plan.fraction[s]==f){ slot = s; break; }
        if (slot < 0){
            if (plan.distinct == 4) return false;
            slot = plan.distinct++;
            plan.fraction[slot] = f;
        }
        plan.same[p] = slot;
    }
    return true;
}

void smallClusterDump(QubitArray* qa, BathArray* cluster, Config* cnf, Pulse* pls, const MatrixXcd* result){
    static const char* path = getenv("CCEX_SMALL_CLUSTER_DUMP");
    if (!path || rank != 0 || strcasecmp(cnf->method,"cce")) return;
    static long left = getenv("CCEX_SMALL_CLUSTER_DUMP_N") ? atol(getenv("CCEX_SMALL_CLUSTER_DUMP_N")) : 400;
    if (left <= 0) return;
    const int d = BathArray_dim(cluster);
    if (d != 2 && d != 4) return;
    Plan plan; bool ensemble, recursive;
    if (!smallClusterPlan(cnf,pls,plan,ensemble,recursive)) return;
    static FILE* fp = fopen(path,"ab");
    if (!fp){ perror(path); ccexFatal(); }
    MatrixXcd Ha, Hb;
    cceClusterHamiltonians(qa,cluster,cnf,Ha,Hb);
    const int head[6] = {d, ensemble, recursive, plan.intervals, plan.distinct, plan.nstep};
    fwrite(head,sizeof(int),6,fp);
    fwrite(&plan.deltat,sizeof(double),1,fp);
    fwrite(plan.fraction,sizeof(double),plan.distinct,fp);
    fwrite(plan.same,sizeof(int),plan.intervals,fp);
    fwrite(Ha.data(),sizeof(doublec),size_t(d)*d,fp);
    fwrite(Hb.data(),sizeof(doublec),size_t(d)*d,fp);
    if (!ensemble){ const MatrixXcd psi = BathArray_Psi0(cluster); fwrite(psi.data(),sizeof(doublec),d,fp); }
    for (int i=0; i<plan.nstep; ++i){ const doublec v = result[i](0,0); fwrite(&v,sizeof(doublec),1,fp); }
    fflush(fp);
    if (--left == 0){ fclose(fp); fp = nullptr; }
}

#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
namespace {
// The state of one slot of one rank, in the node's shared memory: the rank packs and submits, the feeder runs it, the rank
// reads it back. std::atomic<int> is lock-free here, so it works between the processes of the node.
enum SlotState { slotIdle = 0, slotSubmitted = 1, slotRunning = 2, slotDone = 3 };
struct alignas(64) SlotHeader { std::atomic<int> state; long count[2]; };
struct alignas(64) FeederHeader { std::atomic<int> state; char message[248]; };   // in node rank 0's part: 0 starting, 1 ready, -1 failed
static_assert(std::atomic<int>::is_always_lock_free, "the slot flags need lock-free atomics");
constexpr size_t headerBytes = 2*sizeof(SlotHeader)+sizeof(FeederHeader);
size_t roundUp(size_t n, size_t to){ return (n+to-1)/to*to; }

// Window: 2048 clusters, fewer when the coherences of a slot would pass 6 MiB. The shared memory is allocated, page-locked
// and freed once per run, about 0.25 s per GB each on idm4-gpu; windows of 2048, 4096 and 8192 ran P1 CR60 cCCE2 (2 states)
// in about the same time and 16384 and 32768 slower (A26 step 1, jobs 684-688).
int windowFor(int nstep){
    return int(std::max<long>(256,std::min<long>(2048,(6L<<20)/(16L*std::max(1,nstep)))));
}

// The node's shared memory: every rank's part is the headers, then two slots of (in2, in4, out, status).
struct NodeMemory {
    MPI_Comm node = MPI_COMM_NULL; MPI_Win win = MPI_WIN_NULL;
    int nodeRank = 0, nodeSize = 1, windowSize = 0, per2 = 0, per4 = 0, nstep = 0;
    size_t slotBytes = 0, partBytes = 0;
    std::vector<char*> parts;                  // every rank's part (as this process sees it)
    bool allocated() const { return win != MPI_WIN_NULL; }
    bool matches(int n, bool ensemble, int window) const { return allocated() && nstep==n && per2==(ensemble ? 8 : 10) && windowSize==window; }
    void allocate(int n, bool ensemble, int window){       // collective
        nstep = n; windowSize = window;
        const int extra = ensemble ? 0 : 1;
        per2 = 8+2*extra; per4 = 32+4*extra;
        slotBytes = roundUp(sizeof(C)*size_t(windowSize)*(per2+per4+nstep)+size_t(windowSize),64);
        partBytes = roundUp(headerBytes+2*slotBytes,4096);
        MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
        MPI_Comm_rank(node,&nodeRank);
        MPI_Comm_size(node,&nodeSize);
        char* base = nullptr;
        MPI_Win_allocate_shared(MPI_Aint(partBytes),1,MPI_INFO_NULL,node,&base,&win);
        parts.resize(nodeSize);
        for (int r=0; r<nodeSize; ++r){
            MPI_Aint size; int unit; char* p;
            MPI_Win_shared_query(win,r,&size,&unit,&p);
            parts[r] = p;
        }
        char* own = parts[nodeRank];
        for (int s=0; s<2; ++s){ SlotHeader* h = new (own+s*sizeof(SlotHeader)) SlotHeader; h->state.store(slotIdle); h->count[0] = h->count[1] = 0; }
        if (nodeRank==0){ FeederHeader* fh = new (own+2*sizeof(SlotHeader)) FeederHeader; fh->state.store(0); fh->message[0] = 0; }
        MPI_Barrier(node);
        // one contiguous range (MPI_Win_allocate_shared without alloc_shared_noncontig), page-locked in one call
        for (int r=1; r<nodeSize; ++r)
            if (parts[r] != parts[0]+r*partBytes){
                if (rank==0) fprintf(stderr,"Error: the small-cluster GPU batch: the node's shared memory is not contiguous\n");
                ccexFatal();
            }
    }
    void release(){                                          // collective
        if (!allocated()) return;
        MPI_Win_free(&win);
        MPI_Comm_free(&node);
        parts.clear();
    }
    size_t bytes() const { return partBytes*size_t(nodeSize); }
};

// smallClusterWarmUp: the shared memory allocated early, and on node rank 0 the CUDA context opened and that memory
// page-locked in a thread, while the bath is read and the clusters are made
struct Warm {
    std::thread thread;
    GpuSmallDevice* device = nullptr;
    char description[200] = {};
    double openSeconds = 0, lockSeconds = 0, allocateSeconds = 0;
    bool locked = false;
    NodeMemory memory;
};
Warm warm;
}

struct SmallClusterBatch::Impl {
    QubitArray* qa; BathArray* ba; DefectArray* dfa; Config* cnf;
    Plan plan{};
    bool ensemble = false, recursive = false;
    int maxDim = 4;                            // gpu_options.packed_max_dim
    NodeMemory mem;
    bool fromWarm = false;                     // mem is the one smallClusterWarmUp made
    char* own = nullptr;
    // this rank's windows
    struct Window { int first = 0, last = 0; bool submitted = false, waited = false; std::vector<int> where; std::vector<char> dim; long count[2] = {0,0}; };
    Window window[2];
    int order = 0, ncluster = 0; int** clusters = nullptr;
    MatrixXcd Ha, Hb;
    long long gpuClusters = 0, fallback = 0, cpuStarting = 0;
    double prepareSeconds = 0, waitSeconds = 0;
    // the feeder (node rank 0)
    std::thread feeder;
    std::atomic<bool> stop{false};
    double feederOpen = 0, feederLock = 0, feederWaited = 0; long long feederBatches = 0, feederItems = 0; bool locked = false, early = false;

    SlotHeader* header(char* part, int slot){ return reinterpret_cast<SlotHeader*>(part)+slot; }
    FeederHeader* feederHeader(){ return reinterpret_cast<FeederHeader*>(mem.parts[0]+2*sizeof(SlotHeader)); }
    C* in2(char* part, int slot){ return reinterpret_cast<C*>(part+headerBytes+slot*mem.slotBytes); }
    C* in4(char* part, int slot){ return in2(part,slot)+size_t(mem.windowSize)*mem.per2; }
    C* out(char* part, int slot){ return in4(part,slot)+size_t(mem.windowSize)*mem.per4; }
    unsigned char* status(char* part, int slot){ return reinterpret_cast<unsigned char*>(out(part,slot)+size_t(mem.windowSize)*plan.nstep); }

    void feed();
    void prepare(int slot, int first);
    void wait(int slot);
    // whether the feeder is up (its context opened, the memory page-locked); stops the run if it failed
    bool ready(){
        FeederHeader* fh = feederHeader();
        const int state = fh->state.load(std::memory_order_acquire);
        if (state < 0){
            fprintf(stderr,"Error: the packed GPU kernel (feeder on node rank 0): %s\n",fh->message);
            ccexFatal();
        }
        return state == 1;
    }
    // While the feeder is starting, the next cpuWhileStarting clusters go to the CPU: a window that is never submitted
    void cpuWindow(int slot, int first){
        constexpr int cpuWhileStarting = 256;
        Window& w = window[slot];
        w.first = first; w.last = std::min(first+cpuWhileStarting, ncluster+1);
        w.where.assign(w.last-w.first,-1); w.dim.assign(w.last-w.first,0);
        w.count[0] = w.count[1] = 0;
        w.submitted = w.waited = true;
        cpuStarting += w.last-w.first;
    }
};

// The feeder: runs on node rank 0 beside its main thread. Takes every submitted slot of the node into one batch (up to
// maxClusters clusters), in one of two lanes, and marks the slots done when the lane's batch is back.
void SmallClusterBatch::Impl::feed(){
    FeederHeader* fh = feederHeader();
    GpuSmallDevice* gpu = nullptr;
    try {
        const auto t0 = std::chrono::steady_clock::now();
        if (warm.thread.joinable()) warm.thread.join();      // the context (and the memory's page-lock) of smallClusterWarmUp
        gpu = warm.device; warm.device = nullptr;
        early = gpu != nullptr;
        feederOpen = warm.openSeconds;
        if (fromWarm && warm.locked){ locked = true; warm.locked = false; feederLock = warm.lockSeconds; }
        if (!gpu){
            char description[200] = {};
            gpu = gpuSmallOpen(0,description,sizeof description);
            if (!gpu) throw std::runtime_error(*warm.description ? warm.description : description);
        }
        const auto t1 = std::chrono::steady_clock::now();
        if (!locked){ locked = gpuSmallRegister(mem.parts[0],mem.bytes()); feederLock = std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count(); }
        feederWaited = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        fh->state.store(1,std::memory_order_release);
        const long maxClusters = 1L<<19;
        std::vector<std::pair<int,int>> inLane[2];           // (rank, slot) of each lane's batch
        std::vector<GpuSmallItem> items;
        int idle = 0;
        while (true){
            bool progress = false;
            for (int lane=0; lane<2; ++lane){
                if (!inLane[lane].empty() && gpuSmallReady(gpu,lane)){
                    for (auto& rs : inLane[lane]) header(mem.parts[rs.first],rs.second)->state.store(slotDone,std::memory_order_release);
                    inLane[lane].clear();
                    progress = true;
                }
            }
            for (int lane=0; lane<2; ++lane){
                if (!inLane[lane].empty()) continue;
                items.clear();
                long total = 0;
                for (int r=0; r<mem.nodeSize; ++r) for (int s=0; s<2; ++s){
                    SlotHeader* h = header(mem.parts[r],s);
                    if (h->state.load(std::memory_order_acquire) != slotSubmitted) continue;
                    const long n = h->count[0]+h->count[1];
                    if (!items.empty() && total+n > maxClusters) continue;
                    h->state.store(slotRunning,std::memory_order_relaxed);
                    items.push_back(GpuSmallItem{in2(mem.parts[r],s),in4(mem.parts[r],s),h->count[0],h->count[1],
                                                 out(mem.parts[r],s),status(mem.parts[r],s)});
                    inLane[lane].push_back({r,s});
                    total += n;
                }
                if (items.empty()) break;
                gpuSmallSubmit(gpu,lane,plan,ensemble,recursive,items.data(),int(items.size()));
                ++feederBatches; feederItems += long(items.size());
                progress = true;
            }
            if (progress){ idle = 0; continue; }
            if (stop.load(std::memory_order_acquire) && inLane[0].empty() && inLane[1].empty()) break;
            // idle: 20 us between looks, up to 200 us after a while (it shares node rank 0's core)
            ++idle;
            std::this_thread::sleep_for(std::chrono::microseconds(std::min(200,20*(1+idle/64))));
        }
        if (locked) gpuSmallUnregister(mem.parts[0]);
        gpuSmallClose(gpu);
    } catch (const std::exception& e){
        snprintf(fh->message,sizeof fh->message,"%s",e.what());
        fh->state.store(-1,std::memory_order_release);
        if (gpu) gpuSmallClose(gpu);
    }
}

void SmallClusterBatch::Impl::prepare(int slot, int first){
    const double start = MPI_Wtime();
    Window& w = window[slot];
    w.first = first; w.last = std::min(first+mem.windowSize, ncluster+1);
    w.where.assign(w.last-w.first,-1); w.dim.assign(w.last-w.first,0);
    w.count[0] = w.count[1] = 0;
    // the slot was waited for (or never used): the feeder is done with it
    C* dst[2] = {in2(own,slot), in4(own,slot)};
    const int nqubit = QubitArray_getNqubit(qa);
    for (int ic=w.first; ic<w.last; ++ic){
        BathArray* b = createBathArray(clusters[ic]+1,order,ba,dfa,nqubit);
        const int d = BathArray_dim(b);
        if ((d == 2 || d == 4) && d <= maxDim){
            cceClusterHamiltonians(qa,b,cnf,Ha,Hb);
            if (finiteAndHermitian(Ha) && finiteAndHermitian(Hb)){
                const int k = d==4, per = k ? mem.per4 : mem.per2;
                C* p = dst[k]+w.count[k]*per;
                std::memcpy(p,Ha.data(),sizeof(C)*d*d);
                std::memcpy(p+d*d,Hb.data(),sizeof(C)*d*d);
                if (!ensemble){ const MatrixXcd psi = BathArray_Psi0(b); std::memcpy(p+2*d*d,psi.data(),sizeof(C)*d); }
                w.where[ic-w.first] = int(w.count[k]++);
                w.dim[ic-w.first] = char(d);
            } else ++fallback;
        } else ++fallback;
        BathArray_freeAll(b);
    }
    SlotHeader* h = header(own,slot);
    h->count[0] = w.count[0]; h->count[1] = w.count[1];
    w.submitted = true; w.waited = false;
    h->state.store(w.count[0]+w.count[1] > 0 ? slotSubmitted : slotDone,std::memory_order_release);
    gpuClusters += w.count[0]+w.count[1];
    prepareSeconds += MPI_Wtime()-start;
}

void SmallClusterBatch::Impl::wait(int slot){
    const double start = MPI_Wtime();
    SlotHeader* h = header(own,slot);
    FeederHeader* fh = feederHeader();
    long spins = 0;
    while (h->state.load(std::memory_order_acquire) != slotDone){
        if (fh->state.load(std::memory_order_acquire) < 0){
            fprintf(stderr,"Error: the packed GPU kernel (feeder on node rank 0): %s\n",fh->message);
            ccexFatal();
        }
        if (++spins > 64) std::this_thread::sleep_for(std::chrono::microseconds(10));
    }
    window[slot].waited = true;
    waitSeconds += MPI_Wtime()-start;
}
#else
struct SmallClusterBatch::Impl {};
#endif

void smallClusterWarmUp(const Config* cnf){
    if (cnf->gpu_packed_max_dim <= 0) return;
    // refused here, before the bath is read: the packed kernel is cCCE, and it is in the GPU binary only
    if (strcasecmp(cnf->method,"cce")){
        if (rank==0) fprintf(stderr,"Error: gpu_options.packed_max_dim works with method=CCE (this run: %s)\n",cnf->method);
        ccexFatal();
    }
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (warm.memory.allocated()) return;
    const double start = MPI_Wtime();
    warm.memory.allocate(cnf->nstep,cnf->nstate==0,windowFor(cnf->nstep));
    warm.allocateSeconds = MPI_Wtime()-start;
    if (warm.memory.nodeRank != 0) return;
    // the range by value: SmallClusterBatch takes warm.memory over (and clears it) while this thread may still run
    char* const base = warm.memory.parts[0];
    const size_t bytes = warm.memory.bytes();
    warm.thread = std::thread([base,bytes]{
        const auto t0 = std::chrono::steady_clock::now();
        warm.device = gpuSmallOpen(0,warm.description,sizeof warm.description);
        const auto t1 = std::chrono::steady_clock::now();
        if (warm.device) warm.locked = gpuSmallRegister(base,bytes);
        warm.openSeconds = std::chrono::duration<double>(t1-t0).count();
        warm.lockSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count();
    });
#else
    if (rank==0) fprintf(stderr,"Error: gpu_options.packed_max_dim needs the GPU binary (bin/main_gpu.out, bash do_compile.sh gpu)\n");
    ccexFatal();
#endif
}

void smallClusterWarmDown(){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (warm.thread.joinable()) warm.thread.join();
    if (warm.locked){ gpuSmallUnregister(warm.memory.parts[0]); warm.locked = false; }
    gpuSmallClose(warm.device);
    warm.device = nullptr;
    warm.memory.release();                     // collective: every rank made it in smallClusterWarmUp, or none did
#endif
}

SmallClusterBatch::SmallClusterBatch(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, bool isGCCE){
    if (cnf->gpu_packed_max_dim <= 0) return;
    if (isGCCE || strcasecmp(cnf->method,"cce") || QubitArray_getNqubit(qa)!=1 || strcasecmp(cnf->quantity,"coherence")){
        // given in the input: refused; the default (Config_resolvePackedMaxDim) is simply off here
        if (strcmp(cnf->gpu_packed_source,"given")){
            if (rank==0) printf("[gpu-packed] off: the packed kernel is for method=CCE with one qubit and quantity=coherence\n");
            return;
        }
        if (rank==0) fprintf(stderr,"Error: gpu_options.packed_max_dim supports method=CCE with one qubit and quantity=coherence\n");
        ccexFatal();
    }
    Plan plan; bool ensemble, recursive;
    if (!smallClusterPlan(cnf,pls,plan,ensemble,recursive)){
        if (rank==0) printf("[gpu-packed] more than 4 distinct pulse intervals (or %d pulses): every cluster on the CPU\n",
                            smallcluster::maxIntervals-1);
        return;
    }
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    Impl* m = new Impl;
    m->qa = qa; m->ba = ba; m->dfa = dfa; m->cnf = cnf; m->plan = plan; m->ensemble = ensemble; m->recursive = recursive;
    m->maxDim = cnf->gpu_packed_max_dim;
    const int window = windowFor(plan.nstep);
    // the memory of smallClusterWarmUp when it fits (the same on every rank: it depends on the input only), else a new one
    if (warm.memory.matches(plan.nstep,ensemble,window)){ m->mem = warm.memory; warm.memory = NodeMemory{}; m->fromWarm = true; }
    else m->mem.allocate(plan.nstep,ensemble,window);
    m->own = m->mem.parts[m->mem.nodeRank];
    if (m->mem.nodeRank==0) m->feeder = std::thread([m]{ m->feed(); });
    impl = m;
    if (rank==0) printf("[gpu-packed] on: cCCE clusters of bath dimension %s on the GPU (%s, %s, %d intervals, %d distinct), "
                        "window %d, %.1f MiB shared memory per rank; one feeder per node\n", m->maxDim==4 ? "2 and 4" : "2",
                        ensemble ? "ensemble" : "single state", recursive ? "recursive" : "exact", plan.intervals, plan.distinct,
                        m->mem.windowSize, m->mem.partBytes/1048576.0);
#else
    (void)ba; (void)dfa;
    if (rank==0) fprintf(stderr,"Error: gpu_options.packed_max_dim needs the GPU binary (bin/main_gpu.out, bash do_compile.sh gpu)\n");
    ccexFatal();
#endif
}

SmallClusterBatch::~SmallClusterBatch(){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (impl && impl->feeder.joinable()){ impl->stop.store(true); impl->feeder.join(); }
#endif
    delete impl;
}

void SmallClusterBatch::beginOrder(int n, int** list, int count){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!impl) return;
    for (int s=0; s<2; ++s){
        if (impl->window[s].submitted && !impl->window[s].waited) impl->wait(s);
        impl->window[s] = Impl::Window{};
    }
    impl->order = n; impl->clusters = list; impl->ncluster = count;
#else
    (void)n; (void)list; (void)count;
#endif
}

const std::complex<double>* SmallClusterBatch::get(int ic, int& dim){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!impl) return nullptr;
    Impl& m = *impl;
    int s = -1;
    for (int k=0; k<2; ++k) if (m.window[k].submitted && ic >= m.window[k].first && ic < m.window[k].last) s = k;
    if (s < 0){
        // both windows are behind ic (or none was used in this order): the older slot
        s = !m.window[0].submitted ? 0 : !m.window[1].submitted ? 1 : m.window[0].first < m.window[1].first ? 0 : 1;
        if (m.window[s].submitted && !m.window[s].waited) m.wait(s);
        if (m.ready()) m.prepare(s,ic);
        else m.cpuWindow(s,ic);           // the card is not up yet (the start of the run): these on the CPU meanwhile
    }
    Impl::Window& w = m.window[s];
    if (!w.waited){
        const int o = 1-s;
        if (w.last <= m.ncluster && !(m.window[o].submitted && m.window[o].first == w.last) && m.ready()){
            if (m.window[o].submitted && !m.window[o].waited) m.wait(o);
            m.prepare(o,w.last);                                    // the next window first
        }
        m.wait(s);
    }
    const int where = w.where[ic-w.first];
    if (where < 0) return nullptr;
    dim = w.dim[ic-w.first];
    const long index = dim==2 ? where : w.count[0]+where;
    if (m.status(m.own,s)[index] != smallcluster::statusOk){ ++m.fallback; --m.gpuClusters; return nullptr; }
    return reinterpret_cast<const std::complex<double>*>(m.out(m.own,s)+size_t(index)*m.plan.nstep);
#else
    (void)ic; (void)dim;
    return nullptr;
#endif
}

void SmallClusterBatch::report(){
#if defined(CCEX_USE_CUDA) || defined(CCEX_SMALL_CLUSTER_EMULATE)
    if (!impl) return;
    Impl& m = *impl;
    const double releaseStart = MPI_Wtime();
    MPI_Barrier(m.mem.node);                   // every rank has its last window back
    double feeder[4] = {0,0,0,0}; long long feederCounts[2] = {0,0}; int flags[2] = {0,0};
    if (m.mem.nodeRank==0){
        const double start = MPI_Wtime();
        m.stop.store(true,std::memory_order_release);
        m.feeder.join();
        feeder[0] = m.feederOpen; feeder[1] = m.feederLock; feeder[2] = m.feederWaited; feeder[3] = MPI_Wtime()-start;
        feederCounts[0] = m.feederBatches; feederCounts[1] = m.feederItems;
        flags[0] = m.locked; flags[1] = m.early;
    }
    const double freeStart = MPI_Wtime();
    m.mem.release();
    const double released[3] = {MPI_Wtime()-releaseStart, MPI_Wtime()-freeStart, warm.allocateSeconds};
    long long counts[3] = {m.gpuClusters, m.fallback, m.cpuStarting}, total[3] = {0,0,0}, batches[2] = {0,0};
    double seconds[2] = {m.prepareSeconds, m.waitSeconds}, most[2] = {0,0}, sum[2] = {0,0}, slowest[4] = {0,0,0,0};
    int anyFlags[2] = {0,0};
    MPI_Reduce(counts,total,3,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(feederCounts,batches,2,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(seconds,most,2,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(seconds,sum,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(feeder,slowest,4,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(flags,anyFlags,2,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
    double releasedMost[3] = {0,0,0};
    MPI_Reduce(released,releasedMost,3,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    if (rank==0) printf("[gpu-packed] clusters on the GPU %lld, on the CPU while the GPU started %lld, back to the CPU %lld; "
                        "per rank: prepare %.3f s (max %.3f), "
                        "wait for the GPU %.3f s (max %.3f); feeder: CUDA start %.3f s, page-lock %.3f s (%s)%s, first batch "
                        "after %.3f s, %lld batches of %lld windows, stop %.3f s; shared memory allocated early in %.3f s, release %.3f s (of which freeing it %.3f s)\n", total[0], total[2], total[1], sum[0]/nprocess,
                        most[0], sum[1]/nprocess, most[1], slowest[0], slowest[1], anyFlags[0] ? "done" : "refused, pageable copies",
                        anyFlags[1] ? " early" : "", slowest[2], batches[0], batches[1], slowest[3], releasedMost[2], releasedMost[0], releasedMost[1]);
    delete impl;
    impl = nullptr;
#endif
}
