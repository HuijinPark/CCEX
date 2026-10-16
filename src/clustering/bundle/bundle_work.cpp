/**
 * @file bundle_work.cpp
 * @brief Bundle clustering: the file protocol with the adaptive refinement driver.
 */
#include "../../../include/clustering/bundle/bundle_work.h"
#include "../../../include/clustering/bundle/bundle_text.h"
#include "../../../include/clustering/bundle/bundle_cluster.h"
#include <limits>
#include "../../../include/input/json.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <set>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <chrono>
#include <thread>
#include <unistd.h>
#include "../../../include/calculation/gpu/compute_backend.h"
#include "../../../include/calculation/gpu/gpu_batch_memory.h"

namespace {
using adaptive_bundle::require;
[[noreturn]] void die(const char* message) {
    fprintf(stderr,"Error: adaptive work (rank %d): %s\n",rank,message);
    MPI_Abort(MPI_COMM_WORLD,1);std::abort();
}
std::vector<std::vector<int>> lists(cJSON* root,const char* key) {
    auto* a=cJSON_GetObjectItem(root,key);
    if(!a)return {};
    require(cJSON_IsArray(a),"bundles/probes must be arrays");
    std::vector<std::vector<int>> out;
    for(int j=0;j<cJSON_GetArraySize(a);++j) {
        auto* row=cJSON_GetArrayItem(a,j);
        require(cJSON_IsArray(row) && cJSON_GetArraySize(row)>0,"empty or invalid site list");
        std::vector<int> ids;
        for(int k=0;k<cJSON_GetArraySize(row);++k) {
            auto* v=cJSON_GetArrayItem(row,k);
            require(cJSON_IsNumber(v) && std::isfinite(v->valuedouble) && v->valuedouble>=0 &&
                    v->valuedouble<INT_MAX && std::floor(v->valuedouble)==v->valuedouble,"invalid site ID");
            ids.push_back(int(v->valuedouble));
        }
        require(std::is_sorted(ids.begin(),ids.end()) &&
                std::adjacent_find(ids.begin(),ids.end())==ids.end(),"site lists must be sorted and unique");
        out.push_back(std::move(ids));
    }
    return out;
}
std::ofstream output(const std::string& path) {
    std::ofstream f(path);require(bool(f),"cannot open response file");
    f<<std::setprecision(17);return f;
}
void writeValue(std::ostream& f,int id,int step,std::complex<double> z) {
    f<<id<<'\t'<<step<<'\t'<<z.real()<<'\t'<<z.imag()<<'\n';
}
MatrixXcd* gcceProbe(QubitArray* qa,BathArray* ba,Config* c,Pulse* pulse,Output* op,
                    const std::string& probe,int state) {
    if(strcasecmp(op->savemode,"cluster_density_matrix"))return calCoherenceGcce(qa,ba,c,pulse,op);
    const std::string name=std::string(op->outfile)+".probe-"+probe+".rank-"+
        std::to_string(rank)+".state-"+std::to_string(state);
    size_t tail=7; // Output_save appends _wiDiv or _noDiv, including terminator.
    const int labels=qa->nqubit+(ba ? ba->nspin : 0);
    for(int i=0;i<labels;++i)tail+=1+std::to_string(i).size();
    require(name.size()+tail<=MAX_FILEPATH,"probe cluster_density_matrix output path exceeds MAX_FILEPATH");
    std::vector<char> path(MAX_FILEPATH,0);
    std::copy(name.begin(),name.end(),path.begin());
    Output unique=*op;unique.outfile=path.data();
    return calCoherenceGcce(qa,ba,c,pulse,&unique);
}
}

namespace { AdaptiveWork loadedWork; bool workLoaded=false; }

// Forget the request this process has read, so the next calculation in the same
// process reads its own (plan stage 8). One calculation per process never
// needed this; a refinement loop inside one process does.
void resetAdaptiveWork() { loadedWork=AdaptiveWork(); workLoaded=false; }

const AdaptiveWork& adaptiveWork(const Config* c) {
    AdaptiveWork& work=loadedWork;
    if(workLoaded)return work;
    workLoaded=true;
    if(!*c->adaptive_work_file)return work;
    try {
        require(!strcasecmp(c->cluster_mode,"adaptive"),"work request requires cluster.mode bundle");
        std::ifstream input(c->adaptive_work_file);require(bool(input),"cannot open request file");
        std::stringstream buffer;buffer<<input.rdbuf();
        std::unique_ptr<cJSON,decltype(&cJSON_Delete)> root(cJSON_Parse(buffer.str().c_str()),cJSON_Delete);
        require(root && cJSON_IsObject(root.get()),"invalid request JSON");
        auto* mode=cJSON_GetObjectItem(root.get(),"mode");
        auto* response=cJSON_GetObjectItem(root.get(),"response");
        require(cJSON_IsString(mode) && cJSON_IsString(response) && *response->valuestring,"mode/response required");
        work.mode=mode->valuestring;work.response=response->valuestring;
        auto* physics=cJSON_GetObjectItem(root.get(),"physics");
        require(!physics || cJSON_IsBool(physics),"physics must be boolean");
        work.physics=physics && cJSON_IsTrue(physics);
        require(work.mode=="plan" || work.mode=="probe" || work.mode=="run" || work.mode=="factors" ||
                work.mode=="except","unknown work mode");
        work.bundles=lists(root.get(),"bundles");work.probes=lists(root.get(),"probes");
        work.extra_clusters=lists(root.get(),"extra_clusters");
        work.except=lists(root.get(),"except");
        require(work.extra_clusters.empty() || (!work.bundles.empty() && work.mode!="probe"),
                "extra_clusters require explicit bundles and plan/run mode");
        require(work.mode=="probe" || work.probes.empty(),"probes require probe mode");
        require(work.mode!="probe" || !work.probes.empty(),"probe mode requires probes");
        require(!strcasecmp(c->quantity,"coherence"),"work requires scalar coherence");
        auto* state=cJSON_GetObjectItem(root.get(),"state");
        auto* step=cJSON_GetObjectItem(root.get(),"step");
        require(work.mode=="factors" || (!state && !step),"state and step require factors mode");
        require(work.mode=="except" || work.except.empty(),"except requires except mode");
        if(work.mode=="factors" || work.mode=="except") {
            // The two modes of method D's regroup ((internal notes) and -stage3b.md, 4절 and
            // 1절): defined for cCCE bundle order 2 only - a gCCE factor carries L0 and the t = 0 normalization, and at
            // another order a pair does not define the culprit - over an explicit partition, with nothing else in the
            // product (the power-0 terms that cluster_contributions and the rdip ladder add would mix with P_K).
            require(!strcasecmp(c->method,"cce"),"factors/except work requires method cce");
            require(c->order==2,"factors/except work requires bundle order 2");
            require(!work.bundles.empty(),"factors/except work requires bundles");
            require(work.extra_clusters.empty(),"factors/except work takes no extra_clusters");
            require(!c->cluster_contributions && c->nrdip_ladder==0,
                    "factors/except work cannot be combined with output.cluster_contributions or rdip_ladder");
        }
        if(work.mode=="factors") {
            const auto integer=[](cJSON* v){
                return v && cJSON_IsNumber(v) && std::isfinite(v->valuedouble) &&
                       std::floor(v->valuedouble)==v->valuedouble && std::fabs(v->valuedouble)<INT_MAX;
            };
            require(integer(state) && integer(step),"factors work requires integer state and step");
            work.state=int(state->valuedouble);work.step=int(step->valuedouble);
            // the states calculate() visits: 0 for the ensemble, s_range first..last for sampled states (A94)
            require(c->nstate==0 ? work.state==0 : (work.state>=std::max(c->state_first,1) && work.state<=c->nstate),
                    "factors state is not a state of this run");
            require(work.step>=0 && work.step<c->nstep,"factors step is outside 0..nstep-1");
        }
        if(work.mode=="except") {
            require(!work.except.empty(),"except mode requires a nonempty except list");
            for(size_t i=0;i<work.except.size();++i)
                require(work.exceptIndex.emplace(work.except[i],int(i)).second,"except lists a term twice");
        }
        work.active=true;
    } catch(const std::exception& e) {die(e.what());}
    return work;
}

void adaptiveWritePlan(const Config* c,const AdaptivePlan& plan,const std::vector<int>& dims,
                       const std::vector<AdaptiveEdge>& physical,int qubitDimension,int bundleCap) {
    const auto& work=adaptiveWork(c);if(!work.active || rank!=0)return;
    try {
        auto f=output(work.response+".plan.json");
        f<<"{\"order\":"<<c->order<<",\"expansion_limits\":{\"max_clusters\":"<<c->expansion_max_clusters
         <<",\"max_work\":"<<c->expansion_max_work<<"},\"dimension_limits\":{\"bundle_max_dim\":"<<bundleCap
         <<"},\"merge_cutoff_rad_ms\":"<<c->bundle_merge_cutoff
         <<",\"qubit_dimension\":"<<qubitDimension<<",\"site_dimensions\":";
        adaptive_bundle::jsonArray(f,dims);
        f<<",\"bundles\":[";
        std::vector<int> labels(dims.size(),-1);
        for(size_t g=0;g<plan.bundles.size();++g) {
            f<<(g ? "," : "")<<'[';
            for(size_t i=0;i<plan.bundles[g].size();++i) {
                const int s=plan.bundles[g][i];labels[s]=int(g);f<<(i ? "," : "")<<s;
            }
            f<<']';
        }
        std::map<std::pair<int,int>,long double> strengths;
        for(const auto& e:physical) {
            int a=labels[e.a],b=labels[e.b];if(a==b)continue;if(a>b)std::swap(a,b);
            strengths[{a,b}]+=e.strength;
        }
        f<<"],\"edges\":[";
        for(size_t i=0;i<plan.edges.size();++i) {
            auto e=plan.edges[i];f<<(i ? "," : "")<<'['<<e.first<<','<<e.second<<','<<strengths[e]<<']';
        }
        f<<"],\"candidate_edges\":[";
        size_t candidate=0;
        for(const auto& e:strengths)if(e.second>0) {
            f<<(candidate++ ? "," : "")<<'['<<e.first.first<<','<<e.first.second<<','<<e.second<<']';
        }
        f<<"],\"extra_clusters\":[";
        for(size_t i=0;i<work.extra_clusters.size();++i) {
            f<<(i ? "," : "");adaptive_bundle::jsonArray(f,work.extra_clusters[i]);
        }
        f<<"],\"omitted_score\":"<<plan.omitted_strength*double(c->deltat)*(c->nstep-1)<<"}\n";
        f.flush();require(bool(f),"writing plan failed");
    } catch(const std::exception& e) {die(e.what());}
}

void adaptiveWritePhysics(const Config* c,BathArray* ba,QubitArray* qa,DefectArray* dfa) {
    const auto& w=adaptiveWork(c);if(!w.active || !w.physics || rank!=0)return;
    try {
        auto f=output(w.response+".physics.json");
        auto matrix=[&](const MatrixXcd& m) {
            f<<'[';
            for(int i=0;i<m.rows();++i) {
                f<<(i ? "," : "")<<'[';
                for(int j=0;j<m.cols();++j)f<<(j ? "," : "")<<'['<<m(i,j).real()<<','<<m(i,j).imag()<<']';
                f<<']';
            }
            f<<']';
        };
        auto name=[&](const char* value) {adaptive_bundle::jsonString(f,value);};
        auto component=[&](BathSpin* s) {
            f<<"{\"species\":";name(s->name);
            f<<",\"spin\":"<<s->spin<<",\"gyro\":"<<s->gyro<<",\"xyz\":["
             <<s->xyz[0]<<','<<s->xyz[1]<<','<<s->xyz[2]<<"],\"detuning\":"<<s->detuning
             <<",\"hyperfine\":[";
            for(int j=0;j<qa->nqubit;++j){if(j)f<<',';matrix(s->hypf[j]);}
            f<<"],\"quadrupole\":";matrix(s->quad);f<<'}';
        };
        f<<"{\"units\":\"rad/ms, gauss, angstrom\",\"h_bar_conversion\":"<<H_BAR
         <<",\"bfield\":["<<c->bfield[0]<<','<<c->bfield[1]<<','<<c->bfield[2]<<"],\"qubits\":[";
        for(int i=0;i<qa->nqubit;++i) {
            if(i)f<<',';
            auto* q=qa->qubit[i];
            f<<"{\"spin\":"<<q->spin<<",\"alpha\":";matrix(q->alpha);
            f<<",\"beta\":";matrix(q->beta);f<<'}';
        }
        f<<"],\"sites\":[";
        for(int i=0;i<ba->nspin;++i) {
            auto* s=ba->bath[i];if(i)f<<',';
            f<<"{\"id\":"<<i<<",\"species\":";name(s->name);
            f<<",\"spin\":"<<s->spin<<",\"gyro\":"<<s->gyro<<",\"xyz\":["
             <<s->xyz[0]<<','<<s->xyz[1]<<','<<s->xyz[2]<<"],\"detuning\":"<<s->detuning
             <<",\"attached_spins\":"<<(dfa->naddspins ? dfa->naddspins[i] : 0)<<",\"hyperfine\":[";
            for(int j=0;j<qa->nqubit;++j){if(j)f<<',';matrix(s->hypf[j]);}
            f<<"],\"quadrupole\":";matrix(s->quad);
            // Optional physics metadata only: preserve an indivisible defect's
            // components and the SAME internal tensors used by HamilBath.
            // No propagation, random sampling or default bundleing is changed.
            int site=i;
            std::unique_ptr<BathArray,decltype(&BathArray_freeAll)> local(
                createBathArray(&site,1,ba,dfa,qa->nqubit),BathArray_freeAll);
            f<<",\"components\":[";
            for(int j=0;j<local->nspin;++j){if(j)f<<',';component(local->bath[j]);}
            f<<"],\"internal_interactions\":[";
            bool first=true;
            for(int a=0;a<local->nspin;++a)for(int b=a+1;b<local->nspin;++b) {
                if(!first)f<<',';
                first=false;
                f<<"{\"a\":"<<a<<",\"b\":"<<b<<",\"tensor\":";
                matrix(BathArray_int_i_j(local.get(),a,b));f<<'}';
            }
            f<<"]}";
        }
        f<<"]}\n";f.flush();require(bool(f),"writing physics metadata failed");
    } catch(const std::exception& e) {die(e.what());}
}

namespace {
// A1 (2026-10-02): the probes of a state go to the ranks in turn ((probe + state) % ranks) and the ranks no longer
// wait for each other between states, so the large probes of several states can run at once. What keeps them from
// filling the node is this budget: one counter per node in an MPI shared window holds the estimated memory of the
// probes running on it; a probe adds its estimate before it starts, waits while the sum would pass the budget, and
// takes it off when it is done. A probe larger than the whole budget runs only when nothing else on the node does,
// so nothing waits forever. The estimate is 6 dense complex dim x dim matrices: the GPU peak of a dim-16384 probe
// was 25.3 GB, 5.9 of them (job 623). The budget is probe_memory_gb, or 80 % of the node's RAM and, when this run
// may use the GPU, of its device memory.
// The memory this process may use on its node: the physical RAM, or less if a cgroup limits it - a Slurm job is given
// its share of the node that way (review 4dd5fd R3). cgroup v2 (memory.max) and v1 (memory.limit_in_bytes), the
// smallest limit from this process's cgroup up to the root; "max" or an absurd value means none.
double nodeMemoryLimit() {
    double limit=double(sysconf(_SC_PHYS_PAGES))*double(sysconf(_SC_PAGESIZE));
    std::ifstream groups("/proc/self/cgroup");
    std::string line;
    while (std::getline(groups,line)) {
        const size_t a=line.find(':'), b=line.find(':',a+1);
        if (a==std::string::npos || b==std::string::npos) continue;
        const std::string controllers=line.substr(a+1,b-a-1), path=line.substr(b+1);
        std::string root, file;
        if (controllers.empty()) {root="/sys/fs/cgroup"; file="memory.max";}                     // v2
        else if (("," + controllers + ",").find(",memory,")!=std::string::npos) {root="/sys/fs/cgroup/memory"; file="memory.limit_in_bytes";}
        else continue;
        for (std::string at=path;;) {
            std::ifstream f(root+at+(at=="/" ? "" : "/")+file);
            double value=0;
            if (f >> value && value>0) limit=std::min(limit,value);
            if (at.empty() || at=="/") break;
            const size_t slash=at.find_last_of('/');
            at=slash==0 || slash==std::string::npos ? "/" : at.substr(0,slash);
        }
    }
    return limit;
}

class ProbeBudget {
    MPI_Comm node_=MPI_COMM_NULL; MPI_Win win_=MPI_WIN_NULL; long long* shared_=nullptr;   // [0] held, [1] budget
    double waited_=0;
public:
    bool active() const {return shared_!=nullptr;}
    void start(const Config* c) {          // collective over MPI_COMM_WORLD
        waited_=0;                          // per probe step: an in-process run starts and finishes several (4dd5fd R4)
        MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node_);
        int nodeRank=0;MPI_Comm_rank(node_,&nodeRank);
        long long* base=nullptr;
        MPI_Win_allocate_shared(nodeRank==0 ? 2*sizeof(long long) : 0,sizeof(long long),MPI_INFO_NULL,node_,&base,&win_);
        MPI_Aint size=0;int unit=0;
        MPI_Win_shared_query(win_,0,&size,&unit,&shared_);
        if (nodeRank==0) {
            const double ram=nodeMemoryLimit();
            // the GPU is asked only when no budget is given (review 5f7834 R4)
            const double device=c->probe_memory_gb>0 ? 0 : double(computeDeviceMemory(c));
            const double automatic=0.8*(device>0 ? std::min(ram,device) : ram);
            const double bytes=c->probe_memory_gb>0 ? c->probe_memory_gb*1e9 : automatic;
            __atomic_store_n(&shared_[0],0LL,__ATOMIC_SEQ_CST);
            __atomic_store_n(&shared_[1],(long long)std::min(bytes,9e18),__ATOMIC_SEQ_CST);
            if (rank==0) printf("[probe-budget] %.2f GB per node for the probes running at once%s\n",bytes/1e9,
                                c->probe_memory_gb>0 ? " (cluster.adaptive.probe_memory_gb)" : device>0 ? " (80 % of the memory this job may use and of the GPU)" : " (80 % of the memory this job may use)");
        }
        MPI_Barrier(node_);
    }
    // What one probe holds, on the path it will take (review 5f7834 R1, 4dd5fd R2, 0ef89c R1). The CPU: 6 dense
    // dim x dim complex matrices (the measured setup peak); for gCCE one more per pulse (prepareGccePulses expands them
    // to the whole space) and, with evolution=matrix, the npulse+1 free propagators calPropagatorGcce keeps (Ufrees).
    // The GPU: what the device allocates for the run's time batch, pulses and vector/matrix evolution. savemode cluster_density_matrix
    // keeps every time's full gCCE density matrix on the host (result_tot), and on the GPU also a host read-back buffer
    // of one time batch and, with evolution=vector, a device one. One counter serves both, so the larger.
    static long long estimate(const Config* c, const char* savemode, long long dim, int qubitDim, bool gcce, int pulses,
                              bool manualPulse, bool onGpu) {
        const long double square=16.0L*(long double)dim*(long double)dim;
        const int batch=std::max(1,std::min(c->gpu_batch_steps,c->nstep));
        const bool matrix=gcce && strcasecmp(Config_getEvolution(const_cast<Config*>(c)),"vector");
        const bool allfull=gcce && savemode && !strcasecmp(savemode,"cluster_density_matrix");
        const long double npulse=std::max(0,pulses);
        long double host=6.0L*square;
        if (gcce && !onGpu) host+=(npulse+(matrix ? npulse+1 : 0))*square;
        // CPU cCCE (calCoherenceCce, review 89e70e R1): the two Hamiltonians, the two eigensolvers' vectors, M_up,
        // M_down and M; a sampled state adds M^dagger, the ensemble CPMG block up to 16 more matrices, the
        // ensemble Manual block three and four caches per interval
        if (!gcce && !onGpu) {
            const bool manual=c->nstate==0 && manualPulse;
            host=(7.0L+(c->nstate!=0 ? 1.0L : manual ? 3.0L+4.0L*(npulse+1) : 16.0L))*square;
        }
        if (allfull) host+=(long double)(c->nstep+(onGpu ? batch : 0))*square;
        long double device=0;
        if (onGpu && dim<=std::numeric_limits<int>::max()) {
            const bool vector=gcce ? !matrix : c->nstate!=0;
            device=gpuClusterMemoryEstimate(int(dim),gcce ? qubitDim : 1,batch,int(npulse),gcce,vector);
            if (allfull && vector) device+=(long double)batch*square;
        }
        const long double bytes=std::max(host,device);
        return bytes>9e18L ? (long long)9e18 : (long long)bytes;
    }
    void acquire(long long bytes) {
        const auto begin=std::chrono::steady_clock::now();
        for (;;) {
            long long held=__atomic_load_n(&shared_[0],__ATOMIC_SEQ_CST);
            const long long limit=__atomic_load_n(&shared_[1],__ATOMIC_SEQ_CST);
            if ((held==0 || held+bytes<=limit) &&
                __atomic_compare_exchange_n(&shared_[0],&held,held+bytes,false,__ATOMIC_SEQ_CST,__ATOMIC_SEQ_CST))
                break;
            std::this_thread::sleep_for(std::chrono::milliseconds(20));
        }
        waited_+=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    }
    void release(long long bytes) {__atomic_fetch_sub(&shared_[0],bytes,__ATOMIC_SEQ_CST);}
    void finish() {                         // collective over MPI_COMM_WORLD
        if (!active()) return;
        if (waited_>0.05) printf("[probe-budget] rank=%d waited %.1f s for memory\n",rank,waited_);
        MPI_Win_free(&win_);MPI_Comm_free(&node_);shared_=nullptr;
    }
};
ProbeBudget probeBudget;
constexpr long long kProbeBudgetFloor=64LL*1024*1024;   // smaller probes do not take part

// One probe's curve: its sites as a bath of their own, calculated under the node budget. The probe and factors modes
// share it, so a factor is the probe's value bit for bit (adaptive-regroup-stage3.md 3절).
std::unique_ptr<MatrixXcd[]> probeValues(QubitArray* qa,BathArray* ba,DefectArray* dfa,Config* c,Pulse* pulse,Output* op,
                                         const std::vector<int>& sites,const std::string& label,int state) {
    const bool gcce=!strcasecmp(c->method,"gcce");
    ProbeSubset probe=adaptiveProbeSubset(sites,ba,dfa,qa,gcce);
    BathArray* const subset=probe.bath.get();
    const long long dim=probe.dim;
    // cCCE solves the bath alone, gCCE with the qubit; either way dense dim x dim matrices
    const bool onGpu=dim<=std::numeric_limits<int>::max() && computeUseCuda(c,int(dim));
    const long long need=ProbeBudget::estimate(c,Output_getSavemode(op),dim,gcce ? QubitArray_dim(qa) : 1,gcce,
                                               pulse->npulse,!strcasecmp(pulse->pulsename,"Manual"),onGpu);
    const bool budgeted=probeBudget.active() && need>=kProbeBudgetFloor;
    if(budgeted)probeBudget.acquire(need);
    std::unique_ptr<MatrixXcd[]> values;
    try {
        values.reset(gcce ? gcceProbe(qa,subset,c,pulse,op,label,state) : calCoherenceCce(qa,subset,c,pulse));
    } catch(...) {if(budgeted)probeBudget.release(need);throw;}
    if(budgeted)probeBudget.release(need);
    return values;
}
}  // namespace

void adaptiveProbeFinish() {probeBudget.finish();}

void adaptiveProbe(QubitArray* qa,BathArray* ba,DefectArray* dfa,Config* c,Pulse* pulse,Output* op,int state) {
    try {
        const auto& w=adaptiveWork(c);
        if(!probeBudget.active())probeBudget.start(c);   // collective: every rank reaches its first probe state
        auto f=output(w.response+".rank-"+std::to_string(rank)+".state-"+std::to_string(state)+".tsv");
        f<<"probe\tstep\treal\timag\n";
        const bool gcce=!strcasecmp(c->method,"gcce");
        if(rank==0) {
            std::unique_ptr<MatrixXcd[]> zero(gcce ? gcceProbe(qa,nullptr,c,pulse,op,"zero",state) : nullptr);
            for(int t=0;t<c->nstep;++t)writeValue(f,-1,t,gcce ? zero[t](0,0) : std::complex<double>(1,0));
        }
        // probe i goes to rank (i + state) % ranks: a large probe moves to another rank at each state (A1)
        const size_t first=size_t(((rank-state)%nprocess+nprocess)%nprocess);
        for(size_t i=first;i<w.probes.size();i+=nprocess) {
            std::unique_ptr<MatrixXcd[]> values=probeValues(qa,ba,dfa,c,pulse,op,w.probes[i],std::to_string(i),state);
            for(int t=0;t<c->nstep;++t)writeValue(f,int(i),t,values[t](0,0));
        }
        f.flush();require(bool(f),"writing probes failed");
    } catch(const std::exception& e) {die(e.what());}
}

void adaptiveFactors(QubitArray* qa,BathArray* ba,DefectArray* dfa,Config* c,Pulse* pulse,Output* op,int state) {
    try {
        const auto& w=adaptiveWork(c);
        if(!probeBudget.active())probeBudget.start(c);   // collective: every rank reaches the requested state
        const auto& bundles=adaptiveInstalledBundles();
        const auto& edges=adaptiveInstalledEdges();
        // The clusters: every bundle with a connection (the culprit factor reads L_a and L_b whatever their power
        // in the run, so the singletons the run leaves out at power 0 are here too), then every connected pair -
        // numbered by plan.bundles, the pair as (a, b), a bundle as (a, -1).
        std::vector<char> connected(bundles.size(),0);
        for(const auto& e:edges)connected[size_t(e.first)]=connected[size_t(e.second)]=1;
        std::vector<std::pair<int,int>> rows;
        for(size_t g=0;g<bundles.size();++g)if(connected[g])rows.emplace_back(int(g),-1);
        rows.insert(rows.end(),edges.begin(),edges.end());
        auto f=output(w.response+".rank-"+std::to_string(rank)+".factors.tsv");
        f<<"# state="<<state<<" step="<<w.step<<" bundles="<<bundles.size()<<" edges="<<edges.size()<<"\n";
        f<<"a\tb\treal\timag\n";
        // Each value is one cluster calculated on its own, so the ranks need nothing from each other: row i on rank
        // i % ranks.
        for(size_t i=size_t(rank);i<rows.size();i+=size_t(nprocess)) {
            std::vector<int> sites=bundles[size_t(rows[i].first)];
            if(rows[i].second>=0) {
                const auto& other=bundles[size_t(rows[i].second)];
                sites.insert(sites.end(),other.begin(),other.end());
                std::sort(sites.begin(),sites.end());
            }
            std::unique_ptr<MatrixXcd[]> values=probeValues(qa,ba,dfa,c,pulse,op,sites,"factor-"+std::to_string(i),state);
            const std::complex<double> z=values[w.step](0,0);
            f<<rows[i].first<<'\t'<<rows[i].second<<'\t'<<z.real()<<'\t'<<z.imag()<<'\n';
        }
        f.flush();require(bool(f),"writing factors failed");
    } catch(const std::exception& e) {die(e.what());}
}

ProbeSubset adaptiveProbeSubset(const std::vector<int>& sites,BathArray* ba,DefectArray* dfa,QubitArray* qa,bool gcce) {
    for(int s:sites)require(s>=0 && s<ba->nspin,"probe site outside loaded bath");
    ProbeSubset out;
    out.bath.reset(createBathArray(const_cast<int*>(sites.data()),int(sites.size()),ba,dfa,qa->nqubit));
    out.dim=gcce ? QubitArray_dim(qa) : 1;
    for(int s=0;s<out.bath->nspin;++s) {
        const int d=int(2*out.bath->bath[s]->spin+1);
        require(d>0 && out.dim<=std::numeric_limits<long long>::max()/d,"probe dimension leaves the long range");
        out.dim*=d;
    }
    return out;
}

void adaptiveWriteCurve(const Config* c,const MatrixXcd* values,int state) {
    const auto& w=adaptiveWork(c);if(!w.active || rank!=0)return;
    try {
        auto f=output(w.response+".curve.state-"+std::to_string(state)+".tsv");
        f<<"probe\tstep\treal\timag\n";
        for(int t=0;t<c->nstep;++t)writeValue(f,0,t,values[t](0,0));
        f.flush();require(bool(f),"writing curve failed");
    } catch(const std::exception& e) {die(e.what());}
}

// ---------------------------------------------------------------------------------------------------------------------
// mode "except" ((internal notes) 1절): the run as it is, with the product of the terms not
// in the except list kept beside it (P_K) and the excepted terms' raw curves written out as probes.
// ---------------------------------------------------------------------------------------------------------------------
int adaptiveExceptIndex(const int* sites,int n) {
    const auto& index=loadedWork.exceptIndex;
    const auto found=index.find(std::vector<int>(sites,sites+n));
    return found==index.end() ? -1 : found->second;
}

std::unique_ptr<std::ofstream> adaptiveExceptCurves(const Config* c,int state) {
    std::unique_ptr<std::ofstream> out;
    try {
        const auto& w=adaptiveWork(c);
        out.reset(new std::ofstream(output(w.response+".rank-"+std::to_string(rank)+".state-"+std::to_string(state)+".tsv")));
        *out<<"probe\tstep\treal\timag\n";
        if(rank==0)for(int t=0;t<c->nstep;++t)writeValue(*out,-1,t,std::complex<double>(1,0));   // cCCE's L0
    } catch(const std::exception& e) {die(e.what());}
    return out;
}

void adaptiveExceptValues(std::ostream& f,int index,const MatrixXcd* values,int nstep) {
    for(int t=0;t<nstep;++t)writeValue(f,index,t,values[t](0,0));
}

void adaptiveExceptCheck(long long matched) {
    long long total=0;
    MPI_Allreduce(&matched,&total,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    if(total!=(long long)loadedWork.except.size())
        die("the except terms met in the run are not the except list once each");
}

void adaptiveWriteExcept(const Config* c,const std::vector<ScaledComplex>& kept,int state) {
    const auto& w=adaptiveWork(c);if(!w.active || rank!=0)return;
    try {
        std::ofstream f(w.response+".except.state-"+std::to_string(state)+".tsv");
        require(bool(f),"cannot open response file");
        f<<"step\treal\timag\texponent\n";
        char buffer[160];
        for(int t=0;t<c->nstep;++t) {
            // %La: the long double mantissa exactly, so the driver reads back the product the run holds
            snprintf(buffer,sizeof buffer,"%d\t%La\t%La\t%lld\n",t,kept[size_t(t)].real,kept[size_t(t)].imag,
                     kept[size_t(t)].exponent);
            f<<buffer;
        }
        f.flush();require(bool(f),"writing except product failed");
    } catch(const std::exception& e) {die(e.what());}
}
