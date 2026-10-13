// A40: the warp kernel (warp_cluster_gpu.cu, gpu_options.warp_max_dim) against clusters the CPU computed.
//   test_warp_cluster_cuda.out data/warp_cluster/*.bin
// Each file holds clusters dumped by the CPU path (CCEX_WARP_CLUSTER_DUMP, small_cluster_batch.h): the plan (time plan;
// gCCE pulses and qubit states), the input (Hamiltonian(s) and, for a single state, the start) and the CPU result (the
// coherence, or the qdim x qdim reduced density matrix for gCCE quantity dm). The card computes each file through the
// batch API the run uses (gpuSmallSubmit: one batch of two items, in lane 0 and again in lane 1); every value must agree
// with the CPU within the GPU tolerance 1e-8 + 2e-6 |CPU|, every cluster must converge, and the two lanes must give the same
// bits. A wrong device compilation - the sm_120 NVVM optimizer gave wrong Jacobi eigenvalues in the A26 and A40 step-0
// prototypes - fails here; the per-dimension lines say which size went wrong.
//   test_warp_cluster_cuda.out --bench N dump.bin...   also times, per file and dimension, one batch of N clusters (the
//                                                      file's clusters of that dimension repeated; copies included)
#include "calculation/gpu/small_cluster_gpu.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <map>
#include <vector>

using smallcluster::C;
using warpcluster::Plan;

struct Record { Plan plan; int dim; std::vector<C> in, cpu; };

static bool readRecord(FILE* f, Record& r){
    int head[8];
    if (fread(head,sizeof(int),8,f) != 8) return false;
    r.plan = Plan{};
    Plan& p = r.plan;
    p.mode = head[0]; r.dim = head[1]; p.qdim = head[2]; p.ensemble = head[3]; p.dm = head[4];
    p.time.intervals = head[5]; p.time.distinct = head[6]; p.time.nstep = head[7];
    if ((p.mode != warpcluster::modeCce && p.mode != warpcluster::modeGcce) || r.dim < 2 || r.dim > warpcluster::maxDim ||
        p.qdim < 1 || p.qdim > warpcluster::maxQubitDim || p.time.intervals < 1 || p.time.intervals > warpcluster::maxIntervals ||
        p.time.distinct < 1 || p.time.distinct > 4 || p.time.nstep < 1) return false;
    if (fread(&p.time.deltat,sizeof(double),1,f) != 1) return false;
    if (fread(p.time.fraction,sizeof(double),p.time.distinct,f) != (size_t)p.time.distinct) return false;
    if (fread(p.time.same,sizeof(int),p.time.intervals,f) != (size_t)p.time.intervals) return false;
    if (p.mode == warpcluster::modeGcce){
        const size_t q = size_t(p.qdim);
        for (int k=0; k<p.time.intervals-1; ++k) if (fread(p.pulse[k],sizeof(C),q*q,f) != q*q) return false;
        if (fread(p.psia,sizeof(C),q,f) != q || fread(p.psib,sizeof(C),q,f) != q || fread(p.qpsi0,sizeof(C),q,f) != q) return false;
    }
    r.in.resize(warpcluster::inputSize(p,r.dim));
    r.cpu.resize(size_t(p.time.nstep)*warpcluster::outputPer(p));
    return fread(r.in.data(),sizeof(C),r.in.size(),f) == r.in.size() && fread(r.cpu.data(),sizeof(C),r.cpu.size(),f) == r.cpu.size();
}

static bool samePlan(const Plan& a, const Plan& b){
    if (a.mode != b.mode || a.qdim != b.qdim || a.ensemble != b.ensemble || a.dm != b.dm || a.time.intervals != b.time.intervals ||
        a.time.distinct != b.time.distinct || a.time.nstep != b.time.nstep || a.time.deltat != b.time.deltat ||
        std::memcmp(a.time.fraction,b.time.fraction,sizeof(double)*a.time.distinct) ||
        std::memcmp(a.time.same,b.time.same,sizeof(int)*a.time.intervals)) return false;
    if (a.mode != warpcluster::modeGcce) return true;
    const size_t q = size_t(a.qdim);
    for (int k=0; k<a.time.intervals-1; ++k) if (std::memcmp(a.pulse[k],b.pulse[k],sizeof(C)*q*q)) return false;
    return !std::memcmp(a.psia,b.psia,sizeof(C)*q) && !std::memcmp(a.psib,b.psib,sizeof(C)*q) && !std::memcmp(a.qpsi0,b.qpsi0,sizeof(C)*q);
}

struct Score { double worst = 0; long outside = 0, points = 0, notOk = 0, clusters = 0; };
static void compare(const Record& r, const C* got, int status, Score& s){
    ++s.clusters;
    if (status != smallcluster::statusOk) ++s.notOk;
    for (size_t t=0; t<r.cpu.size(); ++t){
        const double diff = std::hypot(got[t].re-r.cpu[t].re, got[t].im-r.cpu[t].im), ref = std::hypot(r.cpu[t].re, r.cpu[t].im);
        if (!(diff <= 1e-8 + 2e-6*ref)) ++s.outside;     // a NaN counts as outside
        if (!(diff <= s.worst)) s.worst = std::isnan(diff) ? INFINITY : diff;
        ++s.points;
    }
}

int main(int argc, char** argv){
    long bench = 0;
    if (argc > 2 && !std::strcmp(argv[1],"--bench")){ bench = atol(argv[2]); argv += 2; argc -= 2; }
    if (argc < 2){ fprintf(stderr,"usage: %s [--bench N] dump.bin...\n",argv[0]); return 2; }
    char description[256] = {};
    GpuSmallDevice* gpu = gpuSmallOpen(0,description,sizeof description);
    if (!gpu){ printf("FAIL warp-cluster: no GPU (%s)\n",description); return 1; }
    printf("[warp-cluster] %s\n",description);
    bool failed = false;
    for (int a=1; a<argc; ++a){
        FILE* f = fopen(argv[a],"rb");
        if (!f){ printf("FAIL warp-cluster %s: cannot open\n",argv[a]); failed = true; continue; }
        std::vector<Record> recs; Record r;
        while (readRecord(f,r)) recs.push_back(r);
        const bool atEnd = feof(f) || fgetc(f) == EOF;
        fclose(f);
        if (recs.empty() || !atEnd){ printf("FAIL warp-cluster %s: unreadable (%zu records)\n",argv[a],recs.size()); failed = true; continue; }
        const Plan& plan = recs[0].plan;
        bool fine = true;
        for (const Record& x : recs) fine = fine && samePlan(x.plan,plan) && warpcluster::compiled(plan.mode,plan.qdim,x.dim);
        if (!fine){ printf("FAIL warp-cluster %s: records with different plans or uncompiled dimensions\n",argv[a]); failed = true; continue; }
        const size_t per = size_t(plan.time.nstep)*warpcluster::outputPer(plan);
        std::map<int,Score> scores;
        bool sameBits = true;
        try {
            gpuSmallWarpPlan(gpu,plan);
            // the file as one batch of two items (the first half and the rest), in lane 0 and again in lane 1
            std::vector<C> in; std::vector<unsigned char> dims;
            for (const Record& x : recs){ in.insert(in.end(),x.in.begin(),x.in.end()); dims.push_back((unsigned char)x.dim); }
            const long n = long(recs.size()), half = n/2;
            size_t firstIn = 0;
            for (long c=0; c<half; ++c) firstIn += recs[c].in.size();
            std::vector<C> previous;
            for (int lane=0; lane<2; ++lane){
                std::vector<C> out(size_t(n)*per);
                std::vector<unsigned char> status(n, 0xff);
                GpuSmallItem items[2] = {
                    GpuSmallItem{nullptr,nullptr,0,0,out.data(),status.data(),in.data(),dims.data(),half},
                    GpuSmallItem{nullptr,nullptr,0,0,out.data()+size_t(half)*per,status.data()+half,in.data()+firstIn,dims.data()+half,n-half}};
                gpuSmallSubmit(gpu,lane,plan.time,plan.ensemble,false,items,2);
                gpuSmallFinish(gpu,lane);
                if (lane == 0) for (long c=0; c<n; ++c) compare(recs[c],out.data()+size_t(c)*per,status[c],scores[recs[c].dim]);
                if (lane == 1) sameBits = out.size() == previous.size() && !std::memcmp(out.data(),previous.data(),out.size()*sizeof(C));
                previous.swap(out);
            }
        } catch (const std::exception& e){
            printf("FAIL warp-cluster %s: %s\n",argv[a],e.what()); failed = true; continue;
        }
        bool fileOk = sameBits;
        for (const auto& kv : scores){
            const Score& s = kv.second;
            const bool ok = s.outside == 0 && s.notOk == 0;
            fileOk = fileOk && ok;
            printf("%s warp-cluster %s: %s %s%s, %d intervals (%d distinct), dimension %d: %ld clusters, %ld values, max |GPU-CPU| %.3e, "
                   "outside %ld, not converged %ld\n", ok ? "PASS" : "FAIL", argv[a], plan.mode==warpcluster::modeGcce ? "gCCE" : "cCCE",
                   plan.ensemble ? "ensemble" : "single", plan.dm ? " dm" : "", plan.time.intervals, plan.time.distinct, kv.first,
                   s.clusters, s.points, s.worst, s.outside, s.notOk);
        }
        printf("%s warp-cluster %s: two lanes %s\n", sameBits ? "PASS" : "FAIL", argv[a], sameBits ? "same bits" : "DIFFER");
        failed |= !fileOk;
        for (const auto& kv : scores){
            if (bench <= 0) break;
            std::vector<C> in; std::vector<unsigned char> dims;
            std::vector<const Record*> mine;
            for (const Record& x : recs) if (x.dim == kv.first) mine.push_back(&x);
            for (long c=0; c<bench; ++c){ const Record& x = *mine[c % mine.size()]; in.insert(in.end(),x.in.begin(),x.in.end()); dims.push_back((unsigned char)x.dim); }
            std::vector<C> out(size_t(bench)*per); std::vector<unsigned char> status(bench);
            GpuSmallItem item{nullptr,nullptr,0,0,out.data(),status.data(),in.data(),dims.data(),bench};
            double best = 1e30;
            for (int k=0; k<4; ++k){                  // the first one warms up
                const auto t0 = std::chrono::steady_clock::now();
                gpuSmallSubmit(gpu,0,plan.time,plan.ensemble,false,&item,1);
                gpuSmallFinish(gpu,0);
                const double t = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
                if (k) best = std::min(best,t);
            }
            printf("BENCH warp-cluster %s: dimension %d, %ld clusters, %d points: %.4f s, %.0f clusters/s, %.2f us per cluster\n",
                   argv[a], kv.first, bench, plan.time.nstep, best, bench/best, best/bench*1e6);
        }
    }
    gpuSmallClose(gpu);
    return failed ? 1 : 0;
}
