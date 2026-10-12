// A26 step 1: the batched small-cluster cCCE kernel (small_cluster_gpu.cu) against clusters the CPU computed.
//   test_small_cluster_cuda.out data/small_cluster/*.bin
// Each file holds clusters dumped by the CPU path (CCEX_SMALL_CLUSTER_DUMP, small_cluster_batch.h): the time plan,
// Halpha, Hbeta, psi and the CPU coherence. Every point of every cluster must agree within the GPU tolerance
// 1e-8 + 2e-6 |CPU| - for the same code compiled for the host (the arithmetic) and on the card (the device compilation:
// A26 step 0 had a device build that got dimension-4 eigenvalues wrong). The card runs each file twice, in the two
// lanes and each time as two items, and the two must be the same bits. Built with -DSMALL_CLUSTER_HOST_ONLY the card part is left out.
#include "calculation/gpu/small_cluster_math.h"
#ifndef SMALL_CLUSTER_HOST_ONLY
#include "calculation/gpu/small_cluster_gpu.h"
#endif
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace smallcluster;

struct Record { int d, ens, rec; Plan plan; std::vector<C> data, cpu; };

static bool readRecord(FILE* f, Record& r){
    int head[6];
    if (fread(head,sizeof(int),6,f) != 6) return false;
    r.d = head[0]; r.ens = head[1]; r.rec = head[2];
    r.plan = Plan{};
    r.plan.intervals = head[3]; r.plan.distinct = head[4]; r.plan.nstep = head[5];
    if ((r.d != 2 && r.d != 4) || r.plan.intervals < 1 || r.plan.intervals > maxIntervals || r.plan.distinct < 1 ||
        r.plan.distinct > 4 || r.plan.nstep < 1) return false;
    if (fread(&r.plan.deltat,sizeof(double),1,f) != 1) return false;
    if (fread(r.plan.fraction,sizeof(double),r.plan.distinct,f) != (size_t)r.plan.distinct) return false;
    if (fread(r.plan.same,sizeof(int),r.plan.intervals,f) != (size_t)r.plan.intervals) return false;
    r.data.resize(2*r.d*r.d + (r.ens ? 0 : r.d));
    r.cpu.resize(r.plan.nstep);
    return fread(r.data.data(),sizeof(C),r.data.size(),f) == r.data.size() &&
           fread(r.cpu.data(),sizeof(C),r.cpu.size(),f) == r.cpu.size();
}

template<int D,bool ENS,int NF,bool REC> static int hostOne(const Record& r, C* out){ return cluster<D,ENS,NF,REC>(r.plan,r.data.data(),out); }
template<int D,bool ENS,bool REC> static int hostSlots(const Record& r, C* out){
    switch (slotsFor(r.plan.distinct)){
        case 1: return hostOne<D,ENS,1,REC>(r,out);
        case 2: return hostOne<D,ENS,2,REC>(r,out);
        default: return hostOne<D,ENS,4,REC>(r,out);
    }
}
static int host(const Record& r, C* out){
    if (r.d == 2){
        if (r.ens) return r.rec ? hostSlots<2,true,true>(r,out) : hostSlots<2,true,false>(r,out);
        return r.rec ? hostSlots<2,false,true>(r,out) : hostSlots<2,false,false>(r,out);
    }
    if (r.ens) return r.rec ? hostSlots<4,true,true>(r,out) : hostSlots<4,true,false>(r,out);
    return r.rec ? hostSlots<4,false,true>(r,out) : hostSlots<4,false,false>(r,out);
}

struct Score { double worst = 0; long outside = 0, points = 0, notOk = 0; };
static void compare(const Record& r, const C* got, int status, Score& s){
    if (status != statusOk) ++s.notOk;
    for (int t=0; t<r.plan.nstep; ++t){
        const double diff = std::hypot(got[t].re-r.cpu[t].re, got[t].im-r.cpu[t].im), ref = std::hypot(r.cpu[t].re, r.cpu[t].im);
        if (!(diff <= 1e-8 + 2e-6*ref)) ++s.outside;     // a NaN counts as outside
        if (!(diff <= s.worst)) s.worst = std::isnan(diff) ? INFINITY : diff;
        ++s.points;
    }
}

int main(int argc, char** argv){
    if (argc < 2){ fprintf(stderr,"usage: %s dump.bin...\n",argv[0]); return 2; }
#ifndef SMALL_CLUSTER_HOST_ONLY
    char description[256] = {};
    GpuSmallDevice* gpu = gpuSmallOpen(0,description,sizeof description);
    if (!gpu){ printf("FAIL small-cluster: no GPU (%s)\n",description); return 1; }
    printf("[small-cluster] %s\n",description);
#endif
    bool failed = false;
    for (int a=1; a<argc; ++a){
        FILE* f = fopen(argv[a],"rb");
        if (!f){ printf("FAIL small-cluster %s: cannot open\n",argv[a]); failed = true; continue; }
        std::vector<Record> recs; Record r;
        while (readRecord(f,r)) recs.push_back(r);
        const bool atEnd = feof(f) || fgetc(f) == EOF;
        fclose(f);
        if (recs.empty() || !atEnd){ printf("FAIL small-cluster %s: unreadable (%zu records)\n",argv[a],recs.size()); failed = true; continue; }
        const Record& first = recs[0];
        for (const Record& x : recs)
            if (x.ens != first.ens || x.rec != first.rec || x.plan.nstep != first.plan.nstep || x.plan.intervals != first.plan.intervals ||
                std::memcmp(x.plan.same,first.plan.same,sizeof(int)*first.plan.intervals) ||
                std::memcmp(x.plan.fraction,first.plan.fraction,sizeof(double)*first.plan.distinct) || x.plan.deltat != first.plan.deltat){
                printf("FAIL small-cluster %s: records with different plans\n",argv[a]); failed = true; recs.clear(); break;
            }
        if (recs.empty()) continue;
        long dims[2] = {0,0};
        for (const Record& x : recs) ++dims[x.d==4];
        const int nstep = first.plan.nstep;

        Score hs;
        std::vector<C> out(nstep);
        for (const Record& x : recs){ const int st = host(x,out.data()); compare(x,out.data(),st,hs); }
        const bool hostOk = hs.outside == 0 && hs.notOk == 0;
        printf("%s small-cluster host %s: %s %s, %d intervals (%d distinct), %zu clusters (dim 2: %ld, dim 4: %ld), %ld points, "
               "max |host-CPU| %.3e, outside %ld, not converged %ld\n", hostOk ? "PASS" : "FAIL", argv[a],
               first.ens ? "ensemble" : "single", first.rec ? "recursive" : "exact", first.plan.intervals, first.plan.distinct,
               recs.size(), dims[0], dims[1], hs.points, hs.worst, hs.outside, hs.notOk);
        failed |= !hostOk;
#ifndef SMALL_CLUSTER_HOST_ONLY
        Score gs; bool sameBits = true;
        std::vector<C> previous;
        try {
            // the file as one batch of two items (each dimension split between them), in lane 0 and again in lane 1
            const int extra = first.ens ? 0 : 1, per2 = 8+2*extra, per4 = 32+4*extra;
            std::vector<C> in2, in4;
            std::vector<int> where(recs.size());
            long count[2] = {0,0};
            for (size_t i=0; i<recs.size(); ++i){
                const Record& x = recs[i]; const int k = x.d==4;
                std::vector<C>& in = k ? in4 : in2;
                in.insert(in.end(),x.data.begin(),x.data.end());
                where[i] = int(count[k]++);
            }
            const long half2 = count[0]/2, half4 = count[1]/2;
            const long n2[2] = {half2, count[0]-half2}, n4[2] = {half4, count[1]-half4};
            for (int lane=0; lane<2; ++lane){
                std::vector<C> out(size_t(count[0]+count[1])*nstep);
                std::vector<unsigned char> status(count[0]+count[1]);
                // item 0 the first halves, item 1 the rest; an item's output is its dimension-2 clusters, then its dimension-4
                GpuSmallItem items[2];
                for (int it=0; it<2; ++it){
                    const long offset = it ? n2[0]+n4[0] : 0;
                    items[it] = GpuSmallItem{in2.data()+size_t(it ? half2 : 0)*per2, in4.data()+size_t(it ? half4 : 0)*per4,
                                             n2[it], n4[it], out.data()+size_t(offset)*nstep, status.data()+offset};
                }
                gpuSmallSubmit(gpu,lane,first.plan,first.ens,first.rec,items,2);
                gpuSmallFinish(gpu,lane);
                std::vector<C> all;
                for (size_t i=0; i<recs.size(); ++i){
                    const Record& x = recs[i];
                    const int k = x.d==4, w = where[i];
                    const int it = k ? (w < half4 ? 0 : 1) : (w < half2 ? 0 : 1);
                    const long local = k ? n2[it] + (w - (it ? half4 : 0)) : w - (it ? half2 : 0);
                    const long index = (it ? n2[0]+n4[0] : 0) + local;
                    const C* got = out.data()+size_t(index)*nstep;
                    if (lane == 0) compare(x,got,status[index],gs);
                    all.insert(all.end(),got,got+nstep);
                }
                if (lane == 1) sameBits = all.size() == previous.size() && !std::memcmp(all.data(),previous.data(),all.size()*sizeof(C));
                previous.swap(all);
            }
        } catch (const std::exception& e){
            printf("FAIL small-cluster gpu %s: %s\n",argv[a],e.what()); failed = true; continue;
        }
        const bool gpuOk = gs.outside == 0 && gs.notOk == 0 && sameBits;
        printf("%s small-cluster gpu %s: %ld points, max |GPU-CPU| %.3e, outside %ld, not converged %ld, two lanes %s\n",
               gpuOk ? "PASS" : "FAIL", argv[a], gs.points, gs.worst, gs.outside, gs.notOk, sameBits ? "same bits" : "DIFFER");
        failed |= !gpuOk;
#endif
    }
#ifndef SMALL_CLUSTER_HOST_ONLY
    gpuSmallClose(gpu);
#endif
    return failed ? 1 : 0;
}
