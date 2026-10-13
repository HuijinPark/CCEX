// A26 step 1: the batched small-cluster cCCE kernel (small_cluster_gpu.cu) against clusters the CPU computed.
//   test_small_cluster_cuda.out data/small_cluster/*.bin
// Each file holds clusters dumped by the CPU path (CCEX_SMALL_CLUSTER_DUMP, small_cluster_batch.h): the time plan,
// Halpha, Hbeta, psi and the CPU coherence. Every point of every cluster must agree within the GPU tolerance
// 1e-8 + 2e-6 |CPU| - for the same code compiled for the host (the arithmetic) and on the card (the device compilation:
// A26 step 0 had a device build that got dimension-4 eigenvalues wrong). The card runs each file twice, in the two
// lanes and each time as two items, and the two must be the same bits. Built with -DSMALL_CLUSTER_HOST_ONLY the card part is
// left out (test_small_cluster_host in make numeric-test). Also: the dump reader refuses a cut file (the first one is cut in
// memory), and built-in Hamiltonians of 1e+-200 scale must come out right or be handed back to the CPU.
#include "calculation/gpu/small_cluster_math.h"
#ifndef SMALL_CLUSTER_HOST_ONLY
#include "calculation/gpu/small_cluster_gpu.h"
#endif
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <utility>
#include <string>
#include <vector>

using namespace smallcluster;

// expectFail: a built-in case the kernel must hand back to the CPU (status not ok) instead of computing
struct Record { int d, ens, rec; Plan plan; std::vector<C> data, cpu; bool expectFail = false; };

// One record of a dump. readEnd only at a clean end of the file, before the first byte of a header; a header or body cut
// short, a read error or a header that makes no sense is readBroken (review 20261012-233058-452278 P3-06: they were all
// taken for the end, and a cut file passed with the records before the cut).
enum ReadResult { readOk, readEnd, readBroken };
static bool readBytes(FILE* f, void* p, size_t bytes){ return fread(p,1,bytes,f) == bytes; }
static ReadResult readRecord(FILE* f, Record& r){
    int head[6];
    const size_t got = fread(head,1,sizeof head,f);
    if (got == 0 && feof(f) && !ferror(f)) return readEnd;
    if (got != sizeof head) return readBroken;
    r = Record{};
    r.d = head[0]; r.ens = head[1]; r.rec = head[2];
    r.plan.intervals = head[3]; r.plan.distinct = head[4]; r.plan.nstep = head[5];
    if ((r.d != 2 && r.d != 4) || (r.ens != 0 && r.ens != 1) || (r.rec != 0 && r.rec != 1) || r.plan.intervals < 1 ||
        r.plan.intervals > maxIntervals || r.plan.distinct < 1 || r.plan.distinct > 4 || r.plan.nstep < 1) return readBroken;
    if (!readBytes(f,&r.plan.deltat,sizeof(double)) || !readBytes(f,r.plan.fraction,sizeof(double)*r.plan.distinct) ||
        !readBytes(f,r.plan.same,sizeof(int)*r.plan.intervals)) return readBroken;
    r.data.resize(2*r.d*r.d + (r.ens ? 0 : r.d));
    r.cpu.resize(r.plan.nstep);
    if (!readBytes(f,r.data.data(),sizeof(C)*r.data.size()) || !readBytes(f,r.cpu.data(),sizeof(C)*r.cpu.size())) return readBroken;
    return readOk;
}
// every record of the file; false if it is broken anywhere
static bool readAll(FILE* f, std::vector<Record>& recs){
    recs.clear();
    Record r;
    while (true){
        const ReadResult res = readRecord(f,r);
        if (res == readEnd) return !ferror(f);
        if (res == readBroken) return false;
        recs.push_back(r);
    }
}
static size_t recordBytes(const Record& r){
    return 6*sizeof(int) + sizeof(double)*(1+r.plan.distinct) + sizeof(int)*r.plan.intervals + sizeof(C)*(r.data.size()+r.cpu.size());
}
// The reader on cut copies of a dump (in memory): the whole file and the file cut right after its first record must
// read, cut inside the second header, inside a body or one byte short must not.
static bool readerChecks(const char* path){
    FILE* f = fopen(path,"rb");
    if (!f) return false;
    std::vector<char> bytes;
    char buffer[65536]; size_t n;
    while ((n = fread(buffer,1,sizeof buffer,f)) > 0) bytes.insert(bytes.end(),buffer,buffer+n);
    fclose(f);
    std::vector<Record> recs;
    auto readPrefix = [&](size_t size, std::vector<Record>& out){
        FILE* m = fmemopen(bytes.data(),size,"rb");
        if (!m) return false;
        const bool ok = readAll(m,out);
        fclose(m);
        return ok;
    };
    if (!readPrefix(bytes.size(),recs) || recs.size() < 2) return false;
    const size_t first = recordBytes(recs[0]);
    struct Cut { size_t size; bool reads; const char* what; } cuts[] = {
        {first, true, "after the first record"}, {first+3, false, "3 bytes into the second header"},
        {first+16, false, "16 bytes into the second header"}, {first+24+8, false, "inside the second record's plan"},
        {first+recordBytes(recs[1])/2, false, "inside the second body"}, {bytes.size()-1, false, "one byte short"}};
    bool ok = true;
    for (const Cut& c : cuts){
        std::vector<Record> got;
        const bool reads = readPrefix(c.size,got);
        const bool right = reads == c.reads && (!reads || got.size() == 1);
        if (!right) printf("FAIL small-cluster reader %s cut %s (%zu bytes): %s\n",path,c.what,c.size,reads ? "read" : "refused");
        ok &= right;
    }
    std::vector<char> bad(bytes);
    const int three = 3;
    std::memcpy(bad.data(),&three,sizeof three);           // dimension 3 in the first header
    FILE* m = fmemopen(bad.data(),bad.size(),"rb");
    std::vector<Record> got;
    const bool reads = m && readAll(m,got);
    if (m) fclose(m);
    if (reads){ printf("FAIL small-cluster reader %s: a header of dimension 3 was read\n",path); ok = false; }
    printf("%s small-cluster reader %s: whole file, clean cut and %zu broken cuts\n",ok ? "PASS" : "FAIL",path,
           sizeof cuts/sizeof cuts[0]);
    return ok;
}

// Built-in cases for the Jacobi scaling (review 20261012-233058-452278 P1-01/P3-02): ensemble Ramsey, Hbeta = 0,
// Halpha = scale * K with K of eigenvalues +-w, two time points 0 and deltat = 1/scale: L = 1 and cos(w). Finite Halpha
// whose squared magnitudes overflow or underflow; the old Jacobi returned wrong eigenvalues for them as converged.
// The last case has an eigenvalue beyond the double range: the kernel must hand it back (status not ok).
static std::vector<std::pair<std::string,std::vector<Record>>> scaleCases(){
    std::vector<std::pair<std::string,std::vector<Record>>> sets;
    auto make = [&](const char* name, int d, double scale, const std::vector<std::pair<int,C>>& lower, double w, bool fail){
        Record r;
        r.d = d; r.ens = 1; r.rec = 0;
        r.plan = Plan{};
        r.plan.intervals = 1; r.plan.distinct = 1; r.plan.nstep = 2; r.plan.deltat = fail ? 1.0 : 1.0/scale;
        r.plan.fraction[0] = 1.0; r.plan.same[0] = 0;
        r.data.assign(2*d*d, cmake(0,0));
        for (const auto& e : lower) r.data[e.first] = cmake(e.second.re*scale, e.second.im*scale);   // column-major, i >= j
        r.cpu = {cmake(1,0), cmake(std::cos(w),0)};
        r.expectFail = fail;
        sets.push_back({name,{r}});
    };
    const C one = cmake(1,0), i = cmake(0,1);
    // D = 2: sigma_x (w = 1), [[1,1],[1,-1]] (w = sqrt 2), sigma_y (w = 1)
    for (double scale : {1.0, 1e-200, 1e200}) {
        char name[64]; snprintf(name,sizeof name,"<sigma_x * %g>",scale);
        make(name,2,scale,{{1,one}},1.0,false);
    }
    for (double scale : {1e155, 1e-200}) {
        char name[64]; snprintf(name,sizeof name,"<[[1,1],[1,-1]] * %g>",scale);
        make(name,2,scale,{{0,one},{1,one},{3,cmake(-1,0)}},std::sqrt(2.0),false);
    }
    make("<sigma_y * 1e-170>",2,1e-170,{{1,i}},1.0,false);
    // D = 4: sigma_x (x) 1 and sigma_y (x) 1 (lower triangle: rows 2, 3 of columns 0, 1)
    for (double scale : {1e200, 1e-200}) {
        char name[64]; snprintf(name,sizeof name,"<sigma_x (x) 1 * %g>",scale);
        make(name,4,scale,{{2,one},{1*4+3,one}},1.0,false);
    }
    make("<sigma_y (x) 1 * 1e-170>",4,1e-170,{{2,i},{1*4+3,i}},1.0,false);
    make("<[[1,1],[1,1]] * 1e308, eigenvalue 2e308: back to the CPU>",2,1e308,{{0,one},{1,one},{3,one}},0.0,true);
    return sets;
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
    if (r.expectFail){                       // must be handed back: a computed status is the failure
        if (status == statusOk) ++s.outside;
        s.points += r.plan.nstep;
        return;
    }
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
    bool failed = !readerChecks(argv[1]);
    // the dump files, then the built-in scale cases; each set is run as one batch
    std::vector<std::pair<std::string,std::vector<Record>>> sets;
    for (int a=1; a<argc; ++a){
        FILE* f = fopen(argv[a],"rb");
        if (!f){ printf("FAIL small-cluster %s: cannot open\n",argv[a]); failed = true; continue; }
        std::vector<Record> recs;
        const bool whole = readAll(f,recs);
        fclose(f);
        if (recs.empty() || !whole){ printf("FAIL small-cluster %s: unreadable or cut short (%zu records read)\n",argv[a],recs.size()); failed = true; continue; }
        sets.push_back({argv[a],recs});
    }
    for (auto& set : scaleCases()) sets.push_back(set);
    for (const auto& set : sets){
        const char* name = set.first.c_str();
        std::vector<Record> recs = set.second;
        const Record& first = recs[0];
        for (const Record& x : recs)
            if (x.ens != first.ens || x.rec != first.rec || x.plan.nstep != first.plan.nstep || x.plan.intervals != first.plan.intervals ||
                std::memcmp(x.plan.same,first.plan.same,sizeof(int)*first.plan.intervals) ||
                std::memcmp(x.plan.fraction,first.plan.fraction,sizeof(double)*first.plan.distinct) || x.plan.deltat != first.plan.deltat){
                printf("FAIL small-cluster %s: records with different plans\n",name); failed = true; recs.clear(); break;
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
               "max |host-CPU| %.3e, outside %ld, not converged %ld\n", hostOk ? "PASS" : "FAIL", name,
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
            printf("FAIL small-cluster gpu %s: %s\n",name,e.what()); failed = true; continue;
        }
        const bool gpuOk = gs.outside == 0 && gs.notOk == 0 && sameBits;
        printf("%s small-cluster gpu %s: %ld points, max |GPU-CPU| %.3e, outside %ld, not converged %ld, two lanes %s\n",
               gpuOk ? "PASS" : "FAIL", name, gs.points, gs.worst, gs.outside, gs.notOk, sameBits ? "same bits" : "DIFFER");
        failed |= !gpuOk;
#endif
    }
#ifndef SMALL_CLUSTER_HOST_ONLY
    gpuSmallClose(gpu);
#endif
    return failed ? 1 : 0;
}
