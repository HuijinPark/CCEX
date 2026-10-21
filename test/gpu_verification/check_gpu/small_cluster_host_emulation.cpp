// A26 step 1: the card API of small_cluster_gpu.h emulated on the host (the same arithmetic of small_cluster_math.h, run
// at submit time), so that the host side of the small-cluster batch - the node shared memory, the feeder thread, the
// window protocol of small_cluster_batch.cpp - can be run end to end on a machine without CUDA:
//   compile small_cluster_batch.cpp with -DCCEX_SMALL_CLUSTER_EMULATE, link this file instead of small_cluster_gpu.o,
//   and run main.out with gpu_options.packed_max_dim on several ranks; the curves must match the CPU run within the GPU
//   tolerance. check_packed_emulation.py builds it that way and runs its regression cases.
#include "calculation/gpu/small_cluster_gpu.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <thread>

using namespace smallcluster;

struct GpuSmallDevice { int unused; };

namespace {
template<int D,bool ENS,bool REC> int one(const Plan& plan, const C* in, C* out){
    switch (slotsFor(plan.distinct)){
        case 1: return cluster<D,ENS,1,REC>(plan,in,out);
        case 2: return cluster<D,ENS,2,REC>(plan,in,out);
        default: return cluster<D,ENS,4,REC>(plan,in,out);
    }
}
template<int D> int dispatch(const Plan& plan, bool ensemble, bool recursive, const C* in, C* out){
    if (ensemble) return recursive ? one<D,true,true>(plan,in,out) : one<D,true,false>(plan,in,out);
    return recursive ? one<D,false,true>(plan,in,out) : one<D,false,false>(plan,in,out);
}
}

GpuSmallDevice* gpuSmallOpen(int, char* description, int capacity){
    // CCEX_SMALL_EMULATE_OPEN_MS: an opening as slow as a real context (0.2-0.5 s on idm4-gpu), so that the early start
    // overlaps the calculation's own start as it does on the card
    if (getenv("CCEX_SMALL_EMULATE_OPEN_MS")) std::this_thread::sleep_for(std::chrono::milliseconds(atoi(getenv("CCEX_SMALL_EMULATE_OPEN_MS"))));
    // CCEX_SMALL_EMULATE_OPEN_FAIL=1: a card that cannot be opened, as gpuSmallOpen answers with CUDA_VISIBLE_DEVICES=-1
    if (getenv("CCEX_SMALL_EMULATE_OPEN_FAIL") && atoi(getenv("CCEX_SMALL_EMULATE_OPEN_FAIL"))){
        snprintf(description,capacity,"host emulation: no card (CCEX_SMALL_EMULATE_OPEN_FAIL)");
        return nullptr;
    }
    snprintf(description,capacity,"host emulation");
    return new GpuSmallDevice{0};
}
void gpuSmallClose(GpuSmallDevice* gpu){ delete gpu; }
bool gpuSmallRegister(void* p, std::size_t bytes){   // touch both ends, as page-locking does
    volatile char first = static_cast<char*>(p)[0], last = static_cast<char*>(p)[bytes-1]; (void)first; (void)last;
    return true;
}
void gpuSmallUnregister(void*){}
void gpuSmallSubmit(GpuSmallDevice*, int, const Plan& plan, bool ensemble, bool recursive, const GpuSmallItem* items, int count){
    const int extra = ensemble ? 0 : 1, per2 = 8+2*extra, per4 = 32+4*extra;
    for (int i=0; i<count; ++i){
        const GpuSmallItem& it = items[i];
        for (long c=0; c<it.count2; ++c)
            it.status[c] = (unsigned char)dispatch<2>(plan,ensemble,recursive,it.in2+c*per2,it.out+c*plan.nstep);
        for (long c=0; c<it.count4; ++c)
            it.status[it.count2+c] = (unsigned char)dispatch<4>(plan,ensemble,recursive,it.in4+c*per4,it.out+(it.count2+c)*plan.nstep);
    }
}
// A40: the warp kernel is not emulated (its arithmetic is warp-parallel CUDA); a run with gpu_options.warp_max_dim stops here
void gpuSmallWarpPlan(GpuSmallDevice*, const warpcluster::Plan&){
    fprintf(stderr,"Error: the host emulation has no warp kernel (gpu_options.warp_max_dim)\n");
    exit(1);
}
bool gpuSmallReady(GpuSmallDevice*, int){ return true; }
void gpuSmallFinish(GpuSmallDevice*, int){}
