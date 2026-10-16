/**
 * @file small_cluster_gpu.cu
 * @brief The packed GPU kernel: the cCCE of small clusters, one thread per cluster.
 */
// A26 step 1: the batched cCCE of small clusters on the card - one thread per cluster, the whole calculation of a
// cluster (both eigensystems, the time chain of every step) in that thread (small_cluster_math.h). A batch is a list of
// items (windows of packed clusters, one per rank); its inputs are copied into one device buffer per dimension, one kernel
// per dimension computes them all, and the coherences are copied back item by item. Two lanes, so that one batch is
// copied while the other computes.
#include "calculation/gpu/small_cluster_gpu.h"
#include "calculation/gpu/warp_cluster_gpu.h"
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>

using smallcluster::C;
using smallcluster::Plan;

namespace {
void check(cudaError_t s, const char* op){
    if (s != cudaSuccess) throw std::runtime_error(std::string("small-cluster GPU: ")+op+": "+cudaGetErrorString(s));
}
#define SC_CHECK(op) check((op),#op)

template<class T> struct Device {
    T* p=nullptr; size_t n=0;
    void reserve(size_t want){
        if (want<=n) return;
        if (p) SC_CHECK(cudaFree(p));
        p=nullptr; n=0;
        SC_CHECK(cudaMalloc((void**)&p,want*sizeof(T)));
        n=want;
    }
    ~Device(){ if (p) cudaFree(p); }
};
struct Lane {
    Device<C> in2, in4, out; Device<unsigned char> status;
    Device<C> inW; Device<WarpEntry> entries;      // A40: the warp clusters' inputs and their launch lists
    WarpEntry* hostEntries=nullptr; size_t hostEntriesSize=0;   // page-locked, read by the copy until the lane is done
    cudaStream_t stream=nullptr;
    cudaEvent_t done=nullptr;
    bool busy=false;
};

template<int D,bool ENS,int NF,bool REC>
__global__ void smallCce(const Plan plan, const C* in, C* out, unsigned char* status, long n){
    __shared__ Plan shared;
    if (threadIdx.x==0) shared=plan;
    __syncthreads();
    const long c=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if (c>=n) return;
    constexpr int per=2*D*D+(ENS ? 0 : D);
    status[c]=(unsigned char)smallcluster::cluster<D,ENS,NF,REC>(shared,in+c*per,out+c*(long)shared.nstep);
}
template<int D,bool ENS,int NF,bool REC>
void launchOne(const Plan& plan, const C* in, C* out, unsigned char* status, long n, cudaStream_t stream){
    const int threads=128;
    smallCce<D,ENS,NF,REC><<<(unsigned)((n+threads-1)/threads),threads,0,stream>>>(plan,in,out,status,n);
    SC_CHECK(cudaGetLastError());
}
template<int D,bool ENS,bool REC>
void launchSlots(const Plan& plan, const C* in, C* out, unsigned char* status, long n, cudaStream_t stream){
    switch (smallcluster::slotsFor(plan.distinct)){
        case 1: launchOne<D,ENS,1,REC>(plan,in,out,status,n,stream); break;
        case 2: launchOne<D,ENS,2,REC>(plan,in,out,status,n,stream); break;
        case 4: launchOne<D,ENS,4,REC>(plan,in,out,status,n,stream); break;
        default: throw std::runtime_error("small-cluster GPU: more than 4 distinct intervals");
    }
}
template<int D>
void launchDim(const Plan& plan, bool ensemble, bool recursive, const C* in, C* out, unsigned char* status, long n, cudaStream_t stream){
    if (!n) return;
    if (ensemble){ if (recursive) launchSlots<D,true,true>(plan,in,out,status,n,stream);  else launchSlots<D,true,false>(plan,in,out,status,n,stream); }
    else         { if (recursive) launchSlots<D,false,true>(plan,in,out,status,n,stream); else launchSlots<D,false,false>(plan,in,out,status,n,stream); }
}
}

struct GpuSmallDevice { Lane lane[2]; warpcluster::Plan warp; bool warpSet=false; };

GpuSmallDevice* gpuSmallOpen(int device, char* description, int capacity){
    // the guard of gpuInitialize (gpu_backend.cu): an empty CUDA_VISIBLE_DEVICES aborts inside cudaGetDeviceCount
    const char* visible=std::getenv("CUDA_VISIBLE_DEVICES");
    if (visible && (!*visible || !std::strcmp(visible,"-1"))){
        snprintf(description,capacity,"CUDA_VISIBLE_DEVICES hides all GPUs");
        return nullptr;
    }
    GpuSmallDevice* gpu=nullptr;
    try {
        int count=0;
        SC_CHECK(cudaGetDeviceCount(&count));
        if (!count) throw std::runtime_error("no visible CUDA devices");
        device%=count;
        SC_CHECK(cudaSetDevice(device));
        cudaDeviceProp properties;
        SC_CHECK(cudaGetDeviceProperties(&properties,device));
        gpu=new GpuSmallDevice;
        for (Lane& l : gpu->lane){
            SC_CHECK(cudaStreamCreateWithFlags(&l.stream,cudaStreamNonBlocking));
            SC_CHECK(cudaEventCreateWithFlags(&l.done,cudaEventDisableTiming));
        }
        snprintf(description,capacity,"device=%d name=%s visible_devices=%d",device,properties.name,count);
        return gpu;
    } catch (const std::exception& e){
        if (gpu) gpuSmallClose(gpu);
        snprintf(description,capacity,"%s",e.what());
        return nullptr;
    }
}

void gpuSmallClose(GpuSmallDevice* gpu){
    if (!gpu) return;
    for (Lane& l : gpu->lane){
        if (l.stream) cudaStreamSynchronize(l.stream);
        if (l.done) cudaEventDestroy(l.done);
        if (l.stream) cudaStreamDestroy(l.stream);
        if (l.hostEntries) cudaFreeHost(l.hostEntries);
    }
    delete gpu;
}

bool gpuSmallRegister(void* pointer, std::size_t bytes){
    if (cudaHostRegister(pointer,bytes,cudaHostRegisterDefault)==cudaSuccess) return true;
    cudaGetLastError();                 // clear it: the copies then go through pageable memory
    return false;
}
void gpuSmallUnregister(void* pointer){ if (cudaHostUnregister(pointer)!=cudaSuccess) cudaGetLastError(); }

void gpuSmallWarpPlan(GpuSmallDevice* gpu, const warpcluster::Plan& plan){
    warpClusterUpload(plan);
    gpu->warp=plan; gpu->warpSet=true;
}

void gpuSmallSubmit(GpuSmallDevice* gpu, int lane, const Plan& plan, bool ensemble, bool recursive,
                    const GpuSmallItem* items, int count){
    Lane& l=gpu->lane[lane];
    if (l.busy) throw std::runtime_error("small-cluster GPU: submit to a busy lane");
    const int extra=ensemble ? 0 : 1, per2=8+2*extra, per4=32+4*extra;
    const size_t nstep=size_t(plan.nstep);
    long total2=0,total4=0,totalW=0;
    for (int i=0;i<count;++i){ total2+=items[i].count2; total4+=items[i].count4; totalW+=items[i].countW; }
    if (totalW && !gpu->warpSet) throw std::runtime_error("small-cluster GPU: warp clusters before the warp plan");
    const size_t perW=totalW ? nstep*warpcluster::outputPer(gpu->warp) : 0;
    l.in2.reserve(size_t(total2)*per2); l.in4.reserve(size_t(total4)*per4);
    l.out.reserve(size_t(total2+total4)*nstep+size_t(totalW)*perW); l.status.reserve(size_t(total2+total4+totalW));
    C* out4=l.out.p+size_t(total2)*nstep;
    unsigned char* status4=l.status.p+total2;
    C* outW=l.out.p+size_t(total2+total4)*nstep;
    unsigned char* statusW=l.status.p+total2+total4;
    // A40: the warp clusters' inputs one item after the other, and a launch list per total dimension
    long countDim[warpcluster::maxDim+1]={}, firstDim[warpcluster::maxDim+2]={};
    size_t totalWIn=0;
    for (int i=0;i<count;++i)
        for (long c=0;c<items[i].countW;++c){ const int d=items[i].dimW[c]; ++countDim[d]; totalWIn+=warpcluster::inputSize(gpu->warp,d); }
    if (totalW){
        for (int d=0;d<=warpcluster::maxDim;++d) firstDim[d+1]=firstDim[d]+countDim[d];
        if (l.hostEntriesSize<size_t(totalW)){
            if (l.hostEntries) SC_CHECK(cudaFreeHost(l.hostEntries));
            l.hostEntries=nullptr; l.hostEntriesSize=0;
            SC_CHECK(cudaMallocHost((void**)&l.hostEntries,sizeof(WarpEntry)*totalW));
            l.hostEntriesSize=size_t(totalW);
        }
        long at[warpcluster::maxDim+1];
        for (int d=0;d<=warpcluster::maxDim;++d) at[d]=firstDim[d];
        long in=0, index=0;
        for (int i=0;i<count;++i)
            for (long c=0;c<items[i].countW;++c){
                const int d=items[i].dimW[c];
                l.hostEntries[at[d]++]=WarpEntry{in,index++};
                in+=warpcluster::inputSize(gpu->warp,d);
            }
        l.inW.reserve(totalWIn); l.entries.reserve(size_t(totalW));
        SC_CHECK(cudaMemcpyAsync(l.entries.p,l.hostEntries,sizeof(WarpEntry)*totalW,cudaMemcpyHostToDevice,l.stream));
        size_t atIn=0;
        for (int i=0;i<count;++i){
            size_t n=0;
            for (long c=0;c<items[i].countW;++c) n+=warpcluster::inputSize(gpu->warp,items[i].dimW[c]);
            if (n) SC_CHECK(cudaMemcpyAsync(l.inW.p+atIn,items[i].inW,sizeof(C)*n,cudaMemcpyHostToDevice,l.stream));
            atIn+=n;
        }
    }
    long at2=0,at4=0;
    for (int i=0;i<count;++i){
        const GpuSmallItem& it=items[i];
        if (it.count2) SC_CHECK(cudaMemcpyAsync(l.in2.p+size_t(at2)*per2,it.in2,sizeof(C)*per2*it.count2,cudaMemcpyHostToDevice,l.stream));
        if (it.count4) SC_CHECK(cudaMemcpyAsync(l.in4.p+size_t(at4)*per4,it.in4,sizeof(C)*per4*it.count4,cudaMemcpyHostToDevice,l.stream));
        at2+=it.count2; at4+=it.count4;
    }
    launchDim<2>(plan,ensemble,recursive,l.in2.p,l.out.p,l.status.p,total2,l.stream);
    launchDim<4>(plan,ensemble,recursive,l.in4.p,out4,status4,total4,l.stream);
    for (int d=0;d<=warpcluster::maxDim && totalW;++d)
        if (countDim[d]) warpClusterLaunch(gpu->warp,d,l.inW.p,l.entries.p+firstDim[d],countDim[d],outW,statusW,l.stream);
    at2=at4=0;
    long atW=0;
    for (int i=0;i<count;++i){
        const GpuSmallItem& it=items[i];
        if (it.count2){
            SC_CHECK(cudaMemcpyAsync(it.out,l.out.p+size_t(at2)*nstep,sizeof(C)*nstep*it.count2,cudaMemcpyDeviceToHost,l.stream));
            SC_CHECK(cudaMemcpyAsync(it.status,l.status.p+at2,it.count2,cudaMemcpyDeviceToHost,l.stream));
        }
        if (it.count4){
            SC_CHECK(cudaMemcpyAsync(it.out+size_t(it.count2)*nstep,out4+size_t(at4)*nstep,sizeof(C)*nstep*it.count4,
                                     cudaMemcpyDeviceToHost,l.stream));
            SC_CHECK(cudaMemcpyAsync(it.status+it.count2,status4+at4,it.count4,cudaMemcpyDeviceToHost,l.stream));
        }
        if (it.countW){
            SC_CHECK(cudaMemcpyAsync(it.out+size_t(it.count2+it.count4)*nstep,outW+size_t(atW)*perW,sizeof(C)*perW*it.countW,
                                     cudaMemcpyDeviceToHost,l.stream));
            SC_CHECK(cudaMemcpyAsync(it.status+it.count2+it.count4,statusW+atW,it.countW,cudaMemcpyDeviceToHost,l.stream));
        }
        at2+=it.count2; at4+=it.count4; atW+=it.countW;
    }
    SC_CHECK(cudaEventRecord(l.done,l.stream));
    l.busy=true;
}

bool gpuSmallReady(GpuSmallDevice* gpu, int lane){
    Lane& l=gpu->lane[lane];
    if (!l.busy) return true;
    const cudaError_t s=cudaEventQuery(l.done);
    if (s==cudaErrorNotReady) return false;
    SC_CHECK(s);
    l.busy=false;
    return true;
}

void gpuSmallFinish(GpuSmallDevice* gpu, int lane){
    Lane& l=gpu->lane[lane];
    if (!l.busy) return;
    SC_CHECK(cudaEventSynchronize(l.done));
    l.busy=false;
}
