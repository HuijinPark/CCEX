/**
 * @file small_cluster_gpu.h
 * @brief The card side of the packed kernel for small cCCE clusters.
 */
#ifndef CCEX_SMALL_CLUSTER_GPU_H
#define CCEX_SMALL_CLUSTER_GPU_H
// A26 step 1: the card side of the batched small-cluster cCCE (small_cluster_gpu.cu). One process of a node drives the
// card for every rank of the node (small_cluster_batch.cpp): it hands over batches of items - one item is one rank's window
// of packed clusters - in two lanes, so that one batch is copied while the other computes. No Eigen or MPI types cross
// this boundary; complex numbers are smallcluster::C (two doubles).
#include "small_cluster_math.h"
#include "warp_cluster.h"
#include <cstddef>

struct GpuSmallDevice;
// One window of one rank: count2 clusters of dimension 2 at in2 and count4 of dimension 4 at in4 (each Halpha, Hbeta, then
// psi for a single state). The coherences come back to out (the dimension-2 clusters first, nstep values each, then the
// dimension-4 ones) and the statuses to status (same order).
// A40: countW clusters for the warp kernel, their inputs one after the other at inW (warpcluster::inputSize each), their
// total dimensions at dimW; their outputs (nstep * warpcluster::outputPer values each) and statuses follow the packed ones.
struct GpuSmallItem {
    const smallcluster::C* in2; const smallcluster::C* in4;
    long count2, count4;
    smallcluster::C* out; unsigned char* status;
    const smallcluster::C* inW; const unsigned char* dimW; long countW;
};
// cudaSetDevice(device), its context, two lanes. nullptr (and description says why) without a usable card.
GpuSmallDevice* gpuSmallOpen(int device, char* description, int capacity);
void gpuSmallClose(GpuSmallDevice* gpu);
// Page-lock host memory the batches copy from and to (faster copies); false if the driver refuses (the copies then go
// through pageable memory).
bool gpuSmallRegister(void* pointer, std::size_t bytes);
void gpuSmallUnregister(void* pointer);
// Copy the items in, compute, copy the results out, in lane `lane` (0 or 1) - asynchronous. The lane must be idle.
void gpuSmallSubmit(GpuSmallDevice* gpu, int lane, const smallcluster::Plan& plan, bool ensemble, bool recursive,
                    const GpuSmallItem* items, int count);
// A40: the warp kernel's plan, before the first batch with warp clusters (waits for the card to be idle).
void gpuSmallWarpPlan(GpuSmallDevice* gpu, const warpcluster::Plan& plan);
// Whether the lane's batch has finished (results in place); gpuSmallFinish waits for it.
bool gpuSmallReady(GpuSmallDevice* gpu, int lane);
void gpuSmallFinish(GpuSmallDevice* gpu, int lane);
#endif
