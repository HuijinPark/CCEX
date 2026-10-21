/**
 * @file warp_cluster_gpu.h
 * @brief The card side of the warp kernel.
 */
#ifndef CCEX_WARP_CLUSTER_GPU_H
#define CCEX_WARP_CLUSTER_GPU_H
// A40: the card side of the warp kernel (warp_cluster_gpu.cu), called by small_cluster_gpu.cu within a batch. CUDA only.
#include "warp_cluster.h"
#include <cuda_runtime.h>

// One cluster of a launch: its input at in + entry.in (values), its output at out + entry.out * nstep * outputPer, its
// status at status[entry.out].
struct WarpEntry { long in, out; };

// The plan into the kernel's constant memory; no warp kernel may be running.
void warpClusterUpload(const warpcluster::Plan& plan);
// The n clusters of total dimension dim listed at entries (device memory), in the stream. The plan is the one uploaded.
void warpClusterLaunch(const warpcluster::Plan& plan, int dim, const smallcluster::C* in, const WarpEntry* entries, long n,
                       smallcluster::C* out, unsigned char* status, cudaStream_t stream);
#endif
