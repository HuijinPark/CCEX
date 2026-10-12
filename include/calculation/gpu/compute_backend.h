#ifndef CCEX_COMPUTE_BACKEND_H
#define CCEX_COMPUTE_BACKEND_H

#include "../../input/general.h"

// Called after all JSON/CLI overrides, before reading bath/tensor files.
void computeBackendInit(const Config* config);
bool computeUseCuda(const Config* config, int dimension);
void computeBackendFinish();
// the memory of the node's smallest GPU when this run may put clusters on the GPU (backend cuda or auto, a CUDA build
// and a supported method), else 0; for the probe memory budget (A1). One counter serves the node, and the ranks of a
// node may all land on one GPU, so the smallest device bounds what may run at once (review 5f7834 R3).
std::size_t computeDeviceMemory(const Config* config);
// The measured per-cluster times, written out and forgotten. A separate process
// per calculation flushed these on exit; the in-process refinement loop has to
// flush them at the end of each step, or the cost table reads an empty log.
void computeProfileFlush();
void computeProfile(const Config* config, const char* method, const char* backend,
                    int dimension, double build, double setup, double evolution);
void computeClusterBatch(const Config* config, int dimension, int count);

#endif
