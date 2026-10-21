/**
 * @file accelerator_defaults.h
 * @brief Default sizes of the accelerated paths (GPU buffer cache, refinement sweeps, pCCE detail limit).
 */
#ifndef CCEX_ACCELERATOR_DEFAULTS_H
#define CCEX_ACCELERATOR_DEFAULTS_H
namespace ccex {
constexpr int gpuCacheMiB=64;
constexpr int refinementSweeps=32;
constexpr int partitionDetailLimit=0;   // cluster.pcce.print_trials (A94: the number printed, 0 none; print_limit 32 before)
}
#endif
