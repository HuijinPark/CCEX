/**
 * @file bundle_defaults.h
 * @brief Bundle clustering: the default limits and check tolerances.
 */
#ifndef CCEX_ADAPTIVE_DEFAULTS_H
#define CCEX_ADAPTIVE_DEFAULTS_H
#include <cstddef>
namespace ccex {
constexpr int adaptiveMaxClusters=200000;
constexpr int adaptiveMaxWork=20000000;
// Bounded checks used when rdip is omitted. These are diagnostic tolerances
// for selected domains, never a full-bath error guarantee.
constexpr double adaptiveCheckAtol=1e-7;
constexpr double adaptiveCheckRtol=1e-5;
constexpr std::size_t adaptiveCheckDomains=8;
constexpr std::size_t adaptiveCheckProbes=64;
constexpr std::size_t adaptiveCheckBytes=64*1024*1024;
}
#endif
