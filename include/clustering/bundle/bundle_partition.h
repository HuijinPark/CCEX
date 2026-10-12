#ifndef CCEX_ADAPTIVE_PARTITION_H
#define CCEX_ADAPTIVE_PARTITION_H

#include "bundle_defaults.h"
#include <cstddef>
#include <vector>
#include <utility>

struct AdaptiveEdge { int a, b; double strength; bool allowed=true; };
struct AdaptiveTerm { std::vector<int> sites; int power; int bundle_order; };
struct AdaptiveLimits {
    size_t max_clusters=ccex::adaptiveMaxClusters;
    size_t max_work=ccex::adaptiveMaxWork;
};
struct AdaptivePlan {
    std::vector<std::vector<int>> bundles;
    std::vector<int> dimensions;
    std::vector<std::pair<int,int>> edges;
    size_t eligible=0, merged=0, blocked_bundle=0;
    // Auto-only audit. Bounds are summed over all positive inter-bundle pairs,
    // including pairs outside an optional distance restriction.
    size_t positive_pairs=0, omitted_pairs=0;
    size_t radius_excluded_pairs=0;
    long double total_strength=0, omitted_strength=0;
    long double radius_excluded_strength=0;
    std::vector<AdaptiveEdge> selected_connections;
};

// Edges describe ALL retained pair connections. The merge threshold selects
// candidates to merge, never removes inter-bundle connections from the expansion.
AdaptivePlan makeAdaptivePlan(const std::vector<int>& site_dimensions,
    std::vector<AdaptiveEdge> edges, int order, int bundle_cap,
    double merge_cutoff, bool allow_merges=true, AdaptiveLimits limits={});

// Auto mode uses all physical site pairs. allowed flags restrict merge/link
// candidates, but never hide coupling from the omitted-budget audit.
AdaptivePlan makeAdaptiveAutoPlan(const std::vector<int>& site_dimensions,
    std::vector<AdaptiveEdge> physical_edges, int order, int bundle_cap,
    double merge_cutoff, double duration_ms,
    double omitted_budget, bool allow_merges=true, AdaptiveLimits limits={});

// Connected bundle clusters through the requested order, closed under nonempty
// subsets. Disconnected subsets must also participate in inclusion-exclusion.
std::vector<std::vector<int>> adaptiveBundleClusters(const AdaptivePlan&, int order,
                                                   AdaptiveLimits limits={});
// Independent of the partition policy: bundle inclusion-exclusion followed
// by expansion into actual site IDs. No second spin-level tilde calculation.
std::vector<AdaptiveTerm> expandAdaptivePlan(const AdaptivePlan&, int nsites, int order,
                                           AdaptiveLimits limits={},
                                           const std::vector<std::vector<int>>& extra_bundles={});

#endif
