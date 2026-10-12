#ifndef CCEX_BUNDLE_POWERS_H
#define CCEX_BUNDLE_POWERS_H

// The cluster expansion over bundles, once for both halves of adaptive clustering (review of the duplicates,
// 2026-10-08): the calculation half (bundle_expansion.cpp, expandAdaptivePlan) and the algorithm half
// (algo_adaptive/base_bundles.cpp, bundle_powers) used to carry their own copy of adaptive_expansion.py.
// Integer and set code only, no floating point, so it may be compiled with either half's flags; each half turns a
// Failure into its own exception and message.

#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace bundle_expansion {

enum class Failure {
    BadOrder,           // order < 1
    MaxClusters,        // expansion_max_clusters
    MaxWork,            // expansion_max_work
    BadEdge,            // not 0 <= a < b < n
    DuplicateEdge,
    BadExtra,           // an extra set empty, unsorted or repeating a bundle
    ExtraOutside,       // an extra set naming a bundle outside 0..n-1
    ExtraDisconnected,  // an extra set not connected in the retained graph
    PowerTooLarge,      // a nonzero power outside 32-bit storage
    PowerOverflow,      // an intermediate power outside 64 bits
};

struct Error : std::runtime_error {
    explicit Error(Failure what) : std::runtime_error("bundle expansion"), kind(what) {}
    Failure kind;
};

struct Limits {
    long long max_clusters;
    long long max_work;
};

using Edge = std::pair<int, int>;

// Every connected set of up to `order` bundles, plus the extras, closed downward under removing one bundle; ordered
// by size, then lexicographically. The work counter and the cluster count are checked where adaptive_expansion.py
// checks them, so a run over a limit fails at the same point in both halves.
std::vector<std::vector<int>> family(int n, const std::vector<Edge>& edges, int order, const Limits& limits,
                                     const std::vector<std::vector<int>>& extras = {});

// bundle set -> nonzero inclusion-exclusion power. Order 1 or 2 without extras takes the closed form (singles
// 1 - degree, every edge +1); otherwise the family's powers, from the largest set down.
std::map<std::vector<int>, int> powers(int n, const std::vector<Edge>& edges, int order, const Limits& limits,
                                       const std::vector<std::vector<int>>& extras = {});

}  // namespace bundle_expansion

#endif
