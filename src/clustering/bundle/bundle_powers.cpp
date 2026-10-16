/**
 * @file bundle_powers.cpp
 * @brief The cluster expansion over bundles, shared by both halves of adaptive clustering.
 */
#include "../../../include/clustering/bundle/bundle_powers.h"

#include <algorithm>
#include <climits>
#include <set>

namespace bundle_expansion {
namespace {
struct Budget {
    Limits limits;
    long long used = 0;
    void tick() {
        if (used >= limits.max_work) throw Error(Failure::MaxWork);
        ++used;
    }
    void count(long long n) const {
        if (n > limits.max_clusters) throw Error(Failure::MaxClusters);
    }
};

std::vector<std::set<int>> adjacency(int n, const std::vector<Edge>& edges, Budget& budget) {
    std::vector<std::set<int>> out(static_cast<std::size_t>(n < 0 ? 0 : n));
    for (const Edge& e : edges) {
        if (!(0 <= e.first && e.first < e.second && e.second < n)) throw Error(Failure::BadEdge);
        if (!out[static_cast<std::size_t>(e.first)].insert(e.second).second) throw Error(Failure::DuplicateEdge);
        out[static_cast<std::size_t>(e.second)].insert(e.first);
        budget.tick();
    }
    return out;
}

std::vector<std::vector<int>> grow(int n, const std::vector<Edge>& edges, int order, Budget& budget,
                                   const std::vector<std::vector<int>>& extras) {
    if (order < 1) throw Error(Failure::BadOrder);
    budget.count(n);
    if (order >= 2) budget.count(static_cast<long long>(n) + static_cast<long long>(edges.size()));
    const std::vector<std::set<int>> next_to = adjacency(n, edges, budget);
    std::set<std::vector<int>> all, current;
    for (int i = 0; i < n; ++i) all.insert(std::vector<int>{i});
    current = all;
    budget.count(static_cast<long long>(all.size()));
    // Connected sets, level by level (paths and stars, not only cliques); a limit fails, it never truncates.
    for (int size = 2; size <= std::min(order, n) && !current.empty(); ++size) {
        std::set<std::vector<int>> next;
        for (const std::vector<int>& ids : current)
            for (int a : ids)
                for (int b : next_to[static_cast<std::size_t>(a)]) {
                    budget.tick();
                    if (std::binary_search(ids.begin(), ids.end(), b)) continue;
                    std::vector<int> expanded = ids;
                    expanded.insert(std::lower_bound(expanded.begin(), expanded.end(), b), b);
                    next.insert(expanded);
                    all.insert(std::move(expanded));
                    budget.count(static_cast<long long>(all.size()));
                }
        current = std::move(next);
    }
    for (const std::vector<int>& ids : extras) {
        if (ids.empty() || !std::is_sorted(ids.begin(), ids.end()) ||
            std::adjacent_find(ids.begin(), ids.end()) != ids.end())
            throw Error(Failure::BadExtra);
        for (int g : ids)
            if (g < 0 || g >= n) throw Error(Failure::ExtraOutside);
        std::set<int> reached{ids.front()};
        std::vector<int> pending{ids.front()};
        while (!pending.empty()) {
            const int a = pending.back();
            pending.pop_back();
            for (int b : next_to[static_cast<std::size_t>(a)]) {
                budget.tick();
                if (std::binary_search(ids.begin(), ids.end(), b) && reached.insert(b).second) pending.push_back(b);
            }
        }
        if (reached.size() != ids.size()) throw Error(Failure::ExtraDisconnected);
        budget.tick();
        all.insert(ids);
        budget.count(static_cast<long long>(all.size()));
    }
    // Close under subsets BEFORE the powers: in a loop the intersection of two connected clusters can be
    // disconnected and still have a nonzero power. Each level is taken before it is added to.
    std::size_t largest = 0;
    for (const std::vector<int>& ids : all) largest = std::max(largest, ids.size());
    for (std::size_t size = largest; size >= 2; --size) {
        std::vector<std::vector<int>> level;
        for (const std::vector<int>& ids : all)
            if (ids.size() == size) level.push_back(ids);
        for (const std::vector<int>& ids : level)
            for (std::size_t j = 0; j < size; ++j) {
                budget.tick();
                std::vector<int> smaller = ids;
                smaller.erase(smaller.begin() + static_cast<std::ptrdiff_t>(j));
                all.insert(std::move(smaller));
                budget.count(static_cast<long long>(all.size()));
            }
    }
    std::vector<std::vector<int>> out(all.begin(), all.end());
    std::stable_sort(out.begin(), out.end(),
                     [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() < b.size(); });
    return out;
}

std::map<std::vector<int>, int> closed_form(int n, const std::vector<Edge>& edges, int order, Budget& budget) {
    budget.count(static_cast<long long>(n) + (order == 2 ? static_cast<long long>(edges.size()) : 0LL));
    std::map<std::vector<int>, long long> power;
    for (int i = 0; i < n; ++i) power[std::vector<int>{i}] = 1;
    std::set<Edge> seen;
    for (const Edge& e : edges) {
        budget.tick();
        if (!(0 <= e.first && e.first < e.second && e.second < n)) throw Error(Failure::BadEdge);
        if (!seen.insert(e).second) throw Error(Failure::DuplicateEdge);
        if (order == 2) {
            budget.tick();
            budget.tick();
            power[std::vector<int>{e.first}] -= 1;
            power[std::vector<int>{e.second}] -= 1;
            power[std::vector<int>{e.first, e.second}] = 1;
        }
    }
    std::map<std::vector<int>, int> out;
    for (const auto& entry : power)
        if (entry.second != 0) {
            if (entry.second < INT_MIN || entry.second > INT_MAX) throw Error(Failure::PowerTooLarge);
            out[entry.first] = static_cast<int>(entry.second);
        }
    return out;
}
}  // namespace

std::vector<std::vector<int>> family(int n, const std::vector<Edge>& edges, int order, const Limits& limits,
                                     const std::vector<std::vector<int>>& extras) {
    Budget budget{limits};
    return grow(n, edges, order, budget, extras);
}

std::map<std::vector<int>, int> powers(int n, const std::vector<Edge>& edges, int order, const Limits& limits,
                                       const std::vector<std::vector<int>>& extras) {
    Budget budget{limits};
    if ((order == 1 || order == 2) && extras.empty()) return closed_form(n, edges, order, budget);
    const std::vector<std::vector<int>> sets = grow(n, edges, order, budget, extras);
    // Exact 64-bit sums while the powers are built; a power is refused only when it is read and nonzero.
    std::map<std::vector<int>, long long> power;
    for (const std::vector<int>& ids : sets) power[ids] = 1;
    for (auto at = sets.rbegin(); at != sets.rend(); ++at) {
        const long long p = power[*at];
        if (p == 0) continue;
        if (p < INT_MIN || p > INT_MAX) throw Error(Failure::PowerTooLarge);
        const std::size_t size = at->size();
        // Every proper subset, in itertools.combinations order (by size, then lexicographic index tuples).
        for (std::size_t take = 1; take < size; ++take) {
            std::vector<std::size_t> pick(take);
            for (std::size_t i = 0; i < take; ++i) pick[i] = i;
            for (;;) {
                std::vector<int> sub;
                for (std::size_t i : pick) sub.push_back((*at)[i]);
                budget.tick();
                long long& old = power[sub];   // present: the family is closed under subsets
                if ((p > 0 && old < LLONG_MIN + p) || (p < 0 && old > LLONG_MAX + p))
                    throw Error(Failure::PowerOverflow);
                old -= p;
                std::size_t i = take;
                while (i > 0 && pick[i - 1] == size - take + i - 1) --i;
                if (i == 0) break;
                ++pick[i - 1];
                for (std::size_t k = i; k < take; ++k) pick[k] = pick[k - 1] + 1;
            }
        }
    }
    std::map<std::vector<int>, int> out;
    for (const auto& entry : power)
        if (entry.second != 0) out[entry.first] = static_cast<int>(entry.second);
    return out;
}

}  // namespace bundle_expansion
