// Method D's judgement (policy_regroup, (internal notes) stage 2):
// every function on its own, with NaN, inf, zero and ties, and whole rounds on a bath
// whose answer is known in closed form.
//
// The closed form: let every site set T carry a factor f(T) (sites f({i}) = a_i, every
// other set 1 unless set) and let a probe of S be the product of f over the subsets of
// S. Then the order-2 bundle expansion of a partition with every pair connected is the
// product of f(T) over the T that touch at most two bundles, and a direct solution is
// the product over all of them - so what each candidate predicts can be written down
// without the expansion code being tested.
#include "check.h"
#include "mkdirs.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "base_bundles.h"
#include "base_numeric.h"
#include "policy_regroup.h"

namespace {
using namespace adaptive;
using namespace adaptive::regroup;
const double inf = std::numeric_limits<double>::infinity();
const double nan_ = std::numeric_limits<double>::quiet_NaN();

Curves one(std::vector<std::vector<std::pair<double, double> > > rows) { return rows; }

// ---------------------------------------------------------------------------
// the closed-form bath
// ---------------------------------------------------------------------------
struct Bath {
    int sites = 0;
    int states = 1;
    std::vector<double> x;                                  // positions on a line
    std::map<Sites, std::vector<double> > f;                // factor per state; absent = 1

    double factor(const Sites& t, int s) const {
        const auto found = f.find(t);
        return found == f.end() ? 1. : found->second[static_cast<std::size_t>(s)];
    }
    std::vector<Sites> subsets(const Sites& set) const {
        std::vector<Sites> out;
        for (unsigned long mask = 1; mask < (1UL << set.size()); ++mask) {
            Sites t;
            for (std::size_t k = 0; k < set.size(); ++k)
                if ((mask >> k) & 1UL)
                    t.push_back(set[k]);
            out.push_back(t);
        }
        return out;
    }
    double direct(const Sites& set, int s) const {
        double v = 1.;
        for (const Sites& t : subsets(set))
            v *= factor(t, s);
        return v;
    }
    double model(const Partition& p, int s) const {
        Sites all;
        for (const Sites& b : p)
            all.insert(all.end(), b.begin(), b.end());
        std::sort(all.begin(), all.end());
        double v = 1.;
        for (const Sites& t : subsets(all)) {
            std::set<int> touched;
            for (int site : t)
                for (std::size_t g = 0; g < p.size(); ++g)
                    if (std::find(p[g].begin(), p[g].end(), site) != p[g].end())
                        touched.insert(static_cast<int>(g));
            if (touched.size() <= 2)
                v *= factor(t, s);
        }
        return v;
    }
    // A curve per state: 1 at the first point, the value at the second.
    Curves curves(const std::function<double(int)>& value) const {
        Curves out;
        for (int s = 0; s < states; ++s)
            out.push_back({{1., 0.}, {value(s), 0.}});
        return out;
    }
    Geometry geometry() const {
        Geometry out(static_cast<std::size_t>(sites));
        for (int i = 0; i < sites; ++i) {
            out[static_cast<std::size_t>(i)].xyz = {x[static_cast<std::size_t>(i)], 0., 0.};
            out[static_cast<std::size_t>(i)].components.push_back(out[static_cast<std::size_t>(i)].xyz);
        }
        return out;
    }
};

// Every probe of the bath, handed out as the search asks; a missing one is a KeyError.
struct Fake : Source {
    std::map<Sites, Curves> store, visible;
    long long fetches = 0;
    const std::map<Sites, Curves>& raw() const override { return visible; }
    std::vector<long long> dims;                 // per site; empty: 2 each
    long long work(const std::vector<Sites>& clusters) const override {
        long long total = 0;
        for (const Sites& c : clusters) {
            const long long d = dimension(c);
            total += d*d*d;
        }
        return total;
    }
    long long run_work(const std::map<Sites, int>& terms) const override {
        std::vector<Sites> c;
        for (const auto& t : terms)
            c.push_back(t.first);
        return work(c);
    }
    long long dimension(const Sites& sites) const override {
        if (dims.empty())
            return 1LL << sites.size();
        long long d = 1;
        for (int i : sites)
            d *= dims[static_cast<std::size_t>(i)];
        return d;
    }
    bool fetch(const std::vector<Sites>& missing) override {
        ++fetches;
        for (const Sites& s : missing) {
            const auto found = store.find(s);
            if (found == store.end())
                throw KeyError(sites_repr(s));
            visible[s] = found->second;
        }
        return true;
    }
    // the factors step: the point out of the stored curves (adaptive-regroup-stage3.md)
    long long point_calls = 0;
    long long point_work(const std::vector<Sites>& clusters) const override { return work(clusters); }
    bool point_values(const Partition& units, const std::vector<Edge>& edges, const Point& at,
                      PointValues* out) override {
        ++point_calls;
        *out = point_values_of(units, edges, store, at);
        return true;
    }
    // the except step: P_K and the whole curve from the stored curves (adaptive-regroup-stage3b.md)
    const Adjacency* links = nullptr;
    long long kept_calls = 0;
    double whole_scale = 1., product_scale = 1.;     // != 1: a step that does not reproduce the run
    bool kept(const Partition& units, const std::vector<Sites>& except, Kept* out) override {
        ++kept_calls;
        const std::map<Sites, int> now = expansion_terms(units, *links, ExpansionLimits());
        std::map<Sites, int> rest = now;
        for (const Sites& sites : except)
            rest.erase(sites);
        out->product = scaled_product(rest, store);
        for (std::vector<LongScaled>& state : out->product)
            for (LongScaled& p : state) {
                p.real *= product_scale;
                p.imag *= product_scale;
            }
        out->whole = approximate_weights(now, store, false);
        for (Curve& curve : out->whole)
            for (std::pair<double, double>& z : curve)
                z.first *= whole_scale;
        out->added.clear();
        for (const Sites& sites : except)
            if (visible.count(sites) == 0) {
                visible[sites] = store.at(sites);
                out->added.push_back(sites);
            }
        return true;
    }
    std::vector<Sites> forgotten;
    void forget(const std::vector<Sites>& clusters) override {
        for (const Sites& sites : clusters) {
            visible.erase(sites);
            forgotten.push_back(sites);
        }
    }
};

Fake source_of(const Bath& bath) {
    Fake out;
    std::vector<int> all;
    for (int i = 0; i < bath.sites; ++i)
        all.push_back(i);
    for (const Sites& s : bath.subsets(all))
        out.store[s] = bath.curves([&](int st) { return bath.direct(s, st); });
    out.visible[Sites()] = bath.curves([](int) { return 1.; });
    out.store[Sites()] = out.visible[Sites()];
    return out;
}

Settings settings_for(const Bath& bath, long long cap) {
    Settings s;
    s.cap = cap;
    s.dimensions.assign(static_cast<std::size_t>(bath.sites), 2);
    s.max_regions = 1000;
    s.max_candidates = 1000;
    return s;
}

// What a round should decide, worked out from the closed form alone: the culprits, and
// for the first region that has a candidate (or every region when none does) the
// rejection counts and the accepted partition.
struct Expected {
    std::string stop;
    Partition accepted;
    std::map<std::string, long long> rejected;      // summed over the regions examined
};

Expected expect(const Bath& bath, const Partition& units, const Constants& c, long long cap,
                const std::set<Partition>& visited) {
    Expected out;
    const Geometry geometry = bath.geometry();
    // the worst point: largest |L| over states at the second point (every value is finite here)
    int s_star = 0;
    double largest = -1.;
    long long u_old = 0;
    for (int s = 0; s < bath.states; ++s) {
        const double v = std::fabs(bath.model(units, s));
        if (v > largest) { largest = v; s_star = s; }
        u_old += v > c.limit ? 1 : 0;
    }
    const double x_old = std::max(largest-1., 0.);
    // culprits: the cross factors of each pair, positive log only
    std::vector<std::pair<std::pair<double, Sites>, std::pair<int, int> > > ranked;
    for (int i = 0; i < static_cast<int>(units.size()); ++i)
        for (int j = i+1; j < static_cast<int>(units.size()); ++j) {
            Sites both = units[static_cast<std::size_t>(i)];
            both.insert(both.end(), units[static_cast<std::size_t>(j)].begin(), units[static_cast<std::size_t>(j)].end());
            std::sort(both.begin(), both.end());
            const double g = std::log(std::fabs(bath.model({units[static_cast<std::size_t>(i)],
                                                            units[static_cast<std::size_t>(j)]}, s_star) /
                                                (bath.direct(units[static_cast<std::size_t>(i)], s_star) *
                                                 bath.direct(units[static_cast<std::size_t>(j)], s_star))));
            if (g > 0.)
                ranked.push_back({{-g, both}, {i, j}});
        }
    std::sort(ranked.begin(), ranked.end());
    if (ranked.empty()) {
        out.stop = "no_positive_culprit";
        return out;
    }
    std::set<Partition> tried;
    for (const auto& r : ranked) {
        const std::vector<int> region = region_of(units, r.second.first, r.second.second, geometry);
        Partition reg;
        for (int k : region) reg.push_back(units[static_cast<std::size_t>(k)]);
        reg = canonical(reg);
        if (!tried.insert(reg).second)
            continue;
        const std::vector<int> outside = neighbours_of(units, region, geometry, c.neighbours);
        Sites rs;
        for (const Sites& b : reg) rs.insert(rs.end(), b.begin(), b.end());
        std::sort(rs.begin(), rs.end());
        const auto mean_error = [&](const Partition& p, const Sites& sites) {
            // one time point carries the difference; the other is exactly 1 on both sides
            double m = 0., d = 0.;
            for (int s = 0; s < bath.states; ++s) { m += bath.model(p, s); d += bath.direct(sites, s); }
            return std::fabs(m/bath.states-d/bath.states)/2.;
        };
        const double now = mean_error(reg, rs);
        std::vector<double> before;
        std::vector<Sites> with_sites;
        for (int k : outside) {
            Partition w = reg; w.push_back(units[static_cast<std::size_t>(k)]); w = canonical(w);
            Sites ws = rs; ws.insert(ws.end(), units[static_cast<std::size_t>(k)].begin(), units[static_cast<std::size_t>(k)].end());
            std::sort(ws.begin(), ws.end());
            with_sites.push_back(ws);
            before.push_back(mean_error(w, ws));
        }
        bool found = false;
        double best_e = 0.;
        Partition best;
        for (const Partition& cand : candidates(reg, std::vector<long long>(static_cast<std::size_t>(bath.sites), 2), cap)) {
            const Partition whole = replace(units, region, cand);
            if (visited.count(whole)) { ++out.rejected["revisit"]; continue; }
            bool phys = true;
            for (int s = 0; s < bath.states; ++s)
                phys = phys && std::fabs(bath.model(cand, s)) <= c.limit;
            if (!phys) { ++out.rejected["unphysical"]; continue; }
            const double e = mean_error(cand, rs);
            if (!(e < now-c.epsilon)) { ++out.rejected["a"]; continue; }
            bool worse = false;
            for (std::size_t k = 0; k < outside.size(); ++k) {
                Partition w = cand; w.push_back(units[static_cast<std::size_t>(outside[k])]);
                if (mean_error(canonical(w), with_sites[k]) > before[k]+c.epsilon) worse = true;
            }
            if (worse) { ++out.rejected["b"]; continue; }
            long long u = 0;
            for (int s = 0; s < bath.states; ++s) u += std::fabs(bath.model(whole, s)) > c.limit ? 1 : 0;
            const double x = std::max(std::fabs(bath.model(whole, s_star))-1., 0.);
            if (!(u <= u_old && x < x_old)) { ++out.rejected["c"]; continue; }
            if (!found || e < best_e || (e == best_e && cand < best)) { found = true; best_e = e; best = cand; }
        }
        if (found) {
            out.accepted = replace(units, region, best);
            return out;
        }
    }
    out.stop = "no_candidate";
    return out;
}

std::map<std::string, long long> rejected_of(const Value& row) {
    std::map<std::string, long long> out;
    for (std::size_t i = 0; i < row.at("regions").size(); ++i) {
        const Value& region = row.at("regions").at(i);
        if (!region.has("rejected"))
            continue;
        for (const std::string& key : region.at("rejected").keys())
            out[key] += region.at("rejected").at(key).as_int();
    }
    return out;
}

// A bath of five sites on a line, the factors drawn from a fixed generator.
Bath random_bath(unsigned seed, int states) {
    Bath b;
    b.sites = 5;
    b.states = states;
    b.x = {0., 1.1, 2.5, 4.2, 6.0};
    unsigned state = seed;
    const auto draw = [&state]() {
        state = state*1103515245u+12345u;
        return static_cast<double>((state >> 8) & 0xffff)/65536.;
    };
    for (int i = 0; i < b.sites; ++i)
        for (int s = 0; s < states; ++s)
            b.f[{i}].push_back(0.9+0.1*draw());
    for (int i = 0; i < b.sites; ++i)
        for (int j = i+1; j < b.sites; ++j)
            for (int s = 0; s < states; ++s)
                b.f[{i, j}].push_back(0.85+0.35*draw());
    for (int i = 0; i < b.sites; ++i)
        for (int j = i+1; j < b.sites; ++j)
            for (int k = j+1; k < b.sites; ++k)
                for (int s = 0; s < states; ++s)
                    b.f[{i, j, k}].push_back(0.8+0.3*draw());
    return b;
}

// ---------------------------------------------------------------------------
// the functions one by one
// ---------------------------------------------------------------------------
void test_canonical() {
    CHECK_EQ(canonical({{3, 1}, {0}, {2, 4}}), (Partition{{0}, {1, 3}, {2, 4}}));
    CHECK_RAISES(canonical({{0}, {}}), ValueError);
    CHECK_RAISES(canonical({{0, 1}, {1}}), ValueError);
}

void test_curves() {
    // physical: finite and |L| <= limit, at every point; NaN and inf are not physical
    CHECK(physical(one({{{1., 0.}, {1.01, 0.}}}), 1.01));
    CHECK(!physical(one({{{1., 0.}, {1.0100001, 0.}}}), 1.01));
    CHECK(!physical(one({{{1., 0.}, {nan_, 0.}}}), 1.01));
    CHECK(!physical(one({{{1., 0.}, {0., inf}}}), 1.01));
    // the worst point: the first non-finite one in state then time order...
    Point p = worst_point(one({{{2., 0.}, {1., 0.}}, {{1., 0.}, {nan_, 0.}}, {{inf, 0.}, {1., 0.}}}));
    CHECK_EQ(p.state, 1); CHECK_EQ(p.step, 1); CHECK(p.nonfinite);
    // ...else the largest |L|, the first of equals
    p = worst_point(one({{{1., 0.}, {0., 1.5}}, {{1.5, 0.}, {1., 0.}}}));
    CHECK_EQ(p.state, 0); CHECK_EQ(p.step, 1); CHECK(!p.nonfinite);
    CHECK_EQ(unphysical_states(one({{{1., 0.}, {1.2, 0.}}, {{1., 0.}, {1., 0.}}, {{nan_, 0.}, {1., 0.}}}), 1.01), 2LL);
    CHECK_EQ(excess({{1., 0.}, {0.5, 0.}}), 0.);
    CHECK(std::fabs(excess({{1., 0.}, {0., -1.25}})-0.25) < 1e-15);
    CHECK_EQ(excess({{1., 0.}, {inf, 0.}}), inf);
    CHECK_EQ(excess({{1., 0.}, {nan_, 0.}}), inf);
    CHECK_EQ(excess({{0., 0.}}), 0.);
}

void test_factor() {
    Factor f = culprit_factor({0.5, 0.}, {0.5, 0.}, {0.5, 0.});
    CHECK(std::fabs(f.g-std::log(2.)) < 1e-15); CHECK_EQ(f.cause, std::string());
    // exactly independent: 0, not positive
    f = culprit_factor({0.25, 0.}, {0.5, 0.}, {0.5, 0.});
    CHECK_EQ(f.g, 0.);
    f = culprit_factor({nan_, 0.}, {1., 0.}, {1., 0.});
    CHECK_EQ(f.g, inf); CHECK_EQ(f.cause, std::string("nonfinite_probe"));
    f = culprit_factor({1., 0.}, {inf, 0.}, {1., 0.});
    CHECK_EQ(f.g, inf); CHECK_EQ(f.cause, std::string("nonfinite_probe"));
    f = culprit_factor({1., 0.}, {0., 0.}, {1., 0.});
    CHECK_EQ(f.g, inf); CHECK_EQ(f.cause, std::string("zero_bundle"));
    f = culprit_factor({0., 0.}, {0.5, 0.}, {1., 0.});
    CHECK_EQ(f.g, -inf); CHECK_EQ(f.cause, std::string("zero_pair"));
    // L_a L_b underflows as a plain product; the scaled one does not
    f = culprit_factor({1e-300, 0.}, {1e-200, 0.}, {1e-200, 0.});
    CHECK(std::isfinite(f.g)); CHECK(std::fabs(f.g-100.*std::log(10.)) < 1e-9);
    // a factor past a double's range is still an ordinary log
    f = culprit_factor({1e-300, 0.}, {1e300, 0.}, {1e300, 0.});
    CHECK(std::fabs(f.g-(-900.*std::log(10.))) < 1e-9); CHECK_EQ(f.cause, std::string());
}

void test_culprits() {
    // three singleton bundles, every pair connected; two pairs positive, one negative
    const Partition p{{0}, {1}, {2}};
    std::map<Sites, Curves> raw;
    const auto put = [&raw](Sites s, double v) { raw[s] = Curves{{{1., 0.}, {v, 0.}}}; };
    put({0}, 0.5); put({1}, 0.5); put({2}, 0.5);
    put({0, 1}, 0.4);   // log 1.6
    put({0, 2}, 0.2);   // log 0.8 < 0
    put({1, 2}, 0.4);   // log 1.6: a tie with {0, 1}, broken by the union's sites
    std::vector<Culprit> c = culprits(p, {{0, 1}, {0, 2}, {1, 2}}, raw, Point{0, 1, false});
    CHECK_EQ(c.size(), std::size_t(2));
    CHECK_EQ(c[0].a, 0); CHECK_EQ(c[0].b, 1); CHECK_EQ(c[1].a, 1); CHECK_EQ(c[1].b, 2);
    put({1, 2}, 0.45);  // now the larger
    c = culprits(p, {{0, 1}, {0, 2}, {1, 2}}, raw, Point{0, 1, false});
    CHECK_EQ(c[0].a, 1); CHECK_EQ(c[0].b, 2);
    // +inf ranks first; -inf and 0 are not culprits
    put({0, 2}, nan_);
    put({1, 2}, 0.);
    c = culprits(p, {{0, 1}, {0, 2}, {1, 2}}, raw, Point{0, 1, false});
    CHECK_EQ(c.size(), std::size_t(2));
    CHECK_EQ(c[0].a, 0); CHECK_EQ(c[0].b, 2); CHECK_EQ(c[0].factor.g, inf);
    CHECK_EQ(c[1].a, 0); CHECK_EQ(c[1].b, 1);
    // only connected pairs are culprits
    c = culprits(p, {{1, 2}}, raw, Point{0, 1, false});
    CHECK(c.empty());
    // a missing probe is an error, not a pass
    raw.erase({0, 1});
    CHECK_RAISES(culprits(p, {{0, 1}}, raw, Point{0, 1, false}), KeyError);
}

void test_geometry() {
    Geometry g(4);
    g[0].xyz = {0., 0., 0.}; g[1].xyz = {1., 0., 0.}; g[2].xyz = {3., 0., 0.}; g[3].xyz = {6., 0., 0.};
    for (Site& s : g) s.components.push_back(s.xyz);
    CHECK_EQ(distance(g, 0, 2), 3.);
    CHECK_EQ(closeness(g, {0, 1}, {2, 3}), 2.);
    // the rdip rule: strictly between rdipcut and rdip, any component pair
    Adjacency a = distance_adjacency(g, 3.f, 0.f);
    CHECK_EQ(a[0], (std::vector<int>{1}));            // 3 is not < 3
    CHECK_EQ(a[1], (std::vector<int>{0, 2}));
    a = distance_adjacency(g, 10.f, 1.f);
    CHECK_EQ(a[0], (std::vector<int>{2, 3}));         // 1 is not > 1
    g[3].components.push_back({0.5, 0., 0.});        // an attached spin near site 0
    a = distance_adjacency(g, 1.f, 0.f);
    CHECK_EQ(a[0], (std::vector<int>{3}));
    g[3].components.push_back({0., 0., 0.});          // and one on top of it
    CHECK_RAISES(distance_adjacency(g, 1.f, 0.f), ValueError);
    // connections of a partition come from its sites, whatever bundles they were in before
    Adjacency line(4);
    line[0] = {1}; line[1] = {0, 2}; line[2] = {1, 3}; line[3] = {2};
    CHECK((connections({{0, 3}, {1}, {2}}, line) == std::vector<Edge>{{0, 1}, {0, 2}, {1, 2}}));
    CHECK((connections({{0}, {1}, {2, 3}}, line) == std::vector<Edge>{{0, 1}, {1, 2}}));
    CHECK((connections({{0}, {2}}, line) == std::vector<Edge>{}));
    CHECK_EQ(complete_adjacency(3)[1], (std::vector<int>{0, 2}));
}

void test_region() {
    Geometry g(6);
    const double x[6] = {0., 1., 3., 3.5, 9., 20.};
    for (int i = 0; i < 6; ++i) { g[i].xyz = {x[i], 0., 0.}; g[i].components.push_back(g[i].xyz); }
    const Partition p{{0}, {1}, {2}, {3}, {4}, {5}};
    CHECK_EQ(region_of(p, 0, 1, g), (std::vector<int>{0, 1, 2}));
    CHECK_EQ(region_of(p, 2, 3, g), (std::vector<int>{2, 3, 1}));
    CHECK_EQ(neighbours_of(p, {0, 1, 2}, g, 3), (std::vector<int>{3, 4, 5}));
    CHECK_EQ(neighbours_of(p, {0, 1, 2}, g, 2), (std::vector<int>{3, 4}));
    // equal distances: the bundle with the smaller site list
    Geometry h(4);
    const double y[4] = {0., 1., 2., -1.};
    for (int i = 0; i < 4; ++i) { h[i].xyz = {y[i], 0., 0.}; h[i].components.push_back(h[i].xyz); }
    CHECK_EQ(region_of({{0}, {1}, {2}, {3}}, 0, 1, h), (std::vector<int>{0, 1, 2}));   // {2} and {3} both at 1
    // two bundles: the region is the pair
    CHECK_EQ(region_of({{0}, {1}}, 0, 1, h), (std::vector<int>{0, 1}));
    CHECK(neighbours_of({{0}, {1}}, {0, 1}, h, 3).empty());
}

void test_candidates() {
    const std::vector<long long> two(12, 2);
    // a full three-bundle region: 9 candidates at K = 2, 36 at K = 3 (review 7 counted the prototype's)
    CHECK_EQ(candidates({{0, 1}, {2, 3}, {4, 5}}, two, 4).size(), std::size_t(9));
    CHECK_EQ(candidates({{0, 1, 2}, {3, 4, 5}, {6, 7, 8}}, two, 8).size(), std::size_t(36));
    // singletons: the three merges and nothing else (moves are merges, swaps change nothing)
    const std::vector<Partition> merges = candidates({{0}, {1}, {2}}, two, 4);
    CHECK_EQ(merges, (std::vector<Partition>{{{0}, {1, 2}}, {{0, 1}, {2}}, {{0, 2}, {1}}}));
    // no candidate fits: every merge is over the cap
    CHECK(candidates({{0}, {1}}, two, 2).empty());
    // a split of a two-site bundle, and a merge of the rest within the cap
    CHECK_EQ(candidates({{0, 1}, {2}}, two, 4),
             (std::vector<Partition>{{{0}, {1}, {2}}, {{0}, {1, 2}}, {{0, 2}, {1}}}));
    // mixed dimensions: site 2 is three-dimensional, the cap 6. A swap of 0 (2) and 2 (3)
    // would put 1 and 2 (dimension 6) and 0 alone - allowed; moving 2 into {0, 1} gives 12 - refused.
    std::vector<long long> mixed{2, 2, 3};
    const std::vector<Partition> m = candidates({{0, 1}, {2}}, mixed, 6);
    CHECK_EQ(m, (std::vector<Partition>{{{0}, {1}, {2}}, {{0}, {1, 2}}, {{0, 2}, {1}}}));
    for (const Partition& p : m)
        for (const Sites& b : p)
            CHECK(bundle_dimension(b, mixed) <= 6);
    // replace keeps the rest and comes back canonical
    CHECK_EQ(replace({{0}, {1, 4}, {2}, {3}}, {0, 2}, {{0, 2}}), (Partition{{0, 2}, {1, 4}, {3}}));
}

void test_local_error() {
    // mean over states first (complex), then |.|, then the mean over time points
    const Curves model{{{1., 0.}, {0.5, 0.5}}, {{1., 0.}, {0.5, -0.5}}};
    const Curves direct{{{1., 0.}, {0.25, 0.}}, {{1., 0.}, {0.25, 0.}}};
    LocalError e = local_error(model, direct);
    CHECK(e.verifiable); CHECK_EQ(e.value, 0.125);
    CHECK_EQ(local_error(Curves{{{1., 0.}, {inf, 0.}}}, Curves{{{1., 0.}, {1., 0.}}}).value, inf);
    CHECK_EQ(local_error(Curves{{{1., 0.}, {nan_, 0.}}}, Curves{{{1., 0.}, {1., 0.}}}).value, inf);
    CHECK(!local_error(Curves{{{1., 0.}, {1., 0.}}}, Curves{{{1., 0.}, {nan_, 0.}}}).verifiable);
    CHECK(!local_error(Curves{{{1., 0.}, {1., 0.}}}, Curves{{{inf, 0.}, {1., 0.}}}).verifiable);
    CHECK_RAISES(local_error(Curves{}, Curves{}), ValueError);
}

void test_conditions() {
    // (a): strictly more than epsilon; NaN and inf decide nothing in its favour
    CHECK(improves(0.1, 0.3, 1e-4));
    CHECK(!improves(0.2999, 0.3, 1e-4));                     // a fall of 1e-4 is not more than 1e-4
    CHECK(!improves(nan_, 0.3, 1e-4));
    CHECK(!improves(0.1, nan_, 1e-4));
    CHECK(improves(0.1, inf, 1e-4));                          // from a broken model to a finite one
    CHECK(!improves(inf, inf, 1e-4));
    // (b): fails only when an error grows by more than epsilon
    CHECK(neighbours_hold({0.2, 0.1}, {0.2, 0.15}, 1e-4));
    CHECK(neighbours_hold({0.20005}, {0.2}, 1e-4));
    CHECK(!neighbours_hold({0.2002}, {0.2}, 1e-4));
    CHECK(!neighbours_hold({inf}, {0.2}, 1e-4));
    CHECK(neighbours_hold({inf}, {inf}, 1e-4));               // the prototype's any(after > before + eps)
    CHECK(neighbours_hold({}, {}, 1e-4));
    CHECK_RAISES(neighbours_hold({0.1}, {}, 1e-4), ValueError);
    // (c): U may stay, the excess has to fall
    CHECK(whole_improves(2, 0.1, 2, 0.2));
    CHECK(!whole_improves(3, 0.0, 2, 0.2));
    CHECK(!whole_improves(2, 0.2, 2, 0.2));
    CHECK(!whole_improves(1, inf, 2, inf));
    CHECK(whole_improves(1, 0.5, 2, inf));
}

// ---------------------------------------------------------------------------
// whole rounds against the closed form
// ---------------------------------------------------------------------------
void test_rounds() {
    const Constants c;
    std::map<std::string, long long> seen;
    long long accepted = 0, rounds = 0;
    for (unsigned seed = 1; seed <= 400; ++seed) {
        for (int states : {1, 2}) {
            const Bath bath = random_bath(seed, states);
            const Geometry g = bath.geometry();
            const Adjacency links = complete_adjacency(bath.sites);
            for (const Partition& start : {Partition{{0}, {1}, {2}, {3}, {4}}, Partition{{0, 1}, {2}, {3, 4}},
                                           Partition{{0}, {1, 2}, {3}, {4}}}) {
                Fake source = source_of(bath);
                Ledger ledger;
                ledger.total = inf;
                Search search(settings_for(bath, 4), &g, &links, &source, &ledger);
                for (const auto& t : expansion_terms(start, links, ExpansionLimits()))
                    source.fetch({t.first});
                const Curves curves = expansion(start, links, source.visible, ExpansionLimits());
                if (physical(curves, c.limit))
                    continue;
                ++rounds;
                const Expected want = expect(bath, start, c, 4, {start});
                const Outcome got = search.round(start, curves, {start}, {});
                CHECK_EQ(got.stop, want.stop);
                CHECK_EQ(got.accepted, want.accepted);
                // the counts the closed form gives for the regions it went through; the
                // round went through the same ones, so they are the same numbers
                std::map<std::string, long long> counts = rejected_of(got.row);
                CHECK_EQ(counts.size(), want.rejected.size());
                for (const auto& r : want.rejected) {
                    CHECK_EQ(counts[r.first], r.second);
                    seen[r.first] += r.second;
                }
                accepted += got.stop.empty() ? 1 : 0;
                seen[got.stop.empty() ? "accepted" : got.stop] += 1;
            }
        }
    }
    // every condition, and every outcome, happened somewhere in the draw
    for (const char* what : {"unphysical", "a", "b", "c", "accepted", "no_candidate"})
        if (seen[what] == 0)
            check::fail(__FILE__, __LINE__, std::string("never seen in the rounds: ")+what);
    std::cout << "rounds " << rounds << ", accepted " << accepted << ", rejections unphysical " << seen["unphysical"]
              << " a " << seen["a"] << " b " << seen["b"] << " c " << seen["c"] << " revisit " << seen["revisit"] << "\n";
}

// An unphysical expansion with no positive culprit: a bundle probe over 1 (not a
// physical probe, but the rule has to say what happens).
void test_no_positive_culprit() {
    Bath b;
    b.sites = 3; b.x = {0., 1., 2.};
    b.f[{0}] = {1.25}; b.f[{1}] = {0.5}; b.f[{2}] = {1.};
    b.f[{0, 1}] = {1.}; b.f[{1, 2}] = {0.5};           // the only cross factors: 1 and 1/2
    const Geometry g = b.geometry();
    const Adjacency links = complete_adjacency(3);
    Fake source = source_of(b);
    Ledger ledger; ledger.total = inf;
    Search search(settings_for(b, 4), &g, &links, &source, &ledger);
    const Partition start{{0}, {1}, {2}};
    for (const auto& t : expansion_terms(start, links, ExpansionLimits())) source.fetch({t.first});
    Curves curves = expansion(start, links, source.visible, ExpansionLimits());
    curves[0][1].first = 1.25;                         // |L| over the limit
    const Outcome got = search.round(start, curves, {start}, {});
    CHECK_EQ(got.stop, std::string("no_positive_culprit"));
}

// One unphysical triad and a far site: the merges of the triad are found, and the
// budget, the expansion limit, the revisit block and a missing probe each act.
Bath triad() {
    Bath b;
    b.sites = 4; b.x = {0., 1., 2., 10.};
    for (int i = 0; i < 4; ++i) b.f[{i}] = {0.96875};
    b.f[{0, 1}] = {1.125};
    b.f[{1, 2}] = {1.0625};
    b.f[{0, 2}] = {1.03125};
    b.f[{0, 1, 2}] = {0.75};
    return b;
}

void test_budget_and_limits() {
    const Bath b = triad();
    const Geometry g = b.geometry();
    const Adjacency links = complete_adjacency(4);
    const Partition start{{0}, {1}, {2}, {3}};
    const auto setup = [&](Fake* source) {
        for (const auto& t : expansion_terms(start, links, ExpansionLimits())) source->fetch({t.first});
        return expansion(start, links, source->visible, ExpansionLimits());
    };
    {
        Fake source = source_of(b);
        const Curves curves = setup(&source);
        CHECK(!physical(curves, 1.01));
        Ledger ledger; ledger.total = inf;
        const Outcome got = Search(settings_for(b, 4), &g, &links, &source, &ledger).round(start, curves, {start}, {});
        const Expected want = expect(b, start, Constants(), 4, {start});
        CHECK(got.stop.empty());
        CHECK_EQ(got.accepted, want.accepted);
        CHECK(ledger.spent > 0.);
        CHECK_EQ(ledger.spent_by["culprit"]+ledger.spent_by["region_direct"]+ledger.spent_by["neighbour_direct"]+
                 ledger.spent_by["candidate"]+ledger.spent_by["outside"], ledger.spent);
        // the same with the region's direct solution, but not the neighbours', affordable:
        // nothing is checked without its neighbours, so every candidate is refused for budget
        Fake again = source_of(b);
        const Curves c2 = setup(&again);
        Ledger tight; tight.total = static_cast<double>(again.work({{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
                                                                    {0, 1, 2}}));
        const Outcome refused = Search(settings_for(b, 4), &g, &links, &again, &tight).round(start, c2, {start}, {});
        CHECK_EQ(refused.stop, std::string("budget_limited"));
        CHECK(tight.spent <= tight.total);
        // nothing affordable at all: the culprit probes themselves
        Fake none = source_of(b);
        const Curves c3 = setup(&none);
        Ledger zero; zero.total = 0.;
        CHECK_EQ(Search(settings_for(b, 4), &g, &links, &none, &zero).round(start, c3, {start}, {}).stop,
                 std::string("budget_limited"));
        CHECK_EQ(zero.spent, 0.);
    }
    {
        // the candidate the round chooses, blocked as visited: another one, or none
        Fake source = source_of(b);
        const Curves curves = setup(&source);
        Ledger ledger; ledger.total = inf;
        const Partition first = Search(settings_for(b, 4), &g, &links, &source, &ledger)
                                    .round(start, curves, {start}, {}).accepted;
        Fake again = source_of(b);
        setup(&again);
        const Outcome next = Search(settings_for(b, 4), &g, &links, &again, &ledger)
                                 .round(start, curves, {start, first}, {});
        CHECK(next.accepted != first);
        CHECK(rejected_of(next.row)["revisit"] >= 1);
        const Outcome rolled = Search(settings_for(b, 4), &g, &links, &again, &ledger)
                                   .round(start, curves, {start}, {first});
        CHECK(rolled.accepted != first);
        CHECK(rejected_of(rolled.row)["rolled_back"] >= 1);
    }
    {
        // the expansion limit: four singletons hold 4 + 6 terms; a merge leaves 3 + 3
        // and fits under 9, so every merge passes the limit - and with a limit of 5,
        // none does and the region has nothing to offer
        Fake source = source_of(b);
        const Curves curves = setup(&source);
        Ledger ledger; ledger.total = inf;
        Settings s = settings_for(b, 4);
        s.limits.max_clusters = 9;
        const Outcome fits = Search(s, &g, &links, &source, &ledger).round(start, curves, {start}, {});
        CHECK(fits.stop.empty());
        CHECK_EQ(rejected_of(fits.row)["expansion_limit"], 0LL);
        s.limits.max_clusters = 5;
        const Outcome capped = Search(s, &g, &links, &source, &ledger).round(start, curves, {start}, {});
        CHECK_EQ(capped.stop, std::string("no_candidate"));
        CHECK(rejected_of(capped.row)["expansion_limit"] >= 3);
    }
    {
        // a probe the source does not have: an error, not a pass
        Fake source = source_of(b);
        const Curves curves = setup(&source);
        source.store.erase({0, 1, 2, 3});
        Ledger ledger; ledger.total = inf;
        CHECK_RAISES(Search(settings_for(b, 4), &g, &links, &source, &ledger).round(start, curves, {start}, {}),
                     KeyError);
    }
    {
        // a direct solution that is not finite: the region's candidates are all refused
        Fake source = source_of(b);
        const Curves curves = setup(&source);
        source.store[{0, 1, 2}][0][1].first = nan_;
        Ledger ledger; ledger.total = inf;
        const Outcome got = Search(settings_for(b, 4), &g, &links, &source, &ledger).round(start, curves, {start}, {});
        CHECK(rejected_of(got.row)["direct_nonfinite"] >= 1);
        CHECK(got.stop.empty() || got.stop == "no_candidate");
        CHECK((got.accepted != Partition{{0}, {1, 2}, {3}}));
    }
    {
        // the limits on regions and candidates are refusals, reported as budget
        Fake source = source_of(b);
        const Curves curves = setup(&source);
        Ledger ledger; ledger.total = inf;
        Settings s = settings_for(b, 4);
        s.max_candidates = 2;
        CHECK_EQ(Search(s, &g, &links, &source, &ledger).round(start, curves, {start}, {}).stop,
                 std::string("budget_limited"));
    }
}

void test_drive() {
    const Bath b = triad();
    const Geometry g = b.geometry();
    const Adjacency links = complete_adjacency(4);
    const Partition start{{0}, {1}, {2}, {3}};
    Constants c;
    CHECK_EQ(c.max_changes, 8); CHECK_EQ(c.neighbours, 3); CHECK_EQ(c.epsilon, 1e-4); CHECK_EQ(c.limit, 1.01);
    const auto prepared = [&](Fake* source) {
        for (const auto& t : expansion_terms(start, links, ExpansionLimits())) source->fetch({t.first});
        return expansion(start, links, source->visible, ExpansionLimits());
    };
    {
        // the run confirms: one change, then physical
        Fake source = source_of(b);
        const Curves first = prepared(&source);
        Ledger ledger; ledger.total = inf;
        Search search(settings_for(b, 4), &g, &links, &source, &ledger);
        Value rows = Value::array();
        long long kept = 0;
        const Result r = drive(start, first, &search, c,
            [&](const Partition& p) { return expansion(p, links, source.visible, ExpansionLimits()); },
            [&](bool k) { kept += k ? 1 : 0; }, &rows);
        CHECK_EQ(r.stop, std::string("physical"));
        CHECK_EQ(r.changes, 1LL); CHECK_EQ(r.rollbacks, 0LL); CHECK_EQ(kept, 1LL);
        CHECK(rows.at(0).at("verification").at("confirmed").as_bool());
    }
    {
        // the run contradicts the prediction once: undone, blocked, and the next round
        // takes something else
        Fake source = source_of(b);
        const Curves first = prepared(&source);
        Ledger ledger; ledger.total = inf;
        Search search(settings_for(b, 4), &g, &links, &source, &ledger);
        Value rows = Value::array();
        int calls = 0;
        std::vector<bool> verdicts;
        const Result r = drive(start, first, &search, c,
            [&](const Partition& p) {
                Curves out = expansion(p, links, source.visible, ExpansionLimits());
                if (calls++ == 0)
                    out[0][1].first = 2.;                   // worse than before at s*
                return out;
            },
            [&](bool k) { verdicts.push_back(k); }, &rows);
        CHECK_EQ(r.rollbacks, 1LL);
        CHECK_EQ(verdicts.size(), std::size_t(2));
        CHECK(!verdicts[0]);
        CHECK(rows.at(0).at("accepted").kind() != Value::Kind::Null);
        CHECK(rows.at(1).at("accepted").kind() == Value::Kind::Null ||
              !(rows.at(1).at("accepted").size() == rows.at(0).at("accepted").size() &&
                write_json(rows.at(1).at("accepted"), 0) == write_json(rows.at(0).at("accepted"), 0)));
        CHECK(rejected_of(rows.at(1))["rolled_back"] >= 1);
    }
    {
        // runs that never confirm: every candidate of the region is undone in turn, and
        // the round then has nothing left...
        Fake source = source_of(b);
        const Curves first = prepared(&source);
        Ledger ledger; ledger.total = inf;
        Search search(settings_for(b, 4), &g, &links, &source, &ledger);
        Value rows = Value::array();
        const auto never = [&](const Partition& p) {
            Curves out = expansion(p, links, source.visible, ExpansionLimits());
            out[0][1].first = 3.;
            return out;
        };
        const Result r = drive(start, first, &search, c, never, [](bool) {}, &rows);
        CHECK_EQ(r.stop, std::string("no_candidate"));
        CHECK_EQ(r.rollbacks, 3LL);
        CHECK_EQ(r.changes, 0LL);
        CHECK_EQ(r.final, start);
        // ...and past max_changes rollbacks it is a failure, not a search to go on with
        Fake again = source_of(b);
        const Curves f2 = prepared(&again);
        Search other(settings_for(b, 4), &g, &links, &again, &ledger);
        Constants two = c;
        two.max_changes = 2;
        rows = Value::array();
        CHECK_RAISES(drive(start, f2, &other, two,
                           [&](const Partition& p) {
                               Curves out = expansion(p, links, again.visible, ExpansionLimits());
                               out[0][1].first = 3.;
                               return out;
                           },
                           [](bool) {}, &rows),
                     RuntimeError);
    }
    {
        // the change limit: with room for one change, the run stops after it when the
        // curve is still unphysical - and says physical when it is not
        Fake source = source_of(b);
        const Curves first = prepared(&source);
        Ledger ledger; ledger.total = inf;
        Search search(settings_for(b, 4), &g, &links, &source, &ledger);
        Constants one_change = c;
        one_change.max_changes = 1;
        Value rows = Value::array();
        Result r = drive(start, first, &search, one_change,
            [&](const Partition& p) {
                Curves out = expansion(p, links, source.visible, ExpansionLimits());
                out[0][1].first = 1.02;                    // confirmed, still over the limit
                return out;
            },
            [](bool) {}, &rows);
        CHECK_EQ(r.stop, std::string("change_limit"));
        CHECK_EQ(r.changes, 1LL);
        Fake again = source_of(b);
        const Curves f2 = prepared(&again);
        Search other(settings_for(b, 4), &g, &links, &again, &ledger);
        rows = Value::array();
        r = drive(start, f2, &other, one_change,
                  [&](const Partition& p) { return expansion(p, links, again.visible, ExpansionLimits()); },
                  [](bool) {}, &rows);
        CHECK_EQ(r.stop, std::string("physical"));
        CHECK_EQ(r.changes, 1LL);
        // a physical start changes nothing
        rows = Value::array();
        r = drive(start, Curves{{{1., 0.}, {0.5, 0.}}}, &other, c,
                  [](const Partition&) -> Curves { throw RuntimeError("no run expected"); },
                  [](bool) {}, &rows);
        CHECK_EQ(r.stop, std::string("physical"));
        CHECK_EQ(r.changes, 0LL);
        CHECK_EQ(rows.size(), std::size_t(0));
    }
}

// ---------------------------------------------------------------------------
// stage 3: the culprit factors from a factors step (adaptive-regroup-stage3.md)
// ---------------------------------------------------------------------------
// A factors response as CCEX writes it (bundle_work.cpp adaptiveFactors: an ostream at
// precision 17, row i on rank i % ranks).
void write_factors(const std::string& dir, int ranks, const Partition& bundles, const std::vector<Edge>& edges,
                   const PointValues& values, int state, int step) {
    testing::mkdirs(dir);
    std::vector<std::pair<int, int> > rows;
    std::set<int> connected;
    for (const Edge& e : edges) { connected.insert(e.first); connected.insert(e.second); }
    for (int g : connected) rows.emplace_back(g, -1);
    for (const Edge& e : edges) rows.push_back(e);
    for (int r = 0; r < ranks; ++r) {
        std::ofstream f(dir+"/response.rank-"+std::to_string(r)+".factors.tsv");
        f.precision(17);
        f << "# state=" << state << " step=" << step << " bundles=" << bundles.size() << " edges=" << edges.size() << "\n";
        f << "a\tb\treal\timag\n";
        for (std::size_t i = static_cast<std::size_t>(r); i < rows.size(); i += static_cast<std::size_t>(ranks)) {
            const Sites key = rows[i].second < 0 ? bundles[static_cast<std::size_t>(rows[i].first)]
                                                 : union_of(bundles, {rows[i].first, rows[i].second});
            const std::complex<double> z = values.at(key);
            f << rows[i].first << '\t' << rows[i].second << '\t' << z.real() << '\t' << z.imag() << '\n';
        }
    }
}

bool same_ranking(const std::vector<Culprit>& x, const std::vector<Culprit>& y) {
    if (x.size() != y.size()) return false;
    for (std::size_t i = 0; i < x.size(); ++i)
        if (x[i].a != y[i].a || x[i].b != y[i].b || x[i].factor.cause != y[i].factor.cause ||
            !(x[i].factor.g == y[i].factor.g))
            return false;
    return true;
}

void test_point_culprits() {
    // Four singleton bundles, every pair connected; values from a list that holds NaN, inf, 0,
    // a tie and numbers whose factor is exactly 1. The ranking from the curves (the probe path)
    // and from a factors response written and read back has to be the same, g bit for bit.
    const Partition p{{0}, {1}, {2}, {3}};
    const std::vector<Edge> edges{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    const std::vector<std::complex<double> > pool{
        {0.5, 0.}, {0.25, 0.}, {0.4, 0.1}, {0.3, -0.2}, {1e-200, 0.}, {0., 0.}, {nan_, 0.}, {inf, 0.}, {0.5, 0.},
        {-0.7, 0.3}, {1.3, 0.}, {0.125, 0.}, {0.9999999999999999, 1e-17}};
    unsigned state = 7;
    const auto draw = [&]() { state = state*1103515245u+12345u; return (state >> 8) % pool.size(); };
    long long specials = 0, empty = 0;
    for (int trial = 0; trial < 400; ++trial) {
        std::map<Sites, Curves> raw;
        PointValues direct;
        for (const Sites& key : factor_clusters(p, edges)) {
            const std::complex<double> z = pool[draw()];
            raw[key] = Curves{{{1., 0.}, {0.5, 0.}}, {{1., 0.}, {z.real(), z.imag()}}};
            direct[key] = z;
            specials += (!finite(z) || (z.real() == 0. && z.imag() == 0.)) ? 1 : 0;
        }
        const Point at{1, 1, false};
        const std::vector<Culprit> from_curves = culprits(p, edges, raw, at);
        const std::string dir = "out/factors/trial-"+std::to_string(trial);
        write_factors(dir, 1+trial%3, p, edges, direct, 4, 1);
        const PointValues read = read_factors(dir, 1+trial%3, p, edges, 4, 1);
        const std::vector<Culprit> from_step = culprits(p, edges, read);
        CHECK(same_ranking(from_curves, from_step));
        empty += from_curves.empty() ? 1 : 0;
    }
    CHECK(specials > 0);
    // no positive culprit at all: every pair exactly the product of its bundles
    std::map<Sites, Curves> raw;
    PointValues direct;
    for (int i = 0; i < 4; ++i) direct[{i}] = {0.5, 0.};
    for (const Edge& e : edges) direct[union_of(p, {e.first, e.second})] = {0.25, 0.};
    write_factors("out/factors/none", 2, p, edges, direct, 0, 3);
    CHECK(culprits(p, edges, read_factors("out/factors/none", 2, p, edges, 0, 3)).empty());
    // ties: equal g, ordered by the union's sites
    direct[{0, 1}] = {0.5, 0.};
    direct[{2, 3}] = {0.5, 0.};
    write_factors("out/factors/ties", 1, p, edges, direct, 0, 3);
    const std::vector<Culprit> tied = culprits(p, edges, read_factors("out/factors/ties", 1, p, edges, 0, 3));
    CHECK_EQ(tied.size(), std::size_t(2));
    CHECK_EQ(tied[0].a, 0); CHECK_EQ(tied[0].b, 1); CHECK_EQ(tied[1].a, 2); CHECK_EQ(tied[1].b, 3);
    // a missing value is a KeyError, not a pass
    PointValues partial = direct;
    partial.erase({0, 1});
    CHECK_RAISES(culprits(p, edges, partial), KeyError);
    // only connected bundles are asked for: a bundle without a connection has no row
    CHECK_EQ(factor_clusters({{0}, {1}, {2}}, {{0, 1}}), (std::vector<Sites>{{0}, {0, 1}, {1}}));
    CHECK(factor_clusters({{0}, {1}}, {}).empty());
    std::cout << "point culprits: " << specials << " special values, " << empty << " trials without a culprit\n";
}

void test_read_factors() {
    const Partition p{{0, 1}, {2}, {3}};
    const std::vector<Edge> edges{{0, 1}, {1, 2}};
    PointValues v;
    for (const Sites& key : factor_clusters(p, edges)) v[key] = {0.5, 0.25};
    const auto fresh = [&](const std::string& name, int ranks) {
        const std::string dir = "out/factors/"+name;
        write_factors(dir, ranks, p, edges, v, 2, 5);
        return dir;
    };
    // as written: every row once
    CHECK_EQ(read_factors(fresh("ok", 2), 2, p, edges, 2, 5).size(), std::size_t(5));
    // a rank's file too few or too many
    CHECK_RAISES(read_factors(fresh("ranks", 2), 3, p, edges, 2, 5), RuntimeError);
    CHECK_RAISES(read_factors(fresh("ranks1", 2), 1, p, edges, 2, 5), RuntimeError);
    // a header for another state, step or partition
    CHECK_RAISES(read_factors(fresh("state", 1), 1, p, edges, 3, 5), RuntimeError);
    CHECK_RAISES(read_factors(fresh("step", 1), 1, p, edges, 2, 4), RuntimeError);
    CHECK_RAISES(read_factors(fresh("edges", 1), 1, p, {{0, 1}}, 2, 5), RuntimeError);
    CHECK_RAISES(read_factors(fresh("bundles", 1), 1, Partition{{0}, {1}, {2}, {3}}, edges, 2, 5), RuntimeError);
    const auto rewrite = [](const std::string& path, const std::function<std::string(const std::string&)>& edit) {
        std::ifstream in(path);
        std::stringstream text;
        text << in.rdbuf();
        in.close();
        std::ofstream(path) << edit(text.str());
    };
    // no column header
    std::string dir = fresh("columns", 1);
    rewrite(dir+"/response.rank-0.factors.tsv", [](const std::string& t) {
        const std::size_t a = t.find('\n'), b = t.find('\n', a+1);
        return t.substr(0, a+1)+t.substr(b+1);
    });
    CHECK_RAISES(read_factors(dir, 1, p, edges, 2, 5), RuntimeError);
    // a row missing
    dir = fresh("missing", 1);
    rewrite(dir+"/response.rank-0.factors.tsv", [](const std::string& t) { return t.substr(0, t.rfind('\n', t.size()-2)+1); });
    CHECK_RAISES(read_factors(dir, 1, p, edges, 2, 5), RuntimeError);
    // a row twice (here across the two files) and a row not asked for
    dir = fresh("twice", 2);
    rewrite(dir+"/response.rank-1.factors.tsv", [](const std::string& t) { return t+"0\t-1\t0.5\t0.25\n"; });
    CHECK_RAISES(read_factors(dir, 2, p, edges, 2, 5), RuntimeError);
    dir = fresh("extra", 1);
    rewrite(dir+"/response.rank-0.factors.tsv", [](const std::string& t) { return t+"0\t2\t0.5\t0.25\n"; });
    CHECK_RAISES(read_factors(dir, 1, p, edges, 2, 5), RuntimeError);
    dir = fresh("fields", 1);
    rewrite(dir+"/response.rank-0.factors.tsv", [](const std::string& t) { return t+"1\t2\t0.5\t0.25\t9\n"; });
    CHECK_RAISES(read_factors(dir, 1, p, edges, 2, 5), RuntimeError);
    // non-finite values come back as they were written
    v[{2}] = {nan_, 0.};
    v[{0, 1, 2}] = {-inf, inf};
    const PointValues back = read_factors(fresh("nonfinite", 1), 1, p, edges, 2, 5);
    CHECK(std::isnan(back.at({2}).real()));
    CHECK_EQ(back.at({0, 1, 2}).real(), -inf); CHECK_EQ(back.at({0, 1, 2}).imag(), inf);
}

// ---------------------------------------------------------------------------
// stage 3b: condition (c) from a kept product (adaptive-regroup-stage3b.md)
// ---------------------------------------------------------------------------
int power_in(const std::map<Sites, int>& table, const Sites& s) {
    const auto found = table.find(s);
    return found == table.end() ? 0 : found->second;
}

void test_term_difference() {
    // Six sites on a line, each linked to its first and second neighbours (a sparse graph):
    // connection counts change with the partition, so a bundle outside the region linked to
    // both of two merged bundles changes its singleton's power, terms cancel and terms appear.
    Adjacency line(6);
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            if (i != j && std::abs(i-j) <= 2) line[static_cast<std::size_t>(i)].push_back(j);
    const std::vector<long long> two(6, 2);
    long long cancelled = 0, appeared = 0, outside_singleton = 0, tables = 0;
    for (const Partition& units : {Partition{{0}, {1}, {2}, {3}, {4}, {5}}, Partition{{0, 1}, {2}, {3, 4}, {5}},
                                   Partition{{0}, {1, 2}, {3}, {4, 5}}}) {
        const std::map<Sites, int> now = expansion_terms(units, line, ExpansionLimits());
        for (int a = 0; a+1 < static_cast<int>(units.size()); ++a) {
            const std::vector<int> region{a, a+1};
            Partition reg{units[static_cast<std::size_t>(a)], units[static_cast<std::size_t>(a+1)]};
            std::vector<std::map<Sites, int> > held;
            for (const Partition& cand : candidates(canonical(reg), two, 8))
                held.push_back(expansion_terms(replace(units, region, cand), line, ExpansionLimits()));
            std::vector<const std::map<Sites, int>*> listed;
            for (const auto& t : held) listed.push_back(&t);
            const std::vector<Sites> x = changed_terms(now, listed);
            const std::set<Sites> xs(x.begin(), x.end());
            // X is exactly the terms some table changes
            for (const auto& term : now) {
                bool differs = false;
                for (const auto& t : held) differs = differs || power_in(t, term.first) != term.second;
                CHECK_EQ(xs.count(term.first) != 0, differs);
            }
            // and K x N_c is the candidate's table, term for term
            Sites region_sites = reg[0];
            region_sites.insert(region_sites.end(), reg[1].begin(), reg[1].end());
            for (const auto& t : held) {
                ++tables;
                const std::map<Sites, int> added = new_terms(now, x, t);
                std::set<Sites> keys;
                for (const auto& e : now) keys.insert(e.first);
                for (const auto& e : t) keys.insert(e.first);
                for (const Sites& k : keys) {
                    const int want = power_in(t, k);
                    const int got = (now.count(k) != 0 && xs.count(k) == 0) ? now.at(k) : power_in(added, k);
                    CHECK_EQ(got, want);
                    if (now.count(k) != 0 && want == 0) ++cancelled;               // a term the candidate cancels
                    if (now.count(k) == 0 && want != 0) ++appeared;
                    bool inside = false;
                    for (int site : k) inside = inside || std::find(region_sites.begin(), region_sites.end(), site) != region_sites.end();
                    if (k.size() <= 2 && !inside && power_in(now, k) != want) ++outside_singleton;
                }
                for (const auto& e : added) CHECK(e.second != 0);
            }
        }
    }
    CHECK(tables > 0);
    if (cancelled == 0) check::fail(__FILE__, __LINE__, "no cancelled term in the term-difference test");
    if (appeared == 0) check::fail(__FILE__, __LINE__, "no new term in the term-difference test");
    if (outside_singleton == 0) check::fail(__FILE__, __LINE__, "no singleton outside the region changed power");
    std::cout << "term difference: " << tables << " tables, " << cancelled << " cancelled, " << appeared
              << " new, " << outside_singleton << " outside the region changed\n";
    // a kept term whose power changes: the table X was not worked out from
    const std::map<Sites, int> now{{{0}, 1}, {{1}, -1}, {{0, 1}, 1}};
    CHECK_RAISES(new_terms(now, {{0, 1}}, std::map<Sites, int>{{{0}, 1}, {{1}, 0}}), RuntimeError);
    CHECK((new_terms(now, {{1}, {0, 1}}, std::map<Sites, int>{{{0}, 1}, {{1}, 2}, {{2}, -1}}) ==
           std::map<Sites, int>{{{1}, 2}, {{2}, -1}}));
    CHECK(changed_terms(now, {}).empty());
}

// The product L0 * prod (L_S / L0)^power in __float128 (113-bit mantissa) without scaling: the reference the except
// path's long double reassembly is held to. Values are kept moderate where it is used, so nothing over- or underflows.
struct Quad { __float128 re, im; };
Quad quad_mul(Quad a, Quad b) { return {a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re}; }
Quad quad_inverse(Quad a) { const __float128 n = a.re*a.re+a.im*a.im; return {a.re/n, -a.im/n}; }
Quad quad_power(Quad a, int n) {
    if (n < 0) { a = quad_inverse(a); n = -n; }
    Quad out{1, 0};
    while (n) { if (n & 1) out = quad_mul(out, a); n >>= 1; if (n) a = quad_mul(a, a); }
    return out;
}
Curves quad_reference(const std::map<Sites, int>& weights, const std::map<Sites, Curves>& raw) {
    const Curves& zero = raw.at(Sites());
    Curves out(zero.size());
    for (std::size_t s = 0; s < zero.size(); ++s)
        for (std::size_t t = 0; t < zero[s].size(); ++t) {
            const Quad z0{zero[s][t].first, zero[s][t].second};
            Quad p = z0;
            for (const std::pair<const Sites, int>& term : weights) {
                const std::pair<double, double>& l = raw.at(term.first)[s][t];
                p = quad_mul(p, quad_power(quad_mul(Quad{l.first, l.second}, quad_inverse(z0)), term.second));
            }
            out[s].emplace_back(static_cast<double>(p.re), static_cast<double>(p.im));
        }
    return out;
}
// The largest |a - b| / |b| over the points finite in both, and whether they are non-finite at the same points.
std::pair<double, bool> relative_to(const Curves& a, const Curves& b) {
    double largest = 0.;
    bool same = true;
    for (std::size_t s = 0; s < a.size(); ++s)
        for (std::size_t t = 0; t < a[s].size(); ++t) {
            const std::complex<double> x(a[s][t].first, a[s][t].second), y(b[s][t].first, b[s][t].second);
            if (finite(x) != finite(y)) { same = false; continue; }
            if (finite(x)) largest = std::max(largest, std::abs(x-y)/std::abs(y));
        }
    return {largest, same};
}

void test_kept_product() {
    // approximate_from (long double, as CCEX assembles) with L0 as the start agrees with a __float128 reference and
    // with approximate_weights (double) within 1e-13 and at the same non-finite points - NaN, inf, 0 and negative
    // powers included - and P_K times N_c is the whole table within 1e-12.
    unsigned state = 11;
    const auto draw = [&]() { state = state*1103515245u+12345u; return static_cast<double>((state >> 8) & 0xffff)/65536.; };
    std::map<Sites, Curves> raw;
    raw[Sites()] = Curves{{{1., 0.}, {1., 0.}, {1., 0.}}, {{1., 0.}, {1., 0.}, {1., 0.}}};
    std::vector<Sites> keys{{0}, {1}, {2}, {3}, {0, 1}, {1, 2}, {2, 3}, {0, 3}, {0, 2}};
    for (const Sites& k : keys) {
        Curves c(2, Curve(3));
        for (Curve& curve : c)
            for (auto& z : curve) z = {0.2+draw(), draw()-0.5};
        raw[k] = c;
    }
    raw[{2}][1][1] = {0., 0.};          // a zero under a negative power: not finite
    raw[{0, 3}][0][2] = {nan_, 0.};
    raw[{0, 2}][1][2] = {1e-300, 0.};
    const std::map<Sites, int> all{{{0}, -2}, {{1}, -1}, {{2}, -3}, {{3}, 1}, {{0, 1}, 1}, {{1, 2}, 1}, {{2, 3}, 1},
                                   {{0, 3}, 1}, {{0, 2}, 2}};
    const Curves whole = approximate_weights(all, raw, false);
    const Curves from = approximate_from(scaled_product({}, raw), all, raw);
    const Curves reference = quad_reference(all, raw);
    long long nonfinite = 0;
    for (const Curve& curve : from)
        for (const auto& z : curve) nonfinite += !finite(std::complex<double>(z.first, z.second));
    CHECK(nonfinite >= 2);
    const std::pair<double, bool> to_reference = relative_to(from, reference);
    CHECK(to_reference.second);
    CHECK(to_reference.first <= 1e-13);
    const std::pair<double, bool> to_double = relative_to(from, whole);
    CHECK(to_double.second);
    CHECK(to_double.first <= 1e-13);
    // split into a kept part and the rest, every way of choosing the rest among these keys
    long long splits = 0;
    for (unsigned mask = 0; mask < (1u << keys.size()); mask += 7) {
        std::map<Sites, int> kept, rest;
        for (std::size_t i = 0; i < keys.size(); ++i)
            ((mask >> i) & 1u ? rest : kept)[keys[i]] = all.at(keys[i]);
        const Agreement a = agreement(whole, approximate_from(scaled_product(kept, raw), rest, raw));
        CHECK(a.same_nonfinite);
        CHECK(a.largest <= 1e-12);
        ++splits;
    }
    CHECK(splits > 0);
    CHECK_RAISES(approximate_from(ScaledCurves(1), all, raw), ValueError);
    // agreement: relative above 1, absolute below, non-finite positions
    const Agreement a = agreement(Curves{{{2., 0.}, {0.5, 0.}}}, Curves{{{2.000000002, 0.}, {0.5000000001, 0.}}});
    CHECK(a.same_nonfinite);
    CHECK(std::fabs(a.largest-1e-9) < 1e-15);
    CHECK(!agreement(Curves{{{nan_, 0.}}}, Curves{{{1., 0.}}}).same_nonfinite);
    CHECK(agreement(Curves{{{nan_, 0.}}}, Curves{{{inf, 0.}}}).same_nonfinite);
}

void test_long_exponent_reassembly() {
    // Singletons with powers of order -1000 (1 - connections in a large bundle graph, adaptive-regroup-plan.md, except
    // check diagnosis 2026-10-11): double rounds each (L/L0) and the power multiplies that rounding (~1e-13), the long
    // double reassembly stays below 1e-15 of the __float128 product. P_K holds the pairs, the singletons are excepted.
    unsigned state = 7;
    const auto draw = [&]() { state = state*1103515245u+12345u; return static_cast<double>((state >> 8) & 0xffff)/65536.; };
    std::map<Sites, Curves> raw;
    raw[Sites()] = Curves(3, Curve(40, std::make_pair(1., 0.)));
    const std::vector<Sites> singles{{0}, {1}, {2}, {3}};
    const std::vector<Sites> pairs{{0, 1}, {1, 2}, {2, 3}, {0, 3}};
    for (const std::vector<Sites>* keys : {&singles, &pairs})
        for (const Sites& k : *keys) {
            Curves c(3, Curve(40));
            for (Curve& curve : c)      // |L| within a few 1e-4 of 1: the product of thousands of powers stays moderate
                for (auto& z : curve) {
                    const double phase = 0.02*(draw()-0.5), modulus = 1.+6e-4*(draw()-0.5);
                    z = {modulus*std::cos(phase), modulus*std::sin(phase)};
                }
            raw[k] = c;
        }
    const std::map<Sites, int> kept{{{0, 1}, 1}, {{1, 2}, 1}, {{2, 3}, 1}, {{0, 3}, 1}};
    const std::map<Sites, int> excepted{{{0}, -1000}, {{1}, -999}, {{2}, -1023}, {{3}, -871}};
    std::map<Sites, int> all = kept;
    all.insert(excepted.begin(), excepted.end());
    const Curves reference = quad_reference(all, raw);
    const std::pair<double, bool> long_double = relative_to(approximate_from(scaled_product(kept, raw), excepted, raw),
                                                            reference);
    const std::pair<double, bool> in_double = relative_to(approximate_weights(all, raw, false), reference);
    CHECK(long_double.second);
    CHECK(in_double.second);
    if (kLongDoubleX87) {      // where long double is double the except path is refused instead (adaptive_main.cpp)
        CHECK(long_double.first < 1e-15);
        CHECK(in_double.first > long_double.first);     // the test can tell the two apart
    }
    std::cout << "long exponent reassembly: relative to __float128 " << long_double.first << " (long double), "
              << in_double.first << " (double)\n";
}

void test_read_except() {
    const std::string dir = "out/except/ok";
    testing::mkdirs(dir);
    const std::vector<int> states{1, 2};
    const auto write = [&](const std::string& d, int state, const std::string& extra, int skip) {
        std::ofstream f(d+"/response.except.state-"+std::to_string(state)+".tsv");
        f << "step\treal\timag\texponent\n";
        const long double re[3] = {0.75L, -0.5000000000000000001L, 0.L};
        const long double im[3] = {0.L, 0.25L, 0.L};
        for (int t = 0; t < 3; ++t) {
            if (t == skip) continue;
            char line[160];
            std::snprintf(line, sizeof line, "%d\t%La\t%La\t%lld\n", t, re[t], im[t], 1000LL*state-t);
            f << line;
        }
        f << extra;
    };
    write(dir, 1, "", -1);
    write(dir, 2, "", -1);
    const ScaledCurves back = read_except(dir, states, 3);
    CHECK_EQ(back.size(), std::size_t(2));
    CHECK(back[0][0].real == 0.75L); CHECK_EQ(back[0][0].exponent, 1000LL);
    // the long double mantissa as written, not rounded to the double -0.5
    CHECK(back[1][1].real == -0.5000000000000000001L); CHECK(back[1][1].real != -0.5L);
    CHECK(back[1][1].imag == 0.25L); CHECK_EQ(back[1][1].exponent, 1999LL);
    CHECK_EQ(back[0][0].value(), std::complex<double>(std::ldexp(0.75, 1000), 0.));
    CHECK(back[0][2].real == 0.L && back[0][2].imag == 0.L && back[0][2].exponent == 0LL);
    // a step missing, twice, outside the grid; a header that is not the except step's; a state with no file
    const auto bad = [&](const std::string& name, const std::string& extra, int skip, const std::string& header) {
        const std::string d = "out/except/"+name;
        testing::mkdirs(d);
        write(d, 1, extra, skip);
        write(d, 2, "", -1);
        if (!header.empty()) std::ofstream(d+"/response.except.state-2.tsv") << header;
        return d;
    };
    CHECK_RAISES(read_except(bad("missing", "", 1, ""), states, 3), RuntimeError);
    CHECK_RAISES(read_except(bad("twice", "1\t0x1p-1\t0x0p+0\t3\n", -1, ""), states, 3), RuntimeError);
    CHECK_RAISES(read_except(bad("outside", "3\t0x1p-1\t0x0p+0\t3\n", -1, ""), states, 3), RuntimeError);
    CHECK_RAISES(read_except(bad("fields", "2\t0x1p-1\t0x0p+0\t3\t1\n", -1, ""), states, 3), RuntimeError);
    CHECK_RAISES(read_except(bad("number", "", -1, ""), states, 4), RuntimeError);
    CHECK_RAISES(read_except(bad("header", "", -1, "probe\tstep\treal\timag\n"), states, 3), RuntimeError);
    CHECK_RAISES(read_except(dir, {1, 2, 3}, 3), Error);   // no file for state 3 (an OSError)
}

// Whole rounds of the three paths on the closed-form bath: the factors and except paths have
// to decide what the probe path decides, with every pair connected and with a sparse graph.
void test_paths_agree() {
    const Constants c;
    long long rounds = 0, accepted = 0, excepted = 0, freed = 0;
    Adjacency line(5);
    for (int i = 0; i+1 < 5; ++i) { line[static_cast<std::size_t>(i)].push_back(i+1); line[static_cast<std::size_t>(i+1)].push_back(i); }
    line[0].push_back(2); line[2].push_back(0);
    for (unsigned seed = 1; seed <= 150; ++seed)
        for (int states : {1, 2})
            for (int graph = 0; graph < 2; ++graph) {
                const Bath bath = random_bath(seed, states);
                const Geometry g = bath.geometry();
                const Adjacency links = graph == 0 ? complete_adjacency(bath.sites) : line;
                for (const Partition& start : {Partition{{0}, {1}, {2}, {3}, {4}}, Partition{{0, 1}, {2}, {3, 4}}}) {
                    Fake base = source_of(bath);
                    for (const auto& t : expansion_terms(start, links, ExpansionLimits())) base.fetch({t.first});
                    const Curves curves = expansion(start, links, base.visible, ExpansionLimits());
                    if (physical(curves, c.limit))
                        continue;
                    ++rounds;
                    Outcome got[4];
                    Ledger ledgers[4];
                    Fake sources[4] = {base, base, base, base};
                    for (int path = 0; path < 4; ++path) {
                        Settings s = settings_for(bath, 4);
                        s.culprit_factors = (path & 1) != 0;
                        s.whole_except = (path & 2) != 0;
                        sources[path].links = &links;
                        ledgers[path].total = inf;
                        if (s.whole_except) ledgers[path].freed = 0;
                        got[path] = Search(s, &g, &links, &sources[path], &ledgers[path]).round(start, curves, {start}, {});
                    }
                    for (int path = 1; path < 4; ++path) {
                        CHECK_EQ(got[path].stop, got[0].stop);
                        CHECK_EQ(got[path].accepted, got[0].accepted);
                        CHECK_EQ(write_json(got[path].row.at("culprits"), 0), write_json(got[0].row.at("culprits"), 0));
                        CHECK_EQ(write_json(got[path].row.at("factors"), 0), write_json(got[0].row.at("factors"), 0));
                        CHECK(rejected_of(got[path].row) == rejected_of(got[0]. row));
                        if (got[0].stop.empty()) {
                            const Agreement a = agreement(got[0].predicted, got[path].predicted);
                            CHECK(a.same_nonfinite);
                            CHECK(a.largest <= kExceptTolerance);
                        }
                    }
                    CHECK_EQ(sources[1].point_calls, 1LL);
                    CHECK_EQ(sources[0].point_calls, 0LL);
                    // the except path drops the region's candidate and outside curves when it is done:
                    // it holds less than the probe path (a later region may ask for a dropped curve
                    // again, for its own region, and keep it)
                    CHECK(sources[2].visible.size() < sources[0].visible.size() || sources[2].forgotten.empty());
                    CHECK_EQ(ledgers[2].freed, static_cast<long long>(sources[2].forgotten.size()));
                    freed += ledgers[2].freed;
                    excepted += sources[2].kept_calls;
                    accepted += got[0].stop.empty() ? 1 : 0;
                }
            }
    CHECK(excepted > 0);
    CHECK(freed > 0);
    std::cout << "paths: " << rounds << " rounds, " << accepted << " accepted, " << excepted << " except steps, "
              << freed << " curves freed\n";
}

// The except step's free checks: a whole curve that is not the run's, or a P_K that does not
// give it back, ends the round with an error.
void test_except_checks() {
    const Bath b = triad();
    const Geometry g = b.geometry();
    const Adjacency links = complete_adjacency(4);
    const Partition start{{0}, {1}, {2}, {3}};
    for (int which = 0; which < 3; ++which) {
        Fake source = source_of(b);
        source.links = &links;
        for (const auto& t : expansion_terms(start, links, ExpansionLimits())) source.fetch({t.first});
        const Curves curves = expansion(start, links, source.visible, ExpansionLimits());
        if (which == 1) source.whole_scale = 1.+1e-9;
        if (which == 2) source.product_scale = 1.+1e-9;
        Ledger ledger;
        ledger.total = inf;
        Settings s = settings_for(b, 4);
        s.whole_except = true;
        if (which == 0) {
            const Outcome got = Search(s, &g, &links, &source, &ledger).round(start, curves, {start}, {});
            CHECK(got.stop.empty());
            CHECK_EQ(source.kept_calls, 1LL);
            CHECK_EQ(ledger.except_passes, 0LL);   // -1 + 1: counted from where the driver sets it
        } else {
            CHECK_RAISES(Search(s, &g, &links, &source, &ledger).round(start, curves, {start}, {}), RuntimeError);
        }
    }
    {
        // no room for the except step and one rerun: the region's candidates are refused for budget
        Fake source = source_of(b);
        source.links = &links;
        for (const auto& t : expansion_terms(start, links, ExpansionLimits())) source.fetch({t.first});
        const Curves curves = expansion(start, links, source.visible, ExpansionLimits());
        Ledger probe_ledger;
        probe_ledger.total = inf;
        Search(settings_for(b, 4), &g, &links, &source, &probe_ledger).round(start, curves, {start}, {});
        Fake again = source_of(b);
        again.links = &links;
        for (const auto& t : expansion_terms(start, links, ExpansionLimits())) again.fetch({t.first});
        Ledger tight;
        // what the round spends before stage B: its directs and the candidates' own models
        tight.total = probe_ledger.spent_by["culprit"]+probe_ledger.spent_by["region_direct"]+
                      probe_ledger.spent_by["neighbour_direct"]+probe_ledger.spent_by["candidate"];
        Settings s = settings_for(b, 4);
        s.whole_except = true;
        const Outcome refused = Search(s, &g, &links, &again, &tight).round(start, curves, {start}, {});
        CHECK_EQ(refused.stop, std::string("budget_limited"));
        CHECK_EQ(again.kept_calls, 0LL);
        CHECK(tight.spent <= tight.total);
    }
}

// Fake with a disk cache: what preload is asked for comes from `disk` into raw without a fetch.
struct Warm : Fake {
    std::map<Sites, Curves> disk;
    void preload(const std::set<Sites>& wanted) override {
        for (const Sites& sites : wanted) {
            const auto found = disk.find(sites);
            if (found != disk.end() && visible.count(sites) == 0)
                visible[sites] = found->second;
        }
    }
};

// The second review (2026-10-11): a round that accepts leaves room for its rerun at every
// total (2); a warm disk cache frees what a cold one does (3); a settled change drops the except
// steps of every other partition, and a round on the same partition reuses its own (4).
void test_second_review() {
    const Bath b = triad();
    const Geometry g = b.geometry();
    const Adjacency links = complete_adjacency(4);
    const Partition start{{0}, {1}, {2}, {3}};
    const auto setup = [&](Fake* source) {
        source->links = &links;
        for (const auto& t : expansion_terms(start, links, ExpansionLimits())) source->fetch({t.first});
        return expansion(start, links, source->visible, ExpansionLimits());
    };
    // sites of 2 each, and of 2, 3, 4, 2 (cap 16) - the review's case, where the candidates' reruns differ
    for (int kind = 0; kind < 4; ++kind) {
        const bool except = kind % 2 == 1, mixed = kind >= 2;
        const std::vector<long long> dims = mixed ? std::vector<long long>{2, 3, 4, 2} : std::vector<long long>();
        Settings s = settings_for(b, mixed ? 16 : 4);
        if (mixed) s.dimensions = dims;
        s.whole_except = except;
        Fake full = source_of(b);
        full.dims = dims;
        const Curves curves = setup(&full);
        Ledger open; open.total = inf;
        const Outcome free_round = Search(s, &g, &links, &full, &open).round(start, curves, {start}, {});
        CHECK(free_round.stop.empty());
        const double top = open.spent+2.*static_cast<double>(full.run_work(expansion_terms(free_round.accepted, links,
                                                                                             ExpansionLimits())));
        long long accepted = 0, refused = 0;
        for (int k = 0; k <= 400; ++k) {
            Fake source = source_of(b);
            source.dims = dims;
            setup(&source);
            Ledger ledger; ledger.total = top*k/400.;
            const Outcome got = Search(s, &g, &links, &source, &ledger).round(start, curves, {start}, {});
            if (got.stop.empty()) {
                ++accepted;
                CHECK(ledger.spent+static_cast<double>(source.run_work(expansion_terms(got.accepted, links,
                                                                                       ExpansionLimits()))) <= ledger.total);
            } else {
                ++refused;
                CHECK(ledger.spent <= ledger.total);
            }
        }
        CHECK(accepted > 0 && refused > 0);
    }
    {
        Settings s = settings_for(b, 4);
        s.whole_except = true;
        Fake cold = source_of(b);
        const Curves curves = setup(&cold);
        Ledger a; a.total = inf; a.freed = 0;
        Search(s, &g, &links, &cold, &a).round(start, curves, {start}, {});
        Warm warm;
        static_cast<Fake&>(warm) = source_of(b);
        setup(&warm);
        for (const auto& kv : warm.store)
            if (warm.visible.count(kv.first) == 0) warm.disk.insert(kv);
        Ledger w; w.total = inf; w.freed = 0;
        Search(s, &g, &links, &warm, &w).round(start, curves, {start}, {});
        std::vector<Sites> x = cold.forgotten, y = warm.forgotten;
        std::sort(x.begin(), x.end());
        std::sort(y.begin(), y.end());
        CHECK(!x.empty());
        CHECK_EQ(x, y);
        CHECK_EQ(a.freed, w.freed);
    }
    {
        Settings s = settings_for(b, 4);
        s.whole_except = true;
        Fake source = source_of(b);
        const Curves curves = setup(&source);
        Ledger ledger; ledger.total = inf;
        Search search(s, &g, &links, &source, &ledger);
        const Outcome first = search.round(start, curves, {start}, {});
        CHECK(first.stop.empty());
        CHECK_EQ(source.kept_calls, 1LL);
        search.round(start, curves, {start}, {});                     // the same round again: the same step
        CHECK_EQ(source.kept_calls, 1LL);
        search.settled(start);                                        // still this partition: kept
        search.round(start, curves, {start}, {});
        CHECK_EQ(source.kept_calls, 1LL);
        search.settled(first.accepted);                               // moved on: dropped
        search.round(start, curves, {start}, {});
        CHECK_EQ(source.kept_calls, 2LL);
    }
}
}   // namespace

int main() {
    test_canonical();
    test_curves();
    test_factor();
    test_culprits();
    test_geometry();
    test_region();
    test_candidates();
    test_local_error();
    test_conditions();
    test_rounds();
    test_no_positive_culprit();
    test_budget_and_limits();
    test_drive();
    test_point_culprits();
    test_read_factors();
    test_term_difference();
    test_kept_product();
    test_long_exponent_reassembly();
    test_read_except();
    test_paths_agree();
    test_except_checks();
    test_second_review();
    return check::report("test_regroup");
}
