/**
 * @file policy_regroup.h
 * @brief Adaptive refinement: the regroup option of merge, which regroups bundles where the expansion gives |L| > 1.
 */
// Method D of (internal notes), an option of the merge policy
// (--regroup, "regroup": true in the Adaptive block): where the bundle expansion gives
// |L| > 1, find the bundle pair that drives it and regroup the bundles around it. With
// the option off nothing here runs and merge is what it was.
//
// The judgement is pure functions of curves, partitions and site positions (stage 2 of
// the plan), so the same code is driven by the merge Driver (real CCEX runs and probes)
// and by the replay of the stored P1 probes (stage 0/4,
// test/cluster_verification/check_adaptive/cpp/regroup_replay.cpp).
#pragma once
#include <array>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "base_bundles.h"
#include "base_json.h"
#include "base_run.h"

namespace adaptive {
namespace regroup {
// The registered values (plan, "결정한 것"): the starting values that reproduce the
// stored-probe test, not production defaults anyone has validated. The limit is an
// operating bound that lets |L| exceed 1 by 1 %, not a rounding allowance, and it is
// kept apart from atol, which D does not use.
struct Constants {
    double epsilon = 1e-4;     // what (a) and (b) count as a change of E_loc
    double limit = 1.01;       // |L| <= limit is physical
    int neighbours = 3;        // the outside bundles (b) checks
    int max_changes = 8;       // verified changes before change_limit
};
const char* const kAlgorithm = "regroup-D-v1";

// Every bundle sorted, the bundles sorted: the one spelling of a partition that
// candidates are deduplicated, compared and ordered by.
using Partition = std::vector<Sites>;
// Throws ValueError on an empty bundle or a site that appears twice.
Partition canonical(const Partition& bundles);

// ---------------------------------------------------------------------------
// Sites: positions and the connections between them
// ---------------------------------------------------------------------------
struct Site {
    std::array<double, 3> xyz{};                       // the host spin
    std::vector<std::array<double, 3> > components;    // host and attached subspins
};
using Geometry = std::vector<Site>;                    // indexed by site ID
// The "sites" of a response.physics.json.
Geometry geometry_of(const Value& physics);

// The distance D measures closeness with: between the host positions, the way the
// prototype does (math.dist of the sites' xyz), written as CCEX writes a distance.
double distance(const Geometry& sites, int a, int b);
// The smallest distance between a site of one set and a site of the other.
double closeness(const Geometry& sites, const Sites& a, const Sites& b);

// site -> the sites it is linked to, sorted. A bundle pair is connected when any of
// their sites are linked (plan item 6): the connections of a candidate are worked
// out again from the sites, never carried over from the bundles it replaces.
using Adjacency = std::vector<std::vector<int> >;
// CCEX's distance rule (clustering/bundle/bundle_cluster.cpp, bundle_connections
// distance): two sites are linked when any pair of their components lies strictly
// between rdipcut and rdip. Both are floats in CCEX's Config and are compared as such.
Adjacency distance_adjacency(const Geometry& sites, float rdip, float rdipcut);
// Every pair linked: the stored-probe test's graph.
Adjacency complete_adjacency(int sites);
// The connected bundle pairs (a < b, indices into `bundles`), sorted.
std::vector<Edge> connections(const Partition& bundles, const Adjacency& links);

// ---------------------------------------------------------------------------
// Curves: physical, the worst point, the unphysical states and the excess
// ---------------------------------------------------------------------------
// Every point finite and |L| <= limit.
bool physical(const Curves& curves, double limit);
struct Point {
    int state = -1;
    int step = -1;
    bool nonfinite = false;
};
// Plan item 2: the first non-finite point in state then time order; when there is
// none, the point of largest |L| (the first of equals in the same order).
Point worst_point(const Curves& curves);
// The states with a non-finite point or a point over the limit.
long long unphysical_states(const Curves& curves, double limit);
// max_t (|L| - 1)+, and +inf when a point is not finite.
double excess(const Curve& curve);

// ---------------------------------------------------------------------------
// The culprit (plan item 3)
// ---------------------------------------------------------------------------
// g = log|L_ab / (L_a L_b)| at one point (cCCE: L0 = 1). Non-finite or undefined
// is +inf, a zero factor -inf. `cause` names which: "" for an ordinary value,
// "nonfinite_probe" (a curve not finite), "zero_bundle" (L_a or L_b is 0, the
// factor is undefined), "zero_pair" (L_ab is 0, the factor is 0), "nonfinite_factor"
// (finite curves whose factor does not fit a double).
struct Factor {
    double g = 0.;
    std::string cause;
};
Factor culprit_factor(std::complex<double> pair, std::complex<double> a, std::complex<double> b);
struct Culprit {
    int a = 0;
    int b = 0;
    Factor factor;
};
// The connected pairs with g > 0, largest first, equal g by the union of their
// sites, then by the pair. `at` holds the value at (s*, t*) of every bundle with a
// connection and of every connected pair (KeyError when one is missing): the one place
// the ranking and its ties are decided, whichever way the values were obtained
// (adaptive-regroup-stage3.md 1절).
using PointValues = std::map<Sites, std::complex<double> >;
std::vector<Culprit> culprits(const Partition& bundles, const std::vector<Edge>& edges, const PointValues& at);
// The same from curves: `raw` holds every bundle and every connected pair, and their
// values at the point are read out of it (the probe path).
std::vector<Culprit> culprits(const Partition& bundles, const std::vector<Edge>& edges,
                              const std::map<Sites, Curves>& raw, const Point& at);
// The clusters a ranking reads: every bundle with a connection, every connected pair;
// sorted, each once. A factors step calculates exactly these.
std::vector<Sites> factor_clusters(const Partition& bundles, const std::vector<Edge>& edges);
PointValues point_values_of(const Partition& bundles, const std::vector<Edge>& edges,
                            const std::map<Sites, Curves>& raw, const Point& at);
// What a factors step wrote (adaptive-regroup-stage3.md 1절, 4절): one
// response.rank-R.factors.tsv per rank, each headed "# state=S step=T bundles=B
// edges=E" and "a b real imag", one row per cluster - (a, -1) a bundle, (a, b) a pair,
// numbered by `bundles` - and every row of factor_clusters exactly once over the files.
// Anything else (a file too few or too many, a header that differs, a row missing,
// twice or not asked for) is a RuntimeError, as read_probes is strict.
PointValues read_factors(const std::string& directory, long long ranks, const Partition& bundles,
                         const std::vector<Edge>& edges, int state, int step);

// ---------------------------------------------------------------------------
// The region, its neighbours and the candidates (plan items 4 and 5)
// ---------------------------------------------------------------------------
// The pair and the one other bundle with the smallest site distance to it (equal
// distances: the bundle with the smaller site list); just the pair when it is all.
std::vector<int> region_of(const Partition& bundles, int a, int b, const Geometry& sites);
// The `count` bundles outside the region nearest to it, by the same rule.
std::vector<int> neighbours_of(const Partition& bundles, const std::vector<int>& region,
                               const Geometry& sites, int count);
// One move inside the region: merge two bundles, split one in two (every
// bipartition), move one site, swap two sites. A bundle fits when the product of
// its site dimensions is at most `cap` (both sides of every move are checked);
// empty bundles go; canonical, deduplicated and sorted, the region itself left out.
std::vector<Partition> candidates(const Partition& region, const std::vector<long long>& dimensions,
                                  long long cap);
// `units` with the bundles of `region` replaced by `replacement`, canonical.
Partition replace(const Partition& units, const std::vector<int>& region, const Partition& replacement);

// ---------------------------------------------------------------------------
// D's local error (plan item 7) - not merge's local()
// ---------------------------------------------------------------------------
// mean_t |mean_s L_model - mean_s L_direct|: the state mean taken of the complex
// curves, its distance to the direct solution's, averaged over every time point with
// equal weight. A model that is not finite somewhere is +inf. A direct solution that
// is not finite cannot verify anything: `verifiable` is false and the value 0.
struct LocalError {
    bool verifiable = true;
    double value = 0.;
};
LocalError local_error(const Curves& model, const Curves& direct);

// The three conditions of item 8, each as a pure comparison, so a NaN, an inf or an
// equality lands where it is written to: anything a comparison cannot affirm fails.
// (a) the region's error falls by more than epsilon.
bool improves(double after, double before, double epsilon);
// (b) no neighbour's error grows by more than epsilon.
bool neighbours_hold(const std::vector<double>& after, const std::vector<double>& before, double epsilon);
// (c) no more unphysical states, and the excess of s* falls.
bool whole_improves(long long unphysical_after, double excess_after, long long unphysical_before,
                    double excess_before);

// ---------------------------------------------------------------------------
// Condition (c) from a kept product (adaptive-regroup-stage3b.md 1절)
// ---------------------------------------------------------------------------
// X: the terms of `now` whose power differs in at least one of `tables` (a term a
// table lacks has power 0 there). The two tables are compared whole, so a bundle
// whose connection count changes - its singleton's power - is in X, and a term a
// candidate cancels is too.
std::vector<Sites> changed_terms(const std::map<Sites, int>& now,
                                 const std::vector<const std::map<Sites, int>*>& tables);
// N_c: the terms of `table` with a nonzero power that are not kept (K = now \ except).
// A kept term whose power in `table` differs from `now`'s is a RuntimeError: X was not
// worked out from this table, and P_K would be the wrong product.
std::map<Sites, int> new_terms(const std::map<Sites, int>& now, const std::vector<Sites>& except,
                               const std::map<Sites, int>& table);
// What an except step wrote of P_K: response.except.state-S.tsv for every state, headed
// "step real imag exponent", one row per step exactly once (the mantissa in %La, kept in long double).
ScaledCurves read_except(const std::string& directory, const std::vector<int>& states, long long steps);
// Two curve sets on the same grid: whether they are non-finite at the same points, and
// the largest |a - b| / max(1, |a|) at the points finite in both.
struct Agreement {
    bool same_nonfinite = true;
    double largest = 0.;
};
Agreement agreement(const Curves& a, const Curves& b);
// tau of adaptive-regroup-stage3b.md 4절 (relative): what the except step's whole curve
// and the curve reassembled from P_K may differ from the run by before the round fails.
const double kExceptTolerance = 1e-12;
// An except step's answer: P_K, the whole curve it ran, and which of the excepted
// curves it added to raw (the others were there already).
struct Kept {
    ScaledCurves product;
    Curves whole;
    std::vector<Sites> added;
};

// The site sets the order-2 expansion of these bundles reads (exponent not 0).
std::map<Sites, int> expansion_terms(const Partition& bundles, const Adjacency& links,
                                     const ExpansionLimits& limits);
// That expansion, assembled with scaled products from `raw` (which holds the empty set).
Curves expansion(const Partition& bundles, const Adjacency& links, const std::map<Sites, Curves>& raw,
                 const ExpansionLimits& limits);

// ---------------------------------------------------------------------------
// The search: one round of D, and the loop of rounds
// ---------------------------------------------------------------------------
// Where the curves of site sets come from: the driver's probes, or the stored ones.
class Source {
 public:
    virtual ~Source() {}
    virtual const std::map<Sites, Curves>& raw() const = 0;
    // Takes into raw what of `wanted` can be had without calculating (a cache).
    virtual void preload(const std::set<Sites>& wanted) { (void)wanted; }
    // The D^3 * states * steps of calculating these, summed over the clusters (no
    // zero-order term: cCCE's L0 is 1 and is not calculated).
    virtual long long work(const std::vector<Sites>& clusters) const = 0;
    // The work of running this partition: every term of its expansion.
    virtual long long run_work(const std::map<Sites, int>& terms) const = 0;
    // The Hilbert dimension of one site set, for the report.
    virtual long long dimension(const Sites& sites) const = 0;
    // Calculates `missing` into raw. False: refused (probe count or memory), nothing done.
    virtual bool fetch(const std::vector<Sites>& missing) = 0;
    // --regroup-factors factors (adaptive-regroup-stage3.md): the work of one state of
    // these clusters, and their values at one point, from one factors step. Not in raw:
    // one state's point cannot enter a local error. False: refused, nothing done.
    virtual long long point_work(const std::vector<Sites>& clusters) const;
    virtual bool point_values(const Partition& units, const std::vector<Edge>& edges, const Point& at,
                              PointValues* out);
    // --regroup-whole except (adaptive-regroup-stage3b.md): one except step of `units`;
    // the excepted curves go into raw. False: refused (probe count or memory), nothing done.
    virtual bool kept(const Partition& units, const std::vector<Sites>& except, Kept* out);
    // Drops these curves from raw (the cache freeing of adaptive-regroup-stage3b.md 6절 3).
    virtual void forget(const std::vector<Sites>& clusters);
};

// The running total of the work D asked for (plan stage 5). Every request is
// reserved before it is made and subtracted when it is, so the limit holds for the
// whole run and not for each request on its own.
struct Ledger {
    double total = 0.;
    double spent = 0.;
    std::map<std::string, double> spent_by;   // initial_run, culprit, region_direct, neighbour_direct, candidate, outside, rerun
    long long probes = 0;
    long long max_dimension = 0;
    // Counted, and reported, only where the option asked for them (-1: not reported),
    // so the ledger of a run without them is what it was.
    long long factor_passes = -1;      // factors steps (--regroup-factors factors)
    long long except_passes = -1;      // except steps (--regroup-whole except)
    long long freed = -1;              // curves dropped from raw after their region
    bool fits(double more, double reserve) const { return spent+more+reserve <= total; }
    void charge(const std::string& what, double work) { spent += work; spent_by[what] += work; }
    Value value() const;
};

struct Settings {
    Constants constants;
    long long cap = 0;                         // bundle_max_dim
    std::vector<long long> dimensions;         // per site, attached spins included
    ExpansionLimits limits;
    long long max_regions = 64;                // culprit regions examined per round
    long long max_candidates = 256;            // candidates of one region
    bool trace_curves = false;                 // the round's curve in the trace (the replay)
    // --regroup-factors factors: the culprit factors from Source::point_values.
    bool culprit_factors = false;
    // --regroup-whole except: (c) from Source::kept, and the candidate and outside curves
    // of a region dropped when the region is done (adaptive-regroup-stage3b.md 6절 3).
    bool whole_except = false;
};

// One round's answer: a stop reason, or the accepted partition with what it predicts.
struct Outcome {
    std::string stop;                          // empty when a candidate was accepted
    Partition accepted;
    Curves predicted;                          // the whole expansion of `accepted`
    Point at;
    long long unphysical_before = 0;
    double excess_before = 0.;
    Value row;
};

class Search {
 public:
    Search(Settings settings, const Geometry* sites, const Adjacency* links, Source* source, Ledger* ledger)
        : settings_(std::move(settings)), sites_(sites), links_(links), source_(source), ledger_(ledger) {}
    // Plan items 2-9 for a partition whose curve `curves` is not physical. `visited`
    // and `rolled_back` are whole partitions no candidate may lead to.
    Outcome round(const Partition& units, const Curves& curves, const std::set<Partition>& visited,
                  const std::set<Partition>& rolled_back);
    long long examined() const { return examined_; }
    const Settings& settings() const { return settings_; }
    // A change the run confirmed: `units` is the partition now, and no round comes back to an
    // earlier one (visited), so their except steps are dropped (second review 4).
    void settled(const Partition& units);

 private:
    // Gets these into raw within the ledger. False when they do not fit.
    bool acquire(const std::vector<std::pair<Sites, std::string> >& needed, double reserve);
    Settings settings_;
    const Geometry* sites_;
    const Adjacency* links_;
    Source* source_;
    Ledger* ledger_;
    long long rounds_ = 0;
    long long examined_ = 0;
    // What acquire() brought into raw - calculated, or read from the disk cache - with the purpose
    // it was asked for: the region's share is dropped when the region is done (whole_except only).
    std::vector<std::pair<Sites, std::string> > fetched_;
    // Except steps already made, by (partition, X): a round after a rollback asks again. Only P_K
    // and the added terms are kept; the whole curve is used once, by the check.
    std::map<std::pair<Partition, std::vector<Sites> >, Kept> kept_;
};

struct Result {
    std::string stop;
    Partition final;
    long long changes = 0;
    long long rollbacks = 0;
};
// The loop of plan items 1 and 8-10: run, stop when physical, take a round's
// candidate, run it, and keep it only when the run confirms (c) - U not up, the
// excess of s* down - otherwise go back and block it. `run` gives a partition's curve
// over the selection states; `settle(true)` keeps that run, `settle(false)` drops it.
// The rows of every round go to `rows`.
Result drive(const Partition& start, const Curves& first, Search* search, const Constants& constants,
             const std::function<Curves(const Partition&)>& run,
             const std::function<void(bool)>& settle, Value* rows);

// A number for a trace: JSON has no inf or nan, so those are written as Python's repr.
Value trace_number(double x);
Value partition_value(const Partition& units);
}   // namespace regroup
}   // namespace adaptive
