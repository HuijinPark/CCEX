/**
 * @file base_bundles.h
 * @brief Adaptive refinement: bundles and curves as the policy sees them (plan, expansion and its cost, curve comparison).
 */
// What the policy uses about bundles and curves: the plan CCEX reports, scaled
// products, curve comparison, and the cluster expansion and its cost.
#pragma once
#include <complex>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#include "base_numeric.h"
#include "base_json.h"
#include "base_run.h"
#include "../../calculation/scaled_product.h"

// ----------------------------------------------------------------------------
// The bounded subset-closed expansion the policy counts terms with
// (plan stage 1, Task 2; adaptive_expansion.py).
//
// Order 1 and 2 take a closed-form shortcut; every other order builds the
// subset-closed family and inverts it (adaptive_expansion._family), which is
// what a run at order 3 or higher needs.
// ----------------------------------------------------------------------------
namespace adaptive {
struct ExpansionLimits {
    long long max_clusters = 200000;
    long long max_work = 20000000;
};

using Edge = std::pair<int, int>;

// bundle index tuple -> inclusion-exclusion power. `extras` are connected bundle
// sets to add on top of the order-sized family, keeping their subset closure;
// an extra forces the general path even at order 1 or 2.
std::map<std::vector<int>, int> bundle_powers(int n, const std::vector<Edge>& edges, int order,
                                             const ExpansionLimits& limits,
                                             const std::vector<std::vector<int> >& extras =
                                                 std::vector<std::vector<int> >());


struct Plan;

// A bounded beam of connected (order+1)-bundle domains - never global coverage.
// Needs the plan's graph, so it is declared here and defined beside the family.
std::map<std::vector<int>, double> next_order_patches(const Plan& plan, int order,
                                                      long long neighbors, long long width);

Sites union_of(const std::vector<Sites>& bundles, const std::vector<int>& indices);
long long bundle_dimension(const Sites& bundle, const std::vector<long long>& site_dimensions);

// Dimensions, and the D^3 work proxies summed from them, are held at int64's
// ceiling rather than wrapped: each is only ever compared - against a budget, a
// cap, or another candidate - and a wrap turns "too big" into a small or negative
// number that wins the comparison. CPython has the exact figure; past the ceiling
// this port promises only "more than anything it is compared with" (plan stage 7).
// Operands are nonnegative.
long long saturating_mul(long long a, long long b);
long long saturating_cube(long long d);

// The same powers, keyed by the sites they cover rather than by bundle index.
std::map<Sites, int> terms(const std::vector<Sites>& bundles, const std::vector<Edge>& edges,
                           int order, const ExpansionLimits& limits);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// What the native binary answers with: the partition it used and the graph it saw
// (plan stage 1; adaptive_refine.py reads this out of response.plan.json).
// ----------------------------------------------------------------------------
namespace adaptive {
struct WeightedEdge {
    int a = 0;
    int b = 0;
    double weight = 0.;
};

struct Plan {
    std::vector<Sites> bundles;
    std::vector<WeightedEdge> edges;
    std::vector<WeightedEdge> candidate_edges;   // plan.get('candidate_edges', plan['edges'])
    bool has_candidate_edges = false;
    std::vector<long long> site_dimensions;
    long long qubit_dimension = 1;
    ExpansionLimits limits;
    bool has_limits = false;

    const std::vector<WeightedEdge>& weights() const {
        return has_candidate_edges ? candidate_edges : edges;
    }
    std::vector<Edge> plain_edges() const {
        std::vector<Edge> out;
        for (const WeightedEdge& edge : edges)
            out.emplace_back(edge.a, edge.b);
        return out;
    }
};
}   // namespace adaptive

// ----------------------------------------------------------------------------
// The scaled complex product the merge policy evaluates its expansion with
// (plan stage 1, Task 1).
//
// A value is carried as (z, exponent) meaning z * 2^exponent, with |z| kept in
// [0.5, 1). A CCE expansion multiplies a long chain of curve ratios together, and
// without the exponent held apart the product would overflow or flush to zero
// long before the answer is wrong.
// ----------------------------------------------------------------------------
namespace adaptive {
struct Scaled {
    std::complex<double> z;
    // Python's int, so no bound to defend: power() multiplies this by a weight up
    // to 2^31 and the products are summed over every term of the expansion. A
    // 64-bit exponent does not obviously overflow, but "does not obviously" is not
    // a contract, and an overflow here would silently return a wrong curve.
    long long exponent = 0;
};

// |z| into [0.5, 1), the rest into the exponent. Zero and a non-finite scale are
// left exactly as they are: frexp has nothing to say about them.
Scaled normalize(std::complex<double> z, long long exponent = static_cast<long long>(0));

// Python's complex multiplication, which is the plain formula. std::complex is
// NOT used for this: libstdc++ adds NaN recovery that CPython does not do, so
// (nan+0j)*(1+0j) would differ.
Scaled multiply(const Scaled& a, const Scaled& b);

// 1/a. Zero or non-finite gives (nan+0j, 0) - the policy reads that as "no ratio".
Scaled inverse(const Scaled& a);

Scaled power(Scaled a, long long n);   // by squaring; n may be negative

// Back to an ordinary complex number. An exponent too large for the result is
// infinity with the component's sign, which is what CPython's ldexp overflow
// gives after adaptive_refine.value() catches it.
std::complex<double> value(const Scaled& a);

bool finite(const std::complex<double>& z);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// Coherence curves: how an expansion is evaluated and how two curves are compared
// (plan stage 1, Task 2; adaptive_refine.py).
// ----------------------------------------------------------------------------
namespace adaptive {
struct Discrepancy {
    double scaled_error = 0.;
    double max_abs_error = 0.;
    bool has_location = false;         // Python's None until something is worst
    int state = 0;
    int step = 0;
};

// The largest scaled difference between two curve sets, and where it is. Both
// grids must be the same shape and not empty.
// One point's error in units of the tolerance: 1 is exactly at the limit. The
// absolute difference comes back with it, because discrepancy wants both.
std::pair<double, double> point_error(const std::complex<double>& x,
                                      const std::complex<double>& y, double atol, double rtol);
Discrepancy discrepancy(const Curves& actual, const Curves& reference, double atol, double rtol);
// One point of a curve set as a complex number.
std::complex<double> curve_at(const Curves& curves, std::size_t state, std::size_t step);
// json.dumps of a discrepancy, the way every report holds it.
Value discrepancy_value(const Discrepancy& found);

// An explicit subset-closed expansion, evaluated with scaled products. `raw` must
// hold the empty cluster and every term.
Curves approximate_weights(const std::map<Sites, int>& weights,
                           const std::map<Sites, Curves>& raw, bool gcce);

Curves approximate(const std::vector<Sites>& bundles, const std::vector<Edge>& edges, int order,
                   const std::map<Sites, Curves>& raw, bool gcce, const ExpansionLimits& limits);

// A cCCE product per state and step, held scaled: adaptive-regroup-stage3b.md's kept part P_K of a
// whole expansion. CCEX accumulates P_K in long double (ScaledComplex, include/calculation/scaled_product.h)
// and writes the mantissa in %La; the except path keeps it and reassembles with the same arithmetic,
// because in double the rounding of each (L/L0) is multiplied by the term's power (1 - connections,
// hundreds to a thousand) and reached 1e-10 against the run (adaptive-regroup-plan.md, 2026-10-11).
// Only the except path uses it; Scaled and approximate_weights (merge, probes) stay in double.
using LongScaled = ::ScaledComplex;
// The except path needs long double to be the x87 80-bit format (64-bit mantissa), as CCEX's ScaledComplex
// has it on the x86 builds. Elsewhere (Apple Silicon: long double is double) the build goes on and only
// --regroup-whole except is refused (adaptive_main.cpp; user 2026-10-11 - a compile-time stop broke every build there).
constexpr bool kLongDoubleX87 = std::numeric_limits<long double>::digits == 64;
using ScaledCurves = std::vector<std::vector<LongScaled> >;
// L0 times the factors (L_S/L0)^power of `weights` (cCCE: no t = 0 normalization), in map order, in
// ScaledComplex arithmetic.
ScaledCurves scaled_product(const std::map<Sites, int>& weights, const std::map<Sites, Curves>& raw);
// `start` times the factors of `weights`, in ScaledComplex arithmetic as CCEX assembles a product
// ((ScaledComplex(L_S) * ScaledComplex(L0).inverse()).power(power)), as values (cCCE). Not bit for bit
// approximate_weights, which is double: within its rounding, with the same non-finite points.
Curves approximate_from(const ScaledCurves& start, const std::map<Sites, int>& weights,
                        const std::map<Sites, Curves>& raw);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// What a set of clusters costs, and site lists as text. (The local reassignment
// and the initial partition that were declared here went with the policies that
// used them on 2026-10-10.)
// ----------------------------------------------------------------------------
namespace adaptive {
// repr() of a tuple of site indices: (), (3,), (0, 1). The policy prints this,
// and a run's log is part of what a recording holds.
std::string sites_repr(const Sites& sites);

// The D^3 * states * steps proxy, counting repeated zero-order calculations.
long long cluster_work(const Plan& plan, const std::vector<Sites>& clusters, long long states,
                    long long steps, long long calls = 1);

// Site lists as they are written to a request or a report, and read back.
Value sites_value(const Sites& sites);                        // [0, 3, 4]
Value sites_list(const std::vector<Sites>& listed);           // [[0, 3], [4]]
Sites sites_of(const Value& listed);
}   // namespace adaptive
