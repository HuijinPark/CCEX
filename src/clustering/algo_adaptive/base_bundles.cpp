/**
 * @file base_bundles.cpp
 * @brief Adaptive refinement: bundles and curves as the policy sees them.
 */
// What the policy uses about bundles and curves: the plan CCEX reports, scaled
// products, curve comparison, and the cluster expansion and its cost.
#include "../../../include/clustering/bundle/bundle_powers.h"
#include "base_bundles.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include "base_numeric.h"

// ----------------------------------------------------------------------------
// products kept as mantissa and exponent
// ----------------------------------------------------------------------------
namespace adaptive {
bool finite(const std::complex<double>& z) {
    return std::isfinite(z.real()) && std::isfinite(z.imag());
}

// The exponent is exact and int64. scaled.h says a wrap here returns a silently
// wrong curve, so every step that can leave the range is judged instead. A real
// run keeps |exponent| under 10 (merge-grow, measured), so this never fires on
// anything reachable - it is here because the wrong answer would be quiet.
long long checked_add(long long a, long long b) {
    long long out = 0;
    if (__builtin_add_overflow(a, b, &out))
        throw OverflowError("integer overflow");
    return out;
}

long long checked_negate(long long a) {
    long long out = 0;
    if (__builtin_sub_overflow(static_cast<long long>(0), a, &out))
        throw OverflowError("integer overflow");
    return out;
}

Scaled normalize(std::complex<double> z, long long exponent) {
    // Python's max(abs(z.real), abs(z.imag)), not symmetric when a NaN is in it (py_max): max(nan, 1.0) is nan but
    // max(1.0, nan) is 1.0, so a NaN imaginary part beside a finite real does NOT stop normalize from scaling.
    const double scale = py_max(std::fabs(z.real()), std::fabs(z.imag()));
    if (scale == 0. || !std::isfinite(scale))
        return Scaled{z, 0};
    int shift = 0;
    std::frexp(scale, &shift);
    return Scaled{std::complex<double>(std::ldexp(z.real(), -shift), std::ldexp(z.imag(), -shift)),
                  checked_add(exponent, shift)};
}

Scaled multiply(const Scaled& a, const Scaled& b) {
    const double real = a.z.real()*b.z.real()-a.z.imag()*b.z.imag();
    const double imag = a.z.real()*b.z.imag()+a.z.imag()*b.z.real();
    return normalize(std::complex<double>(real, imag),
                     checked_add(a.exponent, b.exponent));
}

Scaled inverse(const Scaled& a) {
    if ((a.z.real() == 0. && a.z.imag() == 0.) || !finite(a.z))
        return Scaled{std::complex<double>(std::nan(""), 0.), 0};
    // z.conjugate()/(re*re + im*im): a complex over a float, which CPython widens
    // to a complex and hands to _Py_c_quot (ported in 0B-2).
    const double norm = a.z.real()*a.z.real()+a.z.imag()*a.z.imag();
    double real = 0., imag = 0.;
    complex_div(a.z.real(), -a.z.imag(), norm, 0., &real, &imag);
    return normalize(std::complex<double>(real, imag), checked_negate(a.exponent));
}

Scaled power(Scaled a, long long exponent) {
    long long n = exponent;
    if (n < 0) {
        a = inverse(a);
        n = checked_negate(n);
    }
    Scaled out{std::complex<double>(1., 0.), 0};
    while (n != 0) {
        if ((n & 1) != 0)
            out = multiply(out, a);
        n /= 2;
        if (n != 0)
            a = multiply(a, a);
    }
    return out;
}

std::complex<double> value(const Scaled& a) {
    // CPython's math.ldexp clamps the exponent to an int first, then reports
    // overflow; adaptive_refine.value() catches that and returns infinity with the
    // component's sign - which is the number ldexp already produced.
    const int shift = a.exponent > static_cast<long long>(INT_MAX) ? INT_MAX
                    : a.exponent < static_cast<long long>(INT_MIN) ? INT_MIN
                    : static_cast<int>(a.exponent);
    return std::complex<double>(std::ldexp(a.z.real(), shift), std::ldexp(a.z.imag(), shift));
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// the cluster expansion over bundles, and what it costs
// ----------------------------------------------------------------------------
namespace adaptive {
// The expansion is clustering/bundle/bundle_powers.cpp, shared with the calculation half (2026-10-08; integer and set
// code only, so the flags it is compiled with do not matter); a failure becomes the exception Python raised.
std::map<std::vector<int>, int> bundle_powers(int n, const std::vector<Edge>& edges, int order,
                                             const ExpansionLimits& limits,
                                             const std::vector<std::vector<int> >& extras) {
    using bundle_expansion::Failure;
    try {
        return bundle_expansion::powers(n, edges, order, {limits.max_clusters, limits.max_work}, extras);
    } catch (const bundle_expansion::Error& e) {
        switch (e.kind) {
            case Failure::MaxClusters:
                throw ExpansionLimit("expansion_max_clusters exceeded before coherence calculation");
            case Failure::MaxWork:
                throw ExpansionLimit("expansion_max_work exceeded before coherence calculation");
            case Failure::BadOrder:
                throw ValueError("adaptive requires positive integer bundle order");
            case Failure::BadEdge:
            case Failure::DuplicateEdge:
                throw ValueError("invalid or duplicate adaptive bundle edge");
            case Failure::BadExtra:
            case Failure::ExtraOutside:
                throw ValueError("invalid selective bundle cluster");
            case Failure::ExtraDisconnected:
                throw ValueError("selective cluster must be connected in the retained graph");
            case Failure::PowerTooLarge:
                throw ValueError("adaptive inclusion-exclusion power exceeds 32-bit storage");
            case Failure::PowerOverflow:
                throw ValueError("adaptive inclusion-exclusion intermediate overflow");
        }
        throw ValueError("adaptive expansion failed");
    }
}

Sites union_of(const std::vector<Sites>& bundles, const std::vector<int>& indices) {
    Sites out;
    for (int index : indices)
        out.insert(out.end(), bundles[index].begin(), bundles[index].end());
    std::sort(out.begin(), out.end());
    return out;
}

long long saturating_mul(long long a, long long b) {
    long long out;
    return __builtin_mul_overflow(a, b, &out) ? std::numeric_limits<long long>::max() : out;
}

long long saturating_cube(long long d) {
    return saturating_mul(saturating_mul(d, d), d);
}

long long bundle_dimension(const Sites& bundle, const std::vector<long long>& site_dimensions) {
    long long out = 1;                      // math.prod of an empty bundle is 1
    for (int site : bundle) {
        // A list subscript, not a raw array read: Python counts a negative index
        // from the end and raises IndexError past either end. Reading out of
        // bounds here was undefined behaviour on a site number the plan does not
        // have (found 2026-09-18 by an out-of-range case in the change family).
        const long long at = site < 0 ? static_cast<long long>(site_dimensions.size())+site
                                      : site;
        if (at < 0 || at >= static_cast<long long>(site_dimensions.size()))
            throw IndexError("list index out of range");
        // Saturated, not wrapped: a whole-bath union is asked for its dimension
        // only to be compared against a budget, and a wrap would turn "too big"
        // into a small number the budget waves through.
        if (__builtin_mul_overflow(out, site_dimensions[static_cast<std::size_t>(at)], &out))
            return std::numeric_limits<long long>::max();
    }
    return out;
}

std::map<Sites, int> terms(const std::vector<Sites>& bundles, const std::vector<Edge>& edges,
                           int order, const ExpansionLimits& limits) {
    std::map<Sites, int> out;
    const std::map<std::vector<int>, int> powers =
        bundle_powers(static_cast<int>(bundles.size()), edges, order, limits);
    for (const std::pair<const std::vector<int>, int>& entry : powers)
        out[union_of(bundles, entry.first)] = entry.second;
    return out;
}

std::map<std::vector<int>, double> next_order_patches(const Plan& plan, int order,
                                                      long long neighbors, long long width) {
    const int n = static_cast<int>(plan.bundles.size());
    std::vector<std::map<int, double> > adj(static_cast<std::size_t>(n));
    for (const WeightedEdge& edge : plan.weights()) {
        adj[static_cast<std::size_t>(edge.a)][edge.b] = edge.weight;
        adj[static_cast<std::size_t>(edge.b)][edge.a] = edge.weight;
    }
    std::map<std::vector<int>, double> beam;
    for (int i = 0; i < n; ++i)
        beam[std::vector<int>{i}] = 0.;
    for (int size = 2; size <= std::min(order+1, n); ++size) {
        // Later writers win: Python assigns into the dict rather than taking a
        // maximum, and walks it in sorted key order, so which one is last is fixed.
        std::map<std::vector<int>, double> candidates;
        for (const std::pair<const std::vector<int>, double>& entry : beam) {
            for (int a : entry.first) {
                std::vector<std::pair<int, double> > near(
                    adj[static_cast<std::size_t>(a)].begin(),
                    adj[static_cast<std::size_t>(a)].end());
                std::stable_sort(near.begin(), near.end(),
                                 [](const std::pair<int, double>& x,
                                    const std::pair<int, double>& y) {
                                     if (x.second != y.second) return x.second > y.second;
                                     return x.first < y.first;
                                 });
                if (static_cast<long long>(near.size()) > neighbors)
                    near.resize(static_cast<std::size_t>(neighbors < 0 ? 0 : neighbors));
                for (const std::pair<int, double>& step : near) {
                    const int b = step.first;
                    if (std::find(entry.first.begin(), entry.first.end(), b) !=
                        entry.first.end())
                        continue;
                    std::vector<int> next = entry.first;
                    next.push_back(b);
                    std::sort(next.begin(), next.end());
                    std::vector<double> parts;
                    for (int i : entry.first) {
                        const std::map<int, double>::const_iterator at =
                            adj[static_cast<std::size_t>(b)].find(i);
                        parts.push_back(at == adj[static_cast<std::size_t>(b)].end() ? 0. : at->second);
                    }
                    candidates[next] = entry.second+fsum(parts);
                }
            }
        }
        std::vector<std::vector<int> > ranked;
        for (const std::pair<const std::vector<int>, double>& entry : candidates)
            ranked.push_back(entry.first);
        std::stable_sort(ranked.begin(), ranked.end(),
                         [&candidates](const std::vector<int>& x, const std::vector<int>& y) {
                             if (candidates.at(x) != candidates.at(y))
                                 return candidates.at(x) > candidates.at(y);
                             return x < y;
                         });
        if (static_cast<long long>(ranked.size()) > width)
            ranked.resize(static_cast<std::size_t>(width < 0 ? 0 : width));
        beam.clear();
        for (const std::vector<int>& key : ranked)
            beam[key] = candidates.at(key);
        if (beam.empty())
            break;
    }
    std::map<std::vector<int>, double> out;
    for (const std::pair<const std::vector<int>, double>& entry : beam) {
        if (static_cast<int>(entry.first.size()) == order+1)
            out.insert(entry);
    }
    return out;
}

}   // namespace adaptive

// ----------------------------------------------------------------------------
// comparing curves
// ----------------------------------------------------------------------------
namespace adaptive {
std::complex<double> curve_at(const Curves& curves, std::size_t state, std::size_t step) {
    const std::pair<double, double>& z = curves[state][step];
    return std::complex<double>(z.first, z.second);
}

Value discrepancy_value(const Discrepancy& found) {
    Value out = Value::object();
    out.set("scaled_error", Value::number(found.scaled_error));
    out.set("max_abs_error", Value::number(found.max_abs_error));
    if (found.has_location) {
        Value where = Value::array();
        where.push(Value::integer(static_cast<long long>(found.state)));
        where.push(Value::integer(static_cast<long long>(found.step)));
        out.set("worst_state_step", where);
    } else {
        out.set("worst_state_step", Value());
    }
    return out;
}

std::pair<double, double> point_error(const std::complex<double>& x,
                                      const std::complex<double>& y, double atol, double rtol) {
    const bool valid = finite(x) && finite(y);
    // abs() of a complex is hypot, which is libm: this is one of the places the
    // verified-numeric-environment contract applies (0B-2).
    const double d = valid ? complex_abs(x.real()-y.real(), x.imag()-y.imag()) : HUGE_VAL;
    const double scale = valid ? atol+rtol*complex_abs(y.real(), y.imag()) : 1.;
    return {d, scale != 0. ? d/scale : (d == 0. ? 0. : HUGE_VAL)};
}

Discrepancy discrepancy(const Curves& actual, const Curves& reference, double atol, double rtol) {
    bool shaped = !actual.empty() && actual.size() == reference.size();
    for (std::size_t i = 0; shaped && i < actual.size(); ++i)
        shaped = !actual[i].empty() && actual[i].size() == reference[i].size();
    if (!shaped)
        throw ValueError("coherence comparison requires matching nonempty state/time grids");
    Discrepancy out;
    for (std::size_t state = 0; state < actual.size(); ++state) {
        for (std::size_t step = 0; step < actual[state].size(); ++step) {
            const std::pair<double, double> told = point_error(
                curve_at(actual, state, step), curve_at(reference, state, step), atol, rtol);
            const double d = told.first;
            const double error = told.second;
            if (error > out.scaled_error) {
                out.scaled_error = error;
                out.has_location = true;
                out.state = static_cast<int>(state);
                out.step = static_cast<int>(step);
            }
            if (d > out.max_abs_error)          // Python's max(absolute, d)
                out.max_abs_error = d;
        }
    }
    return out;
}

Curves approximate_weights(const std::map<Sites, int>& weights,
                           const std::map<Sites, Curves>& raw, bool gcce) {
    const std::map<Sites, Curves>::const_iterator zero = raw.find(Sites());
    if (zero == raw.end())
        throw KeyError(sites_repr(Sites()));      // Python's KeyError carries the tuple's repr: (), (3,), (0, 1)
    Curves result;
    for (std::size_t state = 0; state < zero->second.size(); ++state) {
        std::vector<Scaled> accumulated;
        for (std::size_t step = 0; step < zero->second[state].size(); ++step) {
            const Scaled z0 = normalize(curve_at(zero->second, state, step));
            Scaled p = z0;
            // std::map is already in the order sorted(weights.items()) gives.
            for (const std::pair<const Sites, int>& term : weights) {
                const std::map<Sites, Curves>::const_iterator found = raw.find(term.first);
                if (found == raw.end())
                    throw KeyError(sites_repr(term.first));
                const Scaled q = multiply(normalize(curve_at(found->second, state, step)), inverse(z0));
                p = multiply(p, power(q, static_cast<long long>(term.second)));
            }
            accumulated.push_back(p);
        }
        const Scaled norm = gcce ? inverse(accumulated.at(0))
                                 : Scaled{std::complex<double>(1., 0.), 0};
        Curve curve;
        for (const Scaled& p : accumulated) {
            const std::complex<double> z = value(multiply(p, norm));
            curve.emplace_back(z.real(), z.imag());
        }
        result.push_back(curve);
    }
    return result;
}

Curves approximate(const std::vector<Sites>& bundles, const std::vector<Edge>& edges, int order,
                   const std::map<Sites, Curves>& raw, bool gcce, const ExpansionLimits& limits) {
    return approximate_weights(terms(bundles, edges, order, limits), raw, gcce);
}

namespace {
// p times (L_S / L0)^power for every term, in map order, as CCEX's assembleFactor forms a factor.
LongScaled multiply_terms(LongScaled p, const std::map<Sites, int>& weights, const std::map<Sites, Curves>& raw,
                          const LongScaled& inverse_zero, std::size_t state, std::size_t step) {
    for (const std::pair<const Sites, int>& term : weights) {
        const std::map<Sites, Curves>::const_iterator found = raw.find(term.first);
        if (found == raw.end())
            throw KeyError(sites_repr(term.first));
        p = p*(LongScaled(curve_at(found->second, state, step))*inverse_zero).power(term.second);
    }
    return p;
}
}   // namespace

ScaledCurves scaled_product(const std::map<Sites, int>& weights, const std::map<Sites, Curves>& raw) {
    const std::map<Sites, Curves>::const_iterator zero = raw.find(Sites());
    if (zero == raw.end())
        throw KeyError(sites_repr(Sites()));
    ScaledCurves out;
    for (std::size_t state = 0; state < zero->second.size(); ++state) {
        out.emplace_back();
        for (std::size_t step = 0; step < zero->second[state].size(); ++step) {
            const LongScaled z0(curve_at(zero->second, state, step));
            out.back().push_back(multiply_terms(z0, weights, raw, z0.inverse(), state, step));
        }
    }
    return out;
}

Curves approximate_from(const ScaledCurves& start, const std::map<Sites, int>& weights,
                        const std::map<Sites, Curves>& raw) {
    const std::map<Sites, Curves>::const_iterator zero = raw.find(Sites());
    if (zero == raw.end())
        throw KeyError(sites_repr(Sites()));
    if (start.size() != zero->second.size())
        throw ValueError("the starting product and L0 have different state grids");
    Curves result;
    for (std::size_t state = 0; state < start.size(); ++state) {
        if (start[state].size() != zero->second[state].size())
            throw ValueError("the starting product and L0 have different time grids");
        Curve curve;
        for (std::size_t step = 0; step < start[state].size(); ++step) {
            const LongScaled inverse_zero = LongScaled(curve_at(zero->second, state, step)).inverse();
            const std::complex<double> z =
                multiply_terms(start[state][step], weights, raw, inverse_zero, state, step).value();
            curve.emplace_back(z.real(), z.imag());
        }
        result.push_back(curve);
    }
    return result;
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// the work of a set of clusters, and site lists as text
// (the local reassignment and the initial partition that were here went with the
// policies that used them on 2026-10-10)
// ----------------------------------------------------------------------------
namespace adaptive {
std::string sites_repr(const Sites& sites) {
    std::string out = "(";
    for (std::size_t i = 0; i < sites.size(); ++i)
        out += (i ? ", " : "")+std::to_string(sites[i]);
    if (sites.size() == 1)
        out += ",";
    return out+")";
}

// The most this proxy ever reports: past it the answer is only "more than any
// budget", and the exact figure has no reader.
const long long kWorkCeiling = std::numeric_limits<long long>::max();

long long cluster_work(const Plan& plan, const std::vector<Sites>& clusters, long long states,
                    long long steps, long long calls) {
    const long long q = plan.qubit_dimension;
    // set(map(tuple, clusters)): a cluster asked for twice is calculated once.
    const std::set<Sites> unique(clusters.begin(), clusters.end());
    // Built in 128 bits and saturated rather than wrapped: this number is the
    // one a budget is compared against, and a wrap would turn "too expensive"
    // into a small - or negative - value the budget waves through.
    __int128 total = 0;
    for (const Sites& bundle : unique) {
        const __int128 d = static_cast<__int128>(q)*bundle_dimension(bundle, plan.site_dimensions);
        // 2^42 cubed is the last cube 128 bits hold, and already past the ceiling.
        if (d > (static_cast<__int128>(1) << 42)) return kWorkCeiling;
        total += d*d*d;
        if (total > kWorkCeiling) return kWorkCeiling;
    }
    total += static_cast<__int128>(calls)*q*q*q;
    total *= static_cast<__int128>(states);
    total *= static_cast<__int128>(steps);
    return total > kWorkCeiling ? kWorkCeiling : static_cast<long long>(total);
}

Value sites_value(const Sites& sites) {
    Value out = Value::array();
    for (int site : sites)
        out.push(Value::integer(static_cast<long long>(site)));
    return out;
}

Value sites_list(const std::vector<Sites>& listed) {
    Value out = Value::array();
    for (const Sites& sites : listed)
        out.push(sites_value(sites));
    return out;
}

Sites sites_of(const Value& listed) {
    Sites out;
    for (std::size_t i = 0; i < listed.size(); ++i)
        out.push_back(static_cast<int>(listed.at(i).as_int()));
    return out;
}
}   // namespace adaptive
