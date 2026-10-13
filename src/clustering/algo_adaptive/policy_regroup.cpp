// Method D of (internal notes) (the merge policy's --regroup):
// regroup the bundles where the bundle expansion gives |L| > 1. The item numbers in the
// comments are the plan's "방식 D (구현 명세)".
#include "policy_regroup.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <dirent.h>
#include <fnmatch.h>
#include <limits>
#include <sstream>
#include "base_numeric.h"
#include "policy_merge.h"

namespace adaptive {
namespace regroup {
namespace {
double magnitude(const std::pair<double, double>& z) {
    // numpy's abs of a complex, which is hypot and never raises: a finite value too
    // large for a double is +inf here, which every limit below treats as over.
    return std::hypot(z.first, z.second);
}

bool finite_point(const std::pair<double, double>& z) {
    return std::isfinite(z.first) && std::isfinite(z.second);
}

Value sites_text(const Sites& sites) {
    return sites_value(sites);
}

Value edges_value(const std::vector<Edge>& edges) {
    Value out = Value::array();
    for (const Edge& edge : edges) {
        Value one = Value::array();
        one.push(Value::integer(static_cast<long long>(edge.first)));
        one.push(Value::integer(static_cast<long long>(edge.second)));
        out.push(one);
    }
    return out;
}

Value curves_value(const Curves& curves) {
    Value out = Value::array();
    for (const Curve& curve : curves) {
        Value row = Value::array();
        for (const std::pair<double, double>& z : curve) {
            Value point = Value::array();
            point.push(trace_number(z.first));
            point.push(trace_number(z.second));
            row.push(point);
        }
        out.push(row);
    }
    return out;
}

Partition pick(const Partition& units, const std::vector<int>& indices) {
    Partition out;
    for (int i : indices)
        out.push_back(units[static_cast<std::size_t>(i)]);
    return out;
}

Sites all_sites(const Partition& bundles) {
    Sites out;
    for (const Sites& bundle : bundles)
        out.insert(out.end(), bundle.begin(), bundle.end());
    std::sort(out.begin(), out.end());
    return out;
}

const Curves& held(const std::map<Sites, Curves>& raw, const Sites& sites) {
    const std::map<Sites, Curves>::const_iterator found = raw.find(sites);
    if (found == raw.end())
        throw KeyError(sites_repr(sites));
    return found->second;
}

// The nearest of `others` to `sites`: smallest closeness, then the smaller site list.
bool nearer(const Geometry& geometry, const Sites& sites, const Sites& x, const Sites& y) {
    const double a = closeness(geometry, sites, x), b = closeness(geometry, sites, y);
    if (a != b)
        return a < b;
    return x < y;
}
}   // namespace

Value trace_number(double x) {
    if (std::isnan(x))
        return Value::text("nan");
    if (std::isinf(x))
        return Value::text(x > 0 ? "inf" : "-inf");
    return Value::number(x);
}

Value partition_value(const Partition& units) {
    return sites_list(units);
}

Partition canonical(const Partition& bundles) {
    Partition out;
    std::set<int> seen;
    for (const Sites& bundle : bundles) {
        if (bundle.empty())
            throw ValueError("a partition must hold nonempty disjoint bundles");
        Sites sorted_bundle = bundle;
        std::sort(sorted_bundle.begin(), sorted_bundle.end());
        for (int site : sorted_bundle)
            if (!seen.insert(site).second)
                throw ValueError("a partition must hold nonempty disjoint bundles");
        out.push_back(sorted_bundle);
    }
    std::sort(out.begin(), out.end());
    return out;
}

// ---------------------------------------------------------------------------
// sites
// ---------------------------------------------------------------------------
Geometry geometry_of(const Value& physics) {
    const Value& listed = physics.at("sites");
    Geometry out(listed.size());
    const auto position = [](const Value& xyz) {
        std::array<double, 3> p{};
        for (std::size_t k = 0; k < 3; ++k)
            p[k] = number_of(xyz.at(k));
        return p;
    };
    for (std::size_t i = 0; i < listed.size(); ++i) {
        const Value& site = listed.at(i);
        const long long id = site.at("id").as_int();
        if (id < 0 || id >= static_cast<long long>(listed.size()))
            throw ValueError("physics site IDs must be 0..n-1");
        Site& here = out[static_cast<std::size_t>(id)];
        here.xyz = position(site.at("xyz"));
        const Value& components = site.at("components");
        for (std::size_t c = 0; c < components.size(); ++c)
            here.components.push_back(position(components.at(c).at("xyz")));
        if (here.components.empty())
            here.components.push_back(here.xyz);
    }
    return out;
}

double distance(const Geometry& sites, int a, int b) {
    const Site& x = sites[static_cast<std::size_t>(a)];
    const Site& y = sites[static_cast<std::size_t>(b)];
    return std::hypot(x.xyz[0]-y.xyz[0], x.xyz[1]-y.xyz[1], x.xyz[2]-y.xyz[2]);
}

double closeness(const Geometry& sites, const Sites& a, const Sites& b) {
    double out = HUGE_VAL;
    for (int x : a)
        for (int y : b)
            out = std::min(out, distance(sites, std::min(x, y), std::max(x, y)));
    return out;
}

Adjacency distance_adjacency(const Geometry& sites, float rdip, float rdipcut) {
    const int n = static_cast<int>(sites.size());
    Adjacency out(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        for (int j = i+1; j < n; ++j) {
            bool linked = false;
            for (const std::array<double, 3>& a : sites[static_cast<std::size_t>(i)].components)
                for (const std::array<double, 3>& b : sites[static_cast<std::size_t>(j)].components) {
                    // bundle_cluster.cpp: r = hypot(a - b), connected when r > rdipcut && r < rdip.
                    const double r = std::hypot(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
                    if (!std::isfinite(r) || r <= 0.)
                        throw ValueError("distinct bath sites have coincident or invalid positions");
                    if (r > rdipcut && r < rdip)
                        linked = true;
                }
            if (linked) {
                out[static_cast<std::size_t>(i)].push_back(j);
                out[static_cast<std::size_t>(j)].push_back(i);
            }
        }
    return out;
}

Adjacency complete_adjacency(int sites) {
    Adjacency out(static_cast<std::size_t>(sites));
    for (int i = 0; i < sites; ++i)
        for (int j = 0; j < sites; ++j)
            if (i != j)
                out[static_cast<std::size_t>(i)].push_back(j);
    return out;
}

std::vector<Edge> connections(const Partition& bundles, const Adjacency& links) {
    std::vector<int> owner(links.size(), -1);
    for (std::size_t g = 0; g < bundles.size(); ++g)
        for (int site : bundles[g]) {
            if (site < 0 || site >= static_cast<int>(links.size()))
                throw IndexError("site outside the bath");
            owner[static_cast<std::size_t>(site)] = static_cast<int>(g);
        }
    std::set<Edge> out;
    for (std::size_t g = 0; g < bundles.size(); ++g)
        for (int site : bundles[g])
            for (int other : links[static_cast<std::size_t>(site)]) {
                const int h = owner[static_cast<std::size_t>(other)];
                if (h >= 0 && h != static_cast<int>(g))
                    out.insert(Edge(std::min(static_cast<int>(g), h), std::max(static_cast<int>(g), h)));
            }
    return std::vector<Edge>(out.begin(), out.end());
}

// ---------------------------------------------------------------------------
// curves
// ---------------------------------------------------------------------------
bool physical(const Curves& curves, double limit) {
    for (const Curve& curve : curves)
        for (const std::pair<double, double>& z : curve)
            if (!finite_point(z) || !(magnitude(z) <= limit))
                return false;
    return true;
}

Point worst_point(const Curves& curves) {
    Point out;
    for (std::size_t s = 0; s < curves.size(); ++s)
        for (std::size_t t = 0; t < curves[s].size(); ++t)
            if (!finite_point(curves[s][t])) {
                out.state = static_cast<int>(s);
                out.step = static_cast<int>(t);
                out.nonfinite = true;
                return out;
            }
    double largest = -1.;
    for (std::size_t s = 0; s < curves.size(); ++s)
        for (std::size_t t = 0; t < curves[s].size(); ++t) {
            const double size = magnitude(curves[s][t]);
            if (size > largest) {
                largest = size;
                out.state = static_cast<int>(s);
                out.step = static_cast<int>(t);
            }
        }
    return out;
}

long long unphysical_states(const Curves& curves, double limit) {
    long long out = 0;
    for (const Curve& curve : curves) {
        bool bad = false;
        for (const std::pair<double, double>& z : curve)
            if (!finite_point(z) || !(magnitude(z) <= limit))
                bad = true;
        out += bad ? 1 : 0;
    }
    return out;
}

double excess(const Curve& curve) {
    double largest = 0.;
    for (const std::pair<double, double>& z : curve) {
        if (!finite_point(z))
            return HUGE_VAL;
        largest = std::max(largest, magnitude(z));
    }
    return std::max(largest-1., 0.);
}

// ---------------------------------------------------------------------------
// the culprit
// ---------------------------------------------------------------------------
Factor culprit_factor(std::complex<double> pair, std::complex<double> a, std::complex<double> b) {
    Factor out;
    const auto zero = [](const std::complex<double>& z) { return z.real() == 0. && z.imag() == 0.; };
    if (!finite(pair) || !finite(a) || !finite(b)) {
        out.g = HUGE_VAL;
        out.cause = "nonfinite_probe";
        return out;
    }
    if (zero(a) || zero(b)) {
        // L_ab / 0: the irreducible factor is not defined, which item 3 counts as +inf.
        out.g = HUGE_VAL;
        out.cause = "zero_bundle";
        return out;
    }
    if (zero(pair)) {
        out.g = -HUGE_VAL;
        out.cause = "zero_pair";
        return out;
    }
    // Scaled, as CCEX assembles its expansion: L_a L_b may underflow where the factor
    // itself is an ordinary number (review 10). Divided rather than multiplied by an
    // inverse, so a pair that is exactly the product of its bundles gives exactly 1 and
    // g = 0 - the reciprocal's rounding made such a pair a culprit of g = 2e-16.
    const Scaled pair_scaled = normalize(pair);
    const Scaled both = multiply(normalize(a), normalize(b));
    double re = 0., im = 0.;
    complex_div(pair_scaled.z.real(), pair_scaled.z.imag(), both.z.real(), both.z.imag(), &re, &im);
    const Scaled f = normalize(std::complex<double>(re, im), pair_scaled.exponent-both.exponent);
    const std::complex<double> v = value(f);
    const double size = std::hypot(v.real(), v.imag());
    // A finite v can still have an infinite modulus (both parts near DBL_MAX): that
    // case takes the scaled path too, or distinct large factors tie at +inf.
    if (finite(v) && std::isfinite(size) && size >= std::numeric_limits<double>::min()) {
        out.g = std::log(size);
        return out;
    }
    // Past a double's range the value is not representable but its log is: the
    // mantissa's log and the exponent's.
    const double mantissa = std::hypot(f.z.real(), f.z.imag());
    if (!(mantissa > 0.) || !std::isfinite(mantissa)) {
        out.g = HUGE_VAL;
        out.cause = "nonfinite_factor";
        return out;
    }
    out.g = std::log(mantissa)+static_cast<double>(f.exponent)*0.69314718055994530942;
    return out;
}

std::vector<Culprit> culprits(const Partition& bundles, const std::vector<Edge>& edges, const PointValues& at) {
    const auto value_of = [&at](const Sites& sites) {
        const PointValues::const_iterator found = at.find(sites);
        if (found == at.end())
            throw KeyError(sites_repr(sites));
        return found->second;
    };
    std::vector<std::pair<Sites, Culprit> > ranked;
    for (const Edge& edge : edges) {
        const Sites& a = bundles[static_cast<std::size_t>(edge.first)];
        const Sites& b = bundles[static_cast<std::size_t>(edge.second)];
        Sites both = a;
        both.insert(both.end(), b.begin(), b.end());
        std::sort(both.begin(), both.end());
        Culprit one;
        one.a = edge.first;
        one.b = edge.second;
        one.factor = culprit_factor(value_of(both), value_of(a), value_of(b));
        if (one.factor.g > 0.)
            ranked.emplace_back(both, one);
    }
    std::stable_sort(ranked.begin(), ranked.end(),
                     [](const std::pair<Sites, Culprit>& x, const std::pair<Sites, Culprit>& y) {
                         if (x.second.factor.g != y.second.factor.g)
                             return x.second.factor.g > y.second.factor.g;
                         if (x.first != y.first)
                             return x.first < y.first;
                         return std::make_pair(x.second.a, x.second.b) < std::make_pair(y.second.a, y.second.b);
                     });
    std::vector<Culprit> out;
    for (const std::pair<Sites, Culprit>& one : ranked)
        out.push_back(one.second);
    return out;
}

std::vector<Culprit> culprits(const Partition& bundles, const std::vector<Edge>& edges,
                              const std::map<Sites, Curves>& raw, const Point& at) {
    return culprits(bundles, edges, point_values_of(bundles, edges, raw, at));
}

std::vector<Sites> factor_clusters(const Partition& bundles, const std::vector<Edge>& edges) {
    std::set<Sites> out;
    for (const Edge& edge : edges) {
        out.insert(bundles[static_cast<std::size_t>(edge.first)]);
        out.insert(bundles[static_cast<std::size_t>(edge.second)]);
        out.insert(union_of(bundles, {edge.first, edge.second}));
    }
    return std::vector<Sites>(out.begin(), out.end());
}

PointValues point_values_of(const Partition& bundles, const std::vector<Edge>& edges,
                            const std::map<Sites, Curves>& raw, const Point& at) {
    const std::size_t s = static_cast<std::size_t>(at.state), t = static_cast<std::size_t>(at.step);
    PointValues out;
    for (const Sites& sites : factor_clusters(bundles, edges))
        out[sites] = curve_at(held(raw, sites), s, t);
    return out;
}

PointValues read_factors(const std::string& directory, long long ranks, const Partition& bundles,
                         const std::vector<Edge>& edges, int state, int step) {
    // the rows a factors step owes: (a, -1) for a bundle with a connection, (a, b) for a pair
    std::map<std::pair<int, int>, Sites> owed;
    for (const Edge& edge : edges) {
        owed[{edge.first, -1}] = bundles[static_cast<std::size_t>(edge.first)];
        owed[{edge.second, -1}] = bundles[static_cast<std::size_t>(edge.second)];
        owed[{edge.first, edge.second}] = union_of(bundles, {edge.first, edge.second});
    }
    std::vector<std::string> files;
    if (DIR* listing = opendir(directory.c_str())) {
        while (const dirent* entry = readdir(listing))
            if (fnmatch("response.rank-*.factors.tsv", entry->d_name, 0) == 0)
                files.push_back(directory+"/"+entry->d_name);
        closedir(listing);
    }
    std::sort(files.begin(), files.end());
    if (static_cast<long long>(files.size()) != ranks)
        throw RuntimeError("missing MPI factors response");
    const std::string header = "# state="+std::to_string(state)+" step="+std::to_string(step)+" bundles="+
                               std::to_string(bundles.size())+" edges="+std::to_string(edges.size());
    PointValues out;
    std::set<std::pair<int, int> > seen;
    for (const std::string& path : files) {
        std::istringstream rows(read_file(path));
        std::string line;
        if (!std::getline(rows, line) || line != header)
            throw RuntimeError("factors response "+path+" is not for state "+std::to_string(state)+", step "+
                               std::to_string(step)+" and this partition");
        if (!std::getline(rows, line) || line != "a\tb\treal\timag")
            throw RuntimeError("factors response "+path+" has no column header");
        while (std::getline(rows, line)) {
            if (line.empty())
                continue;
            std::istringstream fields(line);
            std::string a_text, b_text, real_text, imag_text, rest;
            std::getline(fields, a_text, '\t');
            std::getline(fields, b_text, '\t');
            std::getline(fields, real_text, '\t');
            std::getline(fields, imag_text, '\t');
            if (std::getline(fields, rest, '\t'))
                throw RuntimeError("factors response row with more than four fields");
            const std::pair<int, int> key(static_cast<int>(parse_int(a_text)), static_cast<int>(parse_int(b_text)));
            const std::map<std::pair<int, int>, Sites>::const_iterator found = owed.find(key);
            if (found == owed.end())
                throw RuntimeError("factors response row that was not asked for");
            if (!seen.insert(key).second)
                throw RuntimeError("duplicate factors response row");
            out[found->second] = std::complex<double>(parse_float(real_text), parse_float(imag_text));
        }
    }
    if (seen.size() != owed.size())
        throw RuntimeError("incomplete factors response");
    return out;
}

// ---------------------------------------------------------------------------
// region, neighbours, candidates
// ---------------------------------------------------------------------------
std::vector<int> region_of(const Partition& bundles, int a, int b, const Geometry& sites) {
    std::vector<int> out{a, b};
    Sites pair = bundles[static_cast<std::size_t>(a)];
    pair.insert(pair.end(), bundles[static_cast<std::size_t>(b)].begin(),
                bundles[static_cast<std::size_t>(b)].end());
    int best = -1;
    for (int k = 0; k < static_cast<int>(bundles.size()); ++k) {
        if (k == a || k == b)
            continue;
        if (best < 0 || nearer(sites, pair, bundles[static_cast<std::size_t>(k)],
                               bundles[static_cast<std::size_t>(best)]))
            best = k;
    }
    if (best >= 0)
        out.push_back(best);
    return out;
}

std::vector<int> neighbours_of(const Partition& bundles, const std::vector<int>& region,
                               const Geometry& sites, int count) {
    Sites inside;
    for (int i : region)
        inside.insert(inside.end(), bundles[static_cast<std::size_t>(i)].begin(),
                      bundles[static_cast<std::size_t>(i)].end());
    std::vector<int> outside;
    for (int k = 0; k < static_cast<int>(bundles.size()); ++k)
        if (std::find(region.begin(), region.end(), k) == region.end())
            outside.push_back(k);
    std::stable_sort(outside.begin(), outside.end(), [&](int x, int y) {
        return nearer(sites, inside, bundles[static_cast<std::size_t>(x)], bundles[static_cast<std::size_t>(y)]);
    });
    if (static_cast<int>(outside.size()) > count)
        outside.resize(static_cast<std::size_t>(std::max(count, 0)));
    return outside;
}

std::vector<Partition> candidates(const Partition& region, const std::vector<long long>& dimensions,
                                  long long cap) {
    const auto fits = [&](const Sites& bundle) { return bundle_dimension(bundle, dimensions) <= cap; };
    const std::size_t n = region.size();
    std::set<Partition> out;
    const auto add = [&](Partition next) {
        Partition kept;
        for (Sites& bundle : next)
            if (!bundle.empty())
                kept.push_back(bundle);
        for (const Sites& bundle : kept)
            if (!fits(bundle))
                return;
        out.insert(canonical(kept));
    };
    // merge two
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i+1; j < n; ++j) {
            Partition next;
            for (std::size_t k = 0; k < n; ++k)
                if (k != i && k != j)
                    next.push_back(region[k]);
            Sites joined = region[i];
            joined.insert(joined.end(), region[j].begin(), region[j].end());
            next.push_back(joined);
            add(next);
        }
    // split one in two: every bipartition, once (the part holding the bundle's first site)
    for (std::size_t i = 0; i < n; ++i) {
        const Sites& bundle = region[i];
        if (bundle.size() < 2 || bundle.size() > 30)
            continue;
        for (unsigned long mask = 1; mask < (1UL << bundle.size())-1; ++mask) {
            if ((mask & 1UL) == 0)
                continue;
            Sites part, rest;
            for (std::size_t k = 0; k < bundle.size(); ++k)
                ((mask >> k) & 1UL ? part : rest).push_back(bundle[k]);
            Partition next;
            for (std::size_t k = 0; k < n; ++k)
                if (k != i)
                    next.push_back(region[k]);
            next.push_back(part);
            next.push_back(rest);
            add(next);
        }
    }
    // move one site from bundle i to bundle j
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            if (i == j)
                continue;
            for (int site : region[i]) {
                Partition next = region;
                next[i].erase(std::find(next[i].begin(), next[i].end(), site));
                next[j].push_back(site);
                add(next);
            }
        }
    // swap two sites between bundles i and j
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i+1; j < n; ++j)
            for (int s : region[i])
                for (int t : region[j]) {
                    Partition next = region;
                    *std::find(next[i].begin(), next[i].end(), s) = t;
                    *std::find(next[j].begin(), next[j].end(), t) = s;
                    add(next);
                }
    out.erase(canonical(region));
    return std::vector<Partition>(out.begin(), out.end());
}

Partition replace(const Partition& units, const std::vector<int>& region, const Partition& replacement) {
    Partition out;
    for (int k = 0; k < static_cast<int>(units.size()); ++k)
        if (std::find(region.begin(), region.end(), k) == region.end())
            out.push_back(units[static_cast<std::size_t>(k)]);
    out.insert(out.end(), replacement.begin(), replacement.end());
    return canonical(out);
}

// ---------------------------------------------------------------------------
// the local error and the expansion
// ---------------------------------------------------------------------------
LocalError local_error(const Curves& model, const Curves& direct) {
    if (model.empty() || model.size() != direct.size())
        throw ValueError("local error needs matching state grids");
    LocalError out;
    for (const Curve& curve : direct)
        for (const std::pair<double, double>& z : curve)
            if (!finite_point(z)) {
                out.verifiable = false;
                return out;
            }
    for (const Curve& curve : model)
        for (const std::pair<double, double>& z : curve)
            if (!finite_point(z)) {
                out.value = HUGE_VAL;
                return out;
            }
    const std::size_t states = model.size(), steps = model[0].size();
    std::vector<double> differences;
    for (std::size_t t = 0; t < steps; ++t) {
        double mr = 0., mi = 0., dr = 0., di = 0.;
        for (std::size_t s = 0; s < states; ++s) {
            if (model[s].size() != steps || direct[s].size() != steps)
                throw ValueError("local error needs matching time grids");
            mr += model[s][t].first;
            mi += model[s][t].second;
            dr += direct[s][t].first;
            di += direct[s][t].second;
        }
        const double n = static_cast<double>(states);
        differences.push_back(std::hypot(mr/n-dr/n, mi/n-di/n));
    }
    // Exactly summed, so the mean does not depend on the order of the time points.
    out.value = fsum(differences)/static_cast<double>(steps);
    return out;
}

bool improves(double after, double before, double epsilon) {
    return after < before-epsilon;
}

bool neighbours_hold(const std::vector<double>& after, const std::vector<double>& before, double epsilon) {
    if (after.size() != before.size())
        throw ValueError("one error per neighbour, before and after");
    // The prototype's test, any(after > before + eps): inf after inf holds, since
    // inf is not more than inf; a NaN compares false and holds too, but a NaN cannot
    // reach here - a direct solution that is not finite refuses the region first, and
    // a model that is not finite is +inf, not NaN.
    for (std::size_t k = 0; k < after.size(); ++k)
        if (after[k] > before[k]+epsilon)
            return false;
    return true;
}

bool whole_improves(long long unphysical_after, double excess_after, long long unphysical_before,
                    double excess_before) {
    return unphysical_after <= unphysical_before && excess_after < excess_before;
}

// ---------------------------------------------------------------------------
// condition (c) from a kept product
// ---------------------------------------------------------------------------
std::vector<Sites> changed_terms(const std::map<Sites, int>& now,
                                 const std::vector<const std::map<Sites, int>*>& tables) {
    std::vector<Sites> out;
    for (const std::pair<const Sites, int>& term : now)
        for (const std::map<Sites, int>* table : tables) {
            const std::map<Sites, int>::const_iterator found = table->find(term.first);
            if ((found == table->end() ? 0 : found->second) != term.second) {
                out.push_back(term.first);
                break;
            }
        }
    return out;
}

std::map<Sites, int> new_terms(const std::map<Sites, int>& now, const std::vector<Sites>& except,
                               const std::map<Sites, int>& table) {
    const std::set<Sites> excepted(except.begin(), except.end());
    // L_c = P_K x prod over N_c holds only when every kept term keeps its power
    for (const std::pair<const Sites, int>& term : now) {
        if (excepted.count(term.first) != 0)
            continue;
        const std::map<Sites, int>::const_iterator found = table.find(term.first);
        if ((found == table.end() ? 0 : found->second) != term.second)
            throw RuntimeError("regroup: a kept term "+sites_repr(term.first)+" changes its power in a candidate");
    }
    std::map<Sites, int> out;
    for (const std::pair<const Sites, int>& term : table)
        if (term.second != 0 && (now.count(term.first) == 0 || excepted.count(term.first) != 0))
            out.insert(term);
    return out;
}

ScaledCurves read_except(const std::string& directory, const std::vector<int>& states, long long steps) {
    ScaledCurves out;
    for (int state : states) {
        const std::string path = directory+"/response.except.state-"+std::to_string(state)+".tsv";
        std::istringstream rows(read_file(path));
        std::string line;
        if (!std::getline(rows, line) || line != "step\treal\timag\texponent")
            throw RuntimeError("except response "+path+" has no header");
        std::vector<LongScaled> curve(static_cast<std::size_t>(steps));
        std::vector<bool> seen(static_cast<std::size_t>(steps), false);
        const auto mantissa = [](const std::string& text) {
            // %La: a hexadecimal long double, or inf/nan
            char* end = nullptr;
            const long double x = std::strtold(text.c_str(), &end);
            if (text.empty() || end != text.c_str()+text.size())
                throw RuntimeError("except response mantissa is not a number: "+text);
            return x;
        };
        while (std::getline(rows, line)) {
            if (line.empty())
                continue;
            std::istringstream fields(line);
            std::string step_text, real_text, imag_text, exponent_text, rest;
            std::getline(fields, step_text, '\t');
            std::getline(fields, real_text, '\t');
            std::getline(fields, imag_text, '\t');
            std::getline(fields, exponent_text, '\t');
            if (std::getline(fields, rest, '\t'))
                throw RuntimeError("except response row with more than four fields");
            const long long step = parse_int(step_text);
            if (step < 0 || step >= steps)
                throw RuntimeError("except response step outside the grid");
            if (seen[static_cast<std::size_t>(step)])
                throw RuntimeError("duplicate except response step");
            seen[static_cast<std::size_t>(step)] = true;
            // the long double mantissa as written, its exponent kept apart; normalized as ScaledComplex
            // does (a non-finite mantissa resets the exponent)
            LongScaled& p = curve[static_cast<std::size_t>(step)];
            p.real = mantissa(real_text);
            p.imag = mantissa(imag_text);
            p.exponent = parse_int(exponent_text);
            p.normalize();
        }
        for (bool filled : seen)
            if (!filled)
                throw RuntimeError("incomplete except response");
        out.push_back(curve);
    }
    return out;
}

Agreement agreement(const Curves& a, const Curves& b) {
    Agreement out;
    if (a.size() != b.size())
        throw ValueError("curves of different state grids");
    for (std::size_t s = 0; s < a.size(); ++s) {
        if (a[s].size() != b[s].size())
            throw ValueError("curves of different time grids");
        for (std::size_t t = 0; t < a[s].size(); ++t) {
            const bool x = finite_point(a[s][t]), y = finite_point(b[s][t]);
            if (x != y)
                out.same_nonfinite = false;
            else if (x)
                out.largest = std::max(out.largest, std::hypot(a[s][t].first-b[s][t].first,
                                                               a[s][t].second-b[s][t].second) /
                                                        std::max(1., magnitude(a[s][t])));
        }
    }
    return out;
}

std::map<Sites, int> expansion_terms(const Partition& bundles, const Adjacency& links,
                                     const ExpansionLimits& limits) {
    return terms(bundles, connections(bundles, links), 2, limits);
}

Curves expansion(const Partition& bundles, const Adjacency& links, const std::map<Sites, Curves>& raw,
                 const ExpansionLimits& limits) {
    return approximate_weights(expansion_terms(bundles, links, limits), raw, false);
}

// ---------------------------------------------------------------------------
// the ledger
// ---------------------------------------------------------------------------
long long Source::point_work(const std::vector<Sites>&) const {
    throw RuntimeError("this source has no factors step");
}

bool Source::point_values(const Partition&, const std::vector<Edge>&, const Point&, PointValues*) {
    throw RuntimeError("this source has no factors step");
}

bool Source::kept(const Partition&, const std::vector<Sites>&, Kept*) {
    throw RuntimeError("this source has no except step");
}

void Source::forget(const std::vector<Sites>&) {
    throw RuntimeError("this source cannot drop curves");
}

Value Ledger::value() const {
    Value out = Value::object();
    out.set("max_total_work", Value::number(total));
    out.set("spent", Value::number(spent));
    Value by = Value::object();
    for (const char* what : {"initial_run", "culprit", "region_direct", "neighbour_direct", "candidate",
                             "outside", "rerun"}) {
        const std::map<std::string, double>::const_iterator found = spent_by.find(what);
        by.set(what, Value::number(found == spent_by.end() ? 0. : found->second));
    }
    out.set("spent_by", by);
    out.set("probes_calculated", Value::integer(probes));
    out.set("max_probe_dimension", Value::integer(max_dimension));
    if (factor_passes >= 0)
        out.set("factor_passes", Value::integer(factor_passes));
    if (except_passes >= 0)
        out.set("except_passes", Value::integer(except_passes));
    if (freed >= 0)
        out.set("curves_freed", Value::integer(freed));
    return out;
}

// ---------------------------------------------------------------------------
// the search
// ---------------------------------------------------------------------------
bool Search::acquire(const std::vector<std::pair<Sites, std::string> >& needed, double reserve) {
    std::set<Sites> wanted;
    std::map<Sites, std::string> absent;           // not in raw before: the purpose that first asked
    for (const std::pair<Sites, std::string>& one : needed) {
        wanted.insert(one.first);
        if (source_->raw().count(one.first) == 0 && absent.count(one.first) == 0)
            absent[one.first] = one.second;
    }
    source_->preload(wanted);
    // Deduplicated: what is asked for twice is calculated once and paid for once, under
    // the first purpose that asked for it. What the disk cache supplied costs nothing but is
    // held now as much as a calculated curve, so it is recorded for release too (second review
    // review 3: a warm cache freed nothing).
    std::map<Sites, std::string> missing;
    for (const std::pair<const Sites, std::string>& one : absent) {
        if (source_->raw().count(one.first) == 0)
            missing.insert(one);
        else
            fetched_.push_back(one);
    }
    if (missing.empty())
        return true;
    std::vector<Sites> listed;
    std::map<std::string, std::vector<Sites> > by;
    for (const std::pair<const Sites, std::string>& one : missing) {
        listed.push_back(one.first);
        by[one.second].push_back(one.first);
    }
    const double work = static_cast<double>(source_->work(listed));
    if (!ledger_->fits(work, reserve))
        return false;
    if (!source_->fetch(listed))
        return false;
    // The parts are disjoint and the proxy is a sum over clusters, so they add up to `work`.
    for (const std::pair<const std::string, std::vector<Sites> >& part : by)
        ledger_->charge(part.first, static_cast<double>(source_->work(part.second)));
    for (const std::pair<const Sites, std::string>& one : missing)
        fetched_.push_back(one);
    ledger_->probes += static_cast<long long>(listed.size());
    for (const Sites& sites : listed)
        ledger_->max_dimension = std::max(ledger_->max_dimension, source_->dimension(sites));
    return true;
}

void Search::settled(const Partition& units) {
    for (std::map<std::pair<Partition, std::vector<Sites> >, Kept>::iterator it = kept_.begin(); it != kept_.end();)
        it = it->first.first == units ? std::next(it) : kept_.erase(it);
}

Outcome Search::round(const Partition& units, const Curves& curves, const std::set<Partition>& visited,
                      const std::set<Partition>& rolled_back) {
    const Constants& c = settings_.constants;
    const std::map<Sites, Curves>& raw = source_->raw();
    Outcome out;
    out.at = worst_point(curves);
    out.unphysical_before = unphysical_states(curves, c.limit);
    out.excess_before = excess(curves[static_cast<std::size_t>(out.at.state)]);
    Value& row = out.row;
    row = Value::object();
    row.set("round", Value::integer(rounds_++));
    row.set("units", partition_value(units));
    row.set("s", Value::integer(static_cast<long long>(out.at.state)));
    row.set("t", Value::integer(static_cast<long long>(out.at.step)));
    row.set("U", Value::integer(out.unphysical_before));
    row.set("excess", trace_number(out.excess_before));
    if (settings_.trace_curves)
        row.set("curve", curves_value(curves));

    // Item 3: the factors at (s*, t*) need every bundle and every connected pair.
    const std::vector<Edge> edges = connections(units, *links_);
    row.set("edges", edges_value(edges));
    const auto refused = [&]() {
        row.set("culprits", Value::array());
        row.set("regions", Value::array());
        row.set("accepted", Value());
        out.stop = "budget_limited";
        return out;
    };
    PointValues at_point;
    if (settings_.culprit_factors) {
        // adaptive-regroup-stage3.md 1절: one factors step of the worst state, its work
        // reserved as "culprit" for that state alone. No connection, no step: nothing to rank.
        if (!edges.empty()) {
            const double work = static_cast<double>(source_->point_work(factor_clusters(units, edges)));
            if (!ledger_->fits(work, 0.) || !source_->point_values(units, edges, out.at, &at_point))
                return refused();
            ledger_->charge("culprit", work);
            ++ledger_->factor_passes;
        }
    } else {
        std::vector<std::pair<Sites, std::string> > needed;
        for (const Sites& bundle : units)
            needed.emplace_back(bundle, "culprit");
        for (const Edge& edge : edges)
            needed.emplace_back(union_of(units, {edge.first, edge.second}), "culprit");
        if (!acquire(needed, 0.))
            return refused();
        at_point = point_values_of(units, edges, raw, out.at);
    }
    const std::vector<Culprit> found = culprits(units, edges, at_point);
    Value listed = Value::array();
    for (const Culprit& one : found) {
        Value item = Value::array();
        item.push(sites_text(units[static_cast<std::size_t>(one.a)]));
        item.push(sites_text(units[static_cast<std::size_t>(one.b)]));
        listed.push(item);
    }
    row.set("culprits", listed);
    Value factors = Value::array();
    for (const Culprit& one : found) {
        Value item = Value::object();
        item.set("g", trace_number(one.factor.g));
        item.set("cause", one.factor.cause.empty() ? Value() : Value::text(one.factor.cause));
        factors.push(item);
    }
    row.set("factors", factors);
    Value regions = Value::array();
    row.set("accepted", Value());
    if (found.empty()) {
        row.set("regions", regions);
        out.stop = "no_positive_culprit";
        return out;
    }

    const std::size_t s_star = static_cast<std::size_t>(out.at.state);
    // The run's expansion limits are CCEX's refusal of a whole partition; the models of
    // a region and its neighbours hold at most four bundles and are only evaluated here.
    const ExpansionLimits local;
    std::set<Partition> tried;
    bool budget_refused = false;
    long long examined_here = 0;
    for (const Culprit& culprit : found) {
        const std::vector<int> region = region_of(units, culprit.a, culprit.b, *sites_);
        const Partition region_units = canonical(pick(units, region));
        Value one = Value::object();
        Value pair = Value::array();
        pair.push(sites_text(units[static_cast<std::size_t>(culprit.a)]));
        pair.push(sites_text(units[static_cast<std::size_t>(culprit.b)]));
        one.set("culprit", pair);
        one.set("region", partition_value(region_units));
        if (tried.count(region_units) != 0) {
            one.set("skipped", Value::boolean(true));
            regions.push(one);
            continue;
        }
        one.set("skipped", Value::boolean(false));
        if (examined_here == settings_.max_regions) {
            // The limit on regions per round: a region it leaves out is not checked, so
            // the round cannot say there is no candidate.
            one.set("refused", Value::text("region_limit"));
            regions.push(one);
            budget_refused = true;
            break;
        }
        ++examined_here;
        ++examined_;
        tried.insert(region_units);
        const std::size_t fetched_before = fetched_.size();
        // adaptive-regroup-stage3b.md 6절 3: with the except step the curves a region asked
        // for its candidates (and the excepted ones) are not wanted after it - the accepted
        // partition's terms come from the next run and except step - so they leave raw and
        // the probe count and cache limits do not fill up region by region. The work stays
        // spent in the ledger.
        const auto release = [&]() {
            if (!settings_.whole_except)
                return;
            std::vector<Sites> drop;
            for (std::size_t k = fetched_before; k < fetched_.size(); ++k)
                if (fetched_[k].second == "candidate" || fetched_[k].second == "outside")
                    drop.push_back(fetched_[k].first);
            fetched_.resize(fetched_before);
            if (!drop.empty()) {
                source_->forget(drop);
                ledger_->freed += static_cast<long long>(drop.size());
            }
        };
        const std::vector<int> outside = neighbours_of(units, region, *sites_, c.neighbours);
        Value near = Value::array();
        for (int k : outside)
            near.push(sites_text(units[static_cast<std::size_t>(k)]));
        one.set("neighbours", near);

        // The direct solutions (items 7, 8a, 8b): the region's, and the region with each
        // neighbour. Neither depends on the candidate.
        const Sites region_sites = all_sites(region_units);
        std::vector<Sites> with_sites;
        std::vector<std::pair<Sites, std::string> > directs;
        for (const std::pair<const Sites, int>& term : expansion_terms(region_units, *links_, local))
            directs.emplace_back(term.first, "culprit");
        directs.emplace_back(region_sites, "region_direct");
        for (int k : outside) {
            Partition with = region_units;
            with.push_back(units[static_cast<std::size_t>(k)]);
            with = canonical(with);
            with_sites.push_back(all_sites(with));
            for (const std::pair<const Sites, int>& term : expansion_terms(with, *links_, local))
                directs.emplace_back(term.first, "culprit");
            directs.emplace_back(with_sites.back(), "neighbour_direct");
        }
        const std::vector<Partition> listed_candidates = candidates(region_units, settings_.dimensions, settings_.cap);
        one.set("candidates", Value::integer(static_cast<long long>(listed_candidates.size())));
        std::map<std::string, long long> rejected;
        const auto finish_region = [&](const char* refused) {
            if (refused != nullptr)
                one.set("refused", Value::text(refused));
            Value counts = Value::object();
            for (const std::pair<const std::string, long long>& entry : rejected)
                counts.set(entry.first, Value::integer(entry.second));
            one.set("rejected", counts);
            one.set("accepted", Value());
            regions.push(one);
            release();
        };
        if (static_cast<long long>(listed_candidates.size()) > settings_.max_candidates) {
            rejected["budget_limited"] = static_cast<long long>(listed_candidates.size());
            budget_refused = true;
            finish_region("candidate_limit");
            continue;
        }
        if (!acquire(directs, 0.)) {
            rejected["budget_limited"] = static_cast<long long>(listed_candidates.size());
            budget_refused = true;
            finish_region("budget_limited");
            continue;
        }
        const LocalError now = local_error(expansion(region_units, *links_, raw, local),
                                           held(raw, region_sites));
        std::vector<LocalError> before;
        for (std::size_t k = 0; k < outside.size(); ++k) {
            Partition with = region_units;
            with.push_back(units[static_cast<std::size_t>(outside[k])]);
            before.push_back(local_error(expansion(canonical(with), *links_, raw, local),
                                         held(raw, with_sites[k])));
        }
        bool verifiable = now.verifiable;
        for (const LocalError& e : before)
            verifiable = verifiable && e.verifiable;
        one.set("e_now", trace_number(now.value));
        Value before_value = Value::array();
        for (const LocalError& e : before)
            before_value.push(trace_number(e.value));
        one.set("before", before_value);
        if (!verifiable) {
            // A direct solution that is not finite verifies nothing (item 7): every
            // candidate of the region is refused, not let through by a comparison
            // with inf or NaN.
            rejected["direct_nonfinite"] = static_cast<long long>(listed_candidates.size());
            finish_region("direct_nonfinite");
            continue;
        }

        struct Trial {
            Partition cand;
            Partition whole;
            std::map<Sites, int> whole_terms;
            std::map<Sites, int> added;        // N_c (whole_except)
            double reserve = 0.;
            double e = 0.;
        };
        // Stage A: what the candidate's own model reads, for item 8 (a).
        std::vector<Trial> trials;
        std::vector<std::pair<Sites, std::string> > batch;
        std::set<Sites> batched;
        double batch_work = 0.;
        // The largest rerun reserved by a candidate already in the batch: any of them may be the
        // one chosen, so a later candidate's probes may not spend it (second review 2).
        double batch_reserve = 0.;
        const auto reserve_for = [&](const std::vector<std::pair<Sites, std::string> >& wanted,
                                     double reserve) {
            // The cache first, as acquire does (review 20261012-233058-452278 P4-03): what
            // the disk holds costs nothing, so only the probes still missing after it are
            // charged; whether a curve was in RAM or on disk does not change the choice.
            // preload keeps the cache-memory limit (all or nothing), and what it supplied
            // is recorded for the region's release like acquire's cache hits.
            std::set<Sites> absent;
            for (const std::pair<Sites, std::string>& w : wanted)
                if (raw.count(w.first) == 0 && batched.count(w.first) == 0)
                    absent.insert(w.first);
            if (!absent.empty()) {
                source_->preload(absent);
                for (const std::pair<Sites, std::string>& w : wanted)
                    if (absent.erase(w.first) != 0 && raw.count(w.first) != 0)
                        fetched_.push_back(w);
            }
            std::vector<Sites> extra;
            for (const std::pair<Sites, std::string>& w : wanted)
                if (raw.count(w.first) == 0 && batched.count(w.first) == 0 &&
                    std::find(extra.begin(), extra.end(), w.first) == extra.end())
                    extra.push_back(w.first);
            const double more = extra.empty() ? 0. : static_cast<double>(source_->work(extra));
            if (!ledger_->fits(batch_work+more, std::max(reserve, batch_reserve)))
                return false;
            batch_reserve = std::max(reserve, batch_reserve);
            for (const std::pair<Sites, std::string>& w : wanted)
                if (raw.count(w.first) == 0 && batched.insert(w.first).second)
                    batch.push_back(w);
            batch_work += more;
            return true;
        };
        for (const Partition& cand : listed_candidates) {
            Trial trial;
            trial.cand = cand;
            trial.whole = replace(units, region, cand);
            if (rolled_back.count(trial.whole) != 0) {
                ++rejected["rolled_back"];
                continue;
            }
            if (visited.count(trial.whole) != 0) {
                ++rejected["revisit"];
                continue;
            }
            // CCEX refuses a run whose expansion passes its limits; refused here, before
            // anything is asked, it is one rejected candidate rather than a failed run
            // (in-process, the end of the whole job: review 14).
            try {
                trial.whole_terms = expansion_terms(trial.whole, *links_, settings_.limits);
            } catch (const ExpansionLimit&) {
                ++rejected["expansion_limit"];
                continue;
            }
            trial.reserve = static_cast<double>(source_->run_work(trial.whole_terms));
            std::vector<std::pair<Sites, std::string> > wanted;
            for (const std::pair<const Sites, int>& term : expansion_terms(cand, *links_, local))
                wanted.emplace_back(term.first, "candidate");
            if (!reserve_for(wanted, trial.reserve)) {
                ++rejected["budget_limited"];
                budget_refused = true;
                continue;
            }
            trials.push_back(trial);
        }
        if (!batch.empty() && !acquire(batch, 0.)) {
            rejected["budget_limited"] += static_cast<long long>(trials.size());
            budget_refused = true;
            trials.clear();
        }
        std::vector<Trial> passed;
        for (Trial& trial : trials) {
            const Curves model = expansion(trial.cand, *links_, raw, local);
            if (!physical(model, c.limit)) {
                ++rejected["unphysical"];
                continue;
            }
            trial.e = local_error(model, held(raw, region_sites)).value;
            if (!improves(trial.e, now.value, c.epsilon)) {                  // (a)
                ++rejected["a"];
                continue;
            }
            passed.push_back(trial);
        }
        // Stage B: the candidate with each neighbour (b), and the whole expansion (c).
        // With the except step (adaptive-regroup-stage3b.md 1절, 3절) the whole expansion of a
        // candidate is P_K, the product of the terms no candidate of the region changes, times
        // its own new terms N_c: one except step for the region, then only N_c as probes.
        const Kept* kept = nullptr;
        if (settings_.whole_except && !passed.empty()) {
            const std::map<Sites, int> now_terms = expansion_terms(units, *links_, settings_.limits);
            std::vector<const std::map<Sites, int>*> tables;
            double least_reserve = HUGE_VAL;
            for (const Trial& trial : passed) {
                tables.push_back(&trial.whole_terms);
                least_reserve = std::min(least_reserve, trial.reserve);
            }
            const std::vector<Sites> except = changed_terms(now_terms, tables);
            one.set("except_terms", Value::integer(static_cast<long long>(except.size())));
            const std::pair<Partition, std::vector<Sites> > key(units, except);
            std::map<std::pair<Partition, std::vector<Sites> >, Kept>::const_iterator done = kept_.find(key);
            if (done == kept_.end()) {
                // its work, a run's, reserved as "outside" with room left for one rerun
                const double work = static_cast<double>(source_->run_work(now_terms));
                Kept answer;
                if (ledger_->fits(work, least_reserve) && source_->kept(units, except, &answer)) {
                    ledger_->charge("outside", work);
                    ++ledger_->except_passes;
                    for (const Sites& sites : answer.added)
                        fetched_.emplace_back(sites, "outside");
                    // The free checks of stage3b 1절: the step's whole curve is the run's, and
                    // P_K times the excepted terms at their powers is that curve again - a state,
                    // seed or table that differs shows here, and the round fails.
                    std::map<Sites, int> excepted;
                    for (const Sites& sites : except)
                        excepted[sites] = now_terms.at(sites);
                    const Agreement run_check = agreement(curves, answer.whole);
                    const Agreement rebuilt = agreement(answer.whole, approximate_from(answer.product, excepted, raw));
                    Value check = Value::object();
                    check.set("run", trace_number(run_check.largest));
                    check.set("reassembled", trace_number(rebuilt.largest));
                    one.set("except_check", check);
                    if (!run_check.same_nonfinite || !(run_check.largest <= kExceptTolerance) ||
                        !rebuilt.same_nonfinite || !(rebuilt.largest <= kExceptTolerance))
                        throw RuntimeError("regroup: the except step does not reproduce the run (relative differences "+
                                           write_json(check, 0)+")");
                    answer.whole = Curves();
                    done = kept_.emplace(key, answer).first;
                } else {
                    rejected["budget_limited"] += static_cast<long long>(passed.size());
                    budget_refused = true;
                    passed.clear();
                }
            }
            if (done != kept_.end()) {
                kept = &done->second;
                for (Trial& trial : passed)
                    trial.added = new_terms(now_terms, except, trial.whole_terms);
            }
        }
        batch.clear();
        batched.clear();
        batch_work = 0.;
        batch_reserve = 0.;
        std::vector<Trial> ready;
        long long new_count = 0;
        for (const Trial& trial : passed) {
            std::vector<std::pair<Sites, std::string> > wanted;
            for (int k : outside) {
                Partition with = trial.cand;
                with.push_back(units[static_cast<std::size_t>(k)]);
                for (const std::pair<const Sites, int>& term :
                     expansion_terms(canonical(with), *links_, local))
                    wanted.emplace_back(term.first, "candidate");
            }
            for (const std::pair<const Sites, int>& term : (kept != nullptr ? trial.added : trial.whole_terms))
                wanted.emplace_back(term.first, "outside");
            new_count += static_cast<long long>(trial.added.size());
            if (!reserve_for(wanted, trial.reserve)) {
                ++rejected["budget_limited"];
                budget_refused = true;
                continue;
            }
            ready.push_back(trial);
        }
        if (!batch.empty() && !acquire(batch, 0.)) {
            rejected["budget_limited"] += static_cast<long long>(ready.size());
            budget_refused = true;
            ready.clear();
        }
        if (kept != nullptr)
            one.set("new_terms", Value::integer(new_count));
        const Trial* best = nullptr;
        std::vector<double> best_after;
        Curves best_curve;
        long long passed_all = 0;
        for (const Trial& trial : ready) {
            std::vector<double> after, was;
            for (std::size_t k = 0; k < outside.size(); ++k) {
                Partition with = trial.cand;
                with.push_back(units[static_cast<std::size_t>(outside[k])]);
                after.push_back(local_error(expansion(canonical(with), *links_, raw, local),
                                            held(raw, with_sites[k])).value);
                was.push_back(before[k].value);
            }
            if (!neighbours_hold(after, was, c.epsilon)) {                   // (b)
                ++rejected["b"];
                continue;
            }
            // (c) on the whole expansion, assembled again from its exponent table - every
            // term, the kept ones included - and never as a ratio update of the old curve
            // (item 9; review 5).
            const Curves whole = kept != nullptr ? approximate_from(kept->product, trial.added, raw)
                                                 : approximate_weights(trial.whole_terms, raw, false);
            if (!whole_improves(unphysical_states(whole, c.limit), excess(whole[s_star]),
                                out.unphysical_before, out.excess_before)) {
                ++rejected["c"];
                continue;
            }
            ++passed_all;
            if (best == nullptr || trial.e < best->e || (trial.e == best->e && trial.cand < best->cand)) {
                best = &trial;
                best_after = after;
                best_curve = whole;
            }
        }
        one.set("passed", Value::integer(passed_all));
        if (best == nullptr) {
            finish_region(nullptr);
            continue;
        }
        Value counts = Value::object();
        for (const std::pair<const std::string, long long>& entry : rejected)
            counts.set(entry.first, Value::integer(entry.second));
        one.set("rejected", counts);
        one.set("accepted", partition_value(best->cand));
        one.set("e", trace_number(best->e));
        Value after_value = Value::array();
        for (double e : best_after)
            after_value.push(trace_number(e));
        one.set("after", after_value);
        one.set("U_new", Value::integer(unphysical_states(best_curve, c.limit)));
        one.set("excess_new", trace_number(excess(best_curve[s_star])));
        one.set("edges_new", edges_value(connections(best->whole, *links_)));
        // How much of the exponent table changed (item 9): site sets whose power differs.
        // The current partition is held to the run's limits too; one CCEX did not run
        // (only a unit test hands one in) has no table to compare, and the count is left out.
        try {
            const std::map<Sites, int> was = expansion_terms(units, *links_, settings_.limits);
            long long changed = 0;
            std::set<Sites> keys;
            for (const std::pair<const Sites, int>& term : was) keys.insert(term.first);
            for (const std::pair<const Sites, int>& term : best->whole_terms) keys.insert(term.first);
            for (const Sites& key : keys) {
                const std::map<Sites, int>::const_iterator a = was.find(key), b = best->whole_terms.find(key);
                if ((a == was.end() ? 0 : a->second) != (b == best->whole_terms.end() ? 0 : b->second))
                    ++changed;
            }
            one.set("changed_terms", Value::integer(changed));
        } catch (const ExpansionLimit&) {
            one.set("changed_terms", Value());
        }
        regions.push(one);
        row.set("regions", regions);
        row.set("accepted", partition_value(best->whole));
        out.accepted = best->whole;
        out.predicted = best_curve;
        release();
        return out;
    }
    row.set("regions", regions);
    // Item 10: no candidate anywhere - unless a check was left undone for want of budget,
    // in which case the round cannot say so.
    out.stop = budget_refused ? "budget_limited" : "no_candidate";
    return out;
}

Result drive(const Partition& start, const Curves& first, Search* search, const Constants& constants,
             const std::function<Curves(const Partition&)>& run,
             const std::function<void(bool)>& settle, Value* rows) {
    Result out;
    Partition units = canonical(start);
    Curves curves = first;
    std::set<Partition> visited{units}, rolled_back;
    for (;;) {
        if (physical(curves, constants.limit)) {
            out.stop = "physical";
            break;
        }
        if (out.changes >= constants.max_changes) {
            out.stop = "change_limit";
            break;
        }
        Outcome found = search->round(units, curves, visited, rolled_back);
        if (!found.stop.empty()) {
            rows->push(found.row);
            out.stop = found.stop;
            break;
        }
        // Item 9: the next real run checks what the round predicted. U may not grow and
        // the excess of s* has to fall; a run that does not confirm it is undone and the
        // partition blocked.
        const Curves actual = run(found.accepted);
        const std::size_t s = static_cast<std::size_t>(found.at.state);
        const long long u = unphysical_states(actual, constants.limit);
        const double x = excess(actual[s]);
        const bool confirmed = whole_improves(u, x, found.unphysical_before, found.excess_before);
        double largest = 0.;
        bool same_nonfinite = actual.size() == found.predicted.size();
        for (std::size_t i = 0; same_nonfinite && i < actual.size(); ++i) {
            same_nonfinite = actual[i].size() == found.predicted[i].size();
            for (std::size_t t = 0; same_nonfinite && t < actual[i].size(); ++t) {
                const bool a = finite_point(actual[i][t]), b = finite_point(found.predicted[i][t]);
                if (a != b)
                    same_nonfinite = false;
                else if (a)
                    largest = std::max(largest, std::hypot(actual[i][t].first-found.predicted[i][t].first,
                                                           actual[i][t].second-found.predicted[i][t].second));
            }
        }
        Value verified = Value::object();
        verified.set("U_actual", Value::integer(u));
        verified.set("excess_actual", trace_number(x));
        verified.set("max_abs_difference_from_prediction", trace_number(largest));
        if (search->settings().whole_except)
            // the bound of adaptive-regroup-stage3b.md 4절 is relative
            verified.set("max_relative_difference_from_prediction",
                         trace_number(same_nonfinite ? agreement(found.predicted, actual).largest : HUGE_VAL));
        verified.set("same_nonfinite_points", Value::boolean(same_nonfinite));
        verified.set("confirmed", Value::boolean(confirmed));
        found.row.set("verification", verified);
        rows->push(found.row);
        settle(confirmed);
        if (!confirmed) {
            rolled_back.insert(found.accepted);
            // A prediction the runs keep contradicting is a fault in the assembly, not a
            // search to go on with.
            if (++out.rollbacks > constants.max_changes)
                throw RuntimeError("regroup: more rollbacks than changes allowed; the assembled predictions do not "
                                   "match the runs");
            continue;
        }
        units = found.accepted;
        curves = actual;
        visited.insert(units);
        search->settled(units);
        ++out.changes;
    }
    out.final = units;
    return out;
}
}   // namespace regroup
}   // namespace adaptive

// ----------------------------------------------------------------------------
// the merge Driver with --regroup: real runs and probes
// ----------------------------------------------------------------------------

namespace adaptive {
int Driver::run_regroup() {
    using namespace regroup;
    std::unique_ptr<Executed> step = execute("run", nullptr, nullptr, false);
    if (!has_physics_)
        throw RuntimeError("regroup needs the site positions, and the run wrote no response.physics.json");
    const Geometry geometry = geometry_of(physics_);
    if (geometry.size() != step->plan.site_dimensions.size())
        throw RuntimeError("the physics metadata and the plan disagree on the number of sites");
    // CCEX keeps both radii as floats (Config.rdip, Config.rdipcut) and compares a double
    // distance with them: the same conversion here gives the same connections.
    const float rdip = static_cast<float>(number_of(cfg_.at("rdip")));
    const float rdipcut = cfg_.has("rdipcut") ? static_cast<float>(number_of(cfg_.at("rdipcut"))) : 0.f;
    const Adjacency links = distance_adjacency(geometry, rdip, rdipcut);
    // Item 6 against CCEX itself, at every run: the connections D works out from the
    // sites have to be the ones the binary planned with (adaptive-regroup-plan.md, the
    // table of reused parts, row 4).
    const auto check_connections = [&links](const Executed& done) {
        if (connections(done.plan.bundles, links) != done.plan.plain_edges())
            throw RuntimeError("the connections worked out from the site distances differ from CCEX's plan");
    };
    check_connections(*step);
    if (canonical(step->plan.bundles) != step->plan.bundles)
        throw RuntimeError("CCEX's partition is not in canonical order");

    Settings settings;
    settings.cap = merge_.bundle_max_dim;
    settings.dimensions = step->plan.site_dimensions;
    settings.limits = config_limits();
    settings.max_regions = limits_.regroup_max_regions;
    settings.max_candidates = limits_.regroup_max_candidates;
    settings.culprit_factors = limits_.regroup_factors;
    settings.whole_except = limits_.regroup_except;
    const long long states = static_cast<long long>(states_.size());
    const long long steps = static_cast<long long>(cfg_.at("nstep").as_int());

    // The driver's probes, as the search asks for them. A local class has the member
    // function's access, so the cache can be asked first: what it holds costs nothing.
    struct Probes : Source {
        Driver* self;
        const Executed** current;
        long long states, steps;
        const std::map<Sites, Curves>& raw() const override { return self->raw_; }
        void preload(const std::set<Sites>& wanted) override { self->load_cached(wanted); }
        long long work(const std::vector<Sites>& clusters) const override {
            return cluster_work((*current)->plan, clusters, states, steps, 0);
        }
        long long run_work(const std::map<Sites, int>& terms) const override {
            std::vector<Sites> clusters;
            for (const std::pair<const Sites, int>& term : terms)
                clusters.push_back(term.first);
            return cluster_work((*current)->plan, clusters, states, steps, 0);
        }
        long long dimension(const Sites& sites) const override {
            return saturating_mul((*current)->plan.qubit_dimension,
                                  bundle_dimension(sites, (*current)->plan.site_dimensions));
        }
        bool fetch(const std::vector<Sites>& missing) override {
            return self->ensure(missing, (*current)->plan.bundles);
        }
        // adaptive-regroup-stage3.md: one factors step of the worst state and step. Its
        // connections are checked against the run's like every step's.
        std::function<void(const Executed&)> check;
        long long point_work(const std::vector<Sites>& clusters) const override {
            return cluster_work((*current)->plan, clusters, 1, steps, 0);
        }
        bool point_values(const Partition& units, const std::vector<Edge>& edges, const Point& at,
                          PointValues* out) override {
            if (at.state < 0 || at.state >= static_cast<int>(self->states_.size()) || at.step < 0 || at.step >= steps)
                throw RuntimeError("regroup: the worst point is outside the run's grid");
            const int state = self->states_[static_cast<std::size_t>(at.state)];
            Value more = Value::object();
            more.set("state", Value::integer(static_cast<long long>(state)));
            more.set("step", Value::integer(static_cast<long long>(at.step)));
            const std::unique_ptr<Executed> done = self->execute("factors", &units, nullptr, false, &more);
            check(*done);
            *out = read_factors(done->directory, self->settings_.ranks, done->plan.bundles, edges, state, at.step);
            return true;
        }
        // adaptive-regroup-stage3b.md: one except step. The excepted curves are held in raw like
        // probes, so the probe count and cache limits apply to them first (ensure's checks).
        bool kept(const Partition& units, const std::vector<Sites>& except, Kept* out) override {
            long long more = 0;
            for (const Sites& sites : except)
                more += self->raw_.count(sites) == 0 ? 1 : 0;
            const long long held = static_cast<long long>(self->raw_.size());
            if (held+more > self->settings_.max_probes ||
                (held+more)*states*steps*64 > self->settings_.cache_mb*1024*1024)
                return false;
            Value request = Value::object();
            request.set("except", sites_list(except));
            const std::unique_ptr<Executed> done = self->execute("except", &units, nullptr, false, &request);
            check(*done);
            const std::map<Sites, Curves> curves =
                read_probe_files(done->directory, self->states_, steps, self->settings_.ranks, except);
            out->whole = read_state_curves(done->directory, self->states_);
            out->product = read_except(done->directory, self->states_, steps);
            out->added.clear();
            // A curve already held stays as it was: a probe, which the earlier decisions read.
            // These are not put in the probe cache either - they come from the run, not a probe step.
            for (const Sites& sites : except)
                if (self->raw_.count(sites) == 0) {
                    self->raw_[sites] = curves.at(sites);
                    out->added.push_back(sites);
                }
            return true;
        }
        void forget(const std::vector<Sites>& clusters) override {
            for (const Sites& sites : clusters)
                if (!sites.empty())
                    self->raw_.erase(sites);
        }
    };
    const Executed* current = step.get();
    Probes probes;
    probes.self = this;
    probes.current = &current;
    probes.states = states;
    probes.steps = steps;
    probes.check = check_connections;

    Ledger ledger;
    ledger.total = settings_.max_total_work;
    if (settings.culprit_factors)
        ledger.factor_passes = 0;
    if (settings.whole_except) {
        ledger.except_passes = 0;
        ledger.freed = 0;
    }
    ledger.charge("initial_run", static_cast<double>(probes.run_work(
        expansion_terms(step->plan.bundles, links, settings.limits))));
    Search search(settings, &geometry, &links, &probes, &ledger);
    const Constants constants;
    Value rows = Value::array();
    const auto record = [&](const std::string& stop, const Result* result) {
        Value block = Value::object();
        block.set("algorithm", Value::text(kAlgorithm));
        Value listed = Value::array();
        for (int state : states_)
            listed.push(Value::integer(static_cast<long long>(state)));
        block.set("states", listed);
        block.set("epsilon", Value::number(constants.epsilon));
        block.set("limit", Value::number(constants.limit));
        block.set("neighbours", Value::integer(static_cast<long long>(constants.neighbours)));
        block.set("max_changes", Value::integer(static_cast<long long>(constants.max_changes)));
        block.set("max_regions", Value::integer(settings.max_regions));
        block.set("max_candidates", Value::integer(settings.max_candidates));
        block.set("connections", Value::text("distance"));
        // the two sources, recorded only when not the probes they were before stages 3 and 3b
        if (settings.culprit_factors)
            block.set("culprit_factors", Value::text("factors"));
        if (settings.whole_except)
            block.set("whole_curve", Value::text("except"));
        block.set("stop", Value::text(stop));
        block.set("changes", Value::integer(result != nullptr ? result->changes : 0));
        block.set("rollbacks", Value::integer(result != nullptr ? result->rollbacks : 0));
        block.set("regions_examined", Value::integer(search.examined()));
        block.set("ledger", ledger.value());
        block.set("rounds", rows);
        report_.set("regroup", block);
    };
    Value limitations = Value::array();
    limitations.push(Value::text("Physical (every |L| <= 1.01) is not accurate: a physical curve can still be far "
                                 "from the answer."));
    limitations.push(Value::text("Each change is checked against direct solutions of its region and of three "
                                 "neighbouring regions only, never against the whole bath."));
    limitations.push(Value::text("D was designed on twelve-site P1 subsystems with ten stored states; "
                                 "elsewhere it is untested."));
    report_.set("limitations", limitations);

    std::unique_ptr<Executed> trial;
    const auto run = [&](const Partition& partition) {
        ledger.charge("rerun", static_cast<double>(probes.run_work(
            expansion_terms(partition, links, settings.limits))));
        trial = execute("run", &partition, nullptr, false);
        check_connections(*trial);
        return read_state_curves(trial->directory, states_);
    };
    const auto settle = [&](bool kept) {
        if (kept) {
            step = std::move(trial);
            current = step.get();
        }
        trial.reset();
        record("running", nullptr);
        persist();
    };
    Result result;
    try {
        result = drive(step->plan.bundles, read_state_curves(step->directory, states_), &search, constants,
                       run, settle, &rows);
    } catch (const Error&) {
        record("failed", nullptr);
        throw;
    }
    record(result.stop, &result);
    return finish(result.stop, step->directory, *step);
}
}   // namespace adaptive
