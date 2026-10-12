// The merge policy's judgement, against CPython's own methods (plan stage 1, Task 4).
#include "check.h"
#include "base_bundles.h"
#include <functional>
#include <complex>
#include "numeric_environment.h"
#include "base_json.h"
#include "policy_merge.h"
#include "base_numeric.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
double from_bits(const std::string& text) {
    unsigned long long word = std::stoull(text, nullptr, 16);
    double value = 0.;
    std::memcpy(&value, &word, sizeof value);
    return value;
}

// The same plans the generator built, in the same order.
std::vector<adaptive::Plan> plans() {
    struct Shape {
        int n;
        std::vector<std::vector<double>> edges;
        std::vector<int> dimensions;
        int qubit;
    };
    std::vector<Shape> shapes{
        {2, {{0, 1, .5}}, {2, 2}, 1},
        {3, {{0, 1, .5}, {1, 2, .25}}, {2, 2, 2}, 1},
        {4, {{0, 1, .9}, {1, 2, .8}, {2, 3, .7}, {0, 3, .6}}, {2, 2, 2, 2}, 1},
        {4, {{0, 1, .9}, {1, 2, .8}}, {2, 2, 4, 2}, 2},
        {6, {}, {2, 2, 2, 2, 2, 2}, 1},
        {9, {}, {2, 2, 2, 2, 2, 2, 2, 2, 2}, 1},
    };
    for (int i = 0; i < 5; ++i)
        shapes[4].edges.push_back({double(i), double(i+1), 1./(i+1)});
    for (int i = 0; i < 8; ++i)
        shapes[5].edges.push_back({double(i), double(i+1), 1./(i+1)});
    std::vector<adaptive::Plan> out;
    for (const Shape& shape : shapes) {
        adaptive::Plan plan;
        for (int i = 0; i < shape.n; ++i)
            plan.bundles.push_back(adaptive::Sites{i});
        for (const std::vector<double>& edge : shape.edges)
            plan.edges.push_back({static_cast<int>(edge[0]), static_cast<int>(edge[1]), edge[2]});
        for (int d : shape.dimensions)
            plan.site_dimensions.push_back(static_cast<long long>(d));
        plan.qubit_dimension = static_cast<long long>(shape.qubit);
        out.push_back(plan);
    }
    return out;
}

std::complex<double> made(int k) { return {1./(1+k%5), (k%3)/8.}; }

adaptive::Curves curves(int seed, int states = 2, int steps = 3) {
    adaptive::Curves out;
    for (int s = 0; s < states; ++s) {
        adaptive::Curve row;
        for (int t = 0; t < steps; ++t) {
            const std::complex<double> z = made(seed+s*7+t*13);
            row.emplace_back(z.real(), z.imag());
        }
        out.push_back(row);
    }
    return out;
}

// raw for every subset of the plan's bundles, exactly as the generator built it.
std::map<adaptive::Sites, adaptive::Curves> raw_for(const adaptive::Plan& plan) {
    std::map<adaptive::Sites, adaptive::Curves> raw;
    raw[adaptive::Sites()] = curves(0);
    const std::size_t n = plan.bundles.size();
    for (std::size_t mask = 1; mask < (std::size_t(1) << n); ++mask) {
        std::vector<int> ids;
        int sum = 0, size = 0;
        for (std::size_t i = 0; i < n; ++i)
            if ((mask >> i) & 1) { ids.push_back(static_cast<int>(i)); sum += static_cast<int>(i); ++size; }
        raw[adaptive::union_of(plan.bundles, ids)] = curves(sum+size);
    }
    return raw;
}

std::string shown_patches(const std::pair<std::vector<std::vector<int>>, adaptive::Coverage>& found) {
    adaptive::Value out = adaptive::Value::array();
    adaptive::Value listed_patches = adaptive::Value::array();
    for (const std::vector<int>& patch : found.first) {
        adaptive::Value ids = adaptive::Value::array();
        for (int i : patch)
            ids.push(adaptive::Value::integer(static_cast<long long>(i)));
        listed_patches.push(ids);
    }
    out.push(listed_patches);
    adaptive::Value coverage = adaptive::Value::object();
    coverage.set("candidate_patches", adaptive::Value::integer(static_cast<long long>(found.second.candidate_patches)));
    coverage.set("checked_patches", adaptive::Value::integer(static_cast<long long>(found.second.checked_patches)));
    coverage.set("work_blocked", adaptive::Value::integer(static_cast<long long>(found.second.work_blocked)));
    coverage.set("coverage_limited", adaptive::Value::boolean(found.second.coverage_limited));
    coverage.set("whole_bath_oracle", adaptive::Value::boolean(found.second.whole_bath_oracle));
    coverage.set("max_inverse_log10_amplification",
                 adaptive::Value::number(found.second.max_inverse_log10_amplification));
    out.push(coverage);
    return adaptive::json_dumps(out, false, false, true);
}

std::string listed(const std::vector<int>& values) {
    adaptive::Value out = adaptive::Value::array();
    for (int v : values)
        out.push(adaptive::Value::integer(static_cast<long long>(v)));
    return adaptive::json_dumps(out, false, false, true);
}
}   // namespace

// The driver's budget question, as the corpus poses it: a work ceiling in the
// D^3 * states * steps currency the policies use.
// The rule itself is the driver's (adaptive::affordable_clusters): a copy of it here
// went on charging clusters the fixture already held after the driver stopped.
std::function<bool(const adaptive::Plan&, const std::vector<adaptive::Sites>&)>
affordable_within(double work, const std::map<adaptive::Sites, adaptive::Curves>* held) {
    return [work, held](const adaptive::Plan& plan, const std::vector<adaptive::Sites>& clusters) {
        return adaptive::affordable_clusters(plan, clusters, *held, 2, 3, static_cast<long long>(work));
    };
}

int main(int argc, char** argv) {
    std::ifstream file(argc > 1 ? argv[1] : "../reference/semantics/merge_decisions.tsv");
    const std::vector<adaptive::Plan> all = plans();
    const bool strict = testing::same_numeric_environment();
    std::string line;
    int rows = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream parts(line);
        std::vector<std::string> column;
        std::string word;
        while (std::getline(parts, word, '\t'))
            column.push_back(word);
        std::string got, message;
        // The fuzz family carries the plan and the probes with it; the fixed one
        // refers to a plan by index. They are told apart by the second column.
        const bool carried = column[1].size() > 1 && column[1][0] == '[';
        if (carried) {
            const std::size_t where = column[0] == "patches" ? 7 : 6;
            const std::string expect = column.size() > where ? column[where] : std::string();
            const std::string expect_message = column.size() > where+1 ? column[where+1]
                                                                       : std::string();
            std::string found, trouble;
            try {
                const adaptive::Value shape = adaptive::read_json(column[1]);
                adaptive::Plan plan;
                for (std::size_t i = 0; i < shape.at(0).size(); ++i) {
                    adaptive::Sites bundle;
                    for (std::size_t k = 0; k < shape.at(0).at(i).size(); ++k)
                        bundle.push_back(static_cast<int>(shape.at(0).at(i).at(k).as_int()));
                    plan.bundles.push_back(bundle);
                }
                for (std::size_t i = 0; i < shape.at(1).size(); ++i)
                    plan.edges.push_back({static_cast<int>(shape.at(1).at(i).at(0).as_int()),
                                          static_cast<int>(shape.at(1).at(i).at(1).as_int()),
                                          shape.at(1).at(i).at(2).as_double()});
                for (std::size_t i = 0; i < shape.at(2).size(); ++i)
                    plan.site_dimensions.push_back(shape.at(2).at(i).as_int());
                plan.qubit_dimension = shape.at(3).as_int();
                if (shape.at(4).kind() != adaptive::Value::Kind::Null) {
                    plan.has_candidate_edges = true;
                    for (std::size_t i = 0; i < shape.at(4).size(); ++i)
                        plan.candidate_edges.push_back(
                            {static_cast<int>(shape.at(4).at(i).at(0).as_int()),
                             static_cast<int>(shape.at(4).at(i).at(1).as_int()),
                             shape.at(4).at(i).at(2).as_double()});
                }
                const adaptive::Value listed_raw = adaptive::read_json(column[2]);
                std::map<adaptive::Sites, adaptive::Curves> raw;
                for (const std::string& key : listed_raw.keys()) {
                    adaptive::Sites sites;
                    for (std::size_t k = 0; k < key.size(); ) {
                        const std::size_t comma = key.find(',', k);
                        sites.push_back(static_cast<int>(adaptive::from_decimal(
                            key.substr(k, comma == std::string::npos ? std::string::npos : comma-k))));
                        if (comma == std::string::npos) break;
                        k = comma+1;
                    }
                    adaptive::Curves curves_here;
                    for (std::size_t c = 0; c < listed_raw.at(key).size(); ++c) {
                        adaptive::Curve curve;
                        for (std::size_t t = 0; t < listed_raw.at(key).at(c).size(); ++t)
                            curve.emplace_back(
                                from_bits(listed_raw.at(key).at(c).at(t).at(0).as_string()),
                                from_bits(listed_raw.at(key).at(c).at(t).at(1).as_string()));
                        curves_here.push_back(curve);
                    }
                    raw[sites] = curves_here;
                }
                adaptive::MergeContext context;
                context.raw = &raw;
                const adaptive::Value cfg =
                    adaptive::read_json(column[column[0] == "patches" ? 6 : 5]);
                context.bundle_max_dim = cfg.at("bundle_max_dim").as_int();
                context.affordable = affordable_within(8e12, &raw);
                if (column[0] == "patches") {
                    context.order = std::stoi(column[3]);
                    context.neighbors = std::stoll(column[4]);
                    context.max_patches = std::stoll(column[5]);
                    found = shown_patches(adaptive::patches(context, plan));
                } else {
                    context.order = std::stoi(column[3]);
                    const adaptive::Value pair = adaptive::read_json(column[4]);
                    adaptive::Feasibility cache;
                    found = adaptive::can_merge(&cache, context, plan,
                                              {static_cast<int>(pair.at(0).as_int()),
                                               static_cast<int>(pair.at(1).as_int())})
                            ? "true" : "false";
                }
            } catch (const adaptive::Error& error) {
                found = error.type();
                trouble = error.what();
            }
            CHECK_EQ(found, expect);
            CHECK_EQ(trouble, expect_message);
            ++rows;
            continue;
        }
        const std::size_t at = column[0] == "patches" ? 6
                             : column[0] == "can_merge" ? 5
                             : column[0] == "bundle_cap" ? 4
                             : column[0] == "merge_bundles" ? 3
                             : column[0] == "proposal_key" ? 5 : 0;
        try {
            if (column[0] == "patches" || column[0] == "can_merge" || column[0] == "bundle_cap") {
                const adaptive::Plan& plan = all.at(std::stoul(column[1]));
                const std::map<adaptive::Sites, adaptive::Curves> raw = raw_for(plan);
                adaptive::MergeContext context;
                context.order = std::stoi(column[2]);
                context.raw = &raw;
                if (column[0] == "patches") {
                    context.neighbors = std::stoll(column[3]);
                    context.max_patches = std::stoll(column[4]);
                    context.affordable = affordable_within(from_bits(column[5]), &raw);
                    const std::pair<std::vector<std::vector<int>>, adaptive::Coverage> found =
                        adaptive::patches(context, plan);
                    adaptive::Value out = adaptive::Value::array();
                    adaptive::Value listed_patches = adaptive::Value::array();
                    for (const std::vector<int>& patch : found.first) {
                        adaptive::Value ids = adaptive::Value::array();
                        for (int i : patch)
                            ids.push(adaptive::Value::integer(static_cast<long long>(i)));
                        listed_patches.push(ids);
                    }
                    out.push(listed_patches);
                    adaptive::Value coverage = adaptive::Value::object();
                    coverage.set("candidate_patches",
                                 adaptive::Value::integer(static_cast<long long>(found.second.candidate_patches)));
                    coverage.set("checked_patches",
                                 adaptive::Value::integer(static_cast<long long>(found.second.checked_patches)));
                    coverage.set("work_blocked",
                                 adaptive::Value::integer(static_cast<long long>(found.second.work_blocked)));
                    coverage.set("coverage_limited", adaptive::Value::boolean(found.second.coverage_limited));
                    coverage.set("whole_bath_oracle", adaptive::Value::boolean(found.second.whole_bath_oracle));
                    coverage.set("max_inverse_log10_amplification",
                                 adaptive::Value::number(found.second.max_inverse_log10_amplification));
                    out.push(coverage);
                    got = adaptive::json_dumps(out, false, false, true);
                } else {
                    context.bundle_max_dim = adaptive::from_decimal(column[3]);
                    if (column[0] == "bundle_cap") {
                        got = adaptive::decimal(adaptive::bundle_cap(context, plan));
                    } else {
                        const adaptive::Value pair = adaptive::read_json(column[4]);
                        adaptive::Feasibility cache;
                        got = adaptive::can_merge(&cache, context, plan,
                                                {static_cast<int>(pair.at(0).as_int()),
                                                 static_cast<int>(pair.at(1).as_int())})
                              ? "true" : "false";
                    }
                }
            } else if (column[0] == "merge_bundles") {
                const adaptive::Plan& plan = all.at(std::stoul(column[1]));
                const adaptive::Value pair = adaptive::read_json(column[2]);
                const std::pair<std::vector<adaptive::Sites>, std::vector<adaptive::Edge>> found =
                    adaptive::merge_bundles(plan.bundles, plan.plain_edges(),
                                         {static_cast<int>(pair.at(0).as_int()),
                                          static_cast<int>(pair.at(1).as_int())});
                adaptive::Value out = adaptive::Value::array();
                adaptive::Value bundles = adaptive::Value::array();
                for (const adaptive::Sites& bundle : found.first) {
                    adaptive::Value ids = adaptive::Value::array();
                    for (int site : bundle)
                        ids.push(adaptive::Value::integer(static_cast<long long>(site)));
                    bundles.push(ids);
                }
                out.push(bundles);
                adaptive::Value edges = adaptive::Value::array();
                for (const adaptive::Edge& edge : found.second) {
                    adaptive::Value ids = adaptive::Value::array();
                    ids.push(adaptive::Value::integer(static_cast<long long>(edge.first)));
                    ids.push(adaptive::Value::integer(static_cast<long long>(edge.second)));
                    edges.push(ids);
                }
                out.push(edges);
                got = adaptive::json_dumps(out, false, false, true);
            } else if (column[0] == "proposal_key") {
                adaptive::MergeContext context;
                context.merge_priority = column[1];
                const adaptive::Value pair = adaptive::read_json(column[4]);
                const adaptive::ProposalKey key = adaptive::proposal_key(
                    context, from_bits(column[2]), from_bits(column[3]),
                    {static_cast<int>(pair.at(0).as_int()), static_cast<int>(pair.at(1).as_int())});
                adaptive::Value out = adaptive::Value::array();
                if (key.has_worst_error)
                    out.push(adaptive::Value::number(key.worst_error));
                out.push(adaptive::Value::number(key.cost));
                adaptive::Value ids = adaptive::Value::array();
                ids.push(adaptive::Value::integer(static_cast<long long>(key.pair.first)));
                ids.push(adaptive::Value::integer(static_cast<long long>(key.pair.second)));
                out.push(ids);
                got = adaptive::json_dumps(out, false, false, true);
            } else {
                continue;
            }
        } catch (const adaptive::Error& error) {
            got = error.type();
            message = error.what();
        }
        const std::string expected = column.size() > at ? column[at] : std::string();
        const std::string expected_message = column.size() > at+1 ? column[at+1] : std::string();
        // patches() reports an amplification built out of log10, so this comparison
        // is one of the places the numeric environment matters: icpx answers one
        // step away from the glibc the expectations were recorded with.
        if (!got.empty() && got[0] == '[' && !expected.empty() && expected[0] == '[') {
            const std::string found = testing::differs(adaptive::read_json(got),
                                                       adaptive::read_json(expected), strict);
            if (!found.empty())
                check::fail(__FILE__, __LINE__, "patches "+found);
        } else {
            CHECK_EQ(got, expected);
        }
        CHECK_EQ(message, expected_message);
        ++rows;
    }
    CHECK(rows >= (argc > 1 ? 1 : 638));
    std::cout << "compared " << rows << " merge decisions\n";
    return check::report("test_merge");
}
