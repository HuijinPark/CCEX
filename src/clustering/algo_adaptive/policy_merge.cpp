// Policy "merge" (the only one since 2026-10-10): probe pairs of bundles and take, one
// step at a time, the merge that lowers a local error the most.
#include "policy_merge.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <dirent.h>
#include <fnmatch.h>
#include <map>
#include <set>
#include <sstream>
#include <sys/stat.h>
#include "base_bundles.h"
#include "base_numeric.h"
#include "base_json.h"
#include "base_run.h"

// ----------------------------------------------------------------------------
// merge decisions: which patches, which pairs, in what order
// ----------------------------------------------------------------------------
namespace adaptive {
bool affordable_clusters(const Plan& plan, const std::vector<Sites>& clusters,
                         const std::map<Sites, Curves>& held, long long states,
                         long long steps, long long work_left) {
    std::vector<Sites> missing;
    for (const Sites& sites : clusters)
        if (held.count(sites) == 0)
            missing.push_back(sites);
    if (missing.empty())
        return true;
    return cluster_work(plan, missing, states, steps) <= work_left;
}

long long bundle_cap(const MergeContext& context, const Plan&) {
    return context.bundle_max_dim;
}

std::pair<std::vector<Sites>, std::vector<Edge>> merge_bundles(
    const std::vector<Sites>& bundles, const std::vector<Edge>& edges,
    const std::vector<int>& members) {
    const std::set<int> joined(members.begin(), members.end());
    const Sites merged = union_of(bundles, std::vector<int>(joined.begin(), joined.end()));
    std::vector<Sites> out;
    for (std::size_t i = 0; i < bundles.size(); ++i)
        if (joined.count(static_cast<int>(i)) == 0)
            out.push_back(bundles[i]);
    out.push_back(merged);
    std::sort(out.begin(), out.end());
    std::map<int, int> lookup;
    for (std::size_t i = 0; i < out.size(); ++i)
        for (int site : out[i])
            lookup[site] = static_cast<int>(i);
    std::set<Edge> mapped;
    for (const Edge& edge : edges) {
        const int x = lookup.at(bundles[static_cast<std::size_t>(edge.first)].at(0));
        const int y = lookup.at(bundles[static_cast<std::size_t>(edge.second)].at(0));
        if (x != y)
            mapped.insert(Edge(std::min(x, y), std::max(x, y)));
    }
    return {out, std::vector<Edge>(mapped.begin(), mapped.end())};
}

bool can_merge(Feasibility* cache, const MergeContext& context, const Plan& plan,
               const std::vector<int>& members) {
    // Only the bundle cap decides. A cluster is `order` bundles and its
    // dimension follows from that cap and the order, so neither the neighbours
    // nor the dimension query is consulted (terminology 4절, 2026-09-18).
    if (cache->plan != &plan) {
        cache->plan = &plan;
        cache->dimensions.clear();
        for (const Sites& bundle : plan.bundles)
            cache->dimensions.push_back(bundle_dimension(bundle, plan.site_dimensions));
    }
    // Compared one factor at a time rather than multiplied out: the answer is
    // the same, and a product that would not fit never happens.
    const long long cap = bundle_cap(context, plan);
    long long product = 1;
    for (int i : members) {
        const long long dimension = cache->dimensions[static_cast<std::size_t>(i)];
        if (product > cap/dimension)
            return false;
        product *= dimension;
    }
    return product <= cap;
}

std::pair<std::vector<std::vector<int>>, Coverage> patches(const MergeContext& context,
                                                           const Plan& plan) {
    const std::size_t n = plan.bundles.size();
    // dict(((a, b), w)) keeps the first occurrence's position but the last value. The
    // map finds an entry's position: a scan of the list per edge was quadratic, and
    // with 3,019 bundles this function took 99 % of the policy's time (A7, job 627).
    std::vector<std::pair<std::pair<int, int>, double>> weights;
    std::map<std::pair<int, int>, std::size_t> weight_at;
    for (const WeightedEdge& edge : plan.weights()) {
        const std::pair<int, int> key(edge.a, edge.b);
        const std::map<std::pair<int, int>, std::size_t>::const_iterator seen = weight_at.find(key);
        if (seen != weight_at.end()) {
            weights[seen->second].second = edge.weight;
        } else {
            weight_at.emplace(key, weights.size());
            weights.emplace_back(key, edge.weight);
        }
    }
    std::map<int, long long> degree;
    for (const WeightedEdge& edge : plan.edges) {
        degree[edge.a] += 1;
        degree[edge.b] += 1;
    }
    std::vector<double> risks;
    for (std::size_t i = 0; i < n; ++i) {
        double risk = 0.;
        if (context.order == 2 && degree[static_cast<int>(i)] > 1) {
            const std::map<Sites, Curves>::const_iterator bundle =
                context.raw->find(plan.bundles[i]);
            const std::map<Sites, Curves>::const_iterator zero = context.raw->find(Sites());
            if (bundle == context.raw->end() || zero == context.raw->end())
                throw KeyError("raw probe missing for the risk estimate");
            for (std::size_t state = 0; state < bundle->second.size(); ++state) {
                for (std::size_t step = 0; step < bundle->second[state].size(); ++step) {
                    const std::complex<double> z = curve_at(bundle->second, state, step);
                    const std::complex<double> z0 = curve_at(zero->second, state, step);
                    if (!finite(z) || !finite(z0) || z == 0. || z0 == 0.) {
                        risk = HUGE_VAL;
                        continue;
                    }
                    const double amplification =
                        std::log10(complex_abs(z0.real(), z0.imag()))-
                        std::log10(complex_abs(z.real(), z.imag()));
                    risk = py_max(risk, static_cast<double>(degree[static_cast<int>(i)]-1)*
                                        py_max(0., amplification));
                }
            }
        }
        risks.push_back(risk);
    }
    std::vector<std::vector<std::pair<int, double>>> neighbours(n);
    for (const std::pair<std::pair<int, int>, double>& entry : weights) {
        neighbours[static_cast<std::size_t>(entry.first.first)].emplace_back(entry.first.second,
                                                                             entry.second);
        neighbours[static_cast<std::size_t>(entry.first.second)].emplace_back(entry.first.first,
                                                                              entry.second);
    }
    std::vector<std::pair<std::vector<int>, double>> candidates;   // insertion order, like a dict
    const std::vector<std::pair<std::vector<int>, double>>::iterator none = candidates.end();
    (void)none;
    std::map<std::vector<int>, std::size_t> candidate_at;            // key -> position, as above
    auto put = [&candidates, &candidate_at](const std::vector<int>& key, double value, bool take_max) {
        const std::map<std::vector<int>, std::size_t>::const_iterator seen = candidate_at.find(key);
        if (seen != candidate_at.end()) {
            std::pair<std::vector<int>, double>& entry = candidates[seen->second];
            entry.second = take_max ? py_max(entry.second, value) : value;
            return;
        }
        candidate_at.emplace(key, candidates.size());
        candidates.emplace_back(key, value);
    };
    bool restricted = false;
    if (context.order == 1) {
        // Order 1 drops the links from the expansion, so plain pairs are the
        // candidates - all of them when the partition is small enough to afford it.
        if (static_cast<long long>(n) <= context.neighbors) {
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = i+1; j < n; ++j)
                    put({static_cast<int>(i), static_cast<int>(j)}, 0., false);
        } else {
            restricted = true;
            for (const std::pair<std::pair<int, int>, double>& entry : weights)
                put({entry.first.first, entry.first.second}, entry.second, false);
        }
    } else if (context.order > 2) {
        // Above order 2 a triple is not the right shape: the domain that has to
        // be checked is one bundle wider than the expansion, and only a bounded
        // beam of those is affordable. Coverage is partial by construction.
        restricted = true;
        for (const std::pair<const std::vector<int>, double>& entry :
                 next_order_patches(plan, context.order, context.neighbors,
                                    context.max_patches))
            put(entry.first, entry.second, false);
    } else {
        for (std::size_t i = 0; i < n; ++i) {
            std::vector<std::pair<int, double>>& adjacent = neighbours[i];
            // key=(-weight, index), and Python's sort is stable.
            std::stable_sort(adjacent.begin(), adjacent.end(),
                             [](const std::pair<int, double>& x, const std::pair<int, double>& y) {
                                 if (x.second != y.second) return x.second > y.second;
                                 return x.first < y.first;
                             });
            if (static_cast<long long>(adjacent.size()) > context.neighbors)
                restricted = true;
            const std::size_t take = std::min<std::size_t>(
                adjacent.size(), static_cast<std::size_t>(std::max<long long>(context.neighbors, 0)));
            for (std::size_t x = 0; x < take; ++x)
                for (std::size_t y = x+1; y < take; ++y) {
                    std::vector<int> key{static_cast<int>(i), adjacent[x].first, adjacent[y].first};
                    std::sort(key.begin(), key.end());
                    put(key, adjacent[x].second+adjacent[y].second, true);
                }
        }
    }
    std::vector<int> full;
    for (std::size_t i = 0; i < n; ++i)
        full.push_back(static_cast<int>(i));
    const bool full_fits = n <= 8 &&
        context.affordable(plan, std::vector<Sites>{union_of(plan.bundles, full)});
    if (full_fits)
        put(full, HUGE_VAL, false);
    std::vector<std::pair<std::vector<int>, double>> ranked = candidates;
    std::stable_sort(ranked.begin(), ranked.end(),
        [&](const std::pair<std::vector<int>, double>& x,
            const std::pair<std::vector<int>, double>& y) {
            const int rank_x = full_fits && x.first == full ? 0 : 1;
            const int rank_y = full_fits && y.first == full ? 0 : 1;
            if (rank_x != rank_y) return rank_x < rank_y;
            double risk_x = 0., risk_y = 0.;
            for (int i : x.first) risk_x = py_max(risk_x, risks[static_cast<std::size_t>(i)]);
            for (int i : y.first) risk_y = py_max(risk_y, risks[static_cast<std::size_t>(i)]);
            if (risk_x != risk_y) return risk_x > risk_y;
            if (x.second != y.second) return x.second > y.second;
            return x.first < y.first;
        });
    if (static_cast<long long>(ranked.size()) > context.max_patches)
        restricted = true;
    Coverage coverage;
    coverage.candidate_patches = static_cast<long long>(candidates.size());
    std::vector<std::vector<int>> out;
    long long blocked = 0;
    const std::size_t limit = std::min<std::size_t>(
        ranked.size(), static_cast<std::size_t>(std::max<long long>(context.max_patches, 0)));
    for (std::size_t i = 0; i < limit; ++i) {
        if (context.affordable(plan, std::vector<Sites>{union_of(plan.bundles, ranked[i].first)}))
            out.push_back(ranked[i].first);
        else
            ++blocked;
    }
    coverage.checked_patches = static_cast<long long>(out.size());
    coverage.work_blocked = blocked;
    coverage.coverage_limited = restricted || blocked > 0;
    coverage.whole_bath_oracle = full_fits;
    double worst = 0.;
    for (double risk : risks)
        worst = py_max(worst, risk);
    coverage.max_inverse_log10_amplification = risks.empty() ? 0. : worst;
    return {out, coverage};
}

Curves local_model(const MergeContext& context, const Plan& plan, const std::vector<int>& patch,
                   const std::pair<std::vector<Sites>, std::vector<Edge>>* changed) {
    std::vector<Sites> bundles;
    std::map<int, int> labels;
    for (std::size_t i = 0; i < patch.size(); ++i) {
        bundles.push_back(plan.bundles[static_cast<std::size_t>(patch[i])]);
        labels[patch[i]] = static_cast<int>(i);
    }
    std::vector<Edge> edges;
    for (const WeightedEdge& edge : plan.edges)
        if (labels.count(edge.a) != 0 && labels.count(edge.b) != 0)
            edges.emplace_back(labels[edge.a], labels[edge.b]);
    if (changed != nullptr) {
        bundles = changed->first;
        edges = changed->second;
    }
    return approximate(bundles, edges, context.order, *context.raw, context.gcce, context.limits);
}

Discrepancy local(const MergeContext& context, const Plan& plan, const std::vector<int>& patch,
                  const std::pair<std::vector<Sites>, std::vector<Edge>>* changed) {
    const Curves approximated = local_model(context, plan, patch, changed);
    const std::vector<Sites> whole{union_of(plan.bundles, patch)};
    const Curves exact = approximate(whole, {}, 1, *context.raw, context.gcce, context.limits);
    return discrepancy(approximated, exact, context.atol, context.rtol);
}

// The termination check's test, applied to one candidate instead of a finished
// run: every point finite, and for cCCE no larger than 1+atol. gCCE is exempt in
// both places.
bool is_physical(const MergeContext& context, const Curves& curves) {
    for (const Curve& curve : curves)
        for (const std::pair<double, double>& z : curve) {
            if (!std::isfinite(z.first) || !std::isfinite(z.second))
                return false;
            // abs() of a complex is hypot, which is libm - one of the places the
            // verified-numeric-environment contract applies (0B-2).
            if (!context.gcce && !(complex_abs(z.first, z.second) <= 1.+context.atol))
                return false;
        }
    return true;
}

ProposalKey proposal_key(const MergeContext& context, double score, double local_before,
                         std::pair<int, int> pair) {
    ProposalKey key;
    key.cost = -score;
    key.pair = pair;
    if (context.merge_priority == "worst-error") {
        key.has_worst_error = true;
        key.worst_error = -local_before;
    }
    return key;
}

bool before(const ProposalKey& left, const ProposalKey& right) {
    if (left.has_worst_error && left.worst_error != right.worst_error)
        return left.worst_error < right.worst_error;
    if (left.cost != right.cost)
        return left.cost < right.cost;
    return left.pair < right.pair;
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// the Driver: steps, probes, requests and the report
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {
// Two curve sets, compared the way Python compares lists of complex numbers:
// 0.0 == -0.0, and a NaN is equal to nothing, not even itself.
bool same(const Curves& left, const Curves& right) {
    if (left.size() != right.size())
        return false;
    for (std::size_t i = 0; i < left.size(); ++i) {
        if (left[i].size() != right[i].size())
            return false;
        for (std::size_t k = 0; k < left[i].size(); ++k)
            if (!(left[i][k].first == right[i][k].first && left[i][k].second == right[i][k].second))
                return false;
    }
    return true;
}
}   // namespace

namespace {
std::string four_digits(long long value) {
    char text[32];
    std::snprintf(text, sizeof text, "%04lld", value);
    return text;
}

// mkdir -p, for the step directory. An existing directory is accepted at every level, the last one included
// (Python's Path.mkdir() would refuse the last; this does not - the comment said it did until 2026-10-08, review of
// the duplicates: base_run's make_new_directory, removed with the scan policies on 2026-10-10, was the one that refused).
void make_directory(const std::string& path) {
    std::string built;
    for (std::size_t at = 0; at <= path.size(); ) {
        const std::size_t slash = path.find('/', at);
        const std::size_t end = slash == std::string::npos ? path.size() : slash;
        built = path.substr(0, end);
        if (!built.empty() && mkdir(built.c_str(), 0777) != 0 && errno != EEXIST)
            throw os_error(errno, built);
        if (slash == std::string::npos)
            break;
        at = slash+1;
    }
}

bool file_exists(const std::string& path) {
    struct stat info;
    return stat(path.c_str(), &info) == 0 && S_ISREG(info.st_mode);
}

void write_file(const std::string& path, const std::string& text) {
    std::FILE* file = std::fopen(path.c_str(), "wb");
    if (file == nullptr)
        throw os_error(errno, path);
    const bool failed = std::fwrite(text.data(), 1, text.size(), file) != text.size();
    std::fclose(file);
    if (failed)
        throw os_error(EIO, path);
}

}   // namespace

bool Driver::affordable(const Plan& plan, const std::vector<Sites>& clusters) const {
    return affordable_clusters(plan, clusters, raw_, static_cast<long long>(states_.size()),
                               cfg_.at("nstep").as_int(), work_left());
}

long long Driver::work_left() const {
    return static_cast<long long>(static_cast<long long>(settings_.max_total_work));
}

void Driver::prepare(const std::string& root, NativeSettings native, Clock* clock) {
    root_ = root;
    native_ = std::move(native);
    clock_ = clock;
    // The report's keys are laid down by begin(), in Python's order. Creating one
    // here would put it first and every later set() would keep it there.
}

std::unique_ptr<Executed> Driver::execute(const std::string& mode,
                                          const std::vector<Sites>* bundles,
                                          const std::vector<Sites>* probes, bool allow_failure,
                                          const Value* more) {
    ++counter_;
    const std::string here = root_+"/"+four_digits(counter_)+"-"+mode;
    make_directory(here);
    Value request = Value::object();
    request.set("mode", Value::text(mode));
    request.set("response", Value::text(here+"/response"));
    if (native_.physics_metadata && !has_physics_)
        request.set("physics", Value::boolean(true));
    if (bundles != nullptr)
        request.set("bundles", sites_list(*bundles));
    if (probes != nullptr)
        request.set("probes", sites_list(*probes));
    if (more != nullptr)
        for (const std::string& key : more->keys())
            request.set(key, more->at(key));
    const std::string request_text = write_json(request, 2);
    write_file(here+"/request.json", request_text);
    Value config = cfg_;
    config.set("adaptive_work_file", Value::text(here+"/request.json"));
    config.set("outfile", Value::text(here+"/coherence"));
    const std::string input_text = write_json(config, 2);
    write_file(here+"/input.json", input_text);

    std::vector<std::string> command{native_.mpiexec, "-n", std::to_string(settings_.ranks),
                                     native_.binary, "-f", here+"/input.json"};
    const double started = clock_->monotonic();
    const NativeResult answer = run_native(command, settings_.cwd, here+"/run.log", native_.timeout);
    const std::string& response_root = here;
    const double wall = clock_->monotonic()-started;

    Value record = Value::object();
    record.set("directory", Value::text(here));
    record.set("mode", Value::text(mode));
    record.set("wall_s", Value::number(wall));
    record.set("returncode", Value::integer(static_cast<long long>(answer.returncode)));
    record.set("timeout", Value::boolean(answer.timed_out));
    // Value hands out const references, so the list is taken, added to and put
    // back. The report is written once per step; nothing here is hot.
    Value runs = report_.at("runs");
    runs.push(record);
    report_.set("runs", runs);
    if (answer.timed_out || answer.returncode != 0) {
        if (allow_failure && !answer.timed_out)
            return nullptr;
        throw RuntimeError("CCEX "+mode+" failed; inspect "+here+"/run.log");
    }

    std::unique_ptr<Executed> out(new Executed());
    out->directory = response_root;
    out->raw_plan = read_json(read_file(response_root+"/response.plan.json"));
    if (file_exists(response_root+"/response.physics.json")) {
        physics_ = read_json(read_file(response_root+"/response.physics.json"));
        has_physics_ = true;
    }
    if (has_physics_)
        out->raw_plan.set("physics", physics_);
    if (out->raw_plan.has("dimension_limits"))
        report_.set("dimension_limits", out->raw_plan.at("dimension_limits"));
    out->plan = plan_of(out->raw_plan);
    if (bundles != nullptr && out->plan.bundles != *bundles)
        throw RuntimeError("requested partition changed");
    return out;
}

void Driver::use_cache(std::unique_ptr<ProbeCache> cache, Value context) {
    cache_path_ = cache->path();
    cache_key_ = cache->key();
    cache_ = std::move(cache);
    cache_context_ = std::move(context);
    // Python puts this into the report when the cache is opened, which is right
    // after input_file_sha256 - and a key's position is set the first time it is
    // written, so it has to be written here rather than at the first persist().
    Value block = Value::object();
    block.set("path", Value::text(cache_path_));
    block.set("context_sha256", Value::text(cache_key_));
    block.set("loaded_clusters", Value::integer(static_cast<long long>(0)));
    block.set("stored_clusters", Value::integer(static_cast<long long>(0)));
    report_.set("raw_cache", block);
}

bool Driver::load_cached(const std::set<Sites>& wanted) {
    if (cache_ == nullptr)
        return true;
    const long long steps = static_cast<long long>(cfg_.at("nstep").as_int());
    // The cache is only trusted while the inputs it was keyed on are unchanged,
    // and the check happens here, not at start-up.
    if (manifest_text(input_manifest(cfg_, settings_.cwd, states_)) !=
        manifest_text(cache_context_.at("files")))
        throw RuntimeError("an input file changed before raw cache access");
    std::vector<Sites> ask;
    std::set<Sites> with_zero = wanted;
    with_zero.insert(Sites());
    for (const Sites& sites : with_zero)
        if (raw_.count(sites) == 0)
            ask.push_back(sites);
    const std::map<Sites, Curves> found = cache_->get(ask);
    std::set<Sites> total(wanted.begin(), wanted.end());
    for (const std::pair<const Sites, Curves>& entry : raw_)
        total.insert(entry.first);
    for (const std::pair<const Sites, Curves>& entry : found)
        total.insert(entry.first);
    const long long count = static_cast<long long>(total.size());
    if (count > settings_.max_probes ||
        count*static_cast<long long>(states_.size())*steps*64 >
            settings_.cache_mb*1024*1024)
        return false;
    for (const std::pair<const Sites, Curves>& entry : found)
        raw_[entry.first] = entry.second;
    loaded_ += static_cast<long long>(found.size());
    return true;
}

bool Driver::ensure(const std::vector<Sites>& needed, const std::vector<Sites>& bundles) {
    std::set<Sites> wanted(needed.begin(), needed.end());
    const long long steps = static_cast<long long>(cfg_.at("nstep").as_int());
    if (!load_cached(wanted))
        return false;
    std::vector<Sites> missing;
    for (const Sites& sites : wanted)
        if (raw_.count(sites) == 0)
            missing.push_back(sites);
    std::sort(missing.begin(), missing.end());
    if (missing.empty())
        return true;
    const long long held = static_cast<long long>(raw_.size());
    if (held+static_cast<long long>(missing.size()) > settings_.max_probes)
        return false;
    const long long estimated = (held+static_cast<long long>(missing.size()))*
                                static_cast<long long>(states_.size())*steps*64;
    if (estimated > settings_.cache_mb*1024*1024)
        return false;
    const long long batch = settings_.probe_batch_clusters > 0 ? settings_.probe_batch_clusters
                                                               : static_cast<long long>(missing.size());
    for (std::size_t start = 0; start < missing.size();
         start += static_cast<std::size_t>(batch)) {
        const std::size_t end = std::min(missing.size(), start+static_cast<std::size_t>(batch));
        read_probes(std::vector<Sites>(missing.begin()+static_cast<long>(start),
                                       missing.begin()+static_cast<long>(end)), bundles);
    }
    return true;
}

std::unique_ptr<Executed> Driver::step_for(const std::string& mode,
                                           const std::vector<Sites>* bundles,
                                           const std::vector<Sites>* probes) {
    if (execute_)
        return std::unique_ptr<Executed>(new Executed(execute_(mode, bundles, probes)));
    return execute(mode, bundles, probes, false);
}

std::map<Sites, Curves> read_probe_files(const std::string& directory, const std::vector<int>& states,
                                         long long steps, long long ranks, const std::vector<Sites>& missing) {
    std::map<Sites, Curves> pending;
    std::map<Sites, std::vector<std::vector<bool>>> seen;
    std::vector<Sites> keys;
    keys.push_back(Sites());
    keys.insert(keys.end(), missing.begin(), missing.end());
    for (const Sites& key : keys) {
        pending[key] = Curves(states.size(),
                              Curve(static_cast<std::size_t>(steps), std::pair<double, double>(0., 0.)));
        seen[key] = std::vector<std::vector<bool>>(
            states.size(), std::vector<bool>(static_cast<std::size_t>(steps), false));
    }
    for (std::size_t si = 0; si < states.size(); ++si) {
        // sorted(d.glob('response.rank-*.state-{state}.tsv')): whatever is there,
        // in name order. Trying rank numbers up to ranks+16 found the same files
        // only while the launcher numbered them from zero without gaps.
        std::vector<std::string> files;
        const std::string pattern = "response.rank-*.state-"+std::to_string(states[si])+".tsv";
        if (DIR* listing = opendir(directory.c_str())) {
            while (const dirent* entry = readdir(listing))
                if (fnmatch(pattern.c_str(), entry->d_name, 0) == 0)
                    files.push_back(directory+"/"+entry->d_name);
            closedir(listing);
        }
        std::sort(files.begin(), files.end());
        if (static_cast<long long>(files.size()) != ranks)
            throw RuntimeError("missing MPI probe response");
        for (const std::string& path : files) {
            const std::string text = read_file(path);   // a read error is an OSError (ferror), not a short text
            std::istringstream rows(text);
            std::string line;
            std::getline(rows, line);                       // the header
            while (std::getline(rows, line)) {
                if (line.empty())
                    continue;
                std::istringstream fields(line);
                std::string probe_text, step_text, real_text, imag_text;
                std::getline(fields, probe_text, '\t');
                std::getline(fields, step_text, '\t');
                std::getline(fields, real_text, '\t');
                std::getline(fields, imag_text, '\t');
                const long long index = static_cast<long long>(parse_int(probe_text));
                const long long step_index = static_cast<long long>(parse_int(step_text));
                // Both are list subscripts in Python: counted from the end when
                // negative, IndexError past either end. Written straight into the
                // vector, a step the response should not have had was a write
                // outside it.
                const auto subscript = [](long long i, std::size_t size) {
                    const long long count = static_cast<long long>(size);
                    const long long at = i < 0 ? i+count : i;
                    if (at < 0 || at >= count)
                        throw IndexError("list index out of range");
                    return static_cast<std::size_t>(at);
                };
                const Sites key = index == -1 ? Sites() : missing[subscript(index, missing.size())];
                const std::size_t which = subscript(step_index, seen[key][si].size());
                if (seen[key][si][which])
                    throw RuntimeError("duplicate raw probe");
                seen[key][si][which] = true;
                pending[key][si][which] =
                    std::pair<double, double>(parse_float(real_text), parse_float(imag_text));
            }
        }
    }
    for (const std::pair<const Sites, Curves>& entry : pending)
        for (const std::vector<bool>& state : seen.at(entry.first))
            for (bool filled : state)
                if (!filled)
                    throw RuntimeError("incomplete raw probe");
    return pending;
}

void Driver::read_probes(const std::vector<Sites>& missing, const std::vector<Sites>& bundles) {
    const std::unique_ptr<Executed> found = step_for("probe", &bundles, &missing);
    const Executed& step = *found;
    const long long steps = static_cast<long long>(cfg_.at("nstep").as_int());
    const std::map<Sites, Curves> pending = read_probe_files(step.directory, states_, steps, settings_.ranks, missing);
    for (const std::pair<const Sites, Curves>& entry : pending) {
        const std::map<Sites, Curves>::const_iterator held = raw_.find(entry.first);
        if (held != raw_.end() && !same(entry.second, held->second))
            throw RuntimeError("L0 changed between identical-state runs");
    }
    for (const std::pair<const Sites, Curves>& entry : pending)
        raw_[entry.first] = entry.second;
    if (cache_ != nullptr) {
        if (manifest_text(input_manifest(cfg_, settings_.cwd, states_)) !=
            manifest_text(cache_context_.at("files")))
            throw RuntimeError("an input file changed during a raw probe");
        cache_->put(pending);
        stored_ += static_cast<long long>(missing.size());
    }
}
}   // namespace adaptive

namespace adaptive {
namespace {
Value coverage_value(const Coverage& coverage) {
    Value out = Value::object();
    out.set("candidate_patches", Value::integer(static_cast<long long>(coverage.candidate_patches)));
    out.set("checked_patches", Value::integer(static_cast<long long>(coverage.checked_patches)));
    out.set("work_blocked", Value::integer(static_cast<long long>(coverage.work_blocked)));
    out.set("coverage_limited", Value::boolean(coverage.coverage_limited));
    out.set("whole_bath_oracle", Value::boolean(coverage.whole_bath_oracle));
    out.set("max_inverse_log10_amplification",
            Value::number(coverage.max_inverse_log10_amplification));
    return out;
}

Value ids_list(const std::vector<int>& listed) {
    Value out = Value::array();
    for (int id : listed)
        out.push(Value::integer(static_cast<long long>(id)));
    return out;
}
}   // namespace

void Driver::begin(const Value& settings, const std::string& binary_sha256) {
    report_.set("status", Value::text("running"));
    report_.set("input", cfg_);
    report_.set("settings", settings);
    report_.set("binary_sha256", Value::text(binary_sha256));
    report_.set("runs", Value::array());
    report_.set("rounds", Value::array());
    Value limitations = Value::array();
    limitations.push(Value::text("Local oracle errors and successive curves are checks, "
                                 "not rigorous global error bounds."));
    limitations.push(Value::text("A greedy merge search cannot undo a poor initial partition."));
    limitations.push(Value::text("One fixed partition and all requested time points/states "
                                 "are used per input."));
    report_.set("limitations", limitations);
    // The same walk input_manifest does, except a file that is not there is left
    // out rather than recorded as absent.
    const Value files = input_manifest(cfg_, settings_.cwd, states_);
    inputs_ = Value::object();
    for (const std::string& name : files.keys())
        if (files.at(name).kind() != Value::Kind::Null)
            inputs_.set(name, files.at(name));
    report_.set("input_file_sha256", inputs_);
}

void Driver::persist() const {
    Value report = report_;
    // Python writes report['input']=self.cfg - the config object itself, not a
    // copy - so the input block is cfg_ as it is at this write (merge does not
    // change it after begin()). The key keeps the position begin() gave it.
    report.set("input", cfg_);
    Value rounds = Value::array();
    for (const Value& row : rounds_)
        rounds.push(row);
    report.set("rounds", rounds);
    if (!cache_path_.empty()) {
        Value cache = Value::object();
        cache.set("path", Value::text(cache_path_));
        cache.set("context_sha256", Value::text(cache_key_));
        cache.set("loaded_clusters", Value::integer(static_cast<long long>(loaded_)));
        cache.set("stored_clusters", Value::integer(static_cast<long long>(stored_)));
        report.set("raw_cache", cache);
    }
    std::FILE* file = std::fopen((root_+"/summary.json").c_str(), "wb");
    if (file == nullptr)
        throw os_error(errno, root_+"/summary.json");
    const std::string text = write_json(report, 2);
    const bool failed = std::fwrite(text.data(), 1, text.size(), file) != text.size();
    std::fclose(file);
    if (failed)
        throw os_error(EIO, root_+"/summary.json");
}

int Driver::finish(const std::string& status, const std::string& directory, const Executed& step) {
    report_.set("status", Value::text(status));
    report_.set("final_directory", Value::text(directory));
    report_.set("final_bundles", sites_list(step.plan.bundles));
    report_.set("cached_probe_clusters", Value::integer(static_cast<long long>(static_cast<long long>(raw_.size()))));
    long long nonfinite = 0, zeros = 0;
    for (const std::pair<const Sites, Curves>& entry : raw_)
        for (const Curve& curve : entry.second)
            for (const std::pair<double, double>& z : curve) {
                if (!std::isfinite(z.first) || !std::isfinite(z.second))
                    ++nonfinite;
                if (z.first == 0. && z.second == 0.)
                    ++zeros;
            }
    report_.set("raw_nonfinite_values", Value::integer(static_cast<long long>(nonfinite)));
    report_.set("raw_zero_values", Value::integer(static_cast<long long>(zeros)));
    for (const std::string& name : inputs_.keys())
        if (sha256_file(name) != inputs_.at(name).as_string())
            throw RuntimeError("an input file changed during refinement");
    write_file(root_+"/final-bundles.json", write_json(step.raw_plan, 2));
    // A rerun keeps the diagnostic rounds and uses the partition that was chosen.
    Value request = Value::object();
    request.set("mode", Value::text("run"));
    request.set("response", Value::text(root_+"/final-response"));
    request.set("bundles", sites_list(step.plan.bundles));
    write_file(root_+"/final-request.json", write_json(request, 2));
    Value config = cfg_;
    config.set("adaptive_work_file", Value::text(root_+"/final-request.json"));
    config.set("outfile", Value::text(root_+"/final-coherence"));
    write_file(root_+"/final-input.json", write_json(config, 2));
    persist();
    std::printf("%s %s/summary.json\n", status.c_str(), root_.c_str());
    std::fflush(stdout);
    // "physical" is regroup's only success (adaptive-regroup-plan.md item 10), and is not
    // called by either of merge's names: a physical curve is not a validated one.
    return status == "validated_full_bath_model" || status == "passed_convergence_checks" ||
           status == "physical" ? 0 : 2;
}
}   // namespace adaptive

namespace adaptive {
// Every field at once. Python keeps the response as a dict and only fails when
// something reads a field that is not there, so a malformed response stops here
// rather than later - which is the port's own behaviour, not Python's.
Plan plan_of(const Value& raw) {
    Plan plan;
    for (std::size_t i = 0; i < raw.at("bundles").size(); ++i)
        plan.bundles.push_back(sites_of(raw.at("bundles").at(i)));
    for (std::size_t i = 0; i < raw.at("edges").size(); ++i) {
        const Value& edge = raw.at("edges").at(i);
        plan.edges.push_back({static_cast<int>(edge.at(0).as_int()),
                              static_cast<int>(edge.at(1).as_int()),
                              edge.at(2).as_double()});
    }
    if (raw.has("candidate_edges")) {
        plan.has_candidate_edges = true;
        for (std::size_t i = 0; i < raw.at("candidate_edges").size(); ++i) {
            const Value& edge = raw.at("candidate_edges").at(i);
            plan.candidate_edges.push_back({static_cast<int>(edge.at(0).as_int()),
                                            static_cast<int>(edge.at(1).as_int()),
                                            edge.at(2).as_double()});
        }
    }
    for (std::size_t i = 0; i < raw.at("site_dimensions").size(); ++i)
        plan.site_dimensions.push_back(raw.at("site_dimensions").at(i).as_int());
    plan.qubit_dimension = raw.at("qubit_dimension").as_int();
    if (raw.has("expansion_limits")) {
        const Value& limits = raw.at("expansion_limits");
        plan.has_limits = true;
        if (limits.has("max_clusters"))
            plan.limits.max_clusters = limits.at("max_clusters").as_int();
        if (limits.has("max_work"))
            plan.limits.max_work = limits.at("max_work").as_int();
    }
    return plan;
}

Curves read_state_curves(const std::string& directory, const std::vector<int>& states) {
    Curves out;
    for (int state : states) {
        const std::string path = directory+"/response.curve.state-"+std::to_string(state)+".tsv";
        const std::string text = read_file(path);       // a read error is an OSError (ferror), not a short text
        std::istringstream rows(text);
        std::string line;
        std::getline(rows, line);                       // header: which column is which
        std::vector<std::string> header;
        {
            std::istringstream names(line);
            std::string name;
            while (std::getline(names, name, '\t'))
                header.push_back(name);
        }
        std::size_t real_at = 0, imag_at = 1;
        for (std::size_t i = 0; i < header.size(); ++i) {
            if (header[i] == "real") real_at = i;
            if (header[i] == "imag") imag_at = i;
        }
        Curve curve;
        while (std::getline(rows, line)) {
            if (line.empty())
                continue;
            std::vector<std::string> fields;
            std::istringstream parts(line);
            std::string field;
            while (std::getline(parts, field, '\t'))
                fields.push_back(field);
            curve.emplace_back(parse_float(fields.at(real_at)), parse_float(fields.at(imag_at)));
        }
        out.push_back(curve);
    }
    return out;
}

namespace {
// Every non-empty subset of a patch, as the sites it covers.
std::vector<Sites> patch_clusters(const std::vector<Sites>& bundles, const std::vector<int>& patch) {
    std::vector<Sites> out;
    for (std::size_t mask = 1; mask < (std::size_t(1) << patch.size()); ++mask) {
        std::vector<int> ids;
        for (std::size_t i = 0; i < patch.size(); ++i)
            if ((mask >> i) & 1)
                ids.push_back(patch[i]);
        out.push_back(union_of(bundles, ids));
    }
    return out;
}
}   // namespace

ExpansionLimits Driver::config_limits() const {
    ExpansionLimits limits;
    if (cfg_.has("expansion_max_clusters"))
        limits.max_clusters = cfg_.at("expansion_max_clusters").as_int();
    if (cfg_.has("expansion_max_work"))
        limits.max_work = cfg_.at("expansion_max_work").as_int();
    return limits;
}

double improvement(double before, double after) {
    // The gain a merge buys, with the convention an infinite before needs.
    if (std::isfinite(before))
        return before-after;
    return std::isfinite(after) ? HUGE_VAL : 0.;
}

std::pair<std::vector<int>, Value> Driver::hub_candidate(
        const MergeContext& context, const Plan& plan,
        const std::map<std::pair<int, int>, Value>& candidates) {
    // Three things about this are measurements, not taste
    // (docs/archive/adaptive-physical-gate-20260917.md, 6절). It grows past the bundle
    // cap, because the cap is a cost limit and this is an accuracy device - a
    // capped hub stopped at two or three spins and came out worse. It grows
    // until the curve is physical, not until the guilty partners run out - that
    // was the third wrong question, and it left a big bundle that cured nothing.
    // And it is never ranked against the others by cost; the caller adopts it
    // whenever CCEX can plan it.
    //
    // What stops the growth is the probe budget, which is where the measured
    // version stopped too. CCEX still applies its own bundle cap to the plan, so
    // a hub that needs more dimension than the input allows is offered and
    // refused, and the round then goes on as it would have without the flag.
    std::vector<std::pair<std::pair<int, int>, double> > bad;    // pair, gain
    std::map<int, double> best;
    for (const std::pair<const std::pair<int, int>, Value>& entry : candidates) {
        const double gain = improvement(entry.second.at("local_before").as_double(),
                                        entry.second.at("local_after").as_double());
        // Python sorts NaN wherever it lands; an inconsistent comparator is
        // undefined behaviour here, so it is refused instead (0B-2 규약).
        if (std::isnan(gain))
            throw ValueError("a merge gain is not a number");
        best[entry.first.first] = std::max(best.count(entry.first.first) != 0
                                           ? best[entry.first.first] : 0., gain);
        best[entry.first.second] = std::max(best.count(entry.first.second) != 0
                                            ? best[entry.first.second] : 0., gain);
        if (entry.second.has("physical") && !entry.second.at("physical").as_bool())
            bad.push_back({entry.first, gain});
    }
    if (bad.size() < 2)
        return {};
    std::map<int, long long> counts;
    for (const std::pair<std::pair<int, int>, double>& one : bad) {
        ++counts[one.first.first];
        ++counts[one.first.second];
    }
    int hub = 0;
    long long shared = -1;
    for (const std::pair<const int, long long>& entry : counts)
        if (entry.second > shared || (entry.second == shared && entry.first > hub)) {
            hub = entry.first;
            shared = entry.second;
        }
    if (shared < static_cast<long long>(bad.size()))
        return {};
    // The guilty partners first, largest measured gain first; then the rest the
    // same way, because the overshoot does not have to end with them.
    std::stable_sort(bad.begin(), bad.end(),
                     [](const std::pair<std::pair<int, int>, double>& x,
                        const std::pair<std::pair<int, int>, double>& y) {
                         return x.second > y.second;
                     });
    std::vector<int> queue;
    std::set<int> taken{hub};
    for (const std::pair<std::pair<int, int>, double>& one : bad)
        for (int i : {one.first.first, one.first.second})
            if (taken.insert(i).second)
                queue.push_back(i);
    std::vector<int> rest;
    for (std::size_t i = 0; i < plan.bundles.size(); ++i)
        if (taken.count(static_cast<int>(i)) == 0)
            rest.push_back(static_cast<int>(i));
    std::stable_sort(rest.begin(), rest.end(), [&](int x, int y) {
        const double left = best.count(x) != 0 ? best.at(x) : 0.;
        const double right = best.count(y) != 0 ? best.at(y) : 0.;
        return left != right ? left > right : x < y;
    });
    queue.insert(queue.end(), rest.begin(), rest.end());

    std::vector<int> whole;
    for (std::size_t i = 0; i < plan.bundles.size(); ++i)
        whole.push_back(static_cast<int>(i));
    std::vector<int> members{hub};
    std::pair<std::vector<int>, Value> found;
    for (int partner : queue) {
        std::vector<int> trial = members;
        trial.push_back(partner);
        std::sort(trial.begin(), trial.end());
        const std::pair<std::vector<Sites>, std::vector<Edge> > merged =
            merge_bundles(plan.bundles, plan.plain_edges(), trial);
        const Sites joined = union_of(plan.bundles, trial);
        const std::pair<std::vector<Sites>, std::vector<Edge> > alone(
            std::vector<Sites>{joined}, std::vector<Edge>{});
        std::set<Sites> wanted;
        for (const std::pair<const Sites, int>& term :
             terms(merged.first, merged.second, order_, config_limits()))
            wanted.insert(term.first);
        std::vector<Sites> patch_bundles;
        std::map<int, int> label;
        for (std::size_t i = 0; i < trial.size(); ++i) {
            patch_bundles.push_back(plan.bundles[static_cast<std::size_t>(trial[i])]);
            label[trial[i]] = static_cast<int>(i);
        }
        std::vector<Edge> patch_edges;
        for (const WeightedEdge& edge : plan.edges)
            if (label.count(edge.a) != 0 && label.count(edge.b) != 0)
                patch_edges.emplace_back(label[edge.a], label[edge.b]);
        for (const std::pair<const Sites, int>& term :
             terms(patch_bundles, patch_edges, order_, config_limits()))
            wanted.insert(term.first);
        wanted.insert(joined);
        if (!ensure(std::vector<Sites>(wanted.begin(), wanted.end()), plan.bundles))
            break;
        members = trial;
        // Physical is asked of the whole model, which is what the check is
        // about; the local errors describe the merge itself, as they do for an
        // ordinary candidate.
        Value proposal = Value::object();
        Value listed = Value::array();
        for (int i : trial)
            listed.push(Value::integer(static_cast<long long>(i)));
        proposal.set("pair", listed);
        proposal.set("hub", Value::boolean(true));
        proposal.set("verification_merge", Value::boolean(false));
        proposal.set("local_before",
                     Value::number(local(context, plan, trial, nullptr).scaled_error));
        proposal.set("local_after",
                     Value::number(local(context, plan, trial, &alone).scaled_error));
        const bool physical = is_physical(context, local_model(context, plan, whole, &merged));
        proposal.set("physical", Value::boolean(physical));
        found = {trial, proposal};
        if (physical)
            break;
    }
    return found;
}

int Driver::run() {
    if (limits_.regroup)
        return run_regroup();
    std::unique_ptr<Executed> step = execute("run", nullptr, nullptr, false);
    std::string directory = step->directory;
    Curves previous;
    bool has_previous = false;
    for (long long iteration = 0; iteration <= limits_.max_rounds; ++iteration) {
        const Curves curves = read_state_curves(directory, states_);
        Discrepancy change;
        if (has_previous)
            change = discrepancy(curves, previous, limits_.atol, limits_.rtol);
        if (!ensure(step->plan.bundles, step->plan.bundles))
            return finish("probe_budget_limited", directory, *step);
        MergeContext context = merge_;
        context.order = order_;
        context.raw = &raw_;
        context.gcce = gcce_;
        context.atol = limits_.atol;
        context.rtol = limits_.rtol;
        context.limits = config_limits();
        context.affordable = [this](const Plan& p, const std::vector<Sites>& c) {
            return affordable(p, c);
        };
        const std::pair<std::vector<std::vector<int>>, Coverage> found = patches(context, step->plan);
        std::vector<Sites> needed;
        for (const std::vector<int>& patch : found.first)
            for (const Sites& sites : patch_clusters(step->plan.bundles, patch))
                needed.push_back(sites);
        Value row = Value::object();
        row.set("iteration", Value::integer(static_cast<long long>(iteration)));
        row.set("bundles", Value::integer(static_cast<long long>(static_cast<long long>(step->plan.bundles.size()))));
        row.set("coverage", coverage_value(found.second));
        row.set("global_change", has_previous ? discrepancy_value(change) : Value());
        rounds_.push_back(row);
        if (!ensure(needed, step->plan.bundles))
            return finish("probe_budget_limited", directory, *step);

        std::vector<std::pair<std::vector<int>, Discrepancy>> assessed;
        for (const std::vector<int>& patch : found.first)
            assessed.emplace_back(patch, local(context, step->plan, patch, nullptr));
        // key=(-scaled_error, patch), and Python's sort is stable.
        std::stable_sort(assessed.begin(), assessed.end(),
            [](const std::pair<std::vector<int>, Discrepancy>& x,
               const std::pair<std::vector<int>, Discrepancy>& y) {
                if (x.second.scaled_error != y.second.scaled_error)
                    return x.second.scaled_error > y.second.scaled_error;
                return x.first < y.first;
            });
        Value listed = Value::array();
        for (const std::pair<std::vector<int>, Discrepancy>& item : assessed) {
            Value one = Value::object();
            one.set("bundles", ids_list(item.first));
            one.set("scaled_error", Value::number(item.second.scaled_error));
            one.set("max_abs_error", Value::number(item.second.max_abs_error));
            if (item.second.has_location) {
                Value where = Value::array();
                where.push(Value::integer(static_cast<long long>(item.second.state)));
                where.push(Value::integer(static_cast<long long>(item.second.step)));
                one.set("worst_state_step", where);
            } else {
                one.set("worst_state_step", Value());
            }
            listed.push(one);
        }
        row.set("patches", listed);
        // max((e['scaled_error'] for _, e in assessed), default=inf). Seeding with
        // the first entry would make the fold unreachable - the list is already
        // sorted by this very number - so it folds from below instead.
        double worst = HUGE_VAL;
        if (!assessed.empty()) {
            worst = -HUGE_VAL;
            for (const std::pair<std::vector<int>, Discrepancy>& item : assessed)
                if (item.second.scaled_error > worst)
                    worst = item.second.scaled_error;
        }
        row.set("max_local_scaled_error", Value::number(worst));
        bool physical = true;
        for (const Curve& curve : curves)
            for (const std::pair<double, double>& z : curve) {
                const std::complex<double> value(z.first, z.second);
                if (!adaptive::finite(value) ||
                    (!gcce_ && !(complex_abs(z.first, z.second) <= 1.+limits_.atol)))
                    physical = false;
            }
        row.set("finite_and_cce_bound", Value::boolean(physical));
        Curves exact;
        if (found.second.whole_bath_oracle) {
            std::vector<int> all;
            for (std::size_t i = 0; i < step->plan.site_dimensions.size(); ++i)
                all.push_back(static_cast<int>(i));
            std::vector<Sites> whole_bundle{union_of(step->plan.bundles, [&]{
                std::vector<int> ids;
                for (std::size_t i = 0; i < step->plan.bundles.size(); ++i)
                    ids.push_back(static_cast<int>(i));
                return ids;
            }())};
            // approximate([tuple(range(len(site_dimensions)))], [], 1, ...): one
            // cluster holding every site, with no edges.
            whole_bundle[0] = Sites(all.begin(), all.end());
            exact = approximate(whole_bundle, {}, 1, raw_, gcce_, config_limits());
            const Discrepancy whole = discrepancy(curves, exact, limits_.atol, limits_.rtol);
            row.set("whole_bath_error", discrepancy_value(whole));
            rounds_.back() = row;
            if (whole.scaled_error <= 1. && physical)
                return finish("validated_full_bath_model", directory, *step);
        }
        rounds_.back() = row;
        if (worst <= 1. && has_previous && change.scaled_error <= 1. && physical &&
            !found.second.coverage_limited)
            return finish("passed_convergence_checks", directory, *step);
        if (iteration == limits_.max_rounds)
            return finish("round_limit", directory, *step);

        std::map<std::pair<int, int>, Value> candidates;
        std::map<std::pair<int, int>, ProposalKey> keys;
        Feasibility cache;
        for (const std::pair<std::vector<int>, Discrepancy>& item : assessed) {
            const std::vector<int>& patch = item.first;
            if (found.second.whole_bath_oracle && patch.size() != step->plan.bundles.size())
                continue;
            std::vector<Sites> local_bundles;
            std::map<int, int> local_index;
            for (std::size_t i = 0; i < patch.size(); ++i) {
                local_bundles.push_back(step->plan.bundles[static_cast<std::size_t>(patch[i])]);
                local_index[patch[i]] = static_cast<int>(i);
            }
            std::vector<Edge> local_edges;
            for (const WeightedEdge& edge : step->plan.edges)
                if (local_index.count(edge.a) != 0 && local_index.count(edge.b) != 0)
                    local_edges.emplace_back(local_index[edge.a], local_index[edge.b]);
            for (std::size_t x = 0; x < patch.size(); ++x)
                for (std::size_t y = x+1; y < patch.size(); ++y) {
                    const std::pair<int, int> pair(patch[x], patch[y]);
                    if (!can_merge(&cache, context, step->plan, {pair.first, pair.second}))
                        continue;
                    const std::pair<std::vector<Sites>, std::vector<Edge>> changed =
                        merge_bundles(local_bundles, local_edges,
                                     {local_index[pair.first], local_index[pair.second]});
                    const double after = local(context, step->plan, patch, &changed).scaled_error;
                    const double before = item.second.scaled_error;
                    const double gain = improvement(before, after);
                    const bool verification = worst <= 1.;
                    if (!verification && !(gain > 1e-12))
                        continue;
                    // Written out rather than returned from a lambda: d*d*d is a
                    // Boost expression holding references to d, and a deduced
                    // return type carries that expression out of the scope where d
                    // lives. AddressSanitizer caught it as a stack-use-after-return
                    // (2026-09-17); the tests did not, because the freed bytes
                    // usually still held the right number.
                    const long long d = bundle_dimension(
                        union_of(step->plan.bundles, {pair.first, pair.second}),
                        step->plan.site_dimensions);
                    const long long cube = saturating_cube(d);
                    // (1 if verification else gain)/(d**3): one is int/int, which
                    // rounds once, the other float/int, which needs the cube as a
                    // double and raises if it will not fit - the way Python does.
                    const double score = verification ? true_div(static_cast<long long>(1), cube)
                                                      : gain/static_cast<double>(cube);
                    Value proposal = Value::object();
                    Value listed_pair = Value::array();
                    listed_pair.push(Value::integer(static_cast<long long>(pair.first)));
                    listed_pair.push(Value::integer(static_cast<long long>(pair.second)));
                    proposal.set("pair", listed_pair);
                    proposal.set("score", Value::number(score));
                    proposal.set("local_before", Value::number(before));
                    proposal.set("local_after", Value::number(after));
                    proposal.set("verification_merge", Value::boolean(verification));
                    if (limits_.prefer_physical)
                        proposal.set("physical", Value::boolean(is_physical(
                            context, local_model(context, step->plan, patch, &changed))));
                    const ProposalKey key = proposal_key(context, score, before, pair);
                    if (candidates.count(pair) == 0 || ::adaptive::before(key, keys.at(pair))) {
                        candidates[pair] = proposal;
                        keys[pair] = key;
                    }
                }
        }
        std::vector<std::pair<int, int>> order;
        for (const std::pair<const std::pair<int, int>, Value>& entry : candidates)
            order.push_back(entry.first);
        std::stable_sort(order.begin(), order.end(),
            [&](const std::pair<int, int>& x, const std::pair<int, int>& y) {
                return ::adaptive::before(keys.at(x), keys.at(y));
            });
        if (limits_.prefer_physical) {
            // Narrow, never refuse: an unphysical-only round still adopts, and a
            // candidate nobody measured counts as physical. Refusing outright
            // ended 140 of 384 baseline runs with no result
            // (docs/archive/adaptive-physical-gate-handoff-20260918.md, 3절).
            std::vector<std::pair<int, int> > good;
            for (const std::pair<int, int>& pair : order) {
                const Value& one = candidates.at(pair);
                if (!one.has("physical") || one.at("physical").as_bool())
                    good.push_back(pair);
            }
            if (!good.empty() && good.size() < order.size()) {
                Value told = Value::object();
                told.set("candidates", Value::integer(static_cast<long long>(order.size())));
                told.set("physical", Value::integer(static_cast<long long>(good.size())));
                row.set("physical_preference", told);
                order = good;
            }
        }
        // Never ranked, always first: the measured rule adopts the hub whenever
        // it exists, and competing on gain is the first of the three wrong
        // questions (252 offered, 8 adopted).
        std::vector<std::vector<int> > tried;
        std::pair<std::vector<int>, Value> hub;
        if (limits_.hub_candidates) {
            hub = hub_candidate(context, step->plan, candidates);
            if (!hub.first.empty()) {
                Value told = Value::object();
                Value listed = Value::array();
                for (int i : hub.first)
                    listed.push(Value::integer(static_cast<long long>(i)));
                told.set("bundles", listed);
                told.set("sites", sites_value(union_of(step->plan.bundles, hub.first)));
                told.set("grown", Value::integer(static_cast<long long>(hub.first.size())-1));
                long long unphysical = 0;
                for (const std::pair<const std::pair<int, int>, Value>& one : candidates)
                    if (one.second.has("physical") && !one.second.at("physical").as_bool())
                        ++unphysical;
                told.set("from_unphysical", Value::integer(unphysical));
                told.set("physical", hub.second.at("physical"));
                told.set("adopted", Value::boolean(false));
                row.set("hub", told);
                tried.push_back(hub.first);
            }
        }
        for (const std::pair<int, int>& pair : order)
            tried.push_back({pair.first, pair.second});
        bool accepted = false;
        std::vector<Sites> chosen;
        Value accepted_proposal;
        for (std::size_t i = 0; i < tried.size() &&
                                static_cast<long long>(i) < limits_.max_candidates; ++i) {
            const std::vector<int>& members = tried[i];
            const bool is_hub = !hub.first.empty() && i == 0;
            const std::pair<int, int> pair(members.front(), members.back());
            if (!is_hub && !can_merge(&cache, context, step->plan, members))
                continue;                               // feasible() checks again
            const std::pair<std::vector<Sites>, std::vector<Edge>> after =
                merge_bundles(step->plan.bundles, step->plan.plain_edges(), members);
            std::unique_ptr<Executed> trial = execute("plan", &after.first, nullptr, true);
            if (trial == nullptr)
                continue;
            if (found.second.whole_bath_oracle) {
                std::vector<Edge> edges;
                for (const WeightedEdge& edge : trial->plan.edges)
                    edges.emplace_back(edge.a, edge.b);
                const Curves predicted = approximate(trial->plan.bundles, edges, order_, raw_,
                                                     gcce_, config_limits());
                const double error =
                    discrepancy(predicted, exact, limits_.atol, limits_.rtol).scaled_error;
                if (!(error < row.at("whole_bath_error").at("scaled_error").as_double()))
                    continue;
                Value with = is_hub ? hub.second : candidates.at(pair);
                with.set("actual_graph_scaled_error", Value::number(error));
                (is_hub ? hub.second : candidates[pair]) = with;
            }
            accepted = true;
            chosen = after.first;
            accepted_proposal = is_hub ? hub.second : candidates.at(pair);
            if (is_hub) {
                Value told = row.at("hub");
                told.set("adopted", Value::boolean(true));
                row.set("hub", told);
            }
            break;
        }
        if (!accepted)
            return finish("no_feasible_improving_merge", directory, *step);
        Value with_sites = accepted_proposal;
        // "pair" holds two bundles for an ordinary merge and more for a hub.
        std::vector<int> joined;
        for (std::size_t i = 0; i < accepted_proposal.at("pair").size(); ++i)
            joined.push_back(static_cast<int>(accepted_proposal.at("pair").at(i).as_int()));
        Value sites = Value::array();
        for (int site : union_of(step->plan.bundles, joined))
            sites.push(Value::integer(static_cast<long long>(site)));
        with_sites.set("sites", sites);
        row.set("accepted", with_sites);
        rounds_.back() = row;
        persist();
        // print('REFINE', iteration, sites, 'local error', before, '->', after):
        // sites is a tuple, so it shows as (1, 2) and not as JSON.
        std::printf("REFINE %lld %s local error %s -> %s\n", iteration,
                    sites_repr(union_of(step->plan.bundles, joined)).c_str(),
                    float_text(accepted_proposal.at("local_before").as_double()).c_str(),
                    float_text(accepted_proposal.at("local_after").as_double()).c_str());
        std::fflush(stdout);
        previous = curves;
        has_previous = true;
        step = execute("run", &chosen, nullptr, false);
        directory = step->directory;
    }
    // Not reachable: the loop runs while iteration <= max_rounds and returns
    // round_limit on the last one. Here so the function has an answer on every
    // path, and so a future change to the loop bound fails loudly rather than
    // falling out of the end.
    return finish("round_limit", directory, *step);
}
}   // namespace adaptive
