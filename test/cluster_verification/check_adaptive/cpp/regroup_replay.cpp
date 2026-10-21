// The stored-probe replay of method D ((internal notes), stages
// 0 and 4): the C++ decision code of --regroup driven by the stored P1 probes the
// prototype ((internal notes)) was run on.
//
//     regroup_replay [--factors] [--except] CASES.jsonl > TRACES.jsonl
//
// --factors and --except take the paths of --regroup-factors factors and --regroup-whole
// except (adaptive-regroup-stage3.md, -stage3b.md): the culprit values read at the worst
// point of the stored curves, and condition (c) assembled from a kept product P_K made
// here from the stored curves of the unchanged terms (CCEX's except step is not run: the
// replay checks the table difference, the ledger and the decisions on the same probes).
//
// Each line of CASES.jsonl is a line of `rebundle2.py trace` (its own trace is not read):
// the case folder (a probe response of every subset of the twelve sites, written by
// CCEX on the server), the states it holds and which of them choose, the N sites, K and
// the starting partition. The line written back has the same fields and the C++ trace.
//
// What is given, and how (review 10): the probes are read as they were stored, never
// recalculated and never put in a cache of the current binary. L0 = 1 (cCCE), as the
// stored responses hold it. The run curve of a partition is its expansion assembled
// from those probes - the same assembly CCEX's run makes, with nothing else to compare
// it with here. Every pair of bundles is connected (the prototype's graph); the site
// IDs are the stored ones, used as they are; a site is two-dimensional (P1 with the
// 14N scalar, the stored plan's site_dimensions), so K = 2, 3 is bundle.max_dim 4, 8.
// The budget is unlimited, which is what the prototype had.
#include <cmath>
#include <dirent.h>
#include <fnmatch.h>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "base_bundles.h"
#include "base_json.h"
#include "base_numeric.h"
#include "base_run.h"
#include "policy_regroup.h"

namespace {
using namespace adaptive;
using namespace adaptive::regroup;

// One case folder: every stored probe over every state, in the order of the states.
struct Stored {
    std::string folder;
    std::vector<int> states;
    std::map<Sites, Curves> curves;
    std::map<int, std::array<double, 3> > positions;
};

std::string parent(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? "." : path.substr(0, slash);
}

Stored load(const std::string& folder) {
    Stored out;
    out.folder = folder;
    const Value request = read_json(read_file(folder+"/request.json"));
    const Value& listed = request.at("probes");
    std::vector<Sites> ids;
    for (std::size_t i = 0; i < listed.size(); ++i)
        ids.push_back(sites_of(listed.at(i)));
    const long long steps = read_json(read_file(folder+"/input.json")).at("nstep").as_int();
    std::vector<std::pair<int, std::string> > files;
    if (DIR* listing = opendir(folder.c_str())) {
        while (const dirent* entry = readdir(listing))
            if (fnmatch("response.rank-*.state-*.tsv", entry->d_name, 0) == 0) {
                const std::string name = entry->d_name;
                const std::size_t at = name.find(".state-")+7;
                files.emplace_back(std::stoi(name.substr(at, name.size()-4-at)), folder+"/"+name);
            }
        closedir(listing);
    }
    std::set<int> states;
    for (const std::pair<int, std::string>& file : files)
        states.insert(file.first);
    out.states.assign(states.begin(), states.end());
    std::map<int, std::size_t> slot;
    for (std::size_t i = 0; i < out.states.size(); ++i)
        slot[out.states[i]] = i;
    std::map<Sites, std::vector<std::vector<bool> > > seen;
    for (const Sites& sites : ids) {
        out.curves[sites] = Curves(out.states.size(), Curve(static_cast<std::size_t>(steps)));
        seen[sites] = std::vector<std::vector<bool> >(out.states.size(),
                                                      std::vector<bool>(static_cast<std::size_t>(steps), false));
    }
    for (const std::pair<int, std::string>& file : files) {
        std::istringstream rows(read_file(file.second));
        std::string line;
        std::getline(rows, line);
        while (std::getline(rows, line)) {
            if (line.empty())
                continue;
            std::istringstream fields(line);
            std::string probe, step, re, im;
            std::getline(fields, probe, '\t');
            std::getline(fields, step, '\t');
            std::getline(fields, re, '\t');
            std::getline(fields, im, '\t');
            const long long id = parse_int(probe), t = parse_int(step);
            if (id < 0)
                continue;
            const Sites& sites = ids.at(static_cast<std::size_t>(id));
            const std::size_t s = slot.at(file.first);
            if (seen[sites][s].at(static_cast<std::size_t>(t)))
                throw RuntimeError(folder+": a probe point written twice");
            seen[sites][s][static_cast<std::size_t>(t)] = true;
            out.curves[sites][s][static_cast<std::size_t>(t)] = {parse_float(re), parse_float(im)};
        }
    }
    for (const auto& entry : seen)
        for (const std::vector<bool>& state : entry.second)
            for (bool filled : state)
                if (!filled)
                    throw RuntimeError(folder+": an incomplete probe");
    const Value physics = read_json(read_file(parent(folder)+"/plan/response.physics.json"));
    for (std::size_t i = 0; i < physics.at("sites").size(); ++i) {
        const Value& site = physics.at("sites").at(i);
        std::array<double, 3> p{};
        for (std::size_t k = 0; k < 3; ++k)
            p[k] = number_of(site.at("xyz").at(k));
        out.positions[static_cast<int>(site.at("id").as_int())] = p;
    }
    return out;
}

// The stored probes as a Source: `visible` is what the search has asked for so far,
// so the ledger counts what a run would have had to calculate.
struct Replay : Source {
    std::map<Sites, Curves> store;
    std::map<Sites, Curves> visible;
    Plan plan;
    long long states = 0, steps = 0;
    const std::map<Sites, Curves>& raw() const override { return visible; }
    long long work(const std::vector<Sites>& clusters) const override {
        return cluster_work(plan, clusters, states, steps, 0);
    }
    long long run_work(const std::map<Sites, int>& terms) const override {
        std::vector<Sites> clusters;
        for (const auto& term : terms)
            clusters.push_back(term.first);
        return work(clusters);
    }
    long long dimension(const Sites& sites) const override {
        return bundle_dimension(sites, plan.site_dimensions);
    }
    bool fetch(const std::vector<Sites>& missing) override {
        for (const Sites& sites : missing) {
            const auto found = store.find(sites);
            if (found == store.end())
                throw RuntimeError("probe "+sites_repr(sites)+" is not in the stored data");
            visible[sites] = found->second;
        }
        return true;
    }
    // --factors: the point the factors step would give, out of the stored curves.
    long long point_work(const std::vector<Sites>& clusters) const override {
        return cluster_work(plan, clusters, 1, steps, 0);
    }
    bool point_values(const Partition& units, const std::vector<Edge>& edges, const Point& at,
                      PointValues* out) override {
        *out = point_values_of(units, edges, store, at);
        return true;
    }
    // --except: P_K and the whole curve from the stored curves (store holds L0 too).
    const Adjacency* links = nullptr;
    ExpansionLimits limits;
    bool kept(const Partition& units, const std::vector<Sites>& except, Kept* out) override {
        const std::map<Sites, int> now = expansion_terms(units, *links, limits);
        std::map<Sites, int> rest = now;
        for (const Sites& sites : except)
            rest.erase(sites);
        out->product = scaled_product(rest, store);
        out->whole = approximate_weights(now, store, false);
        out->added.clear();
        for (const Sites& sites : except)
            if (visible.count(sites) == 0) {
                visible[sites] = store.at(sites);
                out->added.push_back(sites);
            }
        return true;
    }
    void forget(const std::vector<Sites>& clusters) override {
        for (const Sites& sites : clusters)
            visible.erase(sites);
    }
};
}   // namespace

int main(int argc, char** argv) {
    bool factors = false, except = false;
    std::vector<std::string> given;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--factors")
            factors = true;
        else if (arg == "--except")
            except = true;
        else
            given.push_back(arg);
    }
    const std::string cases = given.size() == 1 ? given[0] : std::string();
    if (cases.empty()) {
        std::cerr << "usage: regroup_replay [--factors] [--except] CASES.jsonl > TRACES.jsonl\n";
        return 2;
    }
    try {
        std::istringstream lines(read_file(cases));
        std::string line;
        Stored stored;
        while (std::getline(lines, line)) {
            if (line.empty())
                continue;
            const Value given = read_json(line);
            const std::string folder = given.at("case").as_string();
            if (stored.folder != folder)
                stored = load(folder);
            std::vector<int> states;
            for (std::size_t i = 0; i < given.at("states").size(); ++i)
                states.push_back(static_cast<int>(given.at("states").at(i).as_int()));
            if (states != stored.states)
                throw RuntimeError(folder+": the states differ from the trace's");
            std::vector<std::size_t> select;
            for (std::size_t i = 0; i < given.at("select").size(); ++i)
                select.push_back(static_cast<std::size_t>(given.at("select").at(i).as_int()));
            const Sites sites = sites_of(given.at("sites"));
            const int largest = sites.back();
            const std::set<int> inside(sites.begin(), sites.end());

            Replay source;
            source.plan.site_dimensions.assign(static_cast<std::size_t>(largest)+1, 2);
            source.plan.qubit_dimension = 1;
            source.states = static_cast<long long>(select.size());
            source.steps = static_cast<long long>(stored.curves.begin()->second.at(0).size());
            for (const auto& entry : stored.curves) {
                bool contained = true;
                for (int s : entry.first)
                    contained = contained && inside.count(s) != 0;
                if (!contained)
                    continue;
                Curves chosen;
                for (std::size_t s : select)
                    chosen.push_back(entry.second.at(s));
                source.store[entry.first] = chosen;
            }
            source.visible[Sites()] = Curves(select.size(), Curve(static_cast<std::size_t>(source.steps),
                                                                  std::pair<double, double>(1., 0.)));
            source.store[Sites()] = source.visible[Sites()];
            Geometry geometry(static_cast<std::size_t>(largest)+1);
            Adjacency links(static_cast<std::size_t>(largest)+1);
            for (int s : sites) {
                geometry[static_cast<std::size_t>(s)].xyz = stored.positions.at(s);
                geometry[static_cast<std::size_t>(s)].components.push_back(stored.positions.at(s));
                for (int t : sites)
                    if (t != s)
                        links[static_cast<std::size_t>(s)].push_back(t);
            }
            Settings settings;
            settings.cap = 1LL << given.at("K").as_int();
            settings.dimensions = source.plan.site_dimensions;
            settings.max_regions = std::numeric_limits<long long>::max();
            settings.max_candidates = std::numeric_limits<long long>::max();
            settings.trace_curves = true;
            settings.culprit_factors = factors;
            settings.whole_except = except;
            source.links = &links;
            source.limits = settings.limits;
            Ledger ledger;
            ledger.total = HUGE_VAL;
            Partition start;
            for (std::size_t i = 0; i < given.at("start").size(); ++i)
                start.push_back(sites_of(given.at("start").at(i)));
            start = canonical(start);
            std::vector<Sites> first;
            for (const auto& term : expansion_terms(start, links, settings.limits))
                if (source.visible.count(term.first) == 0)
                    first.push_back(term.first);
            source.fetch(first);
            Search search(settings, &geometry, &links, &source, &ledger);
            const Constants constants;
            Value rows = Value::array();
            const Result result = drive(
                start, expansion(start, links, source.visible, settings.limits), &search, constants,
                [&](const Partition& p) {
                    // the run of a partition: its terms out of the stored curves (with --except they
                    // were not all asked for, or were dropped after their region)
                    std::vector<Sites> missing;
                    for (const auto& term : expansion_terms(p, links, settings.limits))
                        if (source.visible.count(term.first) == 0)
                            missing.push_back(term.first);
                    source.fetch(missing);
                    return expansion(p, links, source.visible, settings.limits);
                },
                [](bool) {}, &rows);
            Value last = Value::object();
            last.set("stop", Value::text(result.stop));
            last.set("final", partition_value(result.final));
            last.set("changes", Value::integer(result.changes));
            last.set("examined", Value::integer(search.examined()));
            last.set("rollbacks", Value::integer(result.rollbacks));
            last.set("ledger", ledger.value());
            rows.push(last);
            Value out = Value::object();
            for (const char* key : {"N", "K", "field", "CR", "cfg", "case", "states", "select", "sites", "start"})
                out.set(key, given.at(key));
            out.set("trace", rows);
            std::cout << json_dumps(out, true, false, true) << "\n";
        }
    } catch (const Error& error) {
        std::cerr << error.type() << ": " << error.what() << "\n";
        return 1;
    }
    return 0;
}
