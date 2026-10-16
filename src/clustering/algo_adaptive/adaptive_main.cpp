/**
 * @file adaptive_main.cpp
 * @brief Adaptive refinement: the entry point, which reads the options and runs the merge policy.
 */
// bin/ccex-adaptive: the refinement driver, without Python (plan stage 1, Task 7).
//
// Structure C of the design: this program makes the decisions and launches the
// CCEX binary a step at a time. It links no MPI and no MKL - the binary it starts
// does that.
#include "base_run.h"
#include "policy_merge.h"
#include "base_json.h"
#include "options.h"
#include "options_table.h"
#include "base_numeric.h"
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {
// argparse's p.error(): the message, and exit 2 - which is not the 1 a failed run
// gives, so the two are kept apart.
struct UsageError {
    std::string message;
};

// Everything Python's Driver.__init__ checks about one input, and the two things
// it derives from it. Run once at start-up.
struct Prepared {
    adaptive::Value cfg;
    std::string method;
    std::vector<int> states;
};

// The input may carry a cluster.adaptive object (a top-level "adaptive" until A65) - that is how bin/main.out is
// asked to refine without a command line. It is this run's request, not part of the calculation: copied into a step's
// input.json it would ask the step to start a refinement of its own. CCEX reads keys without regard to case
// ("Adaptive" until A47), so any spelling is left out; so is a step's "adaptive": {"work_file"}. Its probe_memory_gb is
// kept: it is the steps' setting (the Config reader's cluster.adaptive.probe_memory_gb).
std::string lowered(std::string key) {
    for (char& c : key)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return key;
}
adaptive::Value without_adaptive_block(const adaptive::Value& given) {
    adaptive::Value out = adaptive::Value::object();
    for (const std::string& key : given.keys()) {
        if (lowered(key) != "cluster" || given.at(key).kind() != adaptive::Value::Kind::Object) {
            out.set(key, given.at(key));
            continue;
        }
        const adaptive::Value& block = given.at(key);
        adaptive::Value kept = adaptive::Value::object();
        for (const std::string& inner : block.keys()) {
            if (lowered(inner) != "adaptive") {
                kept.set(inner, block.at(inner));
                continue;
            }
            if (block.at(inner).kind() != adaptive::Value::Kind::Object)
                continue;
            adaptive::Value probe = adaptive::Value::object();
            for (const std::string& option : block.at(inner).keys())
                if (lowered(option) == "probe_memory_gb")
                    probe.set(option, block.at(inner).at(option));
            if (!probe.keys().empty())
                kept.set(inner, probe);
        }
        out.set(key, kept);
    }
    return out;
}

Prepared prepare_config(const adaptive::Value& given) {
    Prepared out;
    out.cfg = given;
    if (!out.cfg.has("rdip") && !adaptive::find_path(out.cfg, {"cluster", "bundle", "connections"}))
        throw adaptive::ValueError(
            "For rdip-omitted checked connections, run CCEX directly. This separate "
            "refinement driver requires rdip or explicit legacy connection settings.");
    const adaptive::Value* mode = adaptive::find_path(out.cfg, {"cluster", "mode"});   // cluster_mode until A65
    if (!mode || mode->kind() != adaptive::Value::Kind::String || mode->as_string() != "bundle" ||
        (out.cfg.has("quantity") && out.cfg.at("quantity").as_string() != "coherence"))
        throw adaptive::ValueError("input must use cluster.mode=bundle and quantity=coherence");
    // method cCCE or gCCE (A47); inside, the policy names them cce and gcce, as CCEX does
    out.method = out.cfg.has("method") ? out.cfg.at("method").as_string() : "ccce";
    for (char& c : out.method)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (out.method == "ccce")
        out.method = "cce";
    else if (out.method != "gcce")
        throw adaptive::ValueError("refinement supports method cCCE/gCCE");
    // seed and s_range are the configuration block's (A64, A94)
    const adaptive::Value* seed = adaptive::find_path(out.cfg, {"configuration", "seed"});
    if (!seed || seed->kind() != adaptive::Value::Kind::Int)
        throw adaptive::ValueError(
            "an explicit integer configuration.seed is required to preserve bath states between runs");
    const long long steps =
        out.cfg.has("nstep") ? static_cast<long long>(out.cfg.at("nstep").as_int()) : 0;
    const double delta = out.cfg.has("deltat") ? out.cfg.at("deltat").as_double() : 0.;
    if (steps < 2 || !(delta > 0.))
        throw adaptive::ValueError("refinement needs at least two time points and positive deltat");
    // the step inputs save the coherence; a savemode block keeps its nodiv (until 2026-10-13 the top-level save_nodiv,
    // which was copied through)
    adaptive::Value savemode = adaptive::Value::text("coherence");
    if (const adaptive::Value* given = adaptive::find_path(out.cfg, {"output", "savemode"}))   // savemode until A65
        if (given->kind() == adaptive::Value::Kind::Object && given->has("nodiv")) {
            savemode = adaptive::Value::object();
            savemode.set("mode", adaptive::Value::text("coherence"));
            savemode.set("nodiv", given->at("nodiv"));
        }
    // A73: the output block and its savemode as the input spells them (key_of), not a second block CCEX would not read
    const std::string output_key = adaptive::key_of(out.cfg, "output");
    adaptive::Value output = out.cfg.has(output_key) && out.cfg.at(output_key).kind() == adaptive::Value::Kind::Object
                                 ? out.cfg.at(output_key) : adaptive::Value::object();
    output.set(adaptive::key_of(output, "savemode"), savemode);
    out.cfg.set(output_key, output);
    out.states.push_back(0);
    // configuration.s_range [first, last] (A94; nstate N before, the states 1..N): the sampled states the run calculates
    const adaptive::Value* range = adaptive::find_path(out.cfg, {"configuration", "s_range"});
    if (range && range->kind() == adaptive::Value::Kind::Array && range->size() == 2) {
        out.states.clear();
        for (long long s = range->at(std::size_t(0)).as_int(); s <= static_cast<long long>(range->at(std::size_t(1)).as_int()); ++s)
            out.states.push_back(static_cast<int>(s));
    }
    return out;
}

}   // namespace

// The tool's body, so that bin/main.out can run the same policy in its own
// process rather than launching this one per step (plan stage 8). bin/ccex-adaptive
// is main() below and passes its own command line straight through.
int adaptive_main(const std::vector<std::string>& given) {
    try {
        std::vector<std::string> arguments = given;
        std::vector<char> here(4096, '\0');
        if (getcwd(here.data(), here.size()) == nullptr)
            throw adaptive::os_error(errno, ".");
        // A removed policy is told why before the parser calls it an invalid
        // choice. The last --policy is the one the parser would keep.
        std::string asked;
        for (std::size_t i = 0; i < arguments.size(); ++i) {
            if (arguments[i] == "--policy" && i+1 < arguments.size())
                asked = arguments[i+1];
            else if (arguments[i].rfind("--policy=", 0) == 0)
                asked = arguments[i].substr(9);
        }
        if (asked == "error-gated" || asked == "budgeted" || asked == "auto-partition" ||
            asked == "cost-preserving" || asked == "auto-order")
            throw adaptive::RuntimeError(
                "the only policy is merge (error-gated, budgeted, auto-partition and "
                "cost-preserving were removed on 2026-10-10, auto-order on 2026-09-25: "
                "use merge, with prefer_physical and hub_candidates if wanted)");
        adaptive::Value args = adaptive::parse_options(arguments, here.data());
        // Driver.__init__ fills this default into args itself, before the report
        // copies them - so the settings the run records carry the number, not None.
        if (args.at("max_candidates").kind() == adaptive::Value::Kind::Null)
            args.set("max_candidates", adaptive::Value::integer(static_cast<long long>(16)));

        const adaptive::Value given = without_adaptive_block(adaptive::read_config(adaptive::without_comments(adaptive::read_file(args.at("input").as_string()))));
        // cfg.pop('adaptive_work_file', None): Value has no remove, so the one key
        // is left out while the rest keeps its order.
        // The checks main() makes before anything is built.
        // The hub is built out of the per-candidate physical verdicts, which
        // nothing measures unless --prefer-physical asked for them.
        if (args.at("hub_candidates").as_bool() && !args.at("prefer_physical").as_bool())
            throw UsageError{"--hub-candidates requires --prefer-physical: "
                             "the hub is built from its verdicts"};
        // Method D's three options. Off, they leave the settings the report copies, so
        // every run without them records exactly what it recorded before they existed.
        const bool regroup = args.at("regroup").as_bool();
        const long long regroup_regions = args.at("regroup_max_regions").as_int();
        const long long regroup_candidates = args.at("regroup_max_candidates").as_int();
        // Stages 3 and 3b (adaptive-regroup-stage3.md, -stage3b.md): the factors work mode
        // for the culprit, the except work mode for condition (c). At their default,
        // "probes", they leave the settings too, so the regroup runs recorded before them
        // record the same.
        const bool regroup_factors = args.at("regroup_factors").as_string() == "factors";
        const bool regroup_except = args.at("regroup_whole").as_string() == "except";
        {
            adaptive::Value kept = adaptive::Value::object();
            for (const std::string& key : args.keys())
                if (!(key == "regroup_factors" && !regroup_factors) && !(key == "regroup_whole" && !regroup_except))
                    kept.set(key, args.at(key));
            args = kept;
        }
        if (!regroup) {
            if (regroup_regions != 64 || regroup_candidates != 256)
                throw UsageError{"--regroup-max-regions and --regroup-max-candidates need --regroup"};
            if (regroup_factors || regroup_except)
                throw UsageError{"--regroup-factors and --regroup-whole need --regroup"};
            adaptive::Value kept = adaptive::Value::object();
            for (const std::string& key : args.keys())
                if (key != "regroup" && key != "regroup_max_regions" && key != "regroup_max_candidates")
                    kept.set(key, args.at(key));
            args = kept;
        } else {
            // D does not narrow or grow merge's candidates: it replaces them
            // (adaptive-regroup-plan.md, 결정한 것).
            if (args.at("prefer_physical").as_bool() || args.at("hub_candidates").as_bool())
                throw UsageError{"--regroup replaces merge's selection and cannot be combined with "
                                 "--prefer-physical or --hub-candidates"};
            if (regroup_regions < 1 || regroup_candidates < 1)
                throw UsageError{"regroup-max-regions and regroup-max-candidates must be positive"};
            if (regroup_except && !adaptive::kLongDoubleX87)
                throw UsageError{"--regroup-whole except needs long double to be the x87 80-bit format, which this "
                                 "build does not have (it reassembles CCEX's long double product, in double the "
                                 "rounding reached 1e-10): use --regroup-whole probes"};
        }
        if (args.at("probe_batch_clusters").as_int() < 0)
            throw UsageError{"probe-batch-clusters must be nonnegative"};
        const double atol = adaptive::number_of(args.at("atol")), rtol = adaptive::number_of(args.at("rtol"));
        if (!std::isfinite(atol) || atol < 0. || !std::isfinite(rtol) || rtol < 0. ||
            atol+rtol == 0.)
            throw UsageError{"tolerances must be finite/nonnegative and at least one must be positive"};
        const double timeout_value = adaptive::number_of(args.at("timeout"));
        if (args.at("max_rounds").as_int() < 0 ||
            args.at("ranks").as_int() < 1 || args.at("max_patches").as_int() < 1 ||
            args.at("max_probes").as_int() < 1 || args.at("max_candidates").as_int() < 1 ||
            args.at("neighbors").as_int() < 1 || args.at("cache_mb").as_int() < 1 ||
            !std::isfinite(timeout_value) || timeout_value <= 0.)
            throw UsageError{"budgets/ranks/timeout must be positive, max-rounds may be zero"};

        const Prepared prepared = prepare_config(given);   // the adaptive block, work_file with it, is left out
        adaptive::Value cfg = prepared.cfg;
        const std::string method = prepared.method;
        const std::vector<int> states = prepared.states;
        if (!cfg.has("order") || cfg.at("order").kind() != adaptive::Value::Kind::Int ||
            cfg.at("order").as_int() < 1)
            throw adaptive::ValueError("adaptive requires positive integer bundle order");
        const int order = static_cast<int>(cfg.at("order").as_int());
        if (regroup) {
            // The first version's scope (adaptive-regroup-plan.md, 결정한 것; review 1, 4):
            // cCCE coherence, bundle order 2, connections by distance. Refused before
            // anything runs rather than run under rules D was not written for.
            if (method != "cce")
                throw adaptive::ValueError("--regroup supports method cce only (gCCE's physical bound is not decided)");
            if (order != 2)
                throw adaptive::ValueError("--regroup needs bundle order 2: the culprit is a pair of bundles");
            const adaptive::Value* given_connections = adaptive::find_path(cfg, {"cluster", "bundle", "connections"});
            std::string connections = given_connections ? given_connections->as_string() : "distance";
            for (char& c : connections)                     // CCEX reads it with strcasecmp
                c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (connections != "distance")
                throw adaptive::ValueError("--regroup needs cluster.bundle.connections distance: its connections are worked "
                                           "out again from the site distances, which auto and checked do not use");
            if (!cfg.has("rdip"))
                throw adaptive::ValueError("--regroup needs rdip: its connections are the sites within it");
            // What CCEX refuses in the factors and except work modes, refused here first: in-process
            // a refused step is an MPI_Abort of the whole job (review 14; adaptive-regroup-stage3.md 4절).
            // Method and order are checked above; the power-0 terms these two options add would mix
            // with the culprit factors and with P_K.
            if (regroup_factors || regroup_except) {
                // the output.cluster_contributions block is on unless it says "enabled": false (A47; in output since A65)
                const adaptive::Value* block = adaptive::find_path(cfg, {"output", "cluster_contributions"});
                const adaptive::Value* enabled = block ? adaptive::find_path(*block, {"enabled"}) : nullptr;
                const bool contributions = block && block->kind() == adaptive::Value::Kind::Object &&
                                           (!enabled || (enabled->kind() == adaptive::Value::Kind::Bool && enabled->as_bool()));
                const bool ladder = cfg.has("rdip_ladder") && cfg.at("rdip_ladder").kind() == adaptive::Value::Kind::Array &&
                                    cfg.at("rdip_ladder").size() > 0;
                if (contributions || ladder)
                    throw adaptive::ValueError("--regroup-factors factors and --regroup-whole except cannot be combined "
                                               "with output.cluster_contributions or rdip_ladder");
            }
        }

        adaptive::DriverSettings settings;
        settings.ranks = static_cast<long long>(args.at("ranks").as_int());
        settings.max_probes = static_cast<long long>(args.at("max_probes").as_int());
        settings.cache_mb = static_cast<long long>(args.at("cache_mb").as_int());
        settings.probe_batch_clusters =
            static_cast<long long>(args.at("probe_batch_clusters").as_int());
        settings.max_total_work = adaptive::number_of(args.at("max_total_work"));
        settings.cwd = args.at("cwd").as_string();

        adaptive::SystemClock clock;
        adaptive::Driver driver(cfg, states, settings, adaptive::Driver::Execute());
        adaptive::NativeSettings native;
        native.mpiexec = args.at("mpiexec").as_string();
        native.binary = adaptive::resolved(settings.cwd, args.at("binary").as_string());
        native.timeout = args.at("timeout").kind() == adaptive::Value::Kind::Int
                             ? static_cast<double>(args.at("timeout").as_int())
                             : args.at("timeout").as_double();
        native.timeout_is_integer = args.at("timeout").kind() == adaptive::Value::Kind::Int;
        // D measures closeness and connections from the site positions.
        native.physics_metadata = regroup;
        const std::string output = adaptive::resolved(settings.cwd, args.at("output").as_string());
        driver.prepare(output, native, &clock);
        driver.begin(args, adaptive::binary_sha256(native.binary));
        // Python's Driver.__init__ opens the cache right after the report skeleton.
        if (args.at("raw_cache").kind() != adaptive::Value::Kind::Null) {
            const adaptive::Value context = adaptive::identity(
                cfg, settings.cwd, states, adaptive::binary_sha256(native.binary), settings.ranks);
            driver.use_cache(std::unique_ptr<adaptive::ProbeCache>(
                                 new adaptive::ProbeCache(args.at("raw_cache").as_string(), context)),
                             context);
        }

        adaptive::MergeContext merge;
        merge.neighbors = static_cast<long long>(args.at("neighbors").as_int());
        merge.max_patches = static_cast<long long>(args.at("max_patches").as_int());
        const adaptive::Value* max_dim = adaptive::find_path(cfg, {"cluster", "bundle", "max_dim"});
        if (!max_dim)
            throw adaptive::ValueError("cluster.mode bundle requires cluster.bundle.max_dim");
        merge.bundle_max_dim = max_dim->as_int();
        merge.merge_priority = args.at("merge_priority").as_string();
        adaptive::Driver::Limits limits;
        limits.max_rounds = static_cast<long long>(args.at("max_rounds").as_int());
        limits.max_candidates = static_cast<long long>(args.at("max_candidates").as_int());
        limits.atol = args.at("atol").as_double();
        limits.rtol = args.at("rtol").as_double();
        limits.prefer_physical = args.at("prefer_physical").as_bool();
        limits.hub_candidates = args.at("hub_candidates").as_bool();
        limits.regroup = regroup;
        limits.regroup_max_regions = regroup_regions;
        limits.regroup_max_candidates = regroup_candidates;
        limits.regroup_factors = regroup_factors;
        limits.regroup_except = regroup_except;
        driver.set_policy(order, method == "gcce", merge, limits);

        int status = 0;
        try {
            status = driver.run();
        } catch (const adaptive::Error& error) {
            // main() keeps the summary: a failed run still says what it did and
            // why, and leaves with 1.
            driver.report().set("status", adaptive::Value::text("failed"));
            driver.report().set("error", adaptive::Value::text(error.what()));
            driver.persist();
            std::cerr << "Refinement failed: " << error.what() << "\n";
            status = 1;
        }
        driver.close_cache();
        return status;
    } catch (const UsageError& error) {
        std::cerr << "ccex-adaptive: error: " << error.message << "\n";
        return 2;
    } catch (const adaptive::Error& error) {
        std::cerr << error.type() << ": " << error.what() << "\n";
        return 1;
    }
}

#ifndef CCEX_ADAPTIVE_NO_MAIN
// bin/ccex-adaptive. bin/main.out compiles this file with CCEX_ADAPTIVE_NO_MAIN and
// keeps only adaptive_main, because it has a main() of its own (plan stage 8).
int main(int argc, char** argv) {
    return adaptive_main(std::vector<std::string>(argv+1, argv+argc));
}
#endif
