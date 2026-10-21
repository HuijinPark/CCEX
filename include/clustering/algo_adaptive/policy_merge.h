/**
 * @file policy_merge.h
 * @brief Adaptive refinement: the merge policy, which joins the bundle pair that lowers a local error the most.
 */
// Policy "merge" (the only one since 2026-10-10): probe pairs of bundles and take, one
// step at a time, the merge that lowers a local error the most.
#pragma once
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "base_bundles.h"
#include "base_run.h"

// ----------------------------------------------------------------------------
// The merge policy's judgement: which patches to look at, which pairs may join,
// and in what order the proposals are tried (plan stage 1, Task 4).
// ----------------------------------------------------------------------------
namespace adaptive {
// What the driver knows that these decisions need. Passed in rather than reached
// for, so the judgement can be tested without a binary or a bath.
// Can CCEX be asked to build these clusters with `work_left` to spend? The cost
// is D^3 * states * steps, and a cluster already held is not asked for - against
// a running total, charging it again shut a search out of its own whole-bath
// reference (2026-09-19).
bool affordable_clusters(const Plan& plan, const std::vector<Sites>& clusters,
                         const std::map<Sites, Curves>& held, long long states,
                         long long steps, long long work_left);

struct MergeContext {
    int order = 2;
    long long neighbors = 8;
    long long max_patches = 64;
    long long bundle_max_dim = 0;
    // "May this tool ask CCEX to build these clusters?" - the driver's budget
    // question, passed in so the judgement stays testable without one.
    std::function<bool(const Plan&, const std::vector<Sites>&)> affordable;
    std::string merge_priority = "gain-per-cost";
    const std::map<Sites, Curves>* raw = nullptr;
    bool gcce = false;
    double atol = 0.;
    double rtol = 0.;
    // The expansion limits come from the config, not from the plan the binary
    // answered with: Driver.approximate reads expansion_max_clusters and
    // expansion_max_work out of cfg every time.
    ExpansionLimits limits;
};

struct Coverage {
    long long candidate_patches = 0;
    long long checked_patches = 0;
    long long work_blocked = 0;
    bool coverage_limited = false;
    bool whole_bath_oracle = false;
    double max_inverse_log10_amplification = 0.;
};

// A merge changes dimensions only at its two endpoints, so the plan's dimensions
// and adjacency are worked out once and reused. Python caches this on the driver
// against the plan object; here the caller holds it and says which plan it is for.
struct Feasibility {
    const Plan* plan = nullptr;
    std::vector<long long> dimensions;
};

long long bundle_cap(const MergeContext& context, const Plan& plan);
bool can_merge(Feasibility* cache, const MergeContext& context, const Plan& plan,
               const std::vector<int>& members);

// The partition after joining these bundles, and the edges that survive it. Two
// is the ordinary case; a hub candidate joins more than two at once.
std::pair<std::vector<Sites>, std::vector<Edge>> merge_bundles(
    const std::vector<Sites>& bundles, const std::vector<Edge>& edges,
    const std::vector<int>& members);

std::pair<std::vector<std::vector<int>>, Coverage> patches(const MergeContext& context,
                                                           const Plan& plan);

// The order proposals are tried in. Returned as the pieces of Python's tuple so
// the caller can compare them the way Python compares tuples.
struct ProposalKey {
    bool has_worst_error = false;
    double worst_error = 0.;             // -local_before, when merge_priority says so
    double cost = 0.;                    // -score
    std::pair<int, int> pair;
};
// The error a patch shows when it is evaluated at the working order instead of
// exactly. `changed` replaces the patch's own bundles and edges, which is how a
// candidate merge is scored before it is asked for.
Discrepancy local(const MergeContext& context, const Plan& plan, const std::vector<int>& patch,
                  const std::pair<std::vector<Sites>, std::vector<Edge>>* changed);

// The curve that patch's model predicts, which local() compares. Split out so
// the physical check can read it without a second discrepancy.
Curves local_model(const MergeContext& context, const Plan& plan, const std::vector<int>& patch,
                   const std::pair<std::vector<Sites>, std::vector<Edge>>* changed);
bool is_physical(const MergeContext& context, const Curves& curves);

ProposalKey proposal_key(const MergeContext& context, double score, double local_before,
                         std::pair<int, int> pair);
bool before(const ProposalKey& left, const ProposalKey& right);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// The refinement driver's state and the parts of it that fetch probes
// (plan stage 1, Task 5; adaptive_refine.py).
//
// execute() is injected: the real one runs the binary (Task 6), a unit test hands
// in its own, and ensure()/read_probes() do not care which.
// ----------------------------------------------------------------------------
namespace adaptive {
// The curves a run wrote, one file per state.
Curves read_state_curves(const std::string& directory, const std::vector<int>& states);

// The response.rank-*.state-S.tsv files of a probe step (or of an except step, which writes its
// excepted terms the same way): one file per rank and state, every (listed cluster or L0, state,
// step) exactly once. `listed[i]` is probe i, L0 is probe -1 and comes back under the empty set.
std::map<Sites, Curves> read_probe_files(const std::string& directory, const std::vector<int>& states,
                                         long long steps, long long ranks, const std::vector<Sites>& listed);

// The typed plan inside a response.plan.json. execute() uses it on what the
// binary just answered.
Plan plan_of(const Value& raw);

struct DriverSettings {
    long long ranks = 1;
    // The D^3 * states * steps this tool may ask CCEX for.
    double max_total_work = 8e12;
    long long max_probes = 4096;
    long long cache_mb = 256;
    long long probe_batch_clusters = 0;
    std::string cwd;
};

// What a native execution answered with: where its files are, and the plan - both
// as the port uses it and as it was written, because parts of it (physics,
// dimension limits, selective clusters) travel on untouched.
struct Executed {
    std::string directory;
    Plan plan;
    Value raw_plan;
};

struct NativeSettings {
    std::string mpiexec = "mpirun";
    std::string binary;                  // resolved, the way Python resolves it
    double timeout = 600.;
    bool timeout_is_integer = true;      // 600 and 600.0 are different characters
    bool physics_metadata = false;
};

class Driver {
 public:
    // mode, the bundles to use (or none), the probes to ask for (or none).
    using Execute = std::function<Executed(const std::string& mode,
                                           const std::vector<Sites>* bundles,
                                           const std::vector<Sites>* probes)>;

    Driver(Value cfg, std::vector<int> states, DriverSettings settings, Execute execute)
        : cfg_(std::move(cfg)), states_(std::move(states)), settings_(std::move(settings)),
          execute_(std::move(execute)) {}

    // Makes sure `needed` is in raw, fetching what is missing. False means a
    // budget said no - the caller stops rather than carrying on with less.
    bool ensure(const std::vector<Sites>& needed, const std::vector<Sites>& bundles);
    void read_probes(const std::vector<Sites>& missing, const std::vector<Sites>& bundles);

    // Everything execute() needs beyond what ensure() does. Set once, before the
    // first call: the clock is read inside execute, so it has to be the same one
    // the recording was made with when replaying.
    void prepare(const std::string& root, NativeSettings native, Clock* clock);
    // mode, and the optional parts of the request. Empty means the call failed in
    // a way allow_failure permits.
    std::unique_ptr<Executed> execute(const std::string& mode,
                                      const std::vector<Sites>* bundles,
                                      const std::vector<Sites>* probes,
                                      bool allow_failure,
                                      const Value* more = nullptr);   // further request keys (factors, except)

    // Everything Driver.__init__ does that is not checking the config: the report
    // skeleton, the input hashes, and the cache. `settings` is the parsed command
    // line, which the report keeps as it is.
    void begin(const Value& settings, const std::string& binary_sha256);
    // The merge loop. Returns the exit status: 0 when the result is certified.
    // With Limits::regroup it is run_regroup() instead.
    int run();
    // Method D (policy_regroup.cpp, (internal notes)): 0
    // only when the curve ends physical.
    int run_regroup();
    int finish(const std::string& status, const std::string& directory,
               const Executed& step);
    void persist() const;
    void close_cache() { if (cache_) cache_->close(); }
    ExpansionLimits config_limits() const;

    // Can this tool ask CCEX to build these clusters? The question is asked in
    // the budget's own currency - D^3 * states * steps - rather than against a
    // separate dimension cap, so that the cost of a state grid and a time grid
    // counts.
    bool affordable(const Plan& plan, const std::vector<Sites>& clusters) const;
    long long work_left() const;

    Value& report() { return report_; }
    const std::vector<int>& states() const { return states_; }
    const Value& config() const { return cfg_; }

 protected:
    // read_probes goes through here so a test can stand in for the whole step
    // while the real run uses execute() above.
    std::unique_ptr<Executed> step_for(const std::string& mode, const std::vector<Sites>* bundles,
                                       const std::vector<Sites>* probes);
    // The probe cache part of ensure(): with a cache, check its inputs, take what it holds of `wanted` and the empty
    // cluster into raw_, and return false (nothing taken) when the curves held would pass max_probes or cache_mb.
    // Without a cache, true.
    bool load_cached(const std::set<Sites>& wanted);

 public:
    const std::map<Sites, Curves>& raw() const { return raw_; }
    std::map<Sites, Curves>& raw() { return raw_; }
    void use_cache(std::unique_ptr<ProbeCache> cache, Value context);
    long long stored_clusters() const { return stored_; }

 protected:
    Value cfg_;
    std::vector<int> states_;
    DriverSettings settings_;
    Execute execute_;
    std::map<Sites, Curves> raw_;
    std::string root_;
    long long counter_ = 0;
    NativeSettings native_;
    Clock* clock_ = nullptr;
    Value report_ = Value::object();
    Value physics_;
    bool has_physics_ = false;
    Value inputs_ = Value::object();     // path -> sha256, for the end-of-run check
    // Python appends the round to the report and then keeps writing into it, so a
    // persist() in the middle writes a half-filled row. Keeping the rows here and
    // rebuilding the list on every write reproduces that.
    std::vector<Value> rounds_;
    std::string cache_path_;
    std::string cache_key_;
    bool gcce_ = false;
    int order_ = 2;
    MergeContext merge_;
    // The candidate the merge loop never offers: the bundle every unphysical
    // candidate of this round shares, merged with its neighbours until the whole
    // model is physical. No members means the round has no such bundle.
    std::pair<std::vector<int>, Value> hub_candidate(
        const MergeContext& context, const Plan& plan,
        const std::map<std::pair<int, int>, Value>& candidates);

 public:
    // Set once from the command line, before run().
    struct Limits {
        long long max_rounds = 8;
        long long max_candidates = 16;
        double atol = 0.;
        double rtol = 0.;
        // Narrow the adoption to candidates whose predicted curve is physical,
        // when some are and some are not. Off by default.
        bool prefer_physical = false;
        // Needs prefer_physical: offer the bundle every unphysical candidate
        // shares, grown until the model is physical. Off by default.
        bool hub_candidates = false;
        // Method D in place of merge's own selection (--regroup). Off by default,
        // and off, nothing of it runs.
        bool regroup = false;
        long long regroup_max_regions = 64;
        long long regroup_max_candidates = 256;
        // --regroup-factors factors: the culprit factors from one factors step per round
        // (adaptive-regroup-stage3.md). --regroup-whole except: condition (c)'s whole curve
        // from one except step per region (adaptive-regroup-stage3b.md). Off: probes, as before.
        bool regroup_factors = false;
        bool regroup_except = false;
    };
    void set_policy(int order, bool gcce, MergeContext merge, Limits limits) {
        order_ = order; gcce_ = gcce; merge_ = std::move(merge); limits_ = limits;
    }

 protected:
    Limits limits_;
    std::unique_ptr<ProbeCache> cache_;
    Value cache_context_;
    long long loaded_ = 0;
    long long stored_ = 0;
};
}   // namespace adaptive
