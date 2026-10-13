// The probe cache (plan 0B-4, Task 4).
#include "check.h"
#include "base_run.h"
#include "base_json.h"
#include <cstdio>
#include <string>
#include <vector>

namespace {

// The conditions every cache in this test is opened with.
adaptive::Value context_value() {
    adaptive::Value config = adaptive::Value::object();
    config.set("nstep", adaptive::Value::integer(static_cast<long long>(3)));
    config.set("deltat", adaptive::Value::number(0.01));
    config.set("method", adaptive::Value::text("cCCE"));
    config.set("order", adaptive::Value::integer(static_cast<long long>(2)));
    config.set("seed", adaptive::Value::integer(static_cast<long long>(19)));
    config.set("rdip", adaptive::Value::number(30.0));
    adaptive::Value pulse = adaptive::Value::object();
    pulse.set("count", adaptive::Value::integer(static_cast<long long>(1)));
    config.set("pulse", pulse);
    adaptive::Value states = adaptive::Value::array();
    states.push(adaptive::Value::integer(static_cast<long long>(0)));
    states.push(adaptive::Value::integer(static_cast<long long>(2)));
    adaptive::Value files = adaptive::Value::object();
    files.set("/tmp/ccex-roundtrip/bath.xyz", adaptive::Value::text(std::string(64, 'b')));
    adaptive::Value context = adaptive::Value::object();
    context.set("schema", adaptive::Value::integer(static_cast<long long>(1)));
    context.set("config", config);
    context.set("cwd", adaptive::Value::text("/tmp/ccex-roundtrip"));
    context.set("states", states);
    context.set("ranks", adaptive::Value::integer(static_cast<long long>(1)));
    context.set("binary_sha256", adaptive::Value::text(std::string(64, 'a')));
    context.set("files", files);
    return context;
}

const std::vector<adaptive::Sites> kSites = {{0}, {1, 2}, {3, 4, 5}};

adaptive::Curves curves_for(const adaptive::Sites& sites) {
    const double base[5] = {0.1, 1e-300, -0.0, 1e16, 3.141592653589793};
    adaptive::Curves out;
    for (int state : {0, 2}) {
        adaptive::Curve curve;
        for (int step = 0; step < 3; ++step)
            curve.emplace_back(base[(sites.size()+static_cast<std::size_t>(state+step)) % 5],
                               -base[static_cast<std::size_t>(step+1) % 5]);
        out.push_back(curve);
    }
    return out;
}

}   // namespace

int main() {
    const adaptive::Value context = context_value();
    const std::string path = "out/cache/probe.sqlite";
    std::remove(path.c_str());
    {
        adaptive::ProbeCache cache(path, context);      // creates its own parent directory
        CHECK_EQ(cache.get(kSites).size(), std::size_t(0));          // empty: a miss
        std::map<adaptive::Sites, adaptive::Curves> values;
        values[kSites[0]] = curves_for(kSites[0]);
        cache.put(values);
        CHECK_EQ(cache.get(kSites).size(), std::size_t(1));          // and now a hit
        cache.put(values);                                            // writing the same twice is fine
        CHECK_EQ(cache.get(kSites).size(), std::size_t(1));
        // The same conditions must not produce two different results.
        std::map<adaptive::Sites, adaptive::Curves> other;
        other[kSites[0]] = curves_for(kSites[1]);
        CHECK_THROWS(cache.put(other), adaptive::RuntimeError,
                     "raw cache has conflicting results for identical conditions");
    }
    {
        adaptive::ProbeCache reopened(path, context);   // survives a reopen
        CHECK_EQ(reopened.get(kSites).size(), std::size_t(1));
    }
    {
        // A context that differs by one character is a different context.
        adaptive::Value changed = context_value();
        adaptive::Value config = changed.at("config");
        config.set("deltat", adaptive::Value::number(0.010000000000000002));
        changed.set("config", config);
        adaptive::ProbeCache separate(path, changed);
        CHECK(separate.key() != adaptive::sha256(adaptive::json_dumps(context, true, true, false)));
        CHECK_EQ(separate.get(kSites).size(), std::size_t(0));
    }
    std::cout << "probe cache checked\n";
    return check::report("test_probe_cache");
}
