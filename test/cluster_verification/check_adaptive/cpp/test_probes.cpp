// Fetching probes and the budgets that stop it (plan stage 1, Task 5).
// The recorded probe steps are the input: what the port reads out of them has to
// be what Python read.
#include "check.h"
#include "mkdirs.h"
#include "policy_merge.h"
#include "base_run.h"
#include "base_json.h"
#include "base_numeric.h"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {
std::string slurp(const std::string& path) {
    std::ifstream file(path);
    std::stringstream text;
    text << file.rdbuf();
    return text.str();
}

}   // namespace

int main(int argc, char** argv) {
    // A recorded probe step from merge-grow: six singleton bundles, three probes.
    const std::string root = argc > 1 ? argv[1]
                                      : "../reference/probe-step";
    const adaptive::Value request = adaptive::read_json(slurp(root+"/request.json"));
    std::vector<adaptive::Sites> bundles, probes;
    for (std::size_t i = 0; i < request.at("bundles").size(); ++i)
        bundles.push_back(sites_of(request.at("bundles").at(i)));
    for (std::size_t i = 0; i < request.at("probes").size(); ++i)
        probes.push_back(sites_of(request.at("probes").at(i)));
    CHECK(!probes.empty());

    adaptive::Value cfg = adaptive::Value::object();
    cfg.set("nstep", adaptive::Value::integer(static_cast<long long>(9)));
    adaptive::DriverSettings settings;
    settings.ranks = 1;
    settings.cwd = "/tmp";

    int asked = 0;
    adaptive::Driver::Execute execute =
        [&](const std::string& mode, const std::vector<adaptive::Sites>*,
            const std::vector<adaptive::Sites>* wanted) {
            CHECK_EQ(mode, std::string("probe"));
            CHECK(wanted != nullptr);
            ++asked;
            return adaptive::Executed{root, adaptive::Plan()};
        };

    // The recording answers with exactly these probes, so reading it back fills
    // the empty cluster and every one that was asked for.
    adaptive::Driver driver(cfg, {0}, settings, execute);
    driver.read_probes(probes, bundles);
    CHECK_EQ(asked, 1);
    CHECK_EQ(driver.raw().size(), probes.size()+1);
    CHECK(driver.raw().count(adaptive::Sites()) == 1);
    for (const adaptive::Sites& probe : probes)
        CHECK(driver.raw().count(probe) == 1);
    for (const std::pair<const adaptive::Sites, adaptive::Curves>& entry : driver.raw()) {
        CHECK_EQ(entry.second.size(), std::size_t(1));
        CHECK_EQ(entry.second[0].size(), std::size_t(9));
    }
    // Reading the same step again must agree with itself.
    driver.read_probes(probes, bundles);
    CHECK_EQ(asked, 2);

    // A rank count the recording cannot satisfy is a missing response, not a
    // partial read.
    adaptive::DriverSettings two = settings;
    two.ranks = 2;
    adaptive::Driver pair(cfg, {0}, two, execute);
    CHECK_THROWS(pair.read_probes(probes, bundles), adaptive::RuntimeError,
                 "missing MPI probe response");

    // ensure() asks for what it does not have, and the budgets stop it.
    adaptive::Driver budgeted(cfg, {0}, settings, execute);
    asked = 0;
    CHECK(budgeted.ensure(probes, bundles));
    CHECK_EQ(asked, 1);
    CHECK(budgeted.ensure(probes, bundles));           // already held: nothing asked
    CHECK_EQ(asked, 1);

    adaptive::DriverSettings tight = settings;
    tight.max_probes = 2;
    adaptive::Driver few(cfg, {0}, tight, execute);
    asked = 0;
    CHECK(!few.ensure(probes, bundles));               // three probes will not fit
    CHECK_EQ(asked, 0);

    adaptive::DriverSettings small = settings;
    small.cache_mb = 0;
    adaptive::Driver tiny(cfg, {0}, small, execute);
    CHECK(!tiny.ensure(probes, bundles));

    // The rest needs a response for whatever is asked, which one recorded step
    // cannot give: batching asks three times for one probe each.
    testing::mkdirs("out/probes");
    int made = 0;
    auto answering = [&](double value, bool drop_row, bool repeat_row) {
        return [&, value, drop_row, repeat_row](const std::string&,
                                                const std::vector<adaptive::Sites>*,
                                                const std::vector<adaptive::Sites>* wanted) {
            const std::string here = "out/probes/step-"+std::to_string(++made);
            testing::mkdirs(here);
            std::string text = "probe\tstep\treal\timag\n";
            for (long long index = -1; index < static_cast<long long>(wanted->size()); ++index)
                for (int step = 0; step < 9; ++step) {
                    if (drop_row && index == 0 && step == 4)
                        continue;
                    char line[128];
                    std::snprintf(line, sizeof line, "%lld\t%d\t%.17g\t0\n", index, step,
                                  index == -1 ? 1. : value);
                    text += line;
                    if (repeat_row && index == 0 && step == 0)
                        text += line;
                }
            std::FILE* file = std::fopen((here+"/response.rank-0.state-0.tsv").c_str(), "wb");
            CHECK(file != nullptr);
            if (file != nullptr) {
                std::fwrite(text.data(), 1, text.size(), file);
                std::fclose(file);
            }
            ++asked;
            return adaptive::Executed{here, adaptive::Plan()};
        };
    };

    // A batch size splits the asking without changing what is held.
    adaptive::DriverSettings batched = settings;
    batched.probe_batch_clusters = 1;
    adaptive::Driver stepwise(cfg, {0}, batched, answering(.5, false, false));
    asked = 0;
    CHECK(stepwise.ensure(probes, bundles));
    CHECK_EQ(asked, static_cast<int>(probes.size()));
    CHECK_EQ(stepwise.raw().size(), probes.size()+1);

    // A row that never came, and one that came twice.
    adaptive::Driver gap(cfg, {0}, settings, answering(.5, true, false));
    CHECK_THROWS(gap.read_probes(probes, bundles), adaptive::RuntimeError, "incomplete raw probe");
    adaptive::Driver twice(cfg, {0}, settings, answering(.5, false, true));
    CHECK_THROWS(twice.read_probes(probes, bundles), adaptive::RuntimeError, "duplicate raw probe");

    // The same cluster answered differently the second time round.
    adaptive::Driver changed(cfg, {0}, settings, answering(.5, false, false));
    changed.read_probes(probes, bundles);
    adaptive::Driver::Execute other = answering(.75, false, false);
    adaptive::Driver again(cfg, {0}, settings, other);
    again.raw() = changed.raw();
    CHECK_THROWS(again.read_probes(probes, bundles), adaptive::RuntimeError,
                 "L0 changed between identical-state runs");

    // With a cache in hand, the input files are checked again at the two points
    // Python checks them: before the cache is read, and before a probe is stored.
    // This is what makes a changed input a failure rather than a stale answer.
    testing::mkdirs("out/probes/inputs");
    const std::string bath = "out/probes/inputs/bath.xyz";
    {
        std::FILE* file = std::fopen(bath.c_str(), "wb");
        CHECK(file != nullptr);
        if (file != nullptr) { std::fputs("one\n", file); std::fclose(file); }
    }
    adaptive::Value cached_cfg = cfg;
    adaptive::Value configuration = adaptive::Value::object();   // the configuration block (A64)
    configuration.set("bathfile", adaptive::Value::text("inputs/bath.xyz"));
    cached_cfg.set("configuration", configuration);
    adaptive::DriverSettings with_cache = settings;
    with_cache.cwd = "out/probes";
    const adaptive::Value context = adaptive::identity(cached_cfg, with_cache.cwd, {0}, "c", 1);
    std::remove("out/probes/cache.sqlite");
    adaptive::Driver cached(cached_cfg, {0}, with_cache, answering(.5, false, false));
    cached.use_cache(std::unique_ptr<adaptive::ProbeCache>(
                         new adaptive::ProbeCache("out/probes/cache.sqlite", context)), context);
    CHECK(cached.ensure(probes, bundles));
    CHECK_EQ(cached.stored_clusters(), static_cast<long long>(probes.size()));

    {
        std::FILE* file = std::fopen(bath.c_str(), "wb");
        if (file != nullptr) { std::fputs("two\n", file); std::fclose(file); }
    }
    adaptive::Driver after(cached_cfg, {0}, with_cache, answering(.5, false, false));
    after.use_cache(std::unique_ptr<adaptive::ProbeCache>(
                        new adaptive::ProbeCache("out/probes/cache.sqlite", context)), context);
    CHECK_THROWS(after.ensure(probes, bundles), adaptive::RuntimeError,
                 "an input file changed before raw cache access");

    std::cout << "read " << probes.size() << " probes from a recorded step\n";
    return check::report("test_probes");
}
