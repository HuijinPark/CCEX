// Input manifests (plan stage 1, Task 5).
#include "check.h"
#include "mkdirs.h"
#include "base_run.h"
#include <iostream>
#include <string>

int main() {
    // A path comes back absolute, a missing file is recorded as absent, and
    // neither crashes.
    testing::mkdirs("out/manifest");
    const std::string here = adaptive::resolved("out/manifest", "kept.txt");
    CHECK(!here.empty() && here[0] == '/');
    CHECK_EQ(here.substr(here.size()-9), std::string("/kept.txt"));
    adaptive::Value cfg = adaptive::Value::object();
    adaptive::Value configuration = adaptive::Value::object();   // the configuration block (A64)
    configuration.set("bathfile", adaptive::Value::text("absent.xyz"));
    cfg.set("configuration", configuration);
    const adaptive::Value files = adaptive::input_manifest(cfg, "out/manifest", {});
    CHECK_EQ(files.keys().size(), std::size_t(1));
    CHECK(files.at(files.keys()[0]).kind() == adaptive::Value::Kind::Null);
    // configuration.path is the folder of the block's relative names
    configuration.set("path", adaptive::Value::text("sub"));
    cfg.set("configuration", configuration);
    const adaptive::Value inFolder = adaptive::input_manifest(cfg, "out/manifest", {});
    CHECK_EQ(inFolder.keys().size(), std::size_t(1));
    const std::string key = inFolder.keys()[0];
    CHECK(key.size() > 15 && key.substr(key.size()-15) == "/sub/absent.xyz");
    // A73: "" in a per-bath-file list is a bath file without such a file, not a path (it was cwd/sub/1, cwd/sub/2)
    adaptive::Value states = adaptive::Value::array();
    states.push(adaptive::Value::text("st_"));
    states.push(adaptive::Value::text(""));
    configuration.set("bath_state_file", states);
    cfg.set("configuration", configuration);
    const adaptive::Value withEmpty = adaptive::input_manifest(cfg, "out/manifest", {1, 2});
    CHECK_EQ(withEmpty.keys().size(), std::size_t(3));   // the bath file, st_1, st_2
    for (const std::string& name : withEmpty.keys())
        CHECK(name.substr(name.size()-6) != "/sub/1" && name.substr(name.size()-6) != "/sub/2");
    // A73: keys as CCEX reads them, without regard to case: the "Output" block the input spelled is the output block
    adaptive::Value spelled = adaptive::Value::object();
    adaptive::Value output = adaptive::Value::object();
    output.set("File", adaptive::Value::text("user"));
    spelled.set("Output", output);
    CHECK_EQ(adaptive::key_of(spelled, "output"), std::string("Output"));
    CHECK_EQ(adaptive::key_of(spelled, "cluster"), std::string("cluster"));
    const adaptive::Value* file = adaptive::find_path(spelled, {"output", "file"});
    CHECK(file != nullptr && file->as_string() == "user");
    std::cout << "manifest: resolved one path, one absent input, configuration.path in front, \"\" skipped; key_of\n";
    return check::report("test_manifest");
}
