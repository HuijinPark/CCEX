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
    cfg.set("bathfile", adaptive::Value::text("absent.xyz"));
    const adaptive::Value files = adaptive::input_manifest(cfg, "out/manifest", {});
    CHECK_EQ(files.keys().size(), std::size_t(1));
    CHECK(files.at(files.keys()[0]).kind() == adaptive::Value::Kind::Null);
    std::cout << "manifest: resolved one path, one absent input\n";
    return check::report("test_manifest");
}
