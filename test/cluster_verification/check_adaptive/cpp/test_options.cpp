// The 25 command line options, against their recorded parses (plan 0B-4, Task 5).
#include "check.h"
#include "options.h"
#include "base_json.h"
#include "options_table.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

// The parser's own view of a result, in the shape the expectation records it:
// "name=value" pairs, sorted, absent values left out, cwd excluded (it is the
// directory the tool was started in, not something the arguments decide).
std::string shown(const adaptive::Value& parsed) {
    std::vector<std::string> pairs;
    for (const std::string& key : parsed.keys()) {
        if (key == "cwd")
            continue;
        const adaptive::Value& value = parsed.at(key);
        switch (value.kind()) {
            case adaptive::Value::Kind::Null: continue;
            case adaptive::Value::Kind::Bool:
                pairs.push_back(key+"="+(value.as_bool() ? "True" : "False"));
                break;
            case adaptive::Value::Kind::Int:
                pairs.push_back(key+"="+adaptive::decimal(value.as_int()));
                break;
            case adaptive::Value::Kind::Double:
                pairs.push_back(key+"="+adaptive::float_text(value.as_double()));
                break;
            default:
                pairs.push_back(key+"="+value.as_string());
        }
    }
    std::sort(pairs.begin(), pairs.end());
    std::string out;
    for (std::size_t i = 0; i < pairs.size(); ++i)
        out += (i ? " " : "")+pairs[i];
    return out;
}
}   // namespace

int main(int argc, char** argv) {
    CHECK_EQ(adaptive::options::count(), std::size_t(27));
    // A corpus can be handed in, so the same driver takes the fuzz families.
    std::ifstream file(argc > 1 ? argv[1] : "../reference/semantics/options.tsv");
    CHECK(file.good());
    std::string line;
    int rows = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream parts(line);
        std::string op, argv_text, status, expected;
        std::getline(parts, op, '\t');
        std::getline(parts, argv_text, '\t');
        std::getline(parts, status, '\t');
        std::getline(parts, expected);
        const adaptive::Value argv_list = adaptive::read_json(argv_text);
        std::vector<std::string> argv;
        for (std::size_t i = 0; i < argv_list.size(); ++i)
            argv.push_back(argv_list.at(i).as_string());
        try {
            const adaptive::Value parsed = adaptive::parse_options(argv, "/work");
            CHECK_EQ(std::string("ok"), status);
            CHECK_EQ(shown(parsed), expected);
        } catch (const adaptive::ValueError&) {
            CHECK_EQ(std::string("error"), status);
        }
        ++rows;
    }
    CHECK(rows >= (argc > 1 ? 1 : 9));
    // The default the table marks as dynamic is filled in at start-up, not baked in.
    const adaptive::Value parsed = adaptive::parse_options(
        {"--input", "i", "--binary", "b", "--output", "o", "--atol", "1e-8", "--rtol", "2e-6"},
        "/work");
    CHECK_EQ(parsed.at("cwd").as_string(), std::string("/work"));
    std::cout << "compared " << rows << " option lists over "
              << adaptive::options::count() << " options\n";
    return check::report("test_options");
}
