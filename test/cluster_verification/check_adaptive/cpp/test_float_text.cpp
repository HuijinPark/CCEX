// Float text, the JSON writer and SHA-256 against recorded answers (plan 0B-4, Tasks 1-3).
#include "check.h"
#include <cmath>
#include "base_json.h"
#include "base_run.h"
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

namespace {

double from_bits(const std::string& hex) {
    const std::uint64_t raw = std::stoull(hex, nullptr, 16);
    double value;
    std::memcpy(&value, &raw, sizeof value);
    return value;
}

std::string unquote(const std::string& quoted) {   // Python repr of a plain string
    std::string out;
    for (std::size_t i = 1; i+1 < quoted.size(); ++i) {
        if (quoted[i] != '\\') { out.push_back(quoted[i]); continue; }
        const char kind = quoted[++i];
        if (kind == 'n') out.push_back('\n');
        else if (kind == 't') out.push_back('\t');
        else if (kind == 'r') out.push_back('\r');
        else if (kind == 'x') { out.push_back(static_cast<char>(std::stoul(quoted.substr(i+1, 2), nullptr, 16))); i += 2; }
        else out.push_back(kind);
    }
    return out;
}
}   // namespace

std::string to_bits(double value) {
    std::uint64_t word = 0;
    std::memcpy(&word, &value, sizeof word);
    char text[32];
    std::snprintf(text, sizeof text, "%016llx", static_cast<unsigned long long>(word));
    return text;
}

int main(int argc, char* argv[]) {
    const std::string path = argc > 1 ? argv[1] : "../reference/semantics/texts.tsv";
    std::ifstream file(path);
    CHECK(file.good());
    std::string line;
    int rows = 0, longer = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream parts(line);
        std::string op, argument, expected;
        std::getline(parts, op, '\t');
        std::getline(parts, argument, '\t');
        std::getline(parts, expected);
        if (op == "repr") {
            // `expected` is the shortest text there is (CPython's). Ours has to read
            // back as the same bits and may be one digit longer, at a power of two.
            const std::string text = adaptive::float_text(from_bits(argument));
            CHECK_EQ(to_bits(adaptive::parse_float(text)), argument);
            CHECK(text.size() <= expected.size()+1);
            longer += text != expected;
        } else if (op == "dumps_compact") {
            CHECK_EQ(adaptive::json_dumps(adaptive::read_json(argument), true, true, false), expected);
        } else if (op == "dumps_default") {
            CHECK_EQ(adaptive::json_dumps(adaptive::read_json(argument), false, false, true), expected);
        } else if (op == "float_str") {
            // The way back: CPython promises float(repr(x)) is x, and the port has
            // to keep that promise with the same bits, not merely the same value.
            CHECK_EQ(to_bits(adaptive::parse_float(argument)), expected);
        } else if (op == "sha256") {
            const std::string payload =
                adaptive::json_dumps(adaptive::read_json(argument), true, true, false);
            CHECK_EQ(adaptive::sha256(payload), expected);
        } else if (op == "sha256_text") {
            CHECK_EQ(adaptive::sha256(unquote(argument)), expected);
        }
        ++rows;
    }
    CHECK(rows >= (argc > 1 ? 1 : 59));
    // A non-finite float is refused where the cache refuses it.
    adaptive::Value one = adaptive::Value::object();
    one.set("x", adaptive::Value::number(HUGE_VAL));
    CHECK_THROWS(adaptive::json_dumps(one, true, true, false), adaptive::ValueError,
                 "Out of range float values are not JSON compliant");
    std::cout << "compared " << rows << " text cases, " << longer << " floats written a digit longer\n";
    return check::report("test_float_text");
}
