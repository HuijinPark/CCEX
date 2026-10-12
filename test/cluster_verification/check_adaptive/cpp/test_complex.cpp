// math and complex edges against CPython (plan 0B-2, Task 4).
#include "check.h"
#include "base_numeric.h"
#include "base_json.h"
#include <gnu/libc-version.h>
#include <cmath>
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

std::string bits_of(double value) {
    std::uint64_t raw;
    std::memcpy(&raw, &value, sizeof raw);
    if (std::isnan(value))
        raw = 0x7ff8000000000000ull;   // NaN equals NaN; its sign is not part of the contract

    char text[24];
    std::snprintf(text, sizeof text, "%016llx", static_cast<unsigned long long>(raw));
    return text;
}

std::string two(const std::string& op, const std::string& ar, const std::string& ai,
                const std::string& br, const std::string& bi) {
    try {
        if (op == "abs") return bits_of(adaptive::complex_abs(from_bits(ar), from_bits(ai)));
        if (op == "div") {
            double re = 0., im = 0.;
            adaptive::complex_div(from_bits(ar), from_bits(ai), from_bits(br), from_bits(bi), &re, &im);
            return bits_of(re)+"\t"+bits_of(im);
        }
        return "unknown operation";
    } catch (const adaptive::Error& error) {
        return error.type()+"\t"+error.what();
    }
}
}   // namespace

namespace {
// How far apart two doubles are, counted in representable steps.
std::uint64_t ulp_distance(const std::string& a, const std::string& b) {
    auto ordered = [](const std::string& hex) {
        const std::uint64_t raw = std::stoull(hex, nullptr, 16);
        return raw & 0x8000000000000000ull ? 0x8000000000000000ull-(raw & 0x7fffffffffffffffull)
                                           : raw+0x8000000000000000ull;
    };
    const std::uint64_t x = ordered(a), y = ordered(b);
    return x > y ? x-y : y-x;
}

bool is_bits(const std::string& text) {
    return text.size() == 16 && text.find_first_not_of("0123456789abcdef") == std::string::npos;
}
}   // namespace

// Is this the machine the expectations were recorded on? Only there is a bit-exact
// comparison meaningful for the functions that come from libm. Everywhere else the
// comparison is a reference check, and one step apart is reported rather than
// failed - but only for a value against a value, never for an exception, and never
// between a finite number and an infinity (11th review).
static bool same_numeric_environment() {
    std::ifstream file("../reference/semantics/environment.json");
    if (!file.good())
        return false;
    std::ostringstream text;
    text << file.rdbuf();
    try {
        const adaptive::Value recorded = adaptive::read_json(text.str());
        return recorded.at("libc").as_string() == "glibc" &&
               recorded.at("libc_version").as_string() == gnu_get_libc_version();
    } catch (const adaptive::Error&) {
        return false;
    }
}

int main(int argc, char* argv[]) {
    int rows = 0, libm_differences = 0;
    std::string libm_note;
    const bool strict = same_numeric_environment();
    for (const std::string& path : {std::string(argc > 1 ? argv[1] : "../reference/semantics/complex_ops.tsv")}) {
        std::ifstream file(path);
        CHECK(file.good());
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            std::istringstream parts(line);
            std::vector<std::string> field;
            std::string piece;
            while (std::getline(parts, piece, '\t'))
                field.push_back(piece);
            field.resize(8);          // trailing empty columns are not produced by getline
            const std::string& op = field[0];
            // op, ar, ai, br, bi, result..., message
            const std::string want = op == "div"
                ? (field[7].empty() ? field[5]+"\t"+field[6] : field[5]+"\t"+field[7])
                : (field.size() > 6 && !field[6].empty() ? field[5]+"\t"+field[6] : field[5]);
            const std::string found = two(op, field[1], field[2], field[3], field[4]);
            // hypot() comes from the machine's libm, and libm differs
            // between machines: idm and idm4-gpu answer one step away from WSL for
            // some inputs (measured 2026-09-16). CPython on those machines would
            // differ in exactly the same way, because it calls the same libm. So a
            // single step apart is recorded, not failed; anything further fails.
            if (!strict && found != want && op == "abs" &&
                is_bits(found) && is_bits(want) && ulp_distance(found, want) <= 1) {
                ++libm_differences;
                if (libm_note.empty()) libm_note = op;
            } else {
                CHECK_EQ(found, want);
            }
            ++rows;
        }
    }
    CHECK(rows >= (argc > 1 ? 1 : 208));
    std::cout << "compared " << rows << " float and complex cases"
              << (strict ? " (strict: the recording environment)" : " (another environment)");
    if (libm_differences > 0)
        std::cout << " (" << libm_differences << " one step from the recorded value, this machine's"
                  << " libm, first in " << libm_note << ")";
    std::cout << "\n";
    return check::report("test_complex");
}
