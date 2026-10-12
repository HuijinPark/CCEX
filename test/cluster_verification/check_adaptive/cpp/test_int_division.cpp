// Integer arithmetic against CPython's own answers (plan 0B-2, Task 2).
#include "check.h"
#include "base_numeric.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::string bits_of(double value) {
    std::uint64_t raw;
    std::memcpy(&raw, &value, sizeof raw);
    if (std::isnan(value))
        raw = 0x7ff8000000000000ull;   // NaN equals NaN; its sign is not part of the contract

    char text[24];
    std::snprintf(text, sizeof text, "%016llx", static_cast<unsigned long long>(raw));
    return text;
}

// One case: what the operation produced, or "Type\tmessage" when it raised.
std::string outcome(const std::string& op, const std::string& a, const std::string& b) {
    try {
        if (op == "truediv")
            return bits_of(adaptive::true_div(adaptive::from_decimal(a), adaptive::from_decimal(b)));
        if (op == "comb")
            return adaptive::decimal(adaptive::comb(std::stoll(a), std::stoll(b)));
        return "unknown operation";
    } catch (const adaptive::Error& error) {
        return error.type()+"\t"+error.what();
    }
}
}   // namespace

int main(int argc, char* argv[]) {
    // A path argument points at a fuzz file in the same format.
    const std::string path = argc > 1 ? argv[1] : "../reference/semantics/int_ops.tsv";
    std::ifstream file(path);
    CHECK(file.good());
    std::string line;
    int rows = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream parts(line);
        std::string op, a, b, expected, message;
        std::getline(parts, op, '\t');
        std::getline(parts, a, '\t');
        std::getline(parts, b, '\t');
        std::getline(parts, expected, '\t');
        std::getline(parts, message);
        const std::string want = message.empty() ? expected : expected+"\t"+message;
        CHECK_EQ(outcome(op, a, b), want);
        ++rows;
    }
    CHECK(rows >= (argc > 1 ? 1 : 633));

    // The case that made the single rounding visible in the first place.
    const long long cube = 2097151LL*2097151LL*2097151LL;
    CHECK_EQ(bits_of(adaptive::true_div(static_cast<long long>(1), cube)), std::string("3c00000180001800"));

    // Exact comparison with a double: no conversion, either way.
    bool unordered = false;
    CHECK_EQ(adaptive::compare(adaptive::from_decimal("100000000000000001"), 1e17, &unordered), 1);
    CHECK_EQ(adaptive::compare(adaptive::from_decimal("9007199254740993"), 9007199254740992., &unordered), 1);
    CHECK_EQ(adaptive::compare(static_cast<long long>(1), 1., &unordered), 0);
    CHECK(!unordered);
    adaptive::compare(static_cast<long long>(1), std::nan(""), &unordered);
    CHECK(unordered);
    std::cout << "compared " << rows << " integer cases\n";
    return check::report("test_int_division");
}
