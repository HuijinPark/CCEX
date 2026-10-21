// int(str) and float(str) against CPython (plan 0B-2, Task 3).
#include "check.h"
#include "base_json.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

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

// The expectation file holds Python's repr() of the input; this reads it back.
std::string from_repr(const std::string& quoted) {
    if (quoted.size() < 2)
        return quoted;
    const char quote = quoted[0];
    std::string out;
    for (std::size_t i = 1; i+1 < quoted.size(); ++i) {
        if (quoted[i] != '\\') {
            out.push_back(quoted[i]);
            continue;
        }
        const char kind = quoted[++i];
        if (kind == 'n') out.push_back('\n');
        else if (kind == 't') out.push_back('\t');
        else if (kind == 'r') out.push_back('\r');
        else if (kind == '\\' || kind == quote) out.push_back(kind);
        else if (kind == 'x' || kind == 'u' || kind == 'U') {
            const int width = kind == 'x' ? 2 : kind == 'u' ? 4 : 8;
            const unsigned point = std::stoul(quoted.substr(i+1, width), nullptr, 16);
            i += width;
            if (point < 0x80) {
                out.push_back(static_cast<char>(point));
            } else if (point < 0x800) {
                out.push_back(static_cast<char>(0xc0 | (point >> 6)));
                out.push_back(static_cast<char>(0x80 | (point & 0x3f)));
            } else if (point < 0x10000) {
                out.push_back(static_cast<char>(0xe0 | (point >> 12)));
                out.push_back(static_cast<char>(0x80 | ((point >> 6) & 0x3f)));
                out.push_back(static_cast<char>(0x80 | (point & 0x3f)));
            } else {
                out.push_back(static_cast<char>(0xf0 | (point >> 18)));
                out.push_back(static_cast<char>(0x80 | ((point >> 12) & 0x3f)));
                out.push_back(static_cast<char>(0x80 | ((point >> 6) & 0x3f)));
                out.push_back(static_cast<char>(0x80 | (point & 0x3f)));
            }
        } else {
            out.push_back(kind);
        }
    }
    return out;
}

std::string outcome(const std::string& op, const std::string& text) {
    try {
        return op == "int" ? adaptive::decimal(adaptive::parse_int(text))
                           : bits_of(adaptive::parse_float(text));
    } catch (const adaptive::Error& error) {
        return error.type()+"\t"+error.what();
    }
}
}   // namespace

int main(int argc, char* argv[]) {
    const std::string path = argc > 1 ? argv[1] : "../reference/semantics/number_strings.tsv";
    std::ifstream file(path);
    CHECK(file.good());
    std::string line;
    int rows = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream parts(line);
        std::string op, quoted, expected, message;
        std::getline(parts, op, '\t');
        std::getline(parts, quoted, '\t');
        std::getline(parts, expected, '\t');
        std::getline(parts, message);
        const std::string want = message.empty() ? expected : expected+"\t"+message;
        CHECK_EQ(outcome(op, from_repr(quoted)), want);
        ++rows;
    }
    CHECK(rows >= (argc > 1 ? 1 : 90));
    std::cout << "compared " << rows << " number strings\n";
    return check::report("test_number_text");
}
