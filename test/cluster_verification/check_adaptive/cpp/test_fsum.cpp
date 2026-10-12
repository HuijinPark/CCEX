// fsum() against its recorded answers (plan 0B-2, Task 5).
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
        raw = 0x7ff8000000000000ull;
    char text[24];
    std::snprintf(text, sizeof text, "%016llx", static_cast<unsigned long long>(raw));
    return text;
}

std::vector<std::string> split(const std::string& text, char by) {
    std::vector<std::string> out;
    std::istringstream parts(text);
    std::string piece;
    while (std::getline(parts, piece, by))
        out.push_back(piece);
    return out;
}

std::vector<double> items_of(const std::string& text) {      // "f:<bits>|f:<bits>|..."
    std::vector<double> items;
    if (text.empty())
        return items;
    for (const std::string& piece : split(text, '|'))
        items.push_back(from_bits(piece.substr(2)));
    return items;
}

std::string outcome(const std::string& text) {
    try {
        return "f:"+bits_of(adaptive::fsum(items_of(text)));
    } catch (const adaptive::Error& error) {
        return error.type()+"\t"+error.what();
    }
}
}   // namespace

int main(int argc, char* argv[]) {
    const std::string path = argc > 1 ? argv[1] : "../reference/semantics/sums.tsv";
    std::ifstream file(path);
    CHECK(file.good());
    std::string line;
    int rows = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream parts(line);
        std::string op, items, expected, message;
        std::getline(parts, op, '\t');
        std::getline(parts, items, '\t');
        std::getline(parts, expected, '\t');
        std::getline(parts, message);
        const std::string want = message.empty() ? expected : expected+"\t"+message;
        CHECK_EQ(outcome(items), want);
        ++rows;
    }
    CHECK(rows >= (argc > 1 ? 1 : 14));
    std::cout << "compared " << rows << " sums\n";
    return check::report("test_fsum");
}
