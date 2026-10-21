// The scaled complex product, against CPython (plan stage 1, Task 1).
#include "check.h"
#include "base_numeric.h"
#include "base_bundles.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
double from_bits(const std::string& text) {
    unsigned long long word = std::stoull(text, nullptr, 16);
    double value = 0.;
    std::memcpy(&value, &word, sizeof value);
    return value;
}

std::string to_bits(double value) {
    if (value != value)
        return "7ff8000000000000";
    unsigned long long word = 0;
    std::memcpy(&word, &value, sizeof word);
    char text[32];
    std::snprintf(text, sizeof text, "%016llx", word);
    return text;
}

std::string shown(const adaptive::Scaled& a) {
    return to_bits(a.z.real())+":"+to_bits(a.z.imag())+":"+adaptive::decimal(a.exponent);
}

std::string shown(const std::complex<double>& z) {
    return to_bits(z.real())+":"+to_bits(z.imag());
}
}   // namespace

int main(int argc, char** argv) {
    std::ifstream file(argc > 1 ? argv[1] : "../reference/semantics/scaled.tsv");
    std::string line;
    int rows = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream parts(line);
        std::vector<std::string> column;
        std::string word;
        while (std::getline(parts, word, '\t'))
            column.push_back(word);
        const std::string& op = column[0];
        const std::size_t at = op == "power" ? 5 : op == "multiply" ? 7 : 4;
        const std::string expected = column.size() > at ? column[at] : std::string();
        const std::string expected_message = column.size() > at+1 ? column[at+1] : std::string();
        std::string got, message;
        // The same shape as the other families: a value, or the exception's type
        // and message. inverse() really does divide by zero - a subnormal's norm
        // underflows - so the failure is part of the contract, not an accident.
        try {
            if (op == "normalize") {
                got = shown(adaptive::normalize({from_bits(column[1]), from_bits(column[2])},
                                              adaptive::from_decimal(column[3])));
            } else if (op == "inverse") {
                got = shown(adaptive::inverse({{from_bits(column[1]), from_bits(column[2])},
                                             adaptive::from_decimal(column[3])}));
            } else if (op == "value") {
                got = shown(adaptive::value({{from_bits(column[1]), from_bits(column[2])},
                                           adaptive::from_decimal(column[3])}));
            } else if (op == "power") {
                got = shown(adaptive::power({{from_bits(column[1]), from_bits(column[2])},
                                           adaptive::from_decimal(column[3])}, adaptive::from_decimal(column[4])));
            } else if (op == "multiply") {
                got = shown(adaptive::multiply({{from_bits(column[1]), from_bits(column[2])},
                                              adaptive::from_decimal(column[3])},
                                             {{from_bits(column[4]), from_bits(column[5])},
                                              adaptive::from_decimal(column[6])}));
            } else {
                continue;
            }
        } catch (const adaptive::Error& error) {
            got = error.type();
            message = error.what();
        }
        CHECK_EQ(got, expected);
        CHECK_EQ(message, expected_message);
        ++rows;
    }
    CHECK(rows >= (argc > 1 ? 1 : 2648));
    std::cout << "compared " << rows << " scaled cases\n";
    return check::report("test_scaled");
}
