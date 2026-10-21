// The three counts that can outgrow int64, against CPython (plan stage 7).
//
// The reference is CPython's exact answer, which is unbounded. This port's
// integers are int64, so a row whose answer does not fit is one the port
// refuses rather than computes: the expectation for it is the OverflowError,
// not a number. Every count a real run reaches fits with room to spare - the
// widest integer across all 52 recordings is five bits.
#include "check.h"
#include "base_bundles.h"
#include "base_json.h"
#include "base_numeric.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
// Is the row's expectation a number this port can hold? A row expecting an
// error names the type instead, and that is not a number at all.
bool beyond_int64(const std::string& expected) {
    std::size_t at = !expected.empty() && expected[0] == '-' ? 1 : 0;
    if (at == expected.size()) return false;
    for (std::size_t i = at; i < expected.size(); ++i)
        if (expected[i] < '0' || expected[i] > '9') return false;
    try {
        adaptive::from_decimal(expected);
        return false;
    } catch (const adaptive::Error&) {
        return true;
    }
}

// "[2, 2, 2, 2, \"x62\"]" plus the count: the row writes the head and the length
// so a 786-entry list does not have to travel as text.
std::vector<long long> repeated_dims(const adaptive::Value& head, long long total) {
    std::vector<long long> out;
    for (long long i = 0; i < total; ++i)
        out.push_back(head.at(static_cast<std::size_t>(
                              i < 4 ? i : (head.size() > 4 ? 3 : 0))).as_int());
    return out;
}
}   // namespace

int main(int argc, char** argv) {
    std::ifstream file(argc > 1 ? argv[1] : "../reference/semantics/counts.tsv");
    std::string line;
    int rows = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream parts(line);
        std::vector<std::string> column;
        std::string word;
        while (std::getline(parts, word, '\t'))
            column.push_back(word);
        const std::size_t at = column[0] == "bundle_dimension" ? 5
                             : column[0] == "comb_halved" ? 3 : 0;
        if (at == 0) continue;
        std::string got, message;
        try {
            if (column[0] == "bundle_dimension") {
                const std::vector<long long> dims =
                    repeated_dims(adaptive::read_json(column[1]),
                                  adaptive::from_decimal(column[2]));
                const long long size =
                    adaptive::from_decimal(column[4]);
                const adaptive::Value head = adaptive::read_json(column[3]);
                adaptive::Sites bundle;
                for (long long i = 0; i < size; ++i)
                    bundle.push_back(static_cast<int>(
                        i < 4 ? head.at(static_cast<std::size_t>(i)).as_int() : i));
                got = adaptive::decimal(adaptive::bundle_dimension(bundle, dims));
            } else {
                const long long n = adaptive::from_decimal(column[1]);
                const long long k = adaptive::from_decimal(column[2]);
                got = adaptive::decimal(adaptive::comb(n, k)/static_cast<long long>(2*k == n ? 2 : 1));
            }
        } catch (const adaptive::Error& error) {
            got = error.type();
            message = error.what();
        }
        if (beyond_int64(column[at])) {
            // A dimension past int64 saturates, because it is only ever compared
            // against a budget; a count past int64 is refused, because callers
            // spend it. Both are deliberate (bigint.h).
            if (column[0] == "bundle_dimension") {
                CHECK_EQ(got, std::string("9223372036854775807"));
                CHECK_EQ(message, std::string());
            } else {
                CHECK_EQ(got, std::string("OverflowError"));
                CHECK_EQ(message, std::string("integer overflow"));
            }
        } else {
            CHECK_EQ(got, column[at]);
            CHECK_EQ(message, column.size() > at+1 ? column[at+1] : std::string());
        }
        ++rows;
    }
    CHECK(rows >= (argc > 1 ? 1 : 58));
    std::cout << "compared " << rows << " count cases\n";
    return check::report("test_counts");
}
