// Term counting and curve comparison, against CPython (plan stage 1, Task 2).
#include "check.h"
#include "base_bundles.h"
#include "base_json.h"
#include "base_numeric.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

namespace {
double from_bits(const std::string& text) {
    unsigned long long word = std::stoull(text, nullptr, 16);
    double value = 0.;
    std::memcpy(&value, &word, sizeof value);
    return value;
}

std::string bits_of(double value) {
    if (value != value)
        return "7ff8000000000000";
    unsigned long long word = 0;
    std::memcpy(&word, &value, sizeof word);
    char text[32];
    std::snprintf(text, sizeof text, "%016llx", word);
    return text;
}

adaptive::Curves curves_of(const adaptive::Value& listed) {
    adaptive::Curves out;
    for (std::size_t i = 0; i < listed.size(); ++i) {
        adaptive::Curve curve;
        for (std::size_t t = 0; t < listed.at(i).size(); ++t)
            curve.emplace_back(from_bits(listed.at(i).at(t).at(0).as_string()),
                               from_bits(listed.at(i).at(t).at(1).as_string()));
        out.push_back(curve);
    }
    return out;
}

// json.dumps of the dict discrepancy() returns, with default=str for the floats
// Python prints through repr.
std::string shown(const adaptive::Discrepancy& found) {
    adaptive::Value out = adaptive::Value::object();
    // json.dumps writes these as numbers; default=str never sees a float.
    out.set("scaled_error", adaptive::Value::number(found.scaled_error));
    out.set("max_abs_error", adaptive::Value::number(found.max_abs_error));
    if (found.has_location) {
        adaptive::Value where = adaptive::Value::array();
        where.push(adaptive::Value::integer(static_cast<long long>(found.state)));
        where.push(adaptive::Value::integer(static_cast<long long>(found.step)));
        out.set("worst_state_step", where);
    } else {
        out.set("worst_state_step", adaptive::Value());
    }
    return adaptive::json_dumps(out, false, false, true);
}
}   // namespace

int main(int argc, char** argv) {
    std::ifstream file(argc > 1 ? argv[1] : "../reference/semantics/expansion.tsv");
    std::string line;
    int rows = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream parts(line);
        std::vector<std::string> column;
        std::string word;
        while (std::getline(parts, word, '\t'))
            column.push_back(word);
        std::string got, message;
        // Where the expected value sits has to be known before the call: an
        // exception thrown inside it would otherwise leave the index at zero.
        const std::size_t at = column[0] == "bundle_powers" ? 4
                             : column[0] == "bundle_powers_limited" ? 3
                             : column[0] == "bundle_powers_extras" ? 5
                             : column[0] == "bundle_powers_free" ? 6
                             : column[0] == "discrepancy" ? 5
                             : column[0] == "approximate" ? 6 : 0;
        try {
            if (column[0] == "bundle_powers" || column[0] == "bundle_powers_limited" ||
                column[0] == "bundle_powers_extras" || column[0] == "bundle_powers_free") {
                int n = 4, order = 2;
                std::vector<adaptive::Edge> edges{{0, 1}, {1, 2}, {2, 3}, {0, 3}};
                std::vector<std::vector<int> > extras;
                adaptive::ExpansionLimits limits;
                if (column[0] != "bundle_powers_limited") {
                    n = static_cast<int>(adaptive::from_decimal(column[1]));
                    const adaptive::Value listed = adaptive::read_json(column[2]);
                    edges.clear();
                    for (std::size_t i = 0; i < listed.size(); ++i)
                        edges.emplace_back(static_cast<int>(listed.at(i).at(0).as_int()),
                                           static_cast<int>(listed.at(i).at(1).as_int()));
                    order = static_cast<int>(adaptive::from_decimal(column[3]));
                    if (column[0] == "bundle_powers_free") {
                        const adaptive::Value given = adaptive::read_json(column[4]);
                        for (const std::string& key : given.kind() ==
                                 adaptive::Value::Kind::Object ? given.keys()
                                                             : std::vector<std::string>())
                            (key == "max_clusters" ? limits.max_clusters : limits.max_work) =
                                given.at(key).as_int();
                    }
                    if (column[0] == "bundle_powers_extras" ||
                        column[0] == "bundle_powers_free") {
                        const adaptive::Value given =
                            adaptive::read_json(column[column[0] == "bundle_powers_free" ? 5 : 4]);
                        for (std::size_t i = 0; i < given.size(); ++i) {
                            std::vector<int> one;
                            for (std::size_t k = 0; k < given.at(i).size(); ++k)
                                one.push_back(static_cast<int>(given.at(i).at(k).as_int()));
                            extras.push_back(one);
                        }
                    }
                } else {
                    const adaptive::Value given = adaptive::read_json(column[1]);
                    // A limit that is not a positive integer is refused, and the
                    // float 2.0 is refused too: Python checks type, not value.
                    for (const std::string& key : given.keys()) {
                        const adaptive::Value& v = given.at(key);
                        if (v.kind() != adaptive::Value::Kind::Int || v.as_int() < 1)
                            throw adaptive::ValueError("expansion limits must be positive integers");
                        (key == "max_clusters" ? limits.max_clusters : limits.max_work) = v.as_int();
                    }
                    order = static_cast<int>(adaptive::from_decimal(column[2]));
                }
                const std::map<std::vector<int>, int> powers =
                    adaptive::bundle_powers(n, edges, order, limits, extras);
                adaptive::Value listed = adaptive::Value::array();
                for (const std::pair<const std::vector<int>, int>& entry : powers) {
                    adaptive::Value pair = adaptive::Value::array();
                    adaptive::Value ids = adaptive::Value::array();
                    for (int id : entry.first)
                        ids.push(adaptive::Value::integer(static_cast<long long>(id)));
                    pair.push(ids);
                    pair.push(adaptive::Value::integer(static_cast<long long>(entry.second)));
                    listed.push(pair);
                }
                got = adaptive::json_dumps(listed, false, false, true);
            } else if (column[0] == "approximate") {
                const adaptive::Value listed = adaptive::read_json(column[1]);
                std::vector<adaptive::Sites> bundles;
                for (std::size_t i = 0; i < listed.size(); ++i) {
                    adaptive::Sites bundle;
                    for (std::size_t k = 0; k < listed.at(i).size(); ++k)
                        bundle.push_back(static_cast<int>(listed.at(i).at(k).as_int()));
                    bundles.push_back(bundle);
                }
                const adaptive::Value listed_edges = adaptive::read_json(column[2]);
                std::vector<adaptive::Edge> edges;
                for (std::size_t i = 0; i < listed_edges.size(); ++i)
                    edges.emplace_back(static_cast<int>(listed_edges.at(i).at(0).as_int()),
                                       static_cast<int>(listed_edges.at(i).at(1).as_int()));
                const int order = static_cast<int>(adaptive::from_decimal(column[3]));
                const bool gcce = column[4] == "1";
                const adaptive::Value listed_raw = adaptive::read_json(column[5]);
                std::map<adaptive::Sites, adaptive::Curves> raw;
                for (const std::string& key : listed_raw.keys()) {
                    adaptive::Sites sites;
                    for (std::size_t at_key = 0; at_key < key.size(); ) {
                        const std::size_t comma = key.find(',', at_key);
                        sites.push_back(static_cast<int>(
                            adaptive::from_decimal(key.substr(at_key, comma == std::string::npos
                                                                    ? std::string::npos
                                                                    : comma-at_key))));
                        if (comma == std::string::npos) break;
                        at_key = comma+1;
                    }
                    raw[sites] = curves_of(adaptive::read_json(listed_raw.at(key).as_string()));
                }
                const adaptive::Curves found =
                    adaptive::approximate(bundles, edges, order, raw, gcce, adaptive::ExpansionLimits());
                adaptive::Value out = adaptive::Value::array();
                for (const adaptive::Curve& curve : found) {
                    adaptive::Value listed_curve = adaptive::Value::array();
                    for (const std::pair<double, double>& z : curve) {
                        adaptive::Value point = adaptive::Value::array();
                        point.push(adaptive::Value::text(bits_of(z.first)));
                        point.push(adaptive::Value::text(bits_of(z.second)));
                        listed_curve.push(point);
                    }
                    out.push(listed_curve);
                }
                got = adaptive::json_dumps(out, false, false, true);
            } else if (column[0] == "discrepancy") {
                got = shown(adaptive::discrepancy(curves_of(adaptive::read_json(column[1])),
                                                curves_of(adaptive::read_json(column[2])),
                                                from_bits(column[3]), from_bits(column[4])));
            } else {
                continue;
            }
        } catch (const adaptive::Error& error) {
            got = error.type();
            message = error.what();
        }
        const std::string expected = column.size() > at ? column[at] : std::string();
        const std::string expected_message = column.size() > at+1 ? column[at+1] : std::string();
        CHECK_EQ(got, expected);
        CHECK_EQ(message, expected_message);
        ++rows;
    }
    CHECK(rows >= (argc > 1 ? 1 : 307));
    std::cout << "compared " << rows << " expansion cases\n";
    return check::report("test_expansion");
}
