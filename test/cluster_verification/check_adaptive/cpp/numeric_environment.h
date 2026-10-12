// Whether this machine is the one the expectations were recorded on, and how to
// compare numbers when it is not (design 3.3, plan 0B-2).
//
// libm is not the same everywhere: idm's icpx answers log10 one step away from the
// glibc the expectations came from, and hypot differed between glibc 2.28 and 2.39
// before that. The contract is that a comparison is strict inside the recorded
// numeric environment and value-to-value within one step outside it.
#pragma once
#include "base_json.h"
#include "base_numeric.h"
#include <cmath>
#include <fstream>
#include <gnu/libc-version.h>
#include <sstream>
#include <string>

namespace testing {
inline bool same_numeric_environment() {
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

// One step apart, counted in representable doubles rather than in a made-up
// epsilon: 0.4617142810125527 and 0.46171428101255274 are neighbours.
inline bool within_one_step(double a, double b) {
    if (a == b)
        return true;
    if (std::isnan(a) || std::isnan(b) || std::isinf(a) || std::isinf(b))
        return false;
    return std::nextafter(a, b) == b;
}

// Two JSON values, compared field by field. `strict` demands equal doubles;
// otherwise a double may be one step away. Returns the path of the first
// difference, or "" when they agree.
inline std::string differs(const adaptive::Value& ours, const adaptive::Value& theirs, bool strict,
                           const std::string& path = "") {
    if (ours.kind() != theirs.kind())
        return path+": kind";
    switch (ours.kind()) {
        case adaptive::Value::Kind::Double:
            if (strict ? !(ours.as_double() == theirs.as_double() ||
                           (ours.as_double() != ours.as_double() &&
                            theirs.as_double() != theirs.as_double()))
                       : !within_one_step(ours.as_double(), theirs.as_double()))
                return path+": "+std::to_string(ours.as_double())+" vs "+
                       std::to_string(theirs.as_double());
            return "";
        case adaptive::Value::Kind::Array: {
            if (ours.size() != theirs.size())
                return path+": length";
            for (std::size_t i = 0; i < ours.size(); ++i) {
                const std::string found =
                    differs(ours.at(i), theirs.at(i), strict, path+"["+std::to_string(i)+"]");
                if (!found.empty())
                    return found;
            }
            return "";
        }
        case adaptive::Value::Kind::Object: {
            if (ours.keys() != theirs.keys())
                return path+": keys";
            for (const std::string& key : ours.keys()) {
                const std::string found =
                    differs(ours.at(key), theirs.at(key), strict, path+"."+key);
                if (!found.empty())
                    return found;
            }
            return "";
        }
        case adaptive::Value::Kind::Null:
            return "";
        case adaptive::Value::Kind::Bool:
            return ours.as_bool() == theirs.as_bool() ? "" : path+": bool";
        case adaptive::Value::Kind::Int:
            return ours.as_int() == theirs.as_int() ? "" : path+": int";
        default:
            return ours.as_string() == theirs.as_string() ? "" : path+": text";
    }
}
}   // namespace testing
