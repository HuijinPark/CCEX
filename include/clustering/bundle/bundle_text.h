/**
 * @file bundle_text.h
 * @brief Bundle clustering: small input-check and JSON text helpers.
 */
#ifndef CCEX_BUNDLE_TEXT_H
#define CCEX_BUNDLE_TEXT_H

// Small helpers the files of the calculation half of adaptive clustering each carried a copy of (review of the
// duplicates, 2026-10-08): the input check that throws, a JSON array of integers, a JSON string. Integer and text code
// only; this header is for clustering/bundle/ (the algorithm half has its own JSON writer).

#include "../../input/json.h"
#include <cstdlib>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace adaptive_bundle {

inline void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

// [a,b,c] with the stream's own formatting of each element.
template <class T> void jsonArray(std::ostream& out, const std::vector<T>& values) {
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) out << (i ? "," : "") << values[i];
    out << ']';
}

// A JSON string, escaped by cJSON.
inline void jsonString(std::ostream& out, const char* text) {
    std::unique_ptr<cJSON, decltype(&cJSON_Delete)> value(cJSON_CreateString(text), cJSON_Delete);
    char* encoded = value ? cJSON_PrintUnformatted(value.get()) : nullptr;
    if (!encoded) throw std::runtime_error("cannot encode a JSON string");
    out << encoded;
    std::free(encoded);
}

}  // namespace adaptive_bundle

#endif
