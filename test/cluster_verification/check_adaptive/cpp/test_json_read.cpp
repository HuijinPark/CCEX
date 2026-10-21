// Read every case and write what we read, in the tagged form Python used.
// The comparison is json_expect.py: whole trees, not a few keys (plan 0B-1, Task 7).
#include "check.h"
#include "base_json.h"
#include "mkdirs.h"
#include <clocale>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

static std::string slurp(const std::string& path) {
    std::ifstream file(path);
    std::ostringstream out;
    out << file.rdbuf();
    return out.str();
}

static std::string bits(double value) {
    std::uint64_t raw;
    std::memcpy(&raw, &value, sizeof raw);
    char text[24];
    std::snprintf(text, sizeof text, "%016llx", (unsigned long long)raw);
    return text;
}

// The same tagging the JSON cases use, so the trees compare directly.
static adaptive::Value tagged(const adaptive::Value& value) {
    using K = adaptive::Value::Kind;
    adaptive::Value out = adaptive::Value::object();
    switch (value.kind()) {
        case K::Null: return adaptive::Value();
        case K::Bool: return adaptive::Value::boolean(value.as_bool());
        case K::Int: out.set("int", adaptive::Value::text(adaptive::decimal(value.as_int()))); return out;
        case K::Double: out.set("float", adaptive::Value::text(bits(value.as_double()))); return out;
        case K::String: {
            // Bytes as well as text, so a wrong internal representation cannot hide
            // behind a writer that makes the opposite mistake.
            out.set("str", adaptive::Value::text(value.as_string()));
            std::string hex;
            for (unsigned char byte : value.as_string()) {
                char pair[3];
                std::snprintf(pair, sizeof pair, "%02x", byte);
                hex += pair;
            }
            out.set("bytes", adaptive::Value::text(hex));
            return out;
        }
        case K::Array: {
            adaptive::Value items = adaptive::Value::array();
            for (std::size_t i = 0; i < value.size(); ++i) items.push(tagged(value.at(i)));
            return items;
        }
        case K::Object: {
            adaptive::Value pairs = adaptive::Value::array();
            for (const std::string& key : value.keys()) {
                adaptive::Value pair = adaptive::Value::array();
                pair.push(tagged(adaptive::Value::text(key)));   // the key's bytes count too
                pair.push(tagged(value.at(key)));
                pairs.push(pair);
            }
            out.set("object", pairs);
            return out;
        }
    }
    return adaptive::Value();
}

// The same content without the nesting: a walk, not a tree. The tagged tree is
// deeper than what it describes (an object level becomes {"object": [[key, value]]}),
// so a document near the supported depth cannot be written in that form. Skipping
// those left real inputs unverified, which is not acceptable (11th review).
static adaptive::Value flattened(const adaptive::Value& root) {
    using K = adaptive::Value::Kind;
    adaptive::Value out = adaptive::Value::array();
    struct Frame { const adaptive::Value* value; std::size_t index; };
    std::vector<Frame> stack;
    const adaptive::Value* current = &root;
    auto mark = [&out](const char* which, const char* text) {
        adaptive::Value item = adaptive::Value::object();
        item.set(which, adaptive::Value::text(text));
        out.push(item);
    };
    for (;;) {
        if (current != nullptr) {
            if (current->kind() == K::Object) {
                mark("o", "{");
                stack.push_back({current, 0});
            } else if (current->kind() == K::Array) {
                mark("a", "[");
                stack.push_back({current, 0});
            } else {
                out.push(tagged(*current));      // a scalar: tagged() does not recurse
            }
            current = nullptr;
        }
        if (stack.empty())
            return out;
        Frame& top = stack.back();
        const bool object = top.value->kind() == K::Object;
        const std::size_t count = object ? top.value->keys().size() : top.value->size();
        if (top.index >= count) {
            stack.pop_back();
            mark(object ? "o" : "a", object ? "}" : "]");
            continue;
        }
        if (object) {
            const std::string& key = top.value->keys()[top.index];
            adaptive::Value item = adaptive::Value::object();
            item.set("key", tagged(adaptive::Value::text(key)));
            out.push(item);
            current = &top.value->at(key);
        } else {
            current = &top.value->at(top.index);
        }
        ++top.index;
    }
}

// What one file produced: the tagged tree, or the exception with its position.
static adaptive::Value outcome(const std::string& text) {
    adaptive::Value record = adaptive::Value::object();
    try {
        adaptive::Value value = adaptive::read_json(text);
        record.set("ok", adaptive::Value::boolean(true));
        record.set("value", tagged(value));
    } catch (const adaptive::JsonDecodeError& error) {
        adaptive::Value pair = adaptive::Value::array();
        pair.push(adaptive::Value::text(error.type()));
        pair.push(adaptive::Value::text(error.what()));
        record.set("ok", adaptive::Value::boolean(false));
        record.set("error", pair);
        record.set("pos", adaptive::Value::integer(static_cast<long long>((unsigned long long)error.pos())));
    } catch (const adaptive::Error& error) {
        // Every Python-shaped failure is data, not a crash: RecursionError is a
        // Error but not a ValueError, and catching only the latter aborted the
        // whole run on the first too-deep document (2026-09-16).
        adaptive::Value pair = adaptive::Value::array();
        pair.push(adaptive::Value::text(error.type()));
        pair.push(adaptive::Value::text(error.what()));
        record.set("ok", adaptive::Value::boolean(false));
        record.set("error", pair);
    }
    return record;
}

static const char* CASES[] = {
    "plain", "big", "huge", "notation", "zero", "digits_in_string", "escaped_key", "escapes",
    "surrogate_pair",
    "lone_surrogate", "consecutive_surrogates", "surrogate_key", "raw_emoji",
    "nonfinite_literals", "nonfinite_strings", "literals", "order",
    "duplicate",
    "nested", "deep", "top_scalar", "top_string", "digits_at_limit", "digits_over_limit",
    "trailing_comma", "bad_number", "bad_number_nonascii_key", "incomplete_exponent", "extra_data",
    "bad_escape", "raw_control", "unterminated", "unterminated_escape", "unterminated_pair"};

int main(int argc, char* argv[]) {
    if (argc > 2 && std::string(argv[1]) == "--write-back") {
        // Read, then write with our own writer, for Python to compare against the
        // original. The writer is otherwise only exercised on one built-by-hand
        // object; this puts every recorded file and every fuzz case through it.
        int written = 0;
        for (const std::string& path : adaptive::list_json(argv[2])) {
            adaptive::Value value;
            try {
                value = adaptive::read_json(slurp(std::string(argv[2])+"/"+path));
            } catch (const adaptive::ValueError&) {
                continue;                   // Python refuses it too; nothing to round trip
            }
            const std::string out = "out/back/"+path;
            testing::mkdirs(out.substr(0, out.rfind('/')));
            adaptive::save_json(out, value, 2);
            ++written;
        }
        std::cout << written << " files written back\n";
        return check::report("test_json_read");
    }
    if (argc > 2 && std::string(argv[1]) == "--all") {
        // Every recorded file, written under out/all with the same relative path.
        int read = 0, flat = 0;
        for (const std::string& path : adaptive::list_json(argv[2])) {      // helper in json_read.cpp
            const std::string out = "out/all/"+path;
            testing::mkdirs(out.substr(0, out.rfind('/')));   // save_json never creates directories
            const std::string text = slurp(std::string(argv[2])+"/"+path);
            adaptive::Value record = outcome(text);
            try {
                adaptive::save_json(out, record, 1);
            } catch (const adaptive::RecursionError&) {
                // Too deep for the tree form: describe the same content as a walk, so
                // the file is still compared rather than set aside.
                adaptive::Value walk = adaptive::Value::object();
                walk.set("ok", adaptive::Value::boolean(true));
                walk.set("flat", flattened(adaptive::read_json(text)));
                adaptive::save_json(out, walk, 1);
                ++flat;
            }
            ++read;
        }
        std::cout << read << " files read (" << flat << " as a walk)\n";
        return check::report("test_json_read");
    }
    testing::mkdirs("out/read");
    for (const char* name : CASES)
        adaptive::save_json(std::string("out/read/")+name+".json",
                          outcome(slurp(std::string("../reference/json-cases/")+name+".json")), 1);
    // A decoded key has to be usable as a key: an escape kept verbatim would make
    // has("seed") false while the tagged tree still looked right (9th review).
    adaptive::Value keys = adaptive::read_json(slurp("../reference/json-cases/escaped_key.json"));
    CHECK(keys.has("seed"));
    CHECK_EQ(adaptive::decimal(keys.at("seed").as_int()), std::string("19"));
    CHECK(keys.has("a"));
    // Nesting: allowed up to the declared limit, refused past it, never a crash.
    // Reading 200,000 levels used to be SIGSEGV (10th review).
    CHECK(adaptive::read_json(std::string(500, '[')+std::string(500, ']')).size() == 1);
    CHECK_THROWS(adaptive::read_json(std::string(501, '[')+std::string(501, ']')),
                 adaptive::RecursionError, "maximum JSON nesting depth (500) exceeded");
    CHECK_THROWS(adaptive::read_json(std::string(200000, '[')+std::string(200000, ']')),
                 adaptive::RecursionError, "maximum JSON nesting depth (500) exceeded");
    {
        std::string mixed;
        for (int i = 0; i < 200; ++i) mixed += "{\"a\": [";
        for (int i = 0; i < 200; ++i) mixed += "]}";
        CHECK(adaptive::read_json(mixed).has("a"));      // 400 levels, objects and arrays mixed
        std::string over;
        for (int i = 0; i < 300; ++i) over += "{\"a\": [";
        for (int i = 0; i < 300; ++i) over += "]}";
        CHECK_THROWS(adaptive::read_json(over), adaptive::RecursionError,
                     "maximum JSON nesting depth (500) exceeded");   // 600 levels
    }

    // Bytes are decoded before JSON is parsed, exactly as Python decodes a file.
    CHECK_THROWS(adaptive::read_json("\"\xff\""), adaptive::UnicodeDecodeError,
                 "'utf-8' codec can't decode byte 0xff in position 1: invalid start byte");
    CHECK_THROWS(adaptive::read_json("\"\xc0\xaf\""), adaptive::UnicodeDecodeError,
                 "'utf-8' codec can't decode byte 0xc0 in position 1: invalid start byte");
    CHECK_THROWS(adaptive::read_json("\"\xed\xa0\x80\""), adaptive::UnicodeDecodeError,
                 "'utf-8' codec can't decode byte 0xed in position 1: invalid continuation byte");
    // A cut-short sequence is reported as a range, and named "unexpected end of data".
    CHECK_THROWS(adaptive::read_json("\"\xe2\x82"), adaptive::UnicodeDecodeError,
                 "'utf-8' codec can't decode bytes in position 1-2: unexpected end of data");
    // The literal is split: "\x98A" would be read as one hex escape 0x98A, which
    // icpx rejects outright (g++ let it through - measured on idm).
    CHECK_THROWS(adaptive::read_json("\"\xf0\x9f\x98" "A\""), adaptive::UnicodeDecodeError,
                 "'utf-8' codec can't decode bytes in position 1-3: invalid continuation byte");
    CHECK_THROWS(adaptive::read_json("\"\xc2"), adaptive::UnicodeDecodeError,
                 "'utf-8' codec can't decode byte 0xc2 in position 1: unexpected end of data");

    // A locale must not change what a number means. Skipped where the locale is absent.
    if (std::setlocale(LC_NUMERIC, "de_DE.UTF-8") != nullptr) {
        CHECK_EQ(adaptive::read_json("1.25").as_double(), 1.25);
        CHECK_EQ(adaptive::read_json("-0.25").as_double(), -0.25);
        CHECK_EQ(adaptive::read_json("1.25e2").as_double(), 125.);
        std::setlocale(LC_NUMERIC, "C");
    }
    std::cout << sizeof CASES/sizeof *CASES << " cases read\n";
    return check::report("test_json_read");
}
