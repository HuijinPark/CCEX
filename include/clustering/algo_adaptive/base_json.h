// JSON as the refinement reads and writes it: the value tree, the reader, the writer,
// the input-file reader, and the text of floats and numbers.
#pragma once
#include <functional>
#include <map>
#include <string>
#include <vector>
#include "base_numeric.h"

// ----------------------------------------------------------------------------
// A JSON value that keeps what Python kept: exact integers, insertion order, int vs float.
// ----------------------------------------------------------------------------
namespace adaptive {
class Value {
 public:
    enum class Kind { Null, Bool, Int, Double, String, Array, Object };

    Value() : kind_(Kind::Null) {}
    static Value boolean(bool v) { Value x; x.kind_ = Kind::Bool; x.bool_ = v; return x; }
    static Value integer(long long v) { Value x; x.kind_ = Kind::Int; x.int_ = std::move(v); return x; }
    static Value number(double v) { Value x; x.kind_ = Kind::Double; x.double_ = v; return x; }
    static Value text(std::string v) { Value x; x.kind_ = Kind::String; x.text_ = std::move(v); return x; }
    static Value array() { Value x; x.kind_ = Kind::Array; return x; }
    static Value object() { Value x; x.kind_ = Kind::Object; return x; }

    Kind kind() const { return kind_; }
    // Reading a value as the wrong kind is a defect, not a zero: as_int() on a
    // missing key used to return 0 and as_double() on a string 0.0, which is how
    // a port silently decides something different (10th review).
    bool as_bool() const { expect(Kind::Bool, "bool"); return bool_; }
    long long as_int() const { expect(Kind::Int, "int"); return int_; }
    double as_double() const { expect(Kind::Double, "float"); return double_; }
    const std::string& as_string() const { expect(Kind::String, "str"); return text_; }

    void push(Value v) { items_.push_back(std::move(v)); }
    // Duplicate keys: the last value wins and keeps the first one's position, as Python's dict does.
    void set(const std::string& key, Value v);
    std::size_t size() const { return items_.size(); }
    // Access follows Python: a missing key is a KeyError and an index past the end
    // an IndexError. Use find() for a key that is allowed to be absent.
    const Value& at(std::size_t i) const;
    const Value& at(const std::string& key) const;
    const Value* find(const std::string& key) const;
    bool has(const std::string& key) const { return find(key) != nullptr; }
    const std::vector<std::string>& keys() const { return keys_; }

 private:
    void expect(Kind kind, const char* name) const;

    Kind kind_;
    bool bool_ = false;
    double double_ = 0.;
    long long int_;
    std::string text_;
    std::vector<Value> items_;            // array elements, or object values in insertion order
    std::vector<std::string> keys_;       // object keys in insertion order
    std::map<std::string, std::size_t> index_;
};

// The same object with one key's value replaced, or added at the end if it is
// new - what Python's dict assignment does, and what a report needs when a
// nested block has to be written back.
Value with(const Value& object, const std::string& key, const Value& value);

// adaptive_refine.save: indent 2, non-finite doubles as "inf"/"-inf"/"nan", trailing newline.
std::string write_json(const Value& value, int indent);
void save_json(const std::string& path, const Value& value, int indent);

// Data Python wrote: exact integers kept exact, with Python's own limits (Task 7).
// It is not "anything goes": the text must be valid UTF-8, nesting has a declared
// limit, and an integer past 4,300 digits fails exactly where Python fails.
Value read_json(const std::string& text);

// The same reader, told to look at every integer token as it is made (Task 8).
// The path is "a[0].b"; it is assembled only for this overload.
using IntegerCheck = std::function<void(long long value, const std::string& path)>;
Value read_json(const std::string& text, const IntegerCheck& on_integer);

// Every *.json under root, as sorted relative paths (the --all comparison walks these).
std::vector<std::string> list_json(const std::string& root);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// Floats and JSON as text, laid out the way Python's repr() and json.dumps() do
// (plan 0B-4).
//
// The probe cache hashes its own text: the context key is the SHA-256 of the
// compact, key-sorted JSON of the context, and every stored payload carries the
// SHA-256 of its own JSON. What that needs is one text per value, every time -
// and a float's text has to read back as the same float.
// ----------------------------------------------------------------------------
namespace adaptive {

// The fewest correctly rounded digits that read back exactly, laid out as repr(float).
std::string float_text(double value);

// json.dumps(). compact uses (',', ':'), otherwise (', ', ': ') as Python defaults do.
// Non-ASCII is escaped (ensure_ascii), and a non-finite float raises when allow_nan
// is false, as the cache requires.
// Text that is about to appear in a message or a log line, made safe to print.
std::string printable(const std::string& text);

std::string json_dumps(const Value& value, bool compact, bool sort_keys, bool allow_nan);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// int(str) and float(str), with Python's rules (plan 0B-2, design 3.3).
//
// Not the same as the JSON integer grammar in py_integer.h: these accept surrounding
// whitespace, underscores between digits, signs, and the several spellings of
// infinity and NaN. They are what reads the TSV responses, the cache payloads and
// the [profile] lines, so the failure behaviour on damaged input lives here.
// ----------------------------------------------------------------------------
namespace adaptive {
long long parse_int(const std::string& text);      // ValueError with Python's message
double parse_float(const std::string& text);

// Python's repr() of a string, which is what its error messages quote.
std::string quoted(const std::string& text);

// A number that may have been read as an integer or as a float, as a double.
double number_of(const Value& value);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// User-supplied JSON: the one input check that differs from Python (design 3.2).
// Python accepts an integer of any size here; we refuse past int64 rather than
// carry an exact big integer into a field that will index or count something.
// ----------------------------------------------------------------------------
namespace adaptive {
Value read_config(const std::string& text);
std::string read_file(const std::string& path);   // the whole file, or an OSError
std::string lowered(const std::string& text);     // ASCII lower case
// A CCEX input file as its JSON: every line is cut at its first '!' or '#', the same
// rule as cJSON_ReadFccein (src/input/option.cpp), so an input with comments reads
// the same here as in the calculation it asks for.
std::string without_comments(const std::string& text);
}
