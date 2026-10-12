// JSON as the refinement reads and writes it: the value tree, the reader, the writer,
// the input-file reader, and the text of floats and numbers.
#include "base_json.h"

#include <algorithm>
#include <cerrno>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <limits>
#include <locale>
#include <locale.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

// ----------------------------------------------------------------------------
// Object semantics for adaptive::Value: insertion order, last value wins (plan 0B-1, Task 6).
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {
const char* kind_name_of(Value::Kind kind) {
    switch (kind) {
        case Value::Kind::Null: return "NoneType";
        case Value::Kind::Bool: return "bool";
        case Value::Kind::Int: return "int";
        case Value::Kind::Double: return "float";
        case Value::Kind::String: return "str";
        case Value::Kind::Array: return "list";
        case Value::Kind::Object: return "dict";
    }
    return "?";
}
}   // namespace

void Value::expect(Kind kind, const char* name) const {
    if (kind_ != kind)
        throw TypeError(std::string("expected ")+name+", found "+kind_name_of(kind_));
}

void Value::set(const std::string& key, Value value) {
    auto found = index_.find(key);
    if (found != index_.end()) {
        items_[found->second] = std::move(value);   // the key keeps its first position
        return;
    }
    index_.emplace(key, items_.size());
    keys_.push_back(key);
    items_.push_back(std::move(value));
}

// The container kind is checked before the lookup: asking an array for a key used
// to raise KeyError, where Python raises TypeError. That is not a wording
// difference - it decides which catch the caller lands in (11th review).
const Value* Value::find(const std::string& key) const {
    expect(Kind::Object, "dict");
    auto found = index_.find(key);
    return found == index_.end() ? nullptr : &items_[found->second];
}

const Value& Value::at(const std::string& key) const {
    const Value* found = find(key);
    if (found == nullptr)
        throw KeyError(key);
    return *found;
}

const Value& Value::at(std::size_t i) const {
    if (kind_ != Kind::Array)
        throw TypeError(std::string(kind_name_of(kind_))+" indices must be integers or slices");
    if (i >= items_.size())
        throw IndexError("list index out of range");
    return items_[i];
}

Value with(const Value& object, const std::string& key, const Value& value) {
    Value out = Value::object();
    for (const std::string& name : object.keys())
        out.set(name, object.at(name));
    out.set(key, value);
    return out;
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// floats as text, and JSON text
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {

// The fewest correctly rounded digits that read back as `value`, plus where the
// decimal point sits. Just below a power of two a digit string one shorter can
// exist that is not the correctly rounded one (the gap there is half the gap
// above); this writes the longer one, which reads back exactly all the same.
void shortest_digits(double value, std::string* digits, int* decpt) {
    char text[32];
    for (int precision = 0; precision <= 16; ++precision) {       // 17 digits always read back
        std::snprintf(text, sizeof text, "%.*e", precision, value);
        if (std::strtod(text, nullptr) == value)
            break;
    }
    const char* exponent = std::strchr(text, 'e');                // d[.ddd]e[+-]xx
    digits->assign(1, text[0]);
    if (exponent-text > 2)
        digits->append(text+2, static_cast<std::size_t>(exponent-text-2));
    *decpt = std::atoi(exponent+1)+1;
}
}   // namespace

std::string float_text(double value) {
    if (std::isnan(value))
        return "nan";
    if (std::isinf(value))
        return value > 0 ? "inf" : "-inf";
    const std::string sign = std::signbit(value) ? "-" : "";
    if (value == 0.)
        return sign+"0.0";
    std::string digits;
    int decpt = 0;
    shortest_digits(std::fabs(value), &digits, &decpt);
    // Python switches to the exponent form outside this window, and nowhere else.
    if (decpt <= -4 || decpt > 16) {
        std::string out = sign+digits.substr(0, 1);
        if (digits.size() > 1)
            out += "."+digits.substr(1);
        const int exponent = decpt-1;
        char tail[16];
        std::snprintf(tail, sizeof tail, "e%c%02d", exponent < 0 ? '-' : '+',
                      exponent < 0 ? -exponent : exponent);
        return out+tail;
    }
    if (decpt <= 0)
        return sign+"0."+std::string(static_cast<std::size_t>(-decpt), '0')+digits;
    if (static_cast<std::size_t>(decpt) >= digits.size())
        return sign+digits+std::string(static_cast<std::size_t>(decpt)-digits.size(), '0')+".0";
    return sign+digits.substr(0, static_cast<std::size_t>(decpt))+"."+
           digits.substr(static_cast<std::size_t>(decpt));
}

namespace {

void escape_ascii(const std::string& text, std::string* out) {
    static const char* kHex = "0123456789abcdef";
    out->push_back('"');
    for (std::size_t i = 0; i < text.size();) {
        const unsigned char byte = text[i];
        if (byte == '"' || byte == '\\') {
            out->push_back('\\');
            out->push_back(static_cast<char>(byte));
            ++i;
        } else if (byte == '\n') { *out += "\\n"; ++i; }
        else if (byte == '\r') { *out += "\\r"; ++i; }
        else if (byte == '\t') { *out += "\\t"; ++i; }
        else if (byte == '\b') { *out += "\\b"; ++i; }
        else if (byte == '\f') { *out += "\\f"; ++i; }
        // json.encoder escapes everything outside space..tilde, so DEL goes too:
        // its ESCAPE_ASCII pattern is [^ -~], not "control characters" (measured
        // against CPython over all 1,114,112 code points, 2026-09-16).
        else if (byte < 0x20 || byte == 0x7f) {
            *out += "\\u00";
            out->push_back(kHex[(byte >> 4) & 0xf]);
            out->push_back(kHex[byte & 0xf]);
            ++i;
        } else if (byte < 0x80) {
            out->push_back(static_cast<char>(byte));
            ++i;
        } else {
            // ensure_ascii: every non-ASCII character leaves as \uXXXX, and one
            // outside the basic plane leaves as the surrogate pair Python writes.
            std::size_t length = byte >= 0xf0 ? 4 : byte >= 0xe0 ? 3 : 2;
            length = std::min(length, text.size()-i);
            unsigned point = byte & (0xffu >> (length+1));
            for (std::size_t k = 1; k < length; ++k)
                point = (point << 6) | (static_cast<unsigned char>(text[i+k]) & 0x3f);
            auto write = [&out](unsigned unit) {   // kHex has static storage: naming it is ill-formed
                *out += "\\u";
                for (int shift = 12; shift >= 0; shift -= 4)
                    out->push_back(kHex[(unit >> shift) & 0xf]);
            };
            if (point >= 0x10000) {
                const unsigned rest = point-0x10000;
                write(0xd800+(rest >> 10));
                write(0xdc00+(rest & 0x3ff));
            } else {
                write(point);
            }
            i += length;
        }
    }
    out->push_back('"');
}

void dump(const Value& value, bool compact, bool sort_keys, bool allow_nan, std::string* out) {
    const char* comma = compact ? "," : ", ";
    const char* colon = compact ? ":" : ": ";
    switch (value.kind()) {
        case Value::Kind::Null: *out += "null"; return;
        case Value::Kind::Bool: *out += value.as_bool() ? "true" : "false"; return;
        case Value::Kind::Int: *out += decimal(value.as_int()); return;
        case Value::Kind::Double: {
            const double number = value.as_double();
            if (!std::isfinite(number)) {
                if (!allow_nan)
                    throw ValueError("Out of range float values are not JSON compliant");
                *out += std::isnan(number) ? "NaN" : number > 0 ? "Infinity" : "-Infinity";
                return;
            }
            *out += float_text(number);
            return;
        }
        case Value::Kind::String: escape_ascii(value.as_string(), out); return;
        case Value::Kind::Array:
            *out += "[";
            for (std::size_t i = 0; i < value.size(); ++i) {
                if (i) *out += comma;
                dump(value.at(i), compact, sort_keys, allow_nan, out);
            }
            *out += "]";
            return;
        case Value::Kind::Object: {
            std::vector<std::string> keys = value.keys();
            if (sort_keys)
                std::sort(keys.begin(), keys.end());   // by code point, as Python sorts
            *out += "{";
            for (std::size_t i = 0; i < keys.size(); ++i) {
                if (i) *out += comma;
                escape_ascii(keys[i], out);
                *out += colon;
                dump(value.at(keys[i]), compact, sort_keys, allow_nan, out);
            }
            *out += "}";
            return;
        }
    }
}
}   // namespace

// A key can hold a NUL or a lone surrogate, and both ruin a message: what() hands
// back a C string, so everything after a NUL is lost, and WTF-8 bytes make the log
// line invalid UTF-8. Control characters and lone surrogates are escaped the way
// Python shows them (measured on the fuzz corpus, 2026-09-16).
std::string printable(const std::string& text) {
    static const char* kHex = "0123456789abcdef";
    std::string out;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const unsigned char byte = text[i];
        if (byte == 0xed && i+2 < text.size() && (unsigned char)text[i+1] >= 0xa0 &&
            (unsigned char)text[i+1] <= 0xbf && ((unsigned char)text[i+2] & 0xc0) == 0x80) {
            const unsigned point = 0xd000 | (((unsigned char)text[i+1] & 0x3f) << 6) |
                                   ((unsigned char)text[i+2] & 0x3f);
            out += "\\u";
            for (int shift = 12; shift >= 0; shift -= 4)
                out.push_back(kHex[(point >> shift) & 0xf]);
            i += 2;
        } else if (byte == '\\') {
            out += "\\\\";              // or a real NUL and the text "\\x00" would look alike
        } else if (byte < 0x20 || byte == 0x7f) {
            out += "\\x";
            out.push_back(kHex[(byte >> 4) & 0xf]);
            out.push_back(kHex[byte & 0xf]);
        } else {
            out.push_back(text[i]);
        }
    }
    return out;
}

std::string json_dumps(const Value& value, bool compact, bool sort_keys, bool allow_nan) {
    std::string out;
    dump(value, compact, sort_keys, allow_nan, &out);
    return out;
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// reading numbers from text
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {

// Decode UTF-8 (or the WTF-8 the JSON reader may produce) into code points.
// Anything malformed becomes 0xFFFD, which no table matches, so it fails to parse.
std::vector<unsigned> code_points(const std::string& text) {
    std::vector<unsigned> points;
    for (std::size_t i = 0; i < text.size();) {
        const unsigned char lead = text[i];
        std::size_t length = lead < 0x80 ? 1 : lead >= 0xf0 ? 4 : lead >= 0xe0 ? 3 :
                             lead >= 0xc0 ? 2 : 0;
        if (length == 0 || i+length > text.size()) {
            points.push_back(0xfffd);
            ++i;
            continue;
        }
        unsigned point = length == 1 ? lead : lead & (0xffu >> (length+1));
        for (std::size_t k = 1; k < length; ++k)
            point = (point << 6) | (static_cast<unsigned char>(text[i+k]) & 0x3f);
        points.push_back(point);
        i += length;
    }
    return points;
}

// CPython transforms decimal digits and whitespace to ASCII and then parses
// ASCII, so int("１_２") is 12 there. This port reads ASCII only - a declared
// compatibility exception (설계 3.3절, 2026-09-19). Reproducing CPython's
// version meant carrying its Nd and whitespace tables, 847 generated lines
// pinned to one Unicode release, and the only input that could reach them is a
// command line someone typed non-ASCII digits into.
// Python strips whitespace from a numeric literal before parsing it. Only the
// ASCII whitespace is stripped here - reading a literal is ASCII-only, a
// declared exception (설계 3.3절). Nothing transforms the rest of the text: a
// byte past 0x7f is not a digit, so all_digits and strtod refuse it where a
// substitution character used to.
bool ascii_space(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}

std::string trimmed(const std::string& text) {
    std::size_t start = 0;
    while (start < text.size() && ascii_space(text[start]))
        ++start;
    std::size_t end = text.size();
    while (end > start && ascii_space(text[end-1]))
        --end;
    return text.substr(start, end-start);
}

// Underscores are allowed between digits only, never doubled, leading or trailing.
bool strip_underscores(const std::string& text, std::string* out) {
    for (std::size_t i = 0; i < text.size(); ++i) {
        if (text[i] != '_') {
            out->push_back(text[i]);
            continue;
        }
        const bool digit_before = i > 0 && std::isdigit(static_cast<unsigned char>(text[i-1]));
        const bool digit_after = i+1 < text.size() &&
                                 std::isdigit(static_cast<unsigned char>(text[i+1]));
        if (!digit_before || !digit_after)
            return false;
    }
    return true;
}

bool all_digits(const std::string& text) {
    if (text.empty())
        return false;
    for (char c : text)
        if (!std::isdigit(static_cast<unsigned char>(c)))
            return false;
    return true;
}

double to_double_c(const std::string& text, bool* ok) {
    // strtod() follows LC_NUMERIC, so the conversion is pinned to the C locale
    // (json_read.cpp has the measurement).
    static locale_t c_locale = newlocale(LC_ALL_MASK, "C", static_cast<locale_t>(0));
    char* end = nullptr;
    const double value = c_locale != static_cast<locale_t>(0)
        ? strtod_l(text.c_str(), &end, c_locale) : std::strtod(text.c_str(), &end);
    *ok = end == text.c_str()+text.size() && !text.empty();
    return value;
}

}   // namespace

std::string lowered(const std::string& text) {
    std::string out = text;
    for (char& c : out)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return out;
}

std::string quoted(const std::string& text) {
    const std::vector<unsigned> points = code_points(text);
    const bool has_single = text.find('\'') != std::string::npos;
    const bool has_double = text.find('"') != std::string::npos;
    const char quote = has_single && !has_double ? '"' : '\'';
    std::string out(1, quote);
    static const char* kHex = "0123456789abcdef";
    for (unsigned point : points) {
        if (point == static_cast<unsigned>(quote) || point == '\\') {
            out.push_back('\\');
            out.push_back(static_cast<char>(point));
        } else if (point == '\n') {
            out += "\\n";
        } else if (point == '\r') {
            out += "\\r";
        } else if (point == '\t') {
            out += "\\t";
        } else if (point >= 0x20 && point < 0x7f) {
            // ASCII printables pass; everything else is escaped. CPython keeps a
            // printable non-ASCII character as itself, which needs its printable
            // table - the same ASCII-only exception. This repr is only
            // ever seen inside an int()/float() error message.
            out.push_back(static_cast<char>(point));
        } else {
            const char* format = point < 0x100 ? "\\x" : point < 0x10000 ? "\\u" : "\\U";
            const int width = point < 0x100 ? 2 : point < 0x10000 ? 4 : 8;
            out += format;
            for (int shift = (width-1)*4; shift >= 0; shift -= 4)
                out.push_back(kHex[(point >> shift) & 0xf]);
        }
    }
    out.push_back(quote);
    return out;
}

long long parse_int(const std::string& text) {
    const std::string body = trimmed(text);
    std::string clean;
    bool valid = strip_underscores(body, &clean);
    std::string digits = clean;
    bool negative = false;
    if (valid && !digits.empty() && (digits[0] == '+' || digits[0] == '-')) {
        negative = digits[0] == '-';
        digits = digits.substr(1);
    }
    if (!valid || !all_digits(digits))
        throw ValueError("invalid literal for int() with base 10: "+quoted(text));
    if (digits.size() > 4300)
        throw ValueError("Exceeds the limit (4300 digits) for integer string conversion: "
                           "value has "+std::to_string(digits.size())+" digits; "
                           "use sys.set_int_max_str_digits() to increase the limit");
    // Built in 128 bits, like from_decimal: this port's integers are int64
    // (plan stage 7), so a literal CPython would hold and this cannot is one it
    // refuses. Accumulating in int64 wrapped instead, which is a wrong answer
    // and not a refusal (measured by the fuzz corpus, 2026-09-19).
    std::string shown = digits;
    while (shown.size() > 1 && shown[0] == '0')
        shown = shown.substr(1);
    if (negative && shown != "0")
        shown = "-"+shown;
    __int128 value = 0;
    const __int128 ceiling = static_cast<__int128>(1) << 63;   // |LLONG_MIN|
    for (char c : digits) {
        value = value*10+(c-'0');
        if (value > ceiling)
            throw ValueError("integer "+shown+" does not fit in 64 bits");
    }
    if (!negative && value == ceiling)
        throw ValueError("integer "+shown+" does not fit in 64 bits");
    return static_cast<long long>(negative ? -value : value);
}

double parse_float(const std::string& text) {
    const std::string body = trimmed(text);
    std::string clean;
    if (!strip_underscores(body, &clean))
        throw ValueError("could not convert string to float: "+quoted(text));
    std::string rest = clean;
    double sign = 1.;
    if (!rest.empty() && (rest[0] == '+' || rest[0] == '-')) {
        sign = rest[0] == '-' ? -1. : 1.;
        rest = rest.substr(1);
    }
    const std::string word = lowered(rest);
    if (word == "inf" || word == "infinity")
        return sign*HUGE_VAL;
    if (word == "nan")
        return std::copysign(std::nan(""), sign);   // Python keeps the sign: float("-nan") is negative
    // strtod would also take "0x1p3", "inf" and "nan"; only the plain decimal
    // grammar is Python's here, so the text is checked before it is converted.
    bool digit_seen = false, dot_seen = false, exponent_seen = false;
    for (std::size_t i = 0; i < rest.size(); ++i) {
        const char c = rest[i];
        if (std::isdigit(static_cast<unsigned char>(c))) {
            digit_seen = true;
        } else if (c == '.' && !dot_seen && !exponent_seen) {
            dot_seen = true;
        } else if ((c == 'e' || c == 'E') && digit_seen && !exponent_seen) {
            exponent_seen = true;
            digit_seen = false;
            if (i+1 < rest.size() && (rest[i+1] == '+' || rest[i+1] == '-'))
                ++i;
        } else {
            digit_seen = false;
            break;
        }
    }
    bool ok = false;
    const double value = digit_seen ? to_double_c(rest, &ok) : 0.;
    if (!digit_seen || !ok)
        throw ValueError("could not convert string to float: "+quoted(text));
    return sign*value;
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// Reading the JSON Python wrote, with Python's own failures (plan 0B-1, Task 7).
//
// Written rather than taken from a library because nlohmann 3.12.0 refuses
// 2^1024 before its SAX callback ever runs, so recorded integers could not
// survive (design 3.2). What that costs is spelled out here: the structure, the
// strings and the error positions are all ours to get right, and the tests
// compare them against what CPython's json module does, message included.
//
// Positions are CHARACTER indices, as Python reports them: a byte index would
// differ as soon as a key is not ASCII.
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {

const std::size_t kMaxIntStrDigits = 4300;   // Python's sys.get_int_max_str_digits() default

// Our own limit, and a compatibility exception: CPython's json stops at about
// 9,999 because of its recursion limit, which is a property of its stack rather
// than of the data. Recursing that far here overflowed the stack instead of
// raising (SIGSEGV at 200,000 while reading, at 10,000 while writing - measured
// 2026-09-16). The deepest recorded reference file nests 9 levels.
const int kMaxDepth = 500;

// The text must be valid UTF-8 before JSON is even considered: Python decodes the
// file's bytes first and raises UnicodeDecodeError, so 0xff, an overlong 0xc0 0xaf
// and a raw surrogate ED A0 80 never reach its parser (measured). Only the escape
// \udXXX may produce a lone surrogate, and that stays inside as WTF-8.
void check_utf8(const std::string& text) {
    // CPython names the failure three ways and reports a range when more than one
    // byte is involved: "byte 0xc2 in position 0: unexpected end of data" but
    // "bytes in position 0-1: invalid continuation byte". The second byte also
    // carries the overlong and surrogate limits (E0 needs A0-BF, ED 80-9F, F0
    // 90-BF, F4 80-8F), which is why a raw surrogate is reported at its lead byte.
    // Rules read off CPython 3.13.10 directly (2026-09-16).
    auto fail = [&text](std::size_t at, std::size_t count, const char* why) {
        std::string where;
        if (count == 1) {
            char byte[8];
            std::snprintf(byte, sizeof byte, "0x%02x", static_cast<unsigned char>(text[at]));
            where = std::string("byte ")+byte+" in position "+std::to_string(at);
        } else {
            where = "bytes in position "+std::to_string(at)+"-"+std::to_string(at+count-1);
        }
        throw UnicodeDecodeError("'utf-8' codec can't decode "+where+": "+why);
    };
    for (std::size_t i = 0; i < text.size();) {
        const unsigned char lead = text[i];
        if (lead < 0x80) { ++i; continue; }
        const std::size_t length = lead >= 0xf0 && lead <= 0xf4 ? 4 :
                                   lead >= 0xe0 && lead <= 0xef ? 3 :
                                   lead >= 0xc2 && lead <= 0xdf ? 2 : 0;
        if (length == 0)
            fail(i, 1, "invalid start byte");
        unsigned char low = 0x80, high = 0xbf;
        if (lead == 0xe0) low = 0xa0;
        else if (lead == 0xed) high = 0x9f;
        else if (lead == 0xf0) low = 0x90;
        else if (lead == 0xf4) high = 0x8f;
        for (std::size_t k = 1; k < length; ++k) {
            if (i+k >= text.size())
                fail(i, k, "unexpected end of data");
            const unsigned char next = text[i+k];
            if (next < low || next > high)
                fail(i, k, "invalid continuation byte");
            low = 0x80;
            high = 0xbf;
        }
        i += length;
    }
}

// strtod() follows LC_NUMERIC: under a German locale "1.25" stopped at the '.'
// and read as 1.0 (measured, 10th review). JSON always uses '.', so the
// conversion is pinned to the C locale and the whole token must be consumed.
double to_double_c(const std::string& text) {
    static locale_t c_locale = newlocale(LC_ALL_MASK, "C", static_cast<locale_t>(0));
    char* end = nullptr;
    const double value = c_locale != static_cast<locale_t>(0)
        ? strtod_l(text.c_str(), &end, c_locale) : std::strtod(text.c_str(), &end);
    if (end != text.c_str()+text.size())
        throw ValueError("could not convert string to float: '"+text+"'");
    return value;
}

class Reader {
 public:
    // on_integer is optional. The path it receives is built only when it is set:
    // reading 986 recorded files has no reason to assemble "a[0].b" for every token.
    Reader(const std::string& text, const IntegerCheck* on_integer)
        : text_(text), on_integer_(on_integer) {}

    Value parse() {
        skip_space();
        Value value = read_value();
        skip_space();
        if (at_ < text_.size())
            fail("Extra data", at_);
        return value;
    }

 private:
    const std::string& text_;
    const IntegerCheck* on_integer_ = nullptr;
    int depth_ = 0;
    std::vector<std::string> path_;         // only maintained when on_integer_ is set
    std::size_t at_ = 0;

    std::string where_in_document() const {
        std::string out;
        for (const std::string& part : path_)
            out += part;
        return out.empty() || out[0] != '.' ? out : out.substr(1);
    }

    bool done() const { return at_ >= text_.size(); }
    char peek() const { return done() ? '\0' : text_[at_]; }

    void skip_space() {
        while (!done() && (text_[at_] == ' ' || text_[at_] == '\t' ||
                           text_[at_] == '\n' || text_[at_] == '\r'))
            ++at_;
    }

    // line, column and character index, counting characters and not bytes
    [[noreturn]] void fail(const std::string& what, std::size_t position) const {
        std::size_t line = 1, column = 1, character = 0;
        for (std::size_t i = 0; i < position && i < text_.size(); ++i) {
            if ((text_[i] & 0xc0) == 0x80)
                continue;                   // a UTF-8 continuation byte is not a character
            ++character;
            if (text_[i] == '\n') {
                ++line;
                column = 1;
            } else {
                ++column;
            }
        }
        throw JsonDecodeError(what+": line "+std::to_string(line)+" column "+std::to_string(column)+
                              " (char "+std::to_string(character)+")", character);
    }

    bool literal(const char* word) {
        const std::size_t length = std::strlen(word);
        if (text_.compare(at_, length, word) != 0)
            return false;
        at_ += length;
        return true;
    }

    Value read_value() {
        if (done())
            fail("Expecting value", at_);
        switch (peek()) {
            case '{': return read_object();
            case '[': return read_array();
            case '"': return Value::text(read_string());
            default: break;
        }
        if (literal("true")) return Value::boolean(true);
        if (literal("false")) return Value::boolean(false);
        if (literal("null")) return Value();
        // Python's json accepts these three by default, so we do too.
        if (literal("NaN")) return Value::number(std::nan(""));
        if (literal("Infinity")) return Value::number(HUGE_VAL);
        if (literal("-Infinity")) return Value::number(-HUGE_VAL);
        if (peek() == '-' || (peek() >= '0' && peek() <= '9'))
            return read_number();
        fail("Expecting value", at_);
    }

    class Deeper {                          // one level, released however the parse ends
     public:
        explicit Deeper(Reader* reader) : reader_(reader) {
            if (++reader_->depth_ > kMaxDepth)
                throw RecursionError("maximum JSON nesting depth ("+std::to_string(kMaxDepth)+
                                       ") exceeded");
        }
        ~Deeper() { --reader_->depth_; }
     private:
        Reader* reader_;
    };

    Value read_object() {
        Deeper level(this);
        Value out = Value::object();
        ++at_;                              // '{'
        skip_space();
        if (peek() == '}') { ++at_; return out; }
        for (;;) {
            skip_space();
            if (peek() != '"')
                fail("Expecting property name enclosed in double quotes", at_);
            const std::string key = read_string();
            skip_space();
            if (peek() != ':')
                fail("Expecting ':' delimiter", at_);
            ++at_;
            skip_space();
            if (on_integer_ != nullptr)
                path_.push_back("."+key);
            out.set(key, read_value());     // a repeated key keeps its first position
            if (on_integer_ != nullptr)
                path_.pop_back();
            skip_space();
            if (peek() == '}') { ++at_; return out; }
            if (peek() != ',')
                fail("Expecting ',' delimiter", at_);
            const std::size_t comma = at_;
            ++at_;
            skip_space();
            if (peek() == '}')
                fail("Illegal trailing comma before end of object", comma);
        }
    }

    Value read_array() {
        Deeper level(this);
        Value out = Value::array();
        ++at_;                              // '['
        skip_space();
        if (peek() == ']') { ++at_; return out; }
        for (;;) {
            skip_space();
            if (on_integer_ != nullptr)
                path_.push_back("["+std::to_string(out.size())+"]");
            out.push(read_value());
            if (on_integer_ != nullptr)
                path_.pop_back();
            skip_space();
            if (peek() == ']') { ++at_; return out; }
            if (peek() != ',')
                fail("Expecting ',' delimiter", at_);
            const std::size_t comma = at_;
            ++at_;
            skip_space();
            if (peek() == ']')
                fail("Illegal trailing comma before end of array", comma);
        }
    }

    static void append_utf8(unsigned point, std::string* out) {
        if (point < 0x80) {
            out->push_back(static_cast<char>(point));
        } else if (point < 0x800) {
            out->push_back(static_cast<char>(0xc0 | (point >> 6)));
            out->push_back(static_cast<char>(0x80 | (point & 0x3f)));
        } else if (point < 0x10000) {
            // A lone surrogate lands here and becomes WTF-8, because Python keeps it.
            out->push_back(static_cast<char>(0xe0 | (point >> 12)));
            out->push_back(static_cast<char>(0x80 | ((point >> 6) & 0x3f)));
            out->push_back(static_cast<char>(0x80 | (point & 0x3f)));
        } else {
            out->push_back(static_cast<char>(0xf0 | (point >> 18)));
            out->push_back(static_cast<char>(0x80 | ((point >> 12) & 0x3f)));
            out->push_back(static_cast<char>(0x80 | ((point >> 6) & 0x3f)));
            out->push_back(static_cast<char>(0x80 | (point & 0x3f)));
        }
    }

    // `escape` is the backslash; CPython points at the 'u' after it (char 2 for
    // "\u00q0", measured), so every message from here is reported one further on.
    unsigned read_hex4(std::size_t escape) {
        // At least one character must follow the four digits - the closing quote, if
        // nothing else. CPython refuses "\u0041 with "Invalid \uXXXX escape" rather
        // than "Unterminated string", and a differential run over 20,000 random
        // documents found 49 inputs where we said the latter (2026-09-16).
        if (at_+4 >= text_.size())
            fail("Invalid \\uXXXX escape", escape+1);
        unsigned point = 0;
        for (int i = 0; i < 4; ++i) {
            const char c = text_[at_+i];
            const int digit = c >= '0' && c <= '9' ? c-'0' :
                              c >= 'a' && c <= 'f' ? c-'a'+10 :
                              c >= 'A' && c <= 'F' ? c-'A'+10 : -1;
            if (digit < 0)
                fail("Invalid \\uXXXX escape", escape+1);
            point = point*16+static_cast<unsigned>(digit);
        }
        at_ += 4;
        return point;
    }

    std::string read_string() {
        const std::size_t start = at_;
        ++at_;                              // '"'
        std::string out;
        for (;;) {
            if (done())
                fail("Unterminated string starting at", start);
            const unsigned char c = text_[at_];
            if (c == '"') { ++at_; return out; }
            if (c < 0x20)
                fail("Invalid control character at", at_);
            if (c != '\\') { out.push_back(text_[at_++]); continue; }
            const std::size_t escape = at_;
            ++at_;
            if (done())
                fail("Unterminated string starting at", start);
            switch (text_[at_]) {
                case '"': out.push_back('"'); ++at_; break;
                case '\\': out.push_back('\\'); ++at_; break;
                case '/': out.push_back('/'); ++at_; break;
                case 'b': out.push_back('\b'); ++at_; break;
                case 'f': out.push_back('\f'); ++at_; break;
                case 'n': out.push_back('\n'); ++at_; break;
                case 'r': out.push_back('\r'); ++at_; break;
                case 't': out.push_back('\t'); ++at_; break;
                case 'u': {
                    ++at_;
                    unsigned point = read_hex4(escape);
                    if (point >= 0xd800 && point <= 0xdbff && at_+1 < text_.size() &&
                        text_[at_] == '\\' && text_[at_+1] == 'u') {
                        // Look, but consume only for a real pair. Swallowing the second
                        // escape unconditionally cost "\ud800\ud801\udc00" its pair:
                        // CPython reads U+D800 then U+10400, we read three lone
                        // surrogates (measured 2026-09-16, 10th review).
                        const std::size_t mark = at_;
                        at_ += 2;
                        const unsigned low = read_hex4(mark);
                        if (low >= 0xdc00 && low <= 0xdfff)
                            point = 0x10000+((point-0xd800) << 10)+(low-0xdc00);
                        else
                            at_ = mark;     // put it back; the loop reads it next
                    }
                    append_utf8(point, &out);
                    break;
                }
                default:
                    fail("Invalid \\escape", escape);
            }
        }
    }

    Value read_number() {
        const std::size_t start = at_;
        if (peek() == '-') ++at_;
        // JSON grammar: no leading zeros. CPython stops the match here too and then
        // reports the next character as a missing delimiter, which is why "01" says
        // "Expecting ',' delimiter" rather than complaining about the number.
        if (peek() == '0') {
            ++at_;
        } else {
            while (!done() && peek() >= '0' && peek() <= '9') ++at_;
        }
        if (at_ == start || (at_ == start+1 && text_[start] == '-'))
            fail("Expecting value", start);
        bool real = false;
        if (peek() == '.') {
            std::size_t digits = at_+1;
            while (digits < text_.size() && text_[digits] >= '0' && text_[digits] <= '9') ++digits;
            if (digits > at_+1) { at_ = digits; real = true; }
        }
        if (peek() == 'e' || peek() == 'E') {
            std::size_t exponent = at_+1;
            if (exponent < text_.size() && (text_[exponent] == '+' || text_[exponent] == '-'))
                ++exponent;
            std::size_t digits = exponent;
            while (digits < text_.size() && text_[digits] >= '0' && text_[digits] <= '9') ++digits;
            if (digits > exponent) { at_ = digits; real = true; }
        }
        const std::string text = text_.substr(start, at_-start);
        if (!real) {
            const std::size_t count = text.size()-(text[0] == '-' ? 1 : 0);
            if (count > kMaxIntStrDigits)
                throw ValueError("Exceeds the limit ("+std::to_string(kMaxIntStrDigits)+
                                   " digits) for integer string conversion: value has "+
                                   std::to_string(count)+" digits; "
                                   "use sys.set_int_max_str_digits() to increase the limit");
            long long value;
            try {
                value = from_decimal(text);
            } catch (const ValueError&) {
                // Only the reader knows where in the document this token sits,
                // and the configuration error names it.
                throw ValueError("integer "+text+
                                   (on_integer_ != nullptr ? " at "+printable(where_in_document())
                                                           : std::string())+
                                   " does not fit in 64 bits");
            }
            if (on_integer_ != nullptr)
                (*on_integer_)(value, where_in_document());   // every token, not just the survivors
            return Value::integer(value);
        }
        return Value::number(to_double_c(text));
    }
};

void collect(const std::string& root, const std::string& prefix, std::vector<std::string>* out) {
    DIR* directory = opendir((prefix.empty() ? root : root+"/"+prefix).c_str());
    if (directory == nullptr)
        return;
    while (dirent* entry = readdir(directory)) {
        const std::string name = entry->d_name;
        if (name == "." || name == "..")
            continue;
        const std::string relative = prefix.empty() ? name : prefix+"/"+name;
        struct stat info;
        if (stat((root+"/"+relative).c_str(), &info) != 0)
            continue;
        if (S_ISDIR(info.st_mode))
            collect(root, relative, out);
        else if (relative.size() > 5 && relative.compare(relative.size()-5, 5, ".json") == 0)
            out->push_back(relative);
    }
    closedir(directory);
}
}   // namespace

Value read_json(const std::string& text) {
    check_utf8(text);
    return Reader(text, nullptr).parse();
}

Value read_json(const std::string& text, const IntegerCheck& on_integer) {
    check_utf8(text);
    return Reader(text, &on_integer).parse();
}

std::vector<std::string> list_json(const std::string& root) {
    std::vector<std::string> found;
    collect(root, "", &found);
    std::sort(found.begin(), found.end());
    return found;
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// Writing JSON the way adaptive_refine.py:35 writes it (plan 0B-1, Task 6).
//
// Two rules here exist because Python reads the result back and its types matter:
// a finite double always carries a '.' or an 'e' (otherwise 1.0 comes back as
// int 1 and -0.0 as int 0), and an integer is written exactly, never through a
// double. Non-finite doubles become the strings clean() writes.
// ----------------------------------------------------------------------------
namespace adaptive {
namespace {

// The writer stops at the reader's kMaxDepth too: writing a 10,000-deep value
// overflowed the stack (SIGSEGV, measured 2026-09-16), and a reader-only limit is
// not enough, because a Value can also be built in memory.

// Python's 4,300-digit limit on integer text is not reproduced here: these
// integers are int64, so nineteen digits is the most any of them can print.
std::string integer_text(long long value) {
    std::string text = decimal(value);
    return text;
}

std::string double_text(double value) {
    // json.dumps writes floats with repr(), so this is repr() - which makes our
    // output byte-identical to Python's rather than merely equal in value.
    return float_text(value);
}

void escape(const std::string& text, std::string* out) {
    static const char* kHex = "0123456789abcdef";
    out->push_back('"');
    for (std::size_t i = 0; i < text.size(); ++i) {
        unsigned char byte = text[i];
        // A lone surrogate is held as WTF-8 (ED A0 80 .. ED BF BF) and must go out
        // as an escape: those bytes are not valid UTF-8 and Python could not read them.
        if (byte == 0xed && i+2 < text.size() && (unsigned char)text[i+1] >= 0xa0 &&
            (unsigned char)text[i+1] <= 0xbf && ((unsigned char)text[i+2] & 0xc0) == 0x80) {
            unsigned point = 0xd000 | (((unsigned char)text[i+1] & 0x3f) << 6) |
                             ((unsigned char)text[i+2] & 0x3f);
            out->append("\\u");
            for (int shift = 12; shift >= 0; shift -= 4)
                out->push_back(kHex[(point >> shift) & 0xf]);
            i += 2;
            continue;
        }
        switch (byte) {
            case '"': out->append("\\\""); break;
            case '\\': out->append("\\\\"); break;
            case '\b': out->append("\\b"); break;
            case '\f': out->append("\\f"); break;
            case '\n': out->append("\\n"); break;
            case '\r': out->append("\\r"); break;
            case '\t': out->append("\\t"); break;
            default:
                if (byte < 0x20) {
                    out->append("\\u00");
                    out->push_back(kHex[(byte >> 4) & 0xf]);
                    out->push_back(kHex[byte & 0xf]);
                } else {
                    out->push_back(text[i]);   // valid UTF-8 goes through unescaped
                }
        }
    }
    out->push_back('"');
}

void write(const Value& value, int indent, int depth, std::string* out) {
    if (depth > kMaxDepth)
        throw RecursionError("maximum JSON nesting depth ("+std::to_string(kMaxDepth)+
                               ") exceeded");
    const std::string pad(indent > 0 ? indent*(depth+1) : 0, ' ');
    const std::string close(indent > 0 ? indent*depth : 0, ' ');
    const char* newline = indent > 0 ? "\n" : "";
    switch (value.kind()) {
        case Value::Kind::Null: out->append("null"); return;
        case Value::Kind::Bool: out->append(value.as_bool() ? "true" : "false"); return;
        case Value::Kind::Int: out->append(integer_text(value.as_int())); return;
        case Value::Kind::Double:
            if (std::isnan(value.as_double())) { out->append("\"nan\""); return; }
            if (std::isinf(value.as_double())) {
                out->append(value.as_double() > 0 ? "\"inf\"" : "\"-inf\"");
                return;
            }
            out->append(double_text(value.as_double()));
            return;
        case Value::Kind::String: escape(value.as_string(), out); return;
        case Value::Kind::Array:
            if (value.size() == 0) { out->append("[]"); return; }
            out->append("[").append(newline);
            for (std::size_t i = 0; i < value.size(); ++i) {
                out->append(pad);
                write(value.at(i), indent, depth+1, out);
                out->append(i+1 < value.size() ? "," : "").append(newline);
            }
            out->append(close).append("]");
            return;
        case Value::Kind::Object:
            if (value.keys().empty()) { out->append("{}"); return; }
            out->append("{").append(newline);
            for (std::size_t i = 0; i < value.keys().size(); ++i) {
                out->append(pad);
                escape(value.keys()[i], out);
                out->append(indent > 0 ? ": " : ":");
                write(value.at(value.keys()[i]), indent, depth+1, out);
                out->append(i+1 < value.keys().size() ? "," : "").append(newline);
            }
            out->append(close).append("}");
            return;
    }
}
}   // namespace

std::string write_json(const Value& value, int indent) {
    std::string out;
    write(value, indent, 0, &out);
    out.push_back('\n');
    return out;
}

void save_json(const std::string& path, const Value& value, int indent) {
    // No directory is created here: Python's save() raises FileNotFoundError when
    // the parent is missing, and the port reproduces that (measured 2026-09-16).
    const std::string text = write_json(value, indent);
    const std::string temporary = path+".tmp";
    // errno decides the exception's name, as it does in Python: a path through a
    // plain file is NotADirectoryError and a directory in the way IsADirectoryError,
    // not one blanket FileNotFoundError (10th review).
    errno = 0;
    std::FILE* file = std::fopen(temporary.c_str(), "wb");
    if (file == nullptr)
        throw os_error(errno, temporary);
    errno = 0;
    const bool written = std::fwrite(text.data(), 1, text.size(), file) == text.size();
    const int write_errno = errno;
    const bool closed = std::fclose(file) == 0;
    if (!written)
        throw os_error(write_errno != 0 ? write_errno : EIO, temporary);
    if (!closed)
        throw os_error(errno != 0 ? errno : EIO, temporary);
    errno = 0;
    if (std::rename(temporary.c_str(), path.c_str()) != 0)
        throw os_error(errno, path);
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// Reading a user's configuration (plan 0B-1, Task 8).
// ----------------------------------------------------------------------------
namespace adaptive {
Value read_config(const std::string& text) {
    // Every integer TOKEN is checked, not the values that survive: {"a": <too big>, "a": 1}
    // has a rejected integer in it even though the final value is small. The
    // reader raises with the path; this hook is what turns path tracking on.
    return read_json(text, [](long long, const std::string&) {});
}

std::string without_comments(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    std::size_t at = 0;
    while (at < text.size()) {
        std::size_t end = text.find('\n', at);
        if (end == std::string::npos) end = text.size();
        const std::size_t cut = text.find_first_of("!#", at);
        out.append(text, at, (cut != std::string::npos && cut < end ? cut : end) - at);
        if (end < text.size()) out.push_back('\n');
        at = end + 1;
    }
    return out;
}

std::string read_file(const std::string& path) {
    // fopen() succeeds on a directory; only the read fails, and without ferror()
    // that came back as an empty string - a silent success (10th review).
    errno = 0;
    std::FILE* file = std::fopen(path.c_str(), "rb");
    if (file == nullptr)
        throw os_error(errno, path);
    std::string text;
    char buffer[65536];
    for (;;) {
        errno = 0;
        const std::size_t read = std::fread(buffer, 1, sizeof buffer, file);
        if (read > 0)
            text.append(buffer, read);
        if (read < sizeof buffer) {
            const int failure = std::ferror(file) ? (errno != 0 ? errno : EIO) : 0;
            std::fclose(file);
            if (failure != 0)
                throw os_error(failure, path);
            return text;
        }
    }
}

double number_of(const Value& value) {
    return value.kind() == Value::Kind::Int ? static_cast<double>(value.as_int()) : value.as_double();
}
}   // namespace adaptive
