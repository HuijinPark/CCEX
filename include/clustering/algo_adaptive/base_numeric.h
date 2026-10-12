// Numbers the policies can rely on: error types, checked 64-bit integers, one-rounding
// division, exact summation, and complex abs/division that fail instead of overflowing.
#pragma once
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// ----------------------------------------------------------------------------
// Python's exception identity, as far as plan 0B-1 needs it. 0B-2 extends this.
// The type and the message are both part of the decision (design 3.4).
// ----------------------------------------------------------------------------
namespace adaptive {
class Error : public std::runtime_error {
 public:
    Error(std::string type, const std::string& message)
        : std::runtime_error(message), type_(std::move(type)) {}
    const std::string& type() const { return type_; }
 private:
    std::string type_;
};

class ValueError : public Error {
 public:
    explicit ValueError(const std::string& message) : Error("ValueError", message) {}
 protected:
    ValueError(std::string type, const std::string& message) : Error(std::move(type), message) {}
};

// Python's own tree, because which parent catches a failure changes the decision
// (design 3.4): ArithmeticError covers OverflowError and ZeroDivisionError, and
// RecursionError is a RuntimeError, not a ValueError.
class ArithmeticError : public Error {
 public:
    explicit ArithmeticError(const std::string& message) : Error("ArithmeticError", message) {}
 protected:
    ArithmeticError(std::string type, const std::string& message)
        : Error(std::move(type), message) {}
};

struct OverflowError : ArithmeticError {
    explicit OverflowError(const std::string& message)
        : ArithmeticError("OverflowError", message) {}
};

struct ZeroDivisionError : ArithmeticError {
    explicit ZeroDivisionError(const std::string& message)
        : ArithmeticError("ZeroDivisionError", message) {}
};

class RuntimeError : public Error {
 public:
    explicit RuntimeError(const std::string& message) : Error("RuntimeError", message) {}
 protected:
    RuntimeError(std::string type, const std::string& message) : Error(std::move(type), message) {}
};

// The policy exception. ExpansionLimit is a ValueError subclass, so `except
// ValueError` catches it and the search carries on. (DimensionLimit, ScanLimit and
// ChangeUnverified went with the policies that raised them on 2026-10-10.)
struct ExpansionLimit : ValueError {
    explicit ExpansionLimit(const std::string& message) : ValueError("ExpansionLimit", message) {}
};

// str(KeyError(k)) is repr(k), so a string key carries quotes and an integer one
// does not. Both kinds appear: a cache key is text, a site index is a number.
struct KeyError : Error {
    explicit KeyError(const std::string& key) : Error("KeyError", "'"+key+"'") {}
    static KeyError integer(const std::string& decimal) {
        return KeyError(decimal, 0);
    }
 private:
    KeyError(const std::string& decimal, int) : Error("KeyError", decimal) {}
};

struct IndexError : Error {
    explicit IndexError(const std::string& message) : Error("IndexError", message) {}
};

struct TypeError : Error {
    explicit TypeError(const std::string& message) : Error("TypeError", message) {}
};

// Our own limit, not Python's: see kMaxDepth in json_read.cpp (compatibility exception).
// RecursionError is a RuntimeError in Python, and stays one here.
struct RecursionError : RuntimeError {
    explicit RecursionError(const std::string& message)
        : RuntimeError("RecursionError", message) {}
};

// Raised where Python decodes a file's bytes into text, before json sees them.
struct UnicodeDecodeError : ValueError {
    explicit UnicodeDecodeError(const std::string& message)
        : ValueError("UnicodeDecodeError", message) {}
};

// json.JSONDecodeError is a ValueError subclass that keeps its own name and position.
class JsonDecodeError : public ValueError {
 public:
    JsonDecodeError(const std::string& message, std::size_t pos)
        : ValueError("JSONDecodeError", message), pos_(pos) {}
    std::size_t pos() const { return pos_; }
 private:
    std::size_t pos_;
};

// errno, named the way Python names it: the type carries as much meaning as the text
// (NotADirectoryError and IsADirectoryError are separate classes there, not OSError).
inline Error os_error(int code, const std::string& path) {
    const char* name = code == ENOENT ? "FileNotFoundError" :
                       code == ENOTDIR ? "NotADirectoryError" :
                       code == EISDIR ? "IsADirectoryError" :
                       code == EACCES ? "PermissionError" : "OSError";
    return Error(name, "[Errno "+std::to_string(code)+"] "+std::strerror(code)+": '"+path+"'");
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// Python's integer conversions, for int64 (design 3.1, plan stage 7).
//
// Python's integers are unbounded and this started as boost's cpp_int to match
// them. Measuring ended that: replaying all 52 recordings, the widest integer
// any of them needs is five bits, and across 321 recorded artifacts the widest
// is 27 bits. The three quantities that could outgrow int64 at production scale
// are each computed only to be compared against a budget and then discarded, so
// they are compared instead of computed (Driver::affordable).
//
// What is left is long long. This header keeps only the conversions that
// have to match Python exactly.
// ----------------------------------------------------------------------------
namespace adaptive {
inline std::string decimal(long long value) {
    return std::to_string(value);
}

// A JSON integer, and nothing else: optional '-', then "0" or a digit string
// with no leading zero. An integer past int64 is refused rather than wrapped -
// no configuration this reads has one, and a silent wrap would be a wrong plan.
inline long long from_decimal(const std::string& text) {
    const std::size_t start = !text.empty() && text[0] == '-' ? 1 : 0;
    bool digits = text.size() > start;
    for (std::size_t i = start; digits && i < text.size(); ++i)
        digits = text[i] >= '0' && text[i] <= '9';
    if (digits && text.size() > start+1 && text[start] == '0')
        digits = false;                   // no leading zeros, as JSON requires
    if (!digits)
        throw ValueError("not a JSON integer: '"+text+"'");
    // Built in 128 bits so that -9223372036854775808 - which has no positive
    // counterpart - parses like any other value.
    __int128 value = 0;
    const __int128 ceiling = static_cast<__int128>(1) << 63;   // |LLONG_MIN|
    for (std::size_t i = start; i < text.size(); ++i) {
        value = value*10+(text[i]-'0');                        // base ten, explicitly
        if (value > ceiling)
            throw ValueError("integer "+text+" does not fit in 64 bits");
    }
    if (start == 0 && value == ceiling)
        throw ValueError("integer "+text+" does not fit in 64 bits");
    return static_cast<long long>(start ? -value : value);
}
// Python's max(a, b): b wins only if it is greater, so a NaN on the right is ignored and one on the left is kept
// (not symmetric). Inline in this header, which only the algorithm half includes (its own floating-point flags).
inline double py_max(double a, double b) { return b > a ? b : a; }

}   // namespace adaptive

// ----------------------------------------------------------------------------
// Integer arithmetic with Python's meaning (plan 0B-2, design 3.1).
//
// a/b is rounded ONCE even when both sides are huge. Doing it as a/double(b)
// rounds twice and lands a unit in the last place away, which changes decisions
// (measured: 1/(2097160**3) differs from 1/float(2097160**3)).
// ----------------------------------------------------------------------------
namespace adaptive {

double true_div(long long a, long long b);    // one rounding; OverflowError when too large
long long comb(long long n, long long k);

// Exact three-way comparison against a double: no conversion, so (10**17+1) < 1e17
// stays false. Returns -1, 0, 1; a NaN has no order and reports `unordered`.
int compare(long long a, double b, bool* unordered);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// abs() and division of a complex number, with a failure where C would hand back
// an infinity or a NaN (plan 0B-2, design 3.3 and 3.4): a magnitude that overflows
// and a division by zero are errors, not values to carry on with.
// ----------------------------------------------------------------------------
namespace adaptive {
double complex_abs(double re, double im);          // OverflowError("absolute value too large")
void complex_div(double ar, double ai, double br, double bi, double* re, double* im);
}   // namespace adaptive

// ----------------------------------------------------------------------------
// Exact summation (Shewchuk), as Python's math.fsum() does it: the answer is the
// correctly rounded sum whatever order the items come in, so a score never
// depends on the order its terms were visited in.
// ----------------------------------------------------------------------------
namespace adaptive {
double fsum(const std::vector<double>& items);
}   // namespace adaptive
