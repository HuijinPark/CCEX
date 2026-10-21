/**
 * @file base_numeric.cpp
 * @brief Adaptive refinement: checked integers, exact summation and safe complex arithmetic.
 */
// Numbers the policies can rely on: error types, checked 64-bit integers, one-rounding
// division, exact summation, and complex abs/division that fail instead of overflowing.
#include "base_numeric.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// ----------------------------------------------------------------------------
// 64-bit integers: one-rounding division, binomials, exact comparison with a double
// ----------------------------------------------------------------------------
namespace adaptive {
// CPython's long_true_divide, for operands that fit in int64: long division to
// 54 bits with a sticky bit, then round to nearest even.
//
// Both operands are below 2^63, so the exact quotient lies between 2^-63 and
// 2^63 - always a normal double. The overflow and subnormal branches CPython
// needs for unbounded integers cannot be reached here and are not written.
double true_div(long long a, long long b) {
    if (b == 0)
        throw ZeroDivisionError("division by zero");
    const bool negative = (a < 0) != (b < 0);
    if (a == 0)
        return negative ? -0.0 : 0.0;
    // |LLONG_MIN| does not fit in int64, so the magnitudes are taken in 128 bits.
    __int128 x = a < 0 ? -static_cast<__int128>(a) : a;
    const __int128 y = b < 0 ? -static_cast<__int128>(b) : b;
    __int128 quotient = x/y, remainder = x%y;
    long long bits = 0;
    for (__int128 rest = quotient; rest != 0; rest >>= 1)
        ++bits;
    long long exponent = 0;
    bool sticky = false;
    const long long want = 54;             // 53 for the mantissa, one to round on
    if (bits > want) {
        for (long long i = 0; i < bits-want; ++i) {
            sticky = sticky || (quotient & 1) != 0;
            quotient >>= 1;
        }
        exponent = bits-want;
        sticky = sticky || remainder != 0;
    } else {
        if (bits == 0) {                   // x < y: shift until the leading bit appears
            while (quotient == 0) {
                remainder <<= 1;
                --exponent;
                if (remainder >= y) { quotient = 1; remainder -= y; }
            }
            bits = 1;
        }
        for (; bits < want; ++bits) {      // one more quotient bit per step
            quotient <<= 1;
            remainder <<= 1;
            if (remainder >= y) { ++quotient; remainder -= y; }
            --exponent;
        }
        sticky = remainder != 0;
    }
    const bool guard = (quotient & 1) != 0;
    quotient >>= 1;
    ++exponent;
    if (guard && (sticky || (quotient & 1) != 0))
        ++quotient;                        // nearest, ties to even
    const double value = std::ldexp(static_cast<double>(static_cast<long long>(quotient)),
                                    static_cast<int>(exponent));
    return negative ? -value : value;
}

long long comb(long long n, long long k) {
    if (k < 0)
        throw ValueError("k must be a non-negative integer");
    if (n < 0)
        throw ValueError("n must be a non-negative integer");
    if (k > n)
        return 0;
    k = std::min(k, n-k);
    // The running value is C(n, i), so the step multiplies by (n-i) before
    // dividing by (i+1) - and that product is up to (i+1) times the answer.
    // C(63, 31) fits in int64 while C(63, 31)*31 does not, so the step is taken
    // in 128 bits and only the result has to fit.
    __int128 exact = 1;
    for (long long i = 0; i < k; ++i) {
        if (exact > (static_cast<__int128>(1) << 126)/static_cast<__int128>(n-i))
            throw OverflowError("integer overflow");
        exact = exact*static_cast<__int128>(n-i)/static_cast<__int128>(i+1);
    }
    if (exact > static_cast<__int128>(std::numeric_limits<long long>::max()))
        throw OverflowError("integer overflow");
    return static_cast<long long>(static_cast<long long>(exact));
}

int compare(long long a, double b, bool* unordered) {
    *unordered = false;
    if (std::isnan(b)) {
        *unordered = true;
        return 0;
    }
    if (std::isinf(b))
        return b > 0 ? -1 : 1;
    // |a| < 2^63, so a double at or past that settles the comparison by sign
    // alone. Lining the two up instead - shifting a by the exponent - is what
    // the boost version did, and it is how a merge cutoff of 1.797e308 came to
    // ask for a 1024-bit integer (found 2026-09-18 by running a real bath).
    if (b >= 9223372036854775808.)
        return -1;
    if (b < -9223372036854775808.)
        return 1;
    // b's integer part now fits exactly, and so does the fraction it leaves.
    const double whole = std::trunc(b);
    const long long right = static_cast<long long>(whole);
    if (a != right)
        return a < right ? -1 : 1;
    const double fraction = b-whole;
    return fraction > 0. ? -1 : fraction < 0. ? 1 : 0;
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// complex abs and division that fail instead of overflowing
// ----------------------------------------------------------------------------
namespace adaptive {
double complex_abs(double re, double im) {
    const double value = std::hypot(re, im);
    if (std::isinf(value) && std::isfinite(re) && std::isfinite(im))
        throw OverflowError("absolute value too large");
    return value;
}

// CPython's _Py_c_quot: divide by the larger part first, so the intermediate
// values stay in range. The zero denominator is a ZeroDivisionError, not an inf.
void complex_div(double ar, double ai, double br, double bi, double* re, double* im) {
    const double abs_br = std::fabs(br), abs_bi = std::fabs(bi);
    if (abs_br >= abs_bi) {
        if (abs_br == 0.) {
            if (std::isnan(bi)) {
                *re = *im = std::nan("");
                return;
            }
            throw ZeroDivisionError("complex division by zero");
        }
        const double ratio = bi/br;
        const double denominator = br+bi*ratio;
        *re = (ar+ai*ratio)/denominator;
        *im = (ai-ar*ratio)/denominator;
    } else if (abs_bi >= abs_br) {
        const double ratio = br/bi;
        const double denominator = br*ratio+bi;
        *re = (ar*ratio+ai)/denominator;
        *im = (ai*ratio-ar)/denominator;
    } else {                                 // at least one part is a NaN
        *re = *im = std::nan("");
    }
}
}   // namespace adaptive

// ----------------------------------------------------------------------------
// exact summation
// ----------------------------------------------------------------------------
namespace adaptive {
// Shewchuk's exact summation, as CPython's math_fsum: keep the non-overlapping
// partial sums, then round once at the end.
double fsum(const std::vector<double>& items) {
    std::vector<double> partials;
    double special = 0., infinities = 0.;
    for (double value : items) {
        const double original = value;
        std::size_t kept = 0;
        for (std::size_t i = 0; i < partials.size(); ++i) {
            double x = value, y = partials[i];
            if (std::fabs(x) < std::fabs(y))
                std::swap(x, y);
            const double high = x+y;
            const double low = y-(high-x);
            if (low != 0.)
                partials[kept++] = low;
            value = high;
        }
        partials.resize(kept);
        if (value == 0.)
            continue;
        if (!std::isfinite(value)) {
            if (std::isfinite(original))
                throw OverflowError("intermediate overflow in fsum");
            if (std::isinf(original))
                infinities += original;
            special += original;
            partials.clear();
            continue;
        }
        partials.push_back(value);
    }
    if (special != 0. || std::isnan(special)) {
        if (std::isnan(infinities))
            throw ValueError("-inf + inf in fsum");
        return special;
    }
    double high = 0.;
    if (!partials.empty()) {
        std::size_t n = partials.size();
        high = partials[--n];
        double low = 0.;
        while (n > 0) {
            const double x = high, y = partials[--n];
            high = x+y;
            low = y-(high-x);
            if (low != 0.)
                break;
        }
        if (n > 0 && ((low < 0. && partials[n-1] < 0.) || (low > 0. && partials[n-1] > 0.))) {
            const double y = low*2.;
            const double x = high+y;
            if (y == x-high)
                high = x;
        }
    }
    return high;
}
}   // namespace adaptive
