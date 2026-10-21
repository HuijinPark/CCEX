// The floating-point model the whole port depends on (plan 0B-1, Task 9).
// Two ways a build can quietly change arithmetic: a fused multiply-add, and a
// fast-math mode where isinf() is false. Both were real here - icpx defaults to
// -fp-model=fast, and -fp-model=precise still leaves FMA on (10th review).
#include "check.h"
#include <cmath>
#include <limits>

int main() {
    // a*a - 1 with a = 1 + 2^-27. Rounded multiply first: exactly 2^-26.
    // Fused: 2^-26 + 2^-54, a different double.
    volatile double a = 1.0+0x1p-27;
    const double contracted = a*a-1.0;
    CHECK_EQ(contracted, 0x1p-26);

    volatile double infinity = std::numeric_limits<double>::infinity();
    CHECK(std::isinf(infinity));
    CHECK(infinity == std::numeric_limits<double>::infinity());
    volatile double quiet = std::numeric_limits<double>::quiet_NaN();
    CHECK(std::isnan(quiet));
    CHECK(!(quiet == quiet));
    volatile double zero = 0.0;
    CHECK(std::signbit(-zero));
    CHECK(!std::signbit(zero));
    CHECK(std::numeric_limits<double>::is_iec559);
    return check::report("test_numeric");
}
