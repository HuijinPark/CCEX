// Which parent catches a failure is part of the decision (plan 0B-2, Task 1).
#include "check.h"
#include "base_numeric.h"
#include <string>

// Caught as `Parent`? The policy code relies on exactly this: an `except ValueError`
// around a candidate keeps the search going, an `except RuntimeError` does not.
template <class Parent, class Thrown> bool caught_as(const Thrown& error) {
    try {
        throw error;
    } catch (const Parent&) {
        return true;
    } catch (...) {
        return false;
    }
}

int main() {
    CHECK(caught_as<adaptive::ValueError>(adaptive::ExpansionLimit("work budget")));
    CHECK(!caught_as<adaptive::RuntimeError>(adaptive::ExpansionLimit("work budget")));
    CHECK(caught_as<adaptive::ArithmeticError>(adaptive::OverflowError("too large")));
    CHECK(caught_as<adaptive::ArithmeticError>(adaptive::ZeroDivisionError("division by zero")));
    CHECK(caught_as<adaptive::RuntimeError>(adaptive::RecursionError("too deep")));
    CHECK(!caught_as<adaptive::ValueError>(adaptive::OverflowError("too large")));
    CHECK_EQ(adaptive::ExpansionLimit("x").type(), std::string("ExpansionLimit"));
    return check::report("test_errors");
}
