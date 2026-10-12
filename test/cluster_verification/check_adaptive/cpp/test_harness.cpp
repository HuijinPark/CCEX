// The harness must fail when it should, and must be able to print what it compared.
#include "check.h"
#include <stdexcept>

int main() {
    CHECK(1+1 == 2);
    CHECK_EQ(std::string("a")+"b", std::string("ab"));
    CHECK_EQ((std::vector<std::string>{"z", "a"}), (std::vector<std::string>{"z", "a"}));
    CHECK_THROWS(throw std::runtime_error("boom"), std::runtime_error, "boom");
    return check::report("test_harness");
}
