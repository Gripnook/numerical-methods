#include "catch.hpp"

#include <cmath>

#include "float-test.h"
#include "integral.h"

using namespace Numeric;

TEST_CASE("the integral of cos(x) is correctly computed")
{
    double result =
        integral([](double x) { return std::cos(x); }, 0.0, M_PI / 2, 10000);
    REQUIRE(fequal(result, 1.0));
}

TEST_CASE("the integral of log(x) is correctly computed")
{
    double result =
        integral([](double x) { return std::log(x); }, 1.0, 2.0, 10000);
    REQUIRE(fequal(result, 2.0 * std::log(2) - 1.0));
}
