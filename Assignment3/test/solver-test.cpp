#include "catch.hpp"

#include <cmath>

#include "float-test.h"
#include "solver.h"

using namespace Numeric;

TEST_CASE("Newton-Raphson works for a linear function")
{
    auto f = [](double x) { return x + 1; };
    auto fprime = [](double x) { return 1; };

    double result;
    int iterations;
    std::tie(result, iterations) = newtonRaphson(f, fprime, 0);
    REQUIRE(fequal(result, -1));
}

TEST_CASE("Newton-Raphson works for a cubic function")
{
    auto f = [](double x) { return x * x * x - x; };
    auto fprime = [](double x) { return 3 * x * x - 1; };

    double result;
    int iterations;
    std::tie(result, iterations) = newtonRaphson(f, fprime, 1.2);
    REQUIRE(fequal(result, 1));
}

TEST_CASE("Newton-Raphson works for a cosine function")
{
    auto f = [](double x) { return std::cos(x); };
    auto fprime = [](double x) { return -std::sin(x); };

    double result;
    int iterations;
    std::tie(result, iterations) = newtonRaphson(f, fprime, 1);
    REQUIRE(fequal(result, M_PI / 2));
}

TEST_CASE("Successive substitution works for a linear function")
{
    auto f = [](double x) { return -1; };

    double result;
    int iterations;
    std::tie(result, iterations) = successiveSubstitution(f, 0);
    REQUIRE(fequal(result, -1));
}

TEST_CASE("Successive substitution works for a cubic function")
{
    auto f = [](double x) { return std::cbrt(x); };

    double result;
    int iterations;
    std::tie(result, iterations) = successiveSubstitution(f, 0.96);
    REQUIRE(fequal(result, 1));
}

TEST_CASE("Successive substitution works for a cosine function")
{
    auto f = [](double x) { return x + std::cos(x); };

    double result;
    int iterations;
    std::tie(result, iterations) = successiveSubstitution(f, 1);
    REQUIRE(fequal(result, M_PI / 2));
}
