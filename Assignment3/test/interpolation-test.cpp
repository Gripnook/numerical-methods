#include "catch.hpp"

#include "interpolation.h"

using namespace Numeric;

TEST_CASE("Lagrange interpolation is correct for a linear function")
{
    auto f = lagrange({{0.0, 1.0}, {1.0, 0.0}});
    CHECK(f(0.0) == 1.0);
    CHECK(f(0.5) == 0.5);
    CHECK(f(1.0) == 0.0);
}

TEST_CASE("Lagrange interpolation is correct for a quadratic function")
{
    auto f = lagrange({{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}});
    CHECK(f(0.0) == 1.0);
    CHECK(f(0.5) == 0.25);
    CHECK(f(1.0) == 0.0);
    CHECK(f(1.5) == 0.25);
    CHECK(f(2.0) == 1.0);
}

TEST_CASE("Hermite interpolation is correct for a linear function")
{
    auto f = hermite({{0.0, 1.0}, {1.0, 0.0}}, {-1.0, -1.0});
    CHECK(f(0.0) == 1.0);
    CHECK(f(0.5) == 0.5);
    CHECK(f(1.0) == 0.0);
}

TEST_CASE("Hermite interpolation is correct for a quadratic function")
{
    auto f = hermite({{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}, {-2.0, 0.0, 2.0});
    CHECK(f(0.0) == 1.0);
    CHECK(f(0.5) == 0.25);
    CHECK(f(1.0) == 0.0);
    CHECK(f(1.5) == 0.25);
    CHECK(f(2.0) == 1.0);
}

TEST_CASE("piecewise linear interpolation is correct")
{
    auto f = pwl({{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}, {3.0, 3.0}});
    CHECK(f(-0.5) == 1.5);
    CHECK(f(0.0) == 1.0);
    CHECK(f(0.5) == 0.5);
    CHECK(f(1.0) == 0.0);
    CHECK(f(1.5) == 0.5);
    CHECK(f(2.0) == 1.0);
    CHECK(f(2.5) == 2.0);
    CHECK(f(3.0) == 3.0);
    CHECK(f(3.5) == 4.0);
}

TEST_CASE("piecewise linear interpolation derivatives are correct")
{
    auto f = pwlprime({{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}, {3.0, 3.0}});
    CHECK(f(-0.5) == -1.0);
    CHECK(f(0.0) == -1.0);
    CHECK(f(0.5) == -1.0);
    CHECK(f(1.0) == 1.0);
    CHECK(f(1.5) == 1.0);
    CHECK(f(2.0) == 2.0);
    CHECK(f(2.5) == 2.0);
    CHECK(f(3.0) == 2.0);
    CHECK(f(3.5) == 2.0);
}

TEST_CASE("inverse piecewise linear interpolation is correct")
{
    auto f = pwlinverse({{0.0, 0.0}, {2.0, 1.0}, {3.0, 1.5}});
    CHECK(f(-0.5) == -1.0);
    CHECK(f(0.0) == 0.0);
    CHECK(f(0.5) == 1.0);
    CHECK(f(1.0) == 2.0);
    CHECK(f(1.5) == 3.0);
    CHECK(f(2.0) == 4.0);
}

TEST_CASE("slope computations are correct")
{
    CHECK(slope({0.0, 0.0}, {1.0, 2.0}) == 2.0);
    CHECK(slope({0.0, 1.0}, {1.0, 0.0}) == -1.0);
    CHECK(slope({0.0, 1.0}, {1.0, 0.5}) == -0.5);
}
