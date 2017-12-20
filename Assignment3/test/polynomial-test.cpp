#include "catch.hpp"

#include "polynomial.h"

using namespace Numeric;

TEST_CASE("polynomials can be evaluated")
{
    Polynomial p{{-1.0, -1.0, 1.0}};
    CHECK(p(0.0) == -1.0);
    CHECK(p(1.0) == -1.0);
    CHECK(p(2.0) == 1.0);
    CHECK(p(-1.0) == 1.0);
}

TEST_CASE("polynomials can be added")
{
    Polynomial p1{{-2.0, 1.0}};
    Polynomial p2{{2.0, 1.0}};
    Polynomial result{{0.0, 2.0}};
    CHECK(p1 + p2 == result);
    p1 += p2;
    CHECK(p1 == result);
}

TEST_CASE("polynomials can be subtracted")
{
    Polynomial p1{{-2.0, 1.0, 1.0}};
    Polynomial p2{{2.0, 1.0}};
    Polynomial result{{-4.0, 0.0, 1.0}};
    CHECK(p1 - p2 == result);
    p1 -= p2;
    CHECK(p1 == result);
}

TEST_CASE("polynomials can be multiplied")
{
    Polynomial p1{{-2.0, 1.0}};
    Polynomial p2{{2.0, 1.0}};
    Polynomial result{{-4.0, 0.0, 1.0}};
    CHECK(p1 * p2 == result);
    p1 *= p2;
    CHECK(p1 == result);
}

TEST_CASE("polynomials can be divided by scalars")
{
    Polynomial p1{{-2.0, 1.0}};
    double divisor = 2.0;
    Polynomial result{{-1.0, 0.5}};
    CHECK(p1 / divisor == result);
    p1 /= divisor;
    CHECK(p1 == result);
}
