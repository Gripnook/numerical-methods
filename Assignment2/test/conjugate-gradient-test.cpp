#include "catch.hpp"

#include <algorithm>

#include "conjugate-gradient.h"
#include "matrix.h"
#include "matrix-util.h"

using namespace Numeric;

// Truncates the least significant bits of the value for a better floating point
// comparison that ignores some round-off error.
float fround(double value)
{
    return static_cast<float>(value);
}

// Checks if the two floating point numbers are equal
// to within some round-off error.
bool fequal(double lhs, double rhs)
{
    return fround(lhs) == fround(rhs);
}

TEST_CASE("cgsolve succeeds for an identity matrix")
{
    Matrix<double> m = eye<double>(3);
    Matrix<double> b = "[1;2;3]";
    auto x = cgsolve(m, b).first;
    REQUIRE(std::equal(x.begin(), x.end(), b.begin(), fequal));
}

TEST_CASE("cgsolve succeeds for a 2x2 system of equations")
{
    Matrix<double> lower = "["
                           "1,0;"
                           "2,1"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[4;3]";
    Matrix<double> b = m * x;
    auto result = cgsolve(m, b).first;
    REQUIRE(std::equal(result.begin(), result.end(), x.begin(), fequal));
}

TEST_CASE("cgsolve succeeds for a 3x3 system of equations")
{
    Matrix<double> lower = "["
                           "4, 0,0;"
                           "5, 1,0;"
                           "9,-1,2"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[1.0;2.5;-2.0]";
    Matrix<double> b = m * x;
    auto result = cgsolve(m, b).first;
    REQUIRE(std::equal(result.begin(), result.end(), x.begin(), fequal));
}

TEST_CASE("cgsolve succeeds for a 4x4 system of equations")
{
    Matrix<double> lower = "["
                           " 1, 0, 0, 0;"
                           " 2, 3, 0, 0;"
                           " 5, 7,11, 0;"
                           "13,17,19,23"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[1;2;4;8]";
    Matrix<double> b = m * x;
    auto result = cgsolve(m, b).first;
    REQUIRE(std::equal(result.begin(), result.end(), x.begin(), fequal));
}

TEST_CASE("cgsolve succeeds for a 5x5 system of equations")
{
    Matrix<double> lower = "["
                           " 1,0,0,0,0;"
                           " 2,1,0,0,0;"
                           " 4,2,1,0,0;"
                           " 8,4,2,1,0;"
                           "16,8,4,2,1"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[13;-7;19;-11;3]";
    Matrix<double> b = m * x;
    auto result = cgsolve(m, b).first;
    REQUIRE(std::equal(result.begin(), result.end(), x.begin(), fequal));
}

TEST_CASE("cgsolve throws exceptions for invalid inputs")
{
    Matrix<double> m1 = "[1,2;3,4;5,6]";
    Matrix<double> m2 = "[1,2;3,4]";
    Matrix<double> m3 = "[0,0,0;0,0,0;0,0,0]";
    Matrix<double> b = "[1;2;3]";
    CHECK_THROWS_WITH(cgsolve(m1, b), "cgsolve: matrix must be square");
    CHECK_THROWS_WITH(cgsolve(m2, b), "cgsolve: b must be an [nx1] vector");
    CHECK_THROWS_WITH(
        cgsolve(m3, b), "cgsolve: matrix must be positive-definite");
}
