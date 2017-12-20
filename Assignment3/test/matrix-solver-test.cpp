#include "catch.hpp"

#include "matrix-solver.h"
#include "matrix.h"
#include "matrix-util.h"

using namespace Numeric;

TEST_CASE("lsolve succeeds for an identity matrix")
{
    Matrix<double> m = eye<double>(3);
    Matrix<double> b = "[1;2;3]";
    auto x = lsolve(m, b);
    REQUIRE(x == b);
}

TEST_CASE("lsolve succeeds for a 3x3 system of equations")
{
    Matrix<double> m = "[4,0,0;5,1,0;9,-1,2]";
    Matrix<double> x = "[1.0;2.5;-2.0]";
    Matrix<double> b = m * x;
    auto result = lsolve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("lsolve throws exceptions for invalid inputs")
{
    Matrix<double> m1 = "[1,2;3,4;5,6]";
    Matrix<double> m2 = "[1,2;3,4]";
    Matrix<double> m3 = "[0,0,0;0,0,0;0,0,0]";
    Matrix<double> b = "[1;2;3]";
    CHECK_THROWS_WITH(lsolve(m1, b), "lsolve: matrix must be square");
    CHECK_THROWS_WITH(lsolve(m2, b), "lsolve: b must be an [nx1] vector");
    CHECK_THROWS_WITH(lsolve(m3, b), "lsolve: matrix must be non-singular");
}

TEST_CASE("usolve succeeds for an identity matrix")
{
    Matrix<double> m = eye<double>(3);
    Matrix<double> b = "[1;2;3]";
    auto x = usolve(m, b);
    REQUIRE(x == b);
}

TEST_CASE("usolve succeeds for a 3x3 system of equations")
{
    Matrix<double> m = "[4,5,1;0,9,-1;0,0,2]";
    Matrix<double> x = "[1.0;2.5;-2.0]";
    Matrix<double> b = m * x;
    auto result = usolve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("usolve throws exceptions for invalid inputs")
{
    Matrix<double> m1 = "[1,2;3,4;5,6]";
    Matrix<double> m2 = "[1,2;3,4]";
    Matrix<double> m3 = "[0,0,0;0,0,0;0,0,0]";
    Matrix<double> b = "[1;2;3]";
    CHECK_THROWS_WITH(usolve(m1, b), "usolve: matrix must be square");
    CHECK_THROWS_WITH(usolve(m2, b), "usolve: b must be an [nx1] vector");
    CHECK_THROWS_WITH(usolve(m3, b), "usolve: matrix must be non-singular");
}

TEST_CASE("solve succeeds for an identity matrix")
{
    Matrix<double> m = eye<double>(3);
    Matrix<double> b = "[1;2;3]";
    auto x = solve(m, b);
    REQUIRE(x == b);
}

TEST_CASE("solve succeeds for a 2x2 system of equations")
{
    Matrix<double> lower = "["
                           "1,0;"
                           "2,1"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[4;3]";
    Matrix<double> b = m * x;
    auto result = solve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("solve succeeds for a 3x3 system of equations")
{
    Matrix<double> lower = "["
                           "4, 0,0;"
                           "5, 1,0;"
                           "9,-1,2"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[1.0;2.5;-2.0]";
    Matrix<double> b = m * x;
    auto result = solve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("solve succeeds for a 4x4 system of equations")
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
    auto result = solve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("solve succeeds for a 5x5 system of equations")
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
    auto result = solve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("solve throws exceptions for invalid inputs")
{
    Matrix<double> m1 = "[1,2;3,4;5,6]";
    Matrix<double> m2 = "[1,2;3,4]";
    Matrix<double> m3 = "[0,0,0;0,0,0;0,0,0]";
    Matrix<double> b = "[1;2;3]";
    CHECK_THROWS_WITH(solve(m1, b), "solve: matrix must be square");
    CHECK_THROWS_WITH(solve(m2, b), "solve: b must be an [nx1] vector");
    CHECK_THROWS_WITH(solve(m3, b), "solve: matrix must be positive-definite");
}

TEST_CASE("bsolve succeeds for an identity matrix")
{
    Matrix<double> m = eye<double>(3);
    Matrix<double> b = "[1;2;3]";
    auto x = bsolve(m, b);
    REQUIRE(x == b);
}

TEST_CASE("bsolve succeeds for a 2x2 system of equations")
{
    Matrix<double> lower = "["
                           "1,0;"
                           "2,1"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[4;3]";
    Matrix<double> b = m * x;
    auto result = bsolve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("bsolve succeeds for a 3x3 system of equations")
{
    Matrix<double> lower = "["
                           "4, 0,0;"
                           "0, 1,0;"
                           "0,-1,2"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[1.0;2.5;-2.0]";
    Matrix<double> b = m * x;
    auto result = bsolve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("bsolve succeeds for a 4x4 system of equations")
{
    Matrix<double> lower = "["
                           "1, 0, 0, 0;"
                           "0, 3, 0, 0;"
                           "0, 0,11, 0;"
                           "0,17, 0,23"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[1;2;4;8]";
    Matrix<double> b = m * x;
    auto result = bsolve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("bsolve succeeds for a 5x5 system of equations")
{
    Matrix<double> lower = "["
                           "1,0,0,0,0;"
                           "2,1,0,0,0;"
                           "4,0,1,0,0;"
                           "0,0,2,1,0;"
                           "0,0,0,0,1"
                           "]";
    Matrix<double> m = lower * transpose(lower);
    Matrix<double> x = "[13;-7;19;-11;3]";
    Matrix<double> b = m * x;
    auto result = bsolve(m, b);
    REQUIRE(result == x);
}

TEST_CASE("bsolve throws exceptions for invalid inputs")
{
    Matrix<double> m1 = "[1,2;3,4;5,6]";
    Matrix<double> m2 = "[1,2;3,4]";
    Matrix<double> m3 = "[0,0,0;0,0,0;0,0,0]";
    Matrix<double> b = "[1;2;3]";
    CHECK_THROWS_WITH(bsolve(m1, b), "bsolve: matrix must be square");
    CHECK_THROWS_WITH(bsolve(m2, b), "bsolve: b must be an [nx1] vector");
    CHECK_THROWS_WITH(
        bsolve(m3, b),
        "bsolve: matrix must be "
        "positive-definite");
}
