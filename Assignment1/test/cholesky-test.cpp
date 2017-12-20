#include "catch.hpp"

#include "cholesky.h"
#include "matrix.h"
#include "matrix-util.h"

using namespace Numeric;

TEST_CASE("cholesky decomposition succeeds for an identity matrix")
{
    auto m = eye<double>(2);
    auto p = cholesky(m);
    CHECK(p.second == true);
    CHECK(p.first == m);
}

TEST_CASE("cholesky decomposition succeeds for a positive-definite matrix")
{
    Matrix<double> lower = "[4,0,0;-3,2,0;1,-1,5]";
    auto m = lower * transpose(lower);
    auto p = cholesky(m);
    CHECK(p.second == true);
    CHECK(p.first == lower);
}

TEST_CASE("cholesky decomposition fails for a non-singular, "
          "non-positive-definite matrix")
{
    Matrix<double> lower = "[4,0,0;-3,2,0;1,-1,5]";
    Matrix<double> upper = "[1,3,2;0,1,4;0,0,1]";
    auto m = lower * upper;
    auto p = cholesky(m);
    REQUIRE(p.second == false);
}

TEST_CASE("cholesky decomposition fails for a singular matrix")
{
    Matrix<double> m = "[0,1;0,0]";
    auto p = cholesky(m);
    REQUIRE(p.second == false);
}

TEST_CASE("cholesky decomposition throws for non-square matrices")
{
    Matrix<double> m = "[0,1]";
    CHECK_THROWS_WITH(cholesky(m), "cholesky: matrix must be square");
    m = "[0;1]";
    CHECK_THROWS_WITH(cholesky(m), "cholesky: matrix must be square");
}

TEST_CASE("bcholesky decomposition succeeds for an identity matrix")
{
    auto m = eye<double>(2);
    auto p = bcholesky(m);
    CHECK(p.second == true);
    CHECK(p.first == m);
}

TEST_CASE("bcholesky decomposition succeeds for a positive-definite matrix")
{
    Matrix<double> lower = "["
                           " 4, 0,0,0;"
                           "-3, 2,0,0;"
                           " 0,-1,5,0;"
                           " 0, 0,1,1"
                           "]";
    auto m = lower * transpose(lower);
    auto p = bcholesky(m);
    CHECK(p.second == true);
    CHECK(p.first == lower);
}

TEST_CASE("bcholesky decomposition fails for a non-singular, "
          "non-positive-definite matrix")
{
    Matrix<double> lower = "[4,0,0;-3,2,0;1,-1,5]";
    Matrix<double> upper = "[1,3,2;0,1,4;0,0,1]";
    auto m = lower * upper;
    auto p = bcholesky(m);
    REQUIRE(p.second == false);
}

TEST_CASE("bcholesky decomposition fails for a singular matrix")
{
    Matrix<double> m = "[0,1;0,0]";
    auto p = bcholesky(m);
    REQUIRE(p.second == false);
}

TEST_CASE("bcholesky decomposition throws for non-square matrices")
{
    Matrix<double> m = "[0,1]";
    CHECK_THROWS_WITH(bcholesky(m), "bcholesky: matrix must be square");
    m = "[0;1]";
    CHECK_THROWS_WITH(bcholesky(m), "bcholesky: matrix must be square");
}
