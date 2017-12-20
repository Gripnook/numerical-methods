#include "catch.hpp"

#include <sstream>

#include "matrix.h"
#include "matrix-util.h"

using namespace Numeric;

TEST_CASE("matrices can be added")
{
    Matrix<double> a = "[3,4;2,1]";
    Matrix<double> b = "[10,3;-2,4]";
    Matrix<double> result = "[13,7;0,5]";
    CHECK(a + b == result);
    a += b;
    CHECK(a == result);
}

TEST_CASE("adding matrices of different sizes throws an exception")
{
    Matrix<double> a = "[3,4;2,1]";
    Matrix<double> b = "[10;-2]";
    Matrix<double> c = "[1,2]";
    CHECK_THROWS_WITH(
        a + b,
        "add: inconsistent matrix dimensions: [2x2], "
        "[2x1]");
    CHECK_THROWS_WITH(
        a += b,
        "add: inconsistent matrix dimensions: [2x2], "
        "[2x1]");
    CHECK_THROWS_WITH(
        a + c,
        "add: inconsistent matrix dimensions: [2x2], "
        "[1x2]");
    CHECK_THROWS_WITH(
        a += c,
        "add: inconsistent matrix dimensions: [2x2], "
        "[1x2]");
}

TEST_CASE("matrices can be subtracted")
{
    Matrix<double> a = "[3,4;2,1]";
    Matrix<double> b = "[10,3;-2,4]";
    Matrix<double> result = "[-7,1;4,-3]";
    CHECK(a - b == result);
    a -= b;
    CHECK(a == result);
}

TEST_CASE("subtracting matrices of different sizes throws an exception")
{
    Matrix<double> a = "[3,4;2,1]";
    Matrix<double> b = "[10;-2]";
    Matrix<double> c = "[1,2]";
    CHECK_THROWS_WITH(
        a - b,
        "subtract: inconsistent matrix dimensions: [2x2], "
        "[2x1]");
    CHECK_THROWS_WITH(
        a -= b,
        "subtract: inconsistent matrix dimensions: [2x2], "
        "[2x1]");
    CHECK_THROWS_WITH(
        a - c,
        "subtract: inconsistent matrix dimensions: [2x2], "
        "[1x2]");
    CHECK_THROWS_WITH(
        a -= c,
        "subtract: inconsistent matrix dimensions: [2x2], "
        "[1x2]");
}

TEST_CASE("matrices can be multiplied")
{
    Matrix<double> a = "[3,4;2,1;1,4]";
    Matrix<double> b = "[10,3;-2,4]";
    Matrix<double> result = "[22,25;18,10;2,19]";
    CHECK(a * b == result);
    a *= b;
    CHECK(a == result);
}

TEST_CASE("multiplying matrices of inconsistent sizes throws an exception")
{
    Matrix<double> a = "[3,4;2,1]";
    Matrix<double> b = "[10;-2;3]";
    CHECK_THROWS_WITH(
        a * b,
        "multiply: inconsistent matrix dimensions: [2x2], "
        "[3x1]");
    CHECK_THROWS_WITH(
        a *= b,
        "multiply: inconsistent matrix dimensions: [2x2], "
        "[3x1]");
}

TEST_CASE("a matrix can be multiplied by a scalar")
{
    Matrix<double> a = "[1,2;3,4;5,6]";
    double k = 3.0;
    Matrix<double> result = "[3,6;9,12;15,18]";
    CHECK(k * a == result);
    CHECK(a * k == result);
    a *= k;
    CHECK(a == result);
}

TEST_CASE("a matrix can be divided by a non-zero scalar")
{
    Matrix<double> a = "[1,2;3,4;5,6]";
    double k = 2.0;
    Matrix<double> result = "[0.5,1.0;1.5,2.0;2.5,3.0]";
    CHECK(a / k == result);
    a /= k;
    CHECK(a == result);
}

TEST_CASE("division by zero throws an exception")
{
    Matrix<double> a = "[1,2;3,4]";
    double k = 0.0;
    CHECK_THROWS_WITH(a / k, "divide: division by zero");
    CHECK_THROWS_WITH(a /= k, "divide: division by zero");
}

TEST_CASE("a square matrix can be transposed")
{
    Matrix<double> m = "[1,2;4,5]";
    Matrix<double> result = "[1,4;2,5]";
    REQUIRE(transpose(m) == result);
}

TEST_CASE("a non-square matrix can be transposed")
{
    Matrix<double> m = "[1,2,3;4,5,6]";
    Matrix<double> result = "[1,4;2,5;3,6]";
    REQUIRE(transpose(m) == result);
}

TEST_CASE("a one-element matrix can be converted to its element type")
{
    Matrix<double> m = "[7.0]";
    double x = m;
    REQUIRE(x == 7.0);
}

TEST_CASE("converting a multi-element matrix to its element type throws an "
          "exception")
{
    Matrix<double> m = "[1,2;3,4]";
    double x = 0;
    REQUIRE_THROWS_WITH(x = m, "cannot convert matrix to single element");
}

TEST_CASE("the bandwidth of a [3 x 3] matrix is correctly computed")
{
    Matrix<double> m = eye<double>(3);
    auto b = bandwidth(m);
    CHECK(b.first == 1);
    CHECK(b.second == 1);
    m = "["
        "1,2,0;"
        "7,1,3;"
        "3,2,1"
        "]";
    b = bandwidth(m);
    CHECK(b.first == 3);
    CHECK(b.second == 2);
}

TEST_CASE("the bandwidth of a [4 x 4] matrix is correctly computed")
{
    Matrix<double> m = eye<double>(4);
    auto b = bandwidth(m);
    CHECK(b.first == 1);
    CHECK(b.second == 1);
    m = "["
        "1,2,0,0;"
        "7,1,3,0;"
        "3,2,1,0;"
        "0,0,0,1"
        "]";
    b = bandwidth(m);
    CHECK(b.first == 3);
    CHECK(b.second == 2);
}

TEST_CASE("the bandwidth of a [4 x 3] matrix is correctly computed")
{
    Matrix<double> m = "["
                       "1,2,0;"
                       "7,1,3;"
                       "0,2,1;"
                       "0,1,0"
                       "]";
    auto b = bandwidth(m);
    CHECK(b.first == 3);
    CHECK(b.second == 2);
}

TEST_CASE("a matrix can be read from a string")
{
    Matrix<double> m;
    std::stringstream ss{"[1,2,3;4,5,6]"};
    ss >> m;

    CHECK(m.size() == 6);
    CHECK(m.rows() == 2);
    CHECK(m.cols() == 3);

    CHECK(m(0, 0) == 1);
    CHECK(m(0, 1) == 2);
    CHECK(m(0, 2) == 3);
    CHECK(m(1, 0) == 4);
    CHECK(m(1, 1) == 5);
    CHECK(m(1, 2) == 6);

    CHECK(m(0) == 1);
    CHECK(m(1) == 2);
    CHECK(m(2) == 3);
    CHECK(m(3) == 4);
    CHECK(m(4) == 5);
    CHECK(m(5) == 6);
}

TEST_CASE("a row vector can be read from a string")
{
    Matrix<double> m;
    std::stringstream ss{"[1,2,3,4,5,6]"};
    ss >> m;

    CHECK(m.size() == 6);
    CHECK(m.rows() == 1);
    CHECK(m.cols() == 6);

    CHECK(m(0, 0) == 1);
    CHECK(m(0, 1) == 2);
    CHECK(m(0, 2) == 3);
    CHECK(m(0, 3) == 4);
    CHECK(m(0, 4) == 5);
    CHECK(m(0, 5) == 6);

    CHECK(m(0) == 1);
    CHECK(m(1) == 2);
    CHECK(m(2) == 3);
    CHECK(m(3) == 4);
    CHECK(m(4) == 5);
    CHECK(m(5) == 6);
}

TEST_CASE("a column vector can be read from a string")
{
    Matrix<double> m;
    std::stringstream ss{"[1;2;3;4;5;6]"};
    ss >> m;

    CHECK(m.size() == 6);
    CHECK(m.rows() == 6);
    CHECK(m.cols() == 1);

    CHECK(m(0, 0) == 1);
    CHECK(m(1, 0) == 2);
    CHECK(m(2, 0) == 3);
    CHECK(m(3, 0) == 4);
    CHECK(m(4, 0) == 5);
    CHECK(m(5, 0) == 6);

    CHECK(m(0) == 1);
    CHECK(m(1) == 2);
    CHECK(m(2) == 3);
    CHECK(m(3) == 4);
    CHECK(m(4) == 5);
    CHECK(m(5) == 6);
}

TEST_CASE("reading in badly formatted matrix strings throws exceptions")
{
    Matrix<double> m;
    std::stringstream ss;
    ss.clear();
    ss.str("");
    CHECK_THROWS_WITH(ss >> m, "input: missing '['");
    ss.clear();
    ss.str("1,2,3");
    CHECK_THROWS_WITH(ss >> m, "input: missing '['");
    ss.clear();
    ss.str("[1 2 3]");
    CHECK_THROWS_WITH(ss >> m, "input: invalid separator");
    ss.clear();
    ss.str("[1;2,3]");
    CHECK_THROWS_WITH(ss >> m, "input: inconsistent matrix dimensions");
    ss.clear();
    ss.str("1,2,3]");
    CHECK_THROWS_WITH(ss >> m, "input: missing '['");
    ss.clear();
    ss.str("[1,2,3");
    CHECK_THROWS_WITH(ss >> m, "input: missing ']'");
    ss.clear();
    ss.str("[a,b,c]");
    CHECK_THROWS_WITH(ss >> m, "input: invalid element");
    ss.clear();
    ss.str("[1'2'3]");
    CHECK_THROWS_WITH(ss >> m, "input: invalid separator");
}

TEST_CASE("a matrix can be constructed implicitly from a string")
{
    Matrix<double> m = "[1,2,3;4,5,6]";

    CHECK(m.size() == 6);
    CHECK(m.rows() == 2);
    CHECK(m.cols() == 3);

    CHECK(m(0, 0) == 1);
    CHECK(m(0, 1) == 2);
    CHECK(m(0, 2) == 3);
    CHECK(m(1, 0) == 4);
    CHECK(m(1, 1) == 5);
    CHECK(m(1, 2) == 6);

    CHECK(m(0) == 1);
    CHECK(m(1) == 2);
    CHECK(m(2) == 3);
    CHECK(m(3) == 4);
    CHECK(m(4) == 5);
    CHECK(m(5) == 6);
}
