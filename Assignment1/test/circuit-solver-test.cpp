#include "catch.hpp"

#include "circuit-solver.h"
#include "matrix.h"

using namespace Circuits;
using namespace Numeric;

// Truncates the least significant bits of the value for a better floating point
// comparison that ignores some round-off error.
float fround(double value)
{
    return static_cast<float>(value);
}

TEST_CASE("csolve succeeds for circuit 1")
{
    std::stringstream in{"1 0 0 10 0\n"
                         "1 0 0 10 10\n"};
    auto v = csolve(in);
    CHECK(fround(v(0)) == 5.0f);
}

TEST_CASE("csolve succeeds for circuit 2")
{
    std::stringstream in{"1 0 10 10 0\n"
                         "1 0 0 10 0\n"};
    auto v = csolve(in);
    CHECK(fround(v(0)) == 50.0f);
}

TEST_CASE("csolve succeeds for circuit 3")
{
    std::stringstream in{"1 0 0 10 10\n"
                         "1 0 10 10 0\n"};
    auto v = csolve(in);
    CHECK(fround(v(0)) == 55.0f);
}

TEST_CASE("csolve succeeds for circuit 4")
{
    std::stringstream in{"1 0 0 10 10\n"
                         "1 0 0 10 0\n"
                         "1 2 0 5 0\n"
                         "2 0 10 5 0\n"};
    auto v = csolve(in);
    CHECK(fround(v(0)) == 20.0f);
    CHECK(fround(v(1)) == 35.0f);
}

TEST_CASE("csolve succeeds for circuit 5")
{
    std::stringstream in{"1 0 0 20 10\n"
                         "1 2 0 10 0\n"
                         "1 3 0 10 0\n"
                         "2 3 0 30 0\n"
                         "2 0 0 30 0\n"
                         "3 0 0 30 0\n"};
    auto v = csolve(in);
    CHECK(fround(v(0)) == 5.0f);
    CHECK(fround(v(1)) == 3.75f);
    CHECK(fround(v(2)) == 3.75f);
}

TEST_CASE("cbsolve succeeds for circuit 1")
{
    std::stringstream in{"1 0 0 10 0\n"
                         "1 0 0 10 10\n"};
    auto v = cbsolve(in);
    CHECK(fround(v(0)) == 5.0f);
}

TEST_CASE("cbsolve succeeds for circuit 2")
{
    std::stringstream in{"1 0 10 10 0\n"
                         "1 0 0 10 0\n"};
    auto v = cbsolve(in);
    CHECK(fround(v(0)) == 50.0f);
}

TEST_CASE("cbsolve succeeds for circuit 3")
{
    std::stringstream in{"1 0 0 10 10\n"
                         "1 0 10 10 0\n"};
    auto v = cbsolve(in);
    CHECK(fround(v(0)) == 55.0f);
}

TEST_CASE("cbsolve succeeds for circuit 4")
{
    std::stringstream in{"1 0 0 10 10\n"
                         "1 0 0 10 0\n"
                         "1 2 0 5 0\n"
                         "2 0 10 5 0\n"};
    auto v = cbsolve(in);
    CHECK(fround(v(0)) == 20.0f);
    CHECK(fround(v(1)) == 35.0f);
}

TEST_CASE("cbsolve succeeds for circuit 5")
{
    std::stringstream in{"1 0 0 20 10\n"
                         "1 2 0 10 0\n"
                         "1 3 0 10 0\n"
                         "2 3 0 30 0\n"
                         "2 0 0 30 0\n"
                         "3 0 0 30 0\n"};
    auto v = cbsolve(in);
    CHECK(fround(v(0)) == 5.0f);
    CHECK(fround(v(1)) == 3.75f);
    CHECK(fround(v(2)) == 3.75f);
}

TEST_CASE("parse returns the correct matrices for circuit 1")
{
    Matrix<double> A, J, Y, E;
    std::stringstream in{"1 0 0 10 0\n"
                         "1 0 0 10 10\n"};
    std::tie(A, J, Y, E) = parse(in);
    CHECK(A == "[1,1]");
    CHECK(J == "[0;0]");
    CHECK(Y == "[0.1,0;0,0.1]");
    CHECK(E == "[0;10]");
}

TEST_CASE("parse returns the correct matrices for circuit 2")
{
    Matrix<double> A, J, Y, E;
    std::stringstream in{"1 0 10 10 0\n"
                         "1 0 0 10 0\n"};
    std::tie(A, J, Y, E) = parse(in);
    CHECK(A == "[1,1]");
    CHECK(J == "[10;0]");
    CHECK(Y == "[0.1,0;0,0.1]");
    CHECK(E == "[0;0]");
}

TEST_CASE("parse returns the nodes in sorted order")
{
    Matrix<double> A, J, Y, E;
    std::stringstream in{"1 0 10 10 0\n"
                         "4 3 0 10 0\n"
                         "2 1 0 10 0\n"
                         "3 2 0 10 0\n"
                         "4 0 0 10 0\n"};
    std::tie(A, J, Y, E) = parse(in);
    REQUIRE(
        A == "["
             "1, 0,-1, 0,0;"
             "0, 0, 1,-1,0;"
             "0,-1, 0, 1,0;"
             "0, 1, 0, 0,1"
             "]");
}

TEST_CASE("parse ignores comments and empty lines")
{
    Matrix<double> A, J, Y, E;
    std::stringstream in{"# comment\n"
                         "\n"
                         "1 0 0 10 0\n"
                         "\n"
                         "1 0 0 10 10\n"};
    std::tie(A, J, Y, E) = parse(in);
    CHECK(A == "[1,1]");
    CHECK(J == "[0;0]");
    CHECK(Y == "[0.1,0;0,0.1]");
    CHECK(E == "[0;10]");
}

TEST_CASE("parse throws an exception if the ground is missing")
{
    std::stringstream in{"1 2 0 10 0"};
    REQUIRE_THROWS_WITH(parse(in), "parse: missing ground node");
}

TEST_CASE("parse throws an exception for invalid inputs")
{
    std::stringstream in;
    in.clear();
    in.str(" ");
    CHECK_THROWS_WITH(parse(in), "parse: invalid node N+");
    in.clear();
    in.str("x 0 0 10 0");
    CHECK_THROWS_WITH(parse(in), "parse: invalid node N+");
    in.clear();
    in.str("0 x 0 10 0");
    CHECK_THROWS_WITH(parse(in), "parse: invalid node N-");
    in.clear();
    in.str("1 0 x 10 0");
    CHECK_THROWS_WITH(parse(in), "parse: invalid branch current Jk");
    in.clear();
    in.str("1 0 0 x 0");
    CHECK_THROWS_WITH(parse(in), "parse: invalid branch resistance Rk");
    in.clear();
    in.str("1 0 0 0 0");
    CHECK_THROWS_WITH(parse(in), "parse: branch resistance Rk cannot be zero");
    in.clear();
    in.str("1 0 0 10 x");
    CHECK_THROWS_WITH(parse(in), "parse: invalid branch voltage Ek");
}
