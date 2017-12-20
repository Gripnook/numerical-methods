#include "catch.hpp"

#include "matrix.h"
#include "matrix-util.h"

using namespace Numeric;

TEST_CASE("eye generates an identity matrix")
{
    auto m1 = eye<double>(1);
    CHECK(m1 == "[1]");
    auto m2 = eye<double>(2);
    CHECK(m2 == "[1,0;0,1]");
    auto m3 = eye<double>(3);
    CHECK(m3 == "[1,0,0;0,1,0;0,0,1]");
}
