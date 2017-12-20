#include "catch.hpp"

#include <sstream>

#include "mesh.h"

using namespace Circuits;

TEST_CASE("the mesh generated is correct for N = 1")
{
    std::stringstream ss;
    generate(1, ss);
    REQUIRE(
        ss.str() == "1 0 1 1 0\n"
                    "1 2 0 1 0\n"
                    "1 3 0 1 0\n"
                    "2 4 0 1 0\n"
                    "3 4 0 1 0\n"
                    "3 5 0 1 0\n"
                    "4 0 0 1 0\n"
                    "5 0 0 1 0\n");
}

TEST_CASE("the mesh generated is correct for N = 2")
{
    std::stringstream ss;
    generate(2, ss);
    REQUIRE(
        ss.str() == "1 0 1 1 0\n"
                    "1 2 0 1 0\n"
                    "1 4 0 1 0\n"
                    "2 3 0 1 0\n"
                    "2 5 0 1 0\n"
                    "3 6 0 1 0\n"
                    "4 5 0 1 0\n"
                    "4 7 0 1 0\n"
                    "5 6 0 1 0\n"
                    "5 8 0 1 0\n"
                    "6 9 0 1 0\n"
                    "7 8 0 1 0\n"
                    "7 10 0 1 0\n"
                    "8 9 0 1 0\n"
                    "8 11 0 1 0\n"
                    "9 12 0 1 0\n"
                    "10 11 0 1 0\n"
                    "10 13 0 1 0\n"
                    "11 12 0 1 0\n"
                    "11 14 0 1 0\n"
                    "12 0 0 1 0\n"
                    "13 14 0 1 0\n"
                    "14 0 0 1 0\n");
}
