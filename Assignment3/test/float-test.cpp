#include "float-test.h"

float fround(double value)
{
    return static_cast<float>(value);
}

bool fequal(double lhs, double rhs)
{
    return fround(lhs) == fround(rhs);
}
