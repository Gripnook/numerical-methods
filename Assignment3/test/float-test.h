#pragma once

// Truncates the least significant bits of the value for a better floating point
// comparison that ignores some round-off error.
float fround(double value);

// Checks if the two floating point numbers are equal
// to within some round-off error.
bool fequal(double lhs, double rhs);
