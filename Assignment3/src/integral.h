#pragma once

#include <vector>
#include <functional>

namespace Numeric {

// Computes the definite integral of f(x) from x = a to x = b using
// Gauss-Legendre integration with n equal segments.
double integral(std::function<double(double)> f, double a, double b, int n);

// Computes the definite integral of f(x) starting from x = a using
// Gauss-Legendre integration with the given segment lengths.
double integral(
    std::function<double(double)> f,
    double a,
    const std::vector<double>& segments);

} // namespace Numeric
