#pragma once

#include <vector>
#include <utility>
#include <functional>

namespace Numeric {

// Interpolates the given data using a single Lagrange polynomial.
std::function<double(double)>
    lagrange(const std::vector<std::pair<double, double>>& data);

// Interpolates the given data using piecewise cubic Hermite polynomials.
std::function<double(double)> hermite(
    const std::vector<std::pair<double, double>>& data,
    const std::vector<double>& derivatives);

// Interpolates the given data using piecewise linear polynomials.
std::function<double(double)>
    pwl(const std::vector<std::pair<double, double>>& data);

// Computes the derivative of the piecewise linear interpolation of the data.
std::function<double(double)>
    pwlprime(const std::vector<std::pair<double, double>>& data);

// Interpolates the inverse of the given data using piecewise linear
// polynomials. This function assumes that the y coordinates are increasing.
std::function<double(double)>
    pwlinverse(const std::vector<std::pair<double, double>>& data);

// Computes the slope between the given points.
double slope(std::pair<double, double> p1, std::pair<double, double> p2);

// Samples the given function at uniformly spaced points in [x0, x1].
// This function samples with a difference dx between points.
std::vector<std::pair<double, double>>
    sample(std::function<double(double)> f, double x0, double x1, double dx);

} // namespace Numeric
