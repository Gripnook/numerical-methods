#pragma once

#include <utility>
#include <functional>

#include "matrix.h"

namespace Numeric {

// Solves the equation f(x) = 0 using the Newton-Raphson method.
std::pair<double, int> newtonRaphson(
    std::function<double(double)> f,
    std::function<double(double)> fprime,
    double x0);

// Solves the equation x = f(x) using successive substitution.
std::pair<double, int>
    successiveSubstitution(std::function<double(double)> f, double x0);

// Solves the matrix equation f(x) = 0 using the Newton-Raphson method.
std::pair<Matrix<double>, int> newtonRaphson(
    std::function<Matrix<double>(const Matrix<double>&)> f,
    std::function<Matrix<double>(const Matrix<double>&)> jacobian,
    const Matrix<double>& x0,
    std::function<void(const Matrix<double>&, double, int)> callback =
        [](auto x, auto error, auto iteration) {});

} // namespace Numeric
