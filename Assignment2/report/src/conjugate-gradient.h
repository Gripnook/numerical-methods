#pragma once

#include <utility>
#include <numeric>
#include <functional>

#include "matrix.h"
#include "cholesky.h"

namespace Numeric {

// The relative tolerance for the solutions of the conjugate gradient method.
static const double tol = 1e-5;

// Solves the system of equations m * x = b, where m is a positive-definite
// [n x n] matrix, and b is an [n x 1] vector, using the un-preconditioned
// conjugate gradient method.
template <typename T>
std::pair<Matrix<T>, int> cgsolve(
    const Matrix<T>& m,
    const Matrix<T>& b,
    std::function<void(const Matrix<T>&)> callback = [](...) {})
{
    if (m.rows() != m.cols())
        throw std::runtime_error{"cgsolve: matrix must be square"};
    if (b.rows() != m.cols() || b.cols() != 1)
        throw std::runtime_error{"cgsolve: b must be an [nx1] vector"};

    if (!cholesky(m).second)
        throw std::runtime_error{"cgsolve: matrix must be positive-definite"};

    Matrix<T> x{b.rows(), b.cols()};
    int iterations = 0;
    callback(x);

    auto r = b - m * x;
    auto p = r;
    while (true)
    {
        auto denom = (transpose(p) * m * p);

        auto alpha = (transpose(p) * r) / denom;
        x = x + alpha * p;
        ++iterations;
        callback(x);

        r = b - m * x;
        auto twoNorm =
            std::sqrt(std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
        if (twoNorm <= tol)
            break;

        auto beta = -(transpose(p) * m * r) / denom;
        p = r + beta * p;
    }

    return {x, iterations};
}
}
