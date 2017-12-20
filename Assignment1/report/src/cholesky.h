#pragma once

#include <utility>
#include <cmath>

#include "matrix.h"

namespace Numeric {

// Returns a pair consisting of the lower triangular matrix that forms the
// cholesky decomposition of the matrix m, if it exists, and a boolean
// indicating if the function succeeded (the matrix m is positive-definite).
template <typename T>
std::pair<Matrix<T>, bool> cholesky(const Matrix<T>& m)
{
    if (m.rows() != m.cols())
        throw std::runtime_error{"cholesky: matrix must be square"};

    int n = m.rows();
    Matrix<T> lower{n, n};
    for (int j = 0; j < n; ++j)
    {
        auto result = m(j, j);
        for (int i = 0; i < j; ++i)
        {
            result -= lower(j, i) * lower(j, i);
        }
        // If the square of L(j,j) is negative, the square root would fail.
        // If the square of L(j,j) is zero, then the matrix is singular.
        if (result <= 0)
            return {lower, false};
        lower(j, j) = sqrt(result);

        for (int i = j + 1; i < n; ++i)
        {
            auto result = m(i, j);
            for (int k = 0; k < j; ++k)
            {
                result -= lower(i, k) * lower(j, k);
            }
            result /= lower(j, j);
            lower(i, j) = result;
        }
    }
    return {lower, true};
}

// Returns a pair consisting of the lower triangular matrix that forms the
// cholesky decomposition of the matrix m, if it exists, and a boolean
// indicating if the function succeeded (the matrix m is positive-definite).
// This function takes advantage of the sparsity of the matrix to use the
// half-bandwidth in the decomposition.
template <typename T>
std::pair<Matrix<T>, bool> bcholesky(const Matrix<T>& m)
{
    if (m.rows() != m.cols())
        throw std::runtime_error{"bcholesky: matrix must be square"};

    auto b = bandwidth(m);
    if (b.first != b.second)
        return {{}, false}; // Matrix is not symmetric.
    auto bw = b.first;

    int n = m.rows();
    Matrix<T> lower{n, n};
    for (int j = 0; j < n; ++j)
    {
        auto result = m(j, j);
        for (int i = std::max(j - bw + 1, 0); i < j; ++i)
        {
            result -= lower(j, i) * lower(j, i);
        }
        // If the square of L(j,j) is negative, the square root would fail.
        // If the square of L(j,j) is zero, then the matrix is singular.
        if (result <= 0)
            return {lower, false};
        lower(j, j) = sqrt(result);

        auto upperBound = std::min(j + bw, n);
        for (int i = j + 1; i < upperBound; ++i)
        {
            auto result = m(i, j);
            for (int k = std::max(i - bw + 1, 0); k < j; ++k)
            {
                result -= lower(i, k) * lower(j, k);
            }
            result /= lower(j, j);
            lower(i, j) = result;
        }
    }
    return {lower, true};
}
}
