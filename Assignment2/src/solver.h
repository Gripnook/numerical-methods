#pragma once

#include "matrix.h"
#include "cholesky.h"

namespace Numeric {

// Solves the system of equations L * x = b, where L is a non-singular,
// lower-triangular [n x n] matrix, and b is an [n x 1] vector. Note that this
// function assumes that L is lower-triangular and will only use that part of
// the matrix to solve without checking the validity of this assumption.
template <typename T>
Matrix<T> lsolve(const Matrix<T>& lower, const Matrix<T>& b);

// Solves the system of equations U * x = b, where U is a non-singular,
// upper-triangular [n x n] matrix, and b is an [n x 1] vector. Note that this
// function assumes that U is upper-triangular and will only use that part of
// the matrix to solve without checking the validity of this assumption.
template <typename T>
Matrix<T> usolve(const Matrix<T>& upper, const Matrix<T>& b);

// Solves the system of equations m * x = b, where m is a positive-definite
// [n x n] matrix, and b is an [n x 1] vector.
template <typename T>
Matrix<T> solve(const Matrix<T>& m, const Matrix<T>& b)
{
    if (m.rows() != m.cols())
        throw std::runtime_error{"solve: matrix must be square"};
    if (b.rows() != m.cols() || b.cols() != 1)
        throw std::runtime_error{"solve: b must be an [nx1] vector"};

    auto p = cholesky(m);
    if (!p.second)
        throw std::runtime_error{"solve: matrix must be positive-definite"};

    auto y = lsolve(p.first, b);
    auto x = usolve(transpose(p.first), y);
    return x;
}

// Solves the system of equations m * x = b, where m is a positive-definite
// [n x n] matrix, and b is an [n x 1] vector.
// This function takes advantage of the sparsity of the matrix to use the
// half-bandwidth in the cholesky decomposition.
template <typename T>
Matrix<T> bsolve(const Matrix<T>& m, const Matrix<T>& b)
{
    if (m.rows() != m.cols())
        throw std::runtime_error{"bsolve: matrix must be square"};
    if (b.rows() != m.cols() || b.cols() != 1)
        throw std::runtime_error{"bsolve: b must be an [nx1] vector"};

    auto p = bcholesky(m);
    if (!p.second)
        throw std::runtime_error{"bsolve: matrix must be positive-definite"};

    auto y = lsolve(p.first, b);
    auto x = usolve(transpose(p.first), y);
    return x;
}

template <typename T>
Matrix<T> lsolve(const Matrix<T>& lower, const Matrix<T>& b)
{
    if (lower.rows() != lower.cols())
        throw std::runtime_error{"lsolve: matrix must be square"};
    if (b.rows() != lower.cols() || b.cols() != 1)
        throw std::runtime_error{"lsolve: b must be an [nx1] vector"};

    int n = lower.rows();
    Matrix<T> x{n, 1};
    for (int i = 0; i < n; ++i)
    {
        auto result = b(i);
        for (int j = 0; j < i; ++j)
        {
            result -= lower(i, j) * x(j);
        }
        if (lower(i, i) == 0)
            throw std::runtime_error{"lsolve: matrix must be non-singular"};
        result /= lower(i, i);
        x(i) = result;
    }
    return x;
}

template <typename T>
Matrix<T> usolve(const Matrix<T>& upper, const Matrix<T>& b)
{
    if (upper.rows() != upper.cols())
        throw std::runtime_error{"usolve: matrix must be square"};
    if (b.rows() != upper.cols() || b.cols() != 1)
        throw std::runtime_error{"usolve: b must be an [nx1] vector"};

    int n = upper.rows();
    Matrix<T> x{n, 1};
    for (int i = n - 1; i >= 0; --i)
    {
        auto result = b(i);
        for (int j = n - 1; j >= i + 1; --j)
        {
            result -= upper(i, j) * x(j);
        }
        if (upper(i, i) == 0)
            throw std::runtime_error{"usolve: matrix must be non-singular"};
        result /= upper(i, i);
        x(i) = result;
    }
    return x;
}
}
