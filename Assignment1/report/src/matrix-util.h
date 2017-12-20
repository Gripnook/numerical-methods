#pragma once

#include "matrix.h"

namespace Numeric {

// Returns an identity matrix of size [n x n].
template <typename T>
Matrix<T> eye(int n)
{
    Matrix<T> m{n, n};
    for (int i = 0; i < n; ++i)
        m(i, i) = 1;
    return m;
}
}
