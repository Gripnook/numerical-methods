#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "matrix.h"
#include "matrix-util.h"
#include "solver.h"
#include "finite-differences.h"
#include "conjugate-gradient.h"

using namespace Numeric;

void question3();

void print(const Matrix<double>& nodes, const Matrix<double>& index);

int main()
{
    question3();
    return 0;
}

void question3()
{
    std::cout << "========== Question 3 ==========" << std::endl;
    double h = 0.02;
    Matrix<double> A, b, index;
    std::tie(A, b, index) = createGrid(h);

    std::cout << "Testing for Positive Definite" << std::endl;
    auto isPositiveDefinite = bcholesky(A);
    std::cout << "  A is positive definite: " << std::boolalpha
              << isPositiveDefinite.second << std::endl;
    isPositiveDefinite = bcholesky(transpose(A) * A);
    std::cout << "  A' A is positive definite: " << std::boolalpha
              << isPositiveDefinite.second << std::endl;

    std::cout << "Banded Cholesky Solver" << std::endl;
    auto bresult = bsolve(transpose(A) * A, transpose(A) * b);
    print(bresult, index);

    std::cout << "Conjugate Gradient Solver" << std::endl;
    std::cout << "iteration,2-norm,inf-norm" << std::endl;
    int iteration = 0;
    std::function<void(const Matrix<double>&)> callback =
        [&](const Matrix<double>& x) {
            auto r = b - A * x;
            auto twoNorm = std::sqrt(
                std::inner_product(r.begin(), r.end(), r.begin(), 0.0));
            auto infNorm = std::abs(
                *std::max_element(r.begin(), r.end(), [](auto lhs, auto rhs) {
                    return std::abs(lhs) < std::abs(rhs);
                }));
            std::cout << iteration++ << "," << twoNorm << "," << infNorm
                      << std::endl;
        };
    auto cgresult = cgsolve(transpose(A) * A, transpose(A) * b, callback);
    print(cgresult.first, index);
}

void print(const Matrix<double>& nodes, const Matrix<double>& index)
{
    for (int j = index.rows() - 1; j >= 0; --j)
    {
        for (int i = 0; i < index.cols(); ++i)
        {
            if (index(i, j) > 0)
            {
                std::cout << std::setw(10) << nodes(index(i, j) - 1) << " ";
            }
            else
            {
                std::cout << std::setw(10) << std::abs(index(i, j)) << " ";
            }
        }
        std::cout << std::endl;
    }
}
