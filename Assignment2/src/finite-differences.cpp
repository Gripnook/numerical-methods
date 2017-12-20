#include "finite-differences.h"

namespace Numeric {

std::tuple<Matrix<double>, Matrix<double>, Matrix<double>> createGrid(double h)
{
    int rows = static_cast<int>(0.1 / h + 1);
    ++rows; // Free boundary nodes.
    int cols = static_cast<int>(0.1 / h + 1);
    ++cols; // Free boundary nodes.

    Matrix<double> index{cols, rows};
    std::vector<std::pair<int, int>> reverseLookup;

    // Store the index of each node in the grid for lookup.
    for (int i = 0; i < cols; ++i)
    {
        for (int j = 0; j < rows; ++j)
        {
            if (i == 0 || j == 0)
            {
                index(i, j) = 0; // Outer conductor shell.
            }
            else if (i * h >= 0.06 && j * h >= 0.08)
            {
                if (i * h == 0.06 || j * h == 0.08)
                    // Inner conductor shell. Negative means fixed node.
                    index(i, j) = -15;
                else
                    index(i, j) = 0; // Inside inner conductor.
            }
            else
            {
                reverseLookup.emplace_back(i, j);
                index(i, j) = reverseLookup.size();
            }
        }
    }

    int n = reverseLookup.size();
    Matrix<double> A{n, n};
    Matrix<double> b{n, 1};
    for (int i = 0; i < n; ++i)
    {
        auto p = reverseLookup[i];
        if (p.first == cols - 1)
        {
            A(i, i) = -1;
            A(i, index(p.first - 2, p.second) - 1) = 1; // Neumann boundary.
        }
        else if (p.second == rows - 1)
        {
            A(i, i) = -1;
            A(i, index(p.first, p.second - 2) - 1) = 1; // Neumann boundary.
        }
        else
        {
            auto update = [&](int elem) {
                if (elem > 0)
                {
                    A(i, elem - 1) = 0.25; // Free node.
                }
                else
                {
                    b(i) += elem * 0.25; // Fixed node.
                }
            };
            A(i, i) = -1;
            update(index(p.first - 1, p.second));
            update(index(p.first, p.second - 1));
            update(index(p.first + 1, p.second));
            update(index(p.first, p.second + 1));
        }
    }

    return std::make_tuple(A, b, index);
}
}
