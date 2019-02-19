#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>

#include "matrix.h"
#include "matrix-util.h"
#include "solver.h"
#include "circuit-solver.h"
#include "mesh.h"
#include "finite-differences.h"

using namespace Numeric;
using namespace Circuits;
using namespace FiniteDifferences;

void question1();
void question2();
void question3();

void print(const Matrix<Node>& grid);

int main()
{
    question1();
    question2();
    question3();
    return 0;
}

void question1()
{
    std::cout << "========== Question 1 ==========" << std::endl;
    {
        std::cout << "Circuit 1" << std::endl;
        std::ifstream in{"../circuits/circuit-1.txt"};
        auto v = csolve(in);
        std::cout << "v = " << v;
    }
    {
        std::cout << "Circuit 2" << std::endl;
        std::ifstream in{"../circuits/circuit-2.txt"};
        auto v = csolve(in);
        std::cout << "v = " << v;
    }
    {
        std::cout << "Circuit 3" << std::endl;
        std::ifstream in{"../circuits/circuit-3.txt"};
        auto v = csolve(in);
        std::cout << "v = " << v;
    }
    {
        std::cout << "Circuit 4" << std::endl;
        std::ifstream in{"../circuits/circuit-4.txt"};
        auto v = csolve(in);
        std::cout << "v = " << v;
    }
    {
        std::cout << "Circuit 5" << std::endl;
        std::ifstream in{"../circuits/circuit-5.txt"};
        auto v = csolve(in);
        std::cout << "v = " << v;
    }
}

void question2()
{
    std::cout << "========== Question 2 ==========" << std::endl;
    std::cout << "Standard Solver" << std::endl;
    std::cout << "N,R(ohm),t(us)" << std::endl;
    for (int n = 2; n <= 10; ++n)
    {
        std::stringstream ss;
        generate(n, ss);
        Matrix<double> A, J, Y, E;
        std::tie(A, J, Y, E) = parse(ss);
        // Ignore Y and E for this problem since Y is the identity matrix
        // and E is zero.
        auto m = A * transpose(A);
        auto b = A * J;
        auto t1 = std::chrono::high_resolution_clock::now();
        auto v = solve(m, b);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
                      .count();
        std::cout << n << "," << v(0) / (1 - v(0)) << "," << us << std::endl;
    }
    std::cout << "Banded Solver" << std::endl;
    std::cout << "N,R(ohm),t(us)" << std::endl;
    for (int n = 2; n <= 10; ++n)
    {
        std::stringstream ss;
        generate(n, ss);
        Matrix<double> A, J, Y, E;
        std::tie(A, J, Y, E) = parse(ss);
        // Ignore Y and E for this problem since Y is the identity matrix
        // and E is zero.
        auto m = A * transpose(A);
        auto b = A * J;
        auto t1 = std::chrono::high_resolution_clock::now();
        auto v = bsolve(m, b);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
                      .count();
        std::cout << n << "," << v(0) / (1 - v(0)) << "," << us << std::endl;
    }
}

void question3()
{
    std::cout << "========== Question 3 ==========" << std::endl;

    double h = 0.02;
    double w = 1.4;
    Matrix<Node> grid;
    std::vector<std::pair<int, int>> freeNodes;
    int iterations;

    std::cout << "SOR" << std::endl;
    std::tie(grid, freeNodes) = createGrid(h);
    std::tie(grid, iterations) = sor(grid, freeNodes, w);
    std::cout << "Iterations: " << iterations << std::endl;
    print(grid);

    std::cout << "Jacobi" << std::endl;
    std::tie(grid, freeNodes) = createGrid(h);
    std::tie(grid, iterations) = jacobi(grid, freeNodes);
    std::cout << "Iterations: " << iterations << std::endl;
    print(grid);

    std::cout << "Comparing values of w for SOR with h = 0.02" << std::endl;
    h = 0.02;
    std::tie(grid, freeNodes) = createGrid(h);
    std::cout << "w,iterations,potential (V)" << std::endl;
    for (w = 1.0; w < 2.0; w += 0.1)
    {
        auto result = sor(grid, freeNodes, w);
        auto potential =
            result.first(static_cast<int>(0.04 / h), static_cast<int>(0.06 / h))
                .potential;
        std::cout << w << "," << result.second << "," << potential << std::endl;
    }

    std::cout << "Comparing values of h for SOR with w = 1.4" << std::endl;
    w = 1.4;
    std::cout << "h,iterations,potential (V)" << std::endl;
    std::vector<double> steps = {0.02, 0.01, 0.005, 0.004, 0.002, 0.001};
    for (auto h : steps)
    {
        std::tie(grid, freeNodes) = createGrid(h);
        auto result = sor(grid, freeNodes, w);
        auto potential =
            result.first(static_cast<int>(0.04 / h), static_cast<int>(0.06 / h))
                .potential;
        std::cout << h << "," << result.second << "," << potential << std::endl;
    }

    std::cout << "Comparing values of h for Jacobi" << std::endl;
    std::cout << "h,iterations,potential (V)" << std::endl;
    for (auto h : steps)
    {
        std::tie(grid, freeNodes) = createGrid(h);
        auto result = jacobi(grid, freeNodes);
        auto potential =
            result.first(static_cast<int>(0.04 / h), static_cast<int>(0.06 / h))
                .potential;
        std::cout << h << "," << result.second << "," << potential << std::endl;
    }

    std::cout << "Using smaller h around (0.06, 0.04) for SOR with w = 1.4"
              << std::endl;
    w = 1.4;
    std::vector<double> x = {
        0.0, 0.02, 0.04, 0.05, 0.055, 0.06, 0.065, 0.07, 0.08, 0.1};
    std::vector<double> y = {0.0,
                             0.02,
                             0.03,
                             0.04,
                             0.041,
                             0.042,
                             0.045,
                             0.05,
                             0.06,
                             0.07,
                             0.08,
                             0.1};
    std::tie(grid, freeNodes) = createGrid(x, y);
    std::tie(grid, iterations) = sor(grid, freeNodes, w);
    std::cout << "Iterations: " << iterations << std::endl;
    print(grid);
}

void print(const Matrix<Node>& grid)
{
    for (int i = grid.rows() - 1; i >= 0; --i)
    {
        std::cout << std::setw(7) << grid(i, 0).y << "  | ";
        for (int j = 0; j < grid.cols(); ++j)
            std::cout << std::setw(10) << grid(i, j).potential << " ";
        std::cout << std::endl;
    }
    std::cout << "         | ";
    for (int j = 0; j < grid.cols(); ++j)
        std::cout << "___________";
    std::cout << std::endl;
    std::cout << "           ";
    for (int j = 0; j < grid.cols(); ++j)
        std::cout << std::setw(10) << grid(0, j).x << " ";
    std::cout << std::endl;
}
