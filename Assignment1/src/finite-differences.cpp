#include "finite-differences.h"

#include <cmath>
#include <algorithm>

using namespace Numeric;

namespace FiniteDifferences {

std::pair<Matrix<Node>, std::vector<std::pair<int, int>>> createGrid(double h)
{
    int rows = static_cast<int>(0.1 / h + 1);
    ++rows; // Free boundary nodes.
    int cols = static_cast<int>(0.1 / h + 1);
    ++cols; // Free boundary nodes.

    Matrix<Node> grid{rows, cols};
    std::vector<std::pair<int, int>> freeNodes;

    for (int i = 0; i < cols; ++i)
    {
        for (int j = 0; j < rows; ++j)
        {
            grid(j, i).x = i * h;
            grid(j, i).y = j * h;
            if (i * h >= 0.06 && j * h >= 0.08)
                grid(j, i).potential = 15; // Inner conductor.
            else if (i != 0 && j != 0)
                freeNodes.emplace_back(j, i);
        }
    }

    return {grid, freeNodes};
}

std::pair<Matrix<Node>, std::vector<std::pair<int, int>>>
    createGrid(std::vector<double> x, std::vector<double> y)
{
    std::sort(std::begin(x), std::end(x));
    std::sort(std::begin(y), std::end(y));

    // Free boundary nodes.
    x.push_back(2 * x[x.size() - 1] - x[x.size() - 2]);
    y.push_back(2 * y[y.size() - 1] - y[y.size() - 2]);

    int rows = y.size();
    int cols = x.size();

    Matrix<Node> grid{rows, cols};
    std::vector<std::pair<int, int>> freeNodes;

    for (int i = 0; i < cols; ++i)
    {
        for (int j = 0; j < rows; ++j)
        {
            grid(j, i).x = x[i];
            grid(j, i).y = y[j];
            if (x[i] >= 0.06 && y[j] >= 0.08)
                grid(j, i).potential = 15; // Inner conductor.
            else if (i != 0 && j != 0)
                freeNodes.emplace_back(j, i);
        }
    }

    return {grid, freeNodes};
}

// A class storing the factors by which adjacent nodes must be multiplied in the
// finite difference equation.
struct Node_internal
{
    double modx1{}, modx2{}, mody1{}, mody2{};
};

// Precomputes the factors needed by the finite difference formula.
std::vector<Node_internal> precompute(
    const Matrix<Node>& grid,
    const std::vector<std::pair<int, int>>& freeNodes);

// Updates the grid with the next iteration of the solution using successive
// over-relaxation with parameter w.
void updateSor(
    Matrix<Node>& grid,
    const std::vector<std::pair<int, int>>& freeNodes,
    const std::vector<Node_internal>& parameters,
    double w);

// Updates the grid with the next iteration of the solution using simultaneous
// relaxation with the Jacobi method.
void updateJacobi(
    Matrix<Node>& grid,
    const std::vector<std::pair<int, int>>& freeNodes,
    const std::vector<Node_internal>& parameters);

// Computes the residuals of the grid and returns true if all of them are below
// the relative tolerance tol.
bool stop(
    const Matrix<Node>& grid,
    const std::vector<std::pair<int, int>>& freeNodes,
    const std::vector<Node_internal>& parameters);

std::pair<Matrix<Node>, int>
    sor(Matrix<Node> grid,
        const std::vector<std::pair<int, int>>& freeNodes,
        double w)
{
    auto parameters = precompute(grid, freeNodes);
    int iterations = 0;
    do
    {
        ++iterations;
        updateSor(grid, freeNodes, parameters, w);
    } while (!stop(grid, freeNodes, parameters));
    return {grid, iterations};
}

std::pair<Matrix<Node>, int>
    jacobi(Matrix<Node> grid, const std::vector<std::pair<int, int>>& freeNodes)
{
    auto parameters = precompute(grid, freeNodes);
    int iterations = 0;
    do
    {
        ++iterations;
        updateJacobi(grid, freeNodes, parameters);
    } while (!stop(grid, freeNodes, parameters));
    return {grid, iterations};
}

std::vector<Node_internal> precompute(
    const Matrix<Node>& grid, const std::vector<std::pair<int, int>>& freeNodes)
{
    int n = freeNodes.size();
    std::vector<Node_internal> result(n);
    for (int i = 0; i < n; ++i)
    {
        auto p = freeNodes[i];
        auto& param = result[i];
        if (p.first == 0 || p.first == grid.rows() - 1 || p.second == 0 ||
            p.second == grid.cols() - 1)
        {
            // Neumann boundary. Do nothing.
        }
        else
        {
            // Compute the factors by which each adjacent node must be
            // multiplied in the finite difference equation.
            auto hx1 =
                grid(p.first, p.second).x - grid(p.first, p.second - 1).x;
            auto hx2 =
                grid(p.first, p.second + 1).x - grid(p.first, p.second).x;
            auto hy1 =
                grid(p.first, p.second).y - grid(p.first - 1, p.second).y;
            auto hy2 =
                grid(p.first + 1, p.second).y - grid(p.first, p.second).y;
            auto denom = 1 / (hx1 * hx2) + 1 / (hy1 * hy2);
            param.modx1 = 1 / (hx1 * (hx1 + hx2) * denom);
            param.modx2 = 1 / (hx2 * (hx1 + hx2) * denom);
            param.mody1 = 1 / (hy1 * (hy1 + hy2) * denom);
            param.mody2 = 1 / (hy2 * (hy1 + hy2) * denom);
        }
    }
    return result;
}

void updateSor(
    Matrix<Node>& grid,
    const std::vector<std::pair<int, int>>& freeNodes,
    const std::vector<Node_internal>& parameters,
    double w)
{
    int n = freeNodes.size();
    for (int i = 0; i < n; ++i)
    {
        auto p = freeNodes[i];
        auto param = parameters[i];
        if (p.first == 0)
        {
            grid(p.first, p.second).potential =
                grid(p.first + 2, p.second).potential; // Neumann boundary.
        }
        else if (p.first == grid.rows() - 1)
        {
            grid(p.first, p.second).potential =
                grid(p.first - 2, p.second).potential; // Neumann boundary.
        }
        else if (p.second == 0)
        {
            grid(p.first, p.second).potential =
                grid(p.first, p.second + 2).potential; // Neumann boundary.
        }
        else if (p.second == grid.cols() - 1)
        {
            grid(p.first, p.second).potential =
                grid(p.first, p.second - 2).potential; // Neumann boundary.
        }
        else
        {
            grid(p.first, p.second).potential =
                (1 - w) * grid(p.first, p.second).potential +
                w * (param.modx1 * grid(p.first, p.second - 1).potential +
                     param.modx2 * grid(p.first, p.second + 1).potential +
                     param.mody1 * grid(p.first - 1, p.second).potential +
                     param.mody2 * grid(p.first + 1, p.second).potential);
        }
    }
}

void updateJacobi(
    Matrix<Node>& grid,
    const std::vector<std::pair<int, int>>& freeNodes,
    const std::vector<Node_internal>& parameters)
{
    int n = freeNodes.size();
    Matrix<Node> temp = grid;
    for (int i = 0; i < n; ++i)
    {
        auto p = freeNodes[i];
        auto param = parameters[i];
        if (p.first == 0)
        {
            temp(p.first, p.second).potential =
                grid(p.first + 2, p.second).potential; // Neumann boundary.
        }
        else if (p.first == grid.rows() - 1)
        {
            temp(p.first, p.second).potential =
                grid(p.first - 2, p.second).potential; // Neumann boundary.
        }
        else if (p.second == 0)
        {
            temp(p.first, p.second).potential =
                grid(p.first, p.second + 2).potential; // Neumann boundary.
        }
        else if (p.second == grid.cols() - 1)
        {
            temp(p.first, p.second).potential =
                grid(p.first, p.second - 2).potential; // Neumann boundary.
        }
        else
        {
            temp(p.first, p.second).potential =
                param.modx1 * grid(p.first, p.second - 1).potential +
                param.modx2 * grid(p.first, p.second + 1).potential +
                param.mody1 * grid(p.first - 1, p.second).potential +
                param.mody2 * grid(p.first + 1, p.second).potential;
        }
    }
    grid = temp;
}

bool stop(
    const Matrix<Node>& grid,
    const std::vector<std::pair<int, int>>& freeNodes,
    const std::vector<Node_internal>& parameters)
{
    int n = freeNodes.size();
    for (int i = 0; i < n; ++i)
    {
        auto p = freeNodes[i];
        auto param = parameters[i];
        auto potential = grid(p.first, p.second).potential;
        if (p.first == 0)
        {
            // Neumann boundary.
            auto residual = grid(p.first + 2, p.second).potential - potential;
            if (potential != 0)
                residual /= potential;
            if (std::fabs(residual) > tol)
                return false;
        }
        else if (p.first == grid.rows() - 1)
        {
            // Neumann boundary.
            auto residual = grid(p.first - 2, p.second).potential - potential;
            if (potential != 0)
                residual /= potential;
            if (std::fabs(residual) > tol)
                return false;
        }
        else if (p.second == 0)
        {
            // Neumann boundary.
            auto residual = grid(p.first, p.second + 2).potential - potential;
            if (potential != 0)
                residual /= potential;
            if (std::fabs(residual) > tol)
                return false;
        }
        else if (p.second == grid.cols() - 1)
        {
            // Neumann boundary.
            auto residual = grid(p.first, p.second - 2).potential - potential;
            if (potential != 0)
                residual /= potential;
            if (std::fabs(residual) > tol)
                return false;
        }
        else
        {
            auto residual =
                param.modx1 * grid(p.first, p.second - 1).potential +
                param.modx2 * grid(p.first, p.second + 1).potential +
                param.mody1 * grid(p.first - 1, p.second).potential +
                param.mody2 * grid(p.first + 1, p.second).potential - potential;
            if (potential != 0)
                residual /= potential;
            if (std::fabs(residual) > tol)
                return false;
        }
    }
    return true;
}
}
