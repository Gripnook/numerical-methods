#pragma once

#include <utility>
#include <vector>

#include "matrix.h"

namespace FiniteDifferences {

// The relative tolerance for the solutions of the finite difference methods.
static const double tol = 1e-5;

// A class representing a single node within the finite difference grid.
struct Node
{
    double potential{};
    double x{}, y{};
};

// Creates the grid of nodes for the symmetric electric field problem inside a
// square conductor using the uniform spacing h. Returns this grid, as well as a
// list of free node coordinates within the grid.
std::pair<Numeric::Matrix<Node>, std::vector<std::pair<int, int>>>
    createGrid(double h);

// Creates the grid of nodes for the symmetric electric field problem inside a
// square conductor using the x and y coordinates specified. Returns this grid,
// as well as a list of free node coordinates within the grid. Note that this
// method will create the Neumann boundary internally, and this should not be
// specified in the x and y vectors.
std::pair<Numeric::Matrix<Node>, std::vector<std::pair<int, int>>>
    createGrid(std::vector<double> x, std::vector<double> y);

// Performs successive over-relaxation on the grid of free nodes for the given
// value of w. Returns the resulting grid and the number of iterations needed.
std::pair<Numeric::Matrix<Node>, int>
    sor(Numeric::Matrix<Node> grid,
        const std::vector<std::pair<int, int>>& freeNodes,
        double w);

// Performs simultaneous relaxation using the Jacobi method on the grid of free
// nodes. Returns the resulting grid and the number of iterations needed.
std::pair<Numeric::Matrix<Node>, int> jacobi(
    Numeric::Matrix<Node> grid,
    const std::vector<std::pair<int, int>>& freeNodes);
}
