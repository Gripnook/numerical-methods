#pragma once

#include <tuple>

#include "matrix.h"

namespace Numeric {

// Creates the matrices A and b that define a symmetric electric field problem
// inside a square conductor. A third parameter contains the mapping from the
// (x, y) grid to the index into the potential vector. If this mapping is a
// number i > 0, then it corresponds to index i-1. Otherwise, it corresponds to
// a fixed node with potential -i.
std::tuple<Matrix<double>, Matrix<double>, Matrix<double>> createGrid(double h);
}
