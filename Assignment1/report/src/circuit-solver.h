#pragma once

#include <tuple>
#include <iosfwd>

#include "matrix.h"

namespace Circuits {

// Solves the circuit given by a formatted circuit description. Each line in the
// input stream is expected to either be a comment (beginning with #), or a
// branch description of the form (N+ N- Jk Rk Ek). The node voltages are
// returned in sorted order (e.g. a circuit with nodes [0,1,37,4,5] would have
// the following mapping: 1->0, 4->1, 5->2, 37->3).
Numeric::Matrix<double> csolve(std::istream& in);

// Solves the circuit given by a formatted circuit description. Each line in the
// input stream is expected to either be a comment (beginning with #), or a
// branch description of the form (N+ N- Jk Rk Ek). The node voltages are
// returned in sorted order (e.g. a circuit with nodes [0,1,37,4,5] would have
// the following mapping: 1->0, 4->1, 5->2, 37->3).
// This function takes advantage of the sparsity of matrices to use the
// half-bandwidth in the cholesky decomposition.
Numeric::Matrix<double> cbsolve(std::istream& in);

// Parses a formatted circuit description. Each line in the input stream is
// expected to either be a comment (beginning with #), or a branch description
// of the form (N+ N- Jk Rk Ek). Returns a tuple (A, J, Y, E) containing the
// reduced incidence matrix A, the current vector J, the admittance matrix Y,
// and the voltage vector E. The nodes are returned in sorted order in the
// incidence matrix, with ground ignored (e.g. a circuit with nodes [0,1,37,4,5]
// would have the following mapping: 1->0, 4->1, 5->2, 37->3).
std::tuple<
    Numeric::Matrix<double>,
    Numeric::Matrix<double>,
    Numeric::Matrix<double>,
    Numeric::Matrix<double>>
    parse(std::istream& in);
}
