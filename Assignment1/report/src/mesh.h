#pragma once

#include <iosfwd>

namespace Circuits {

// Generates the circuit description for a [N x 2N] mesh of 1 ohm resistors with
// a 1 ampere current source between the bottom-left and top-right corners of
// the mesh. The current source is placed in parallel with a 1 ohm resistor to
// satisfy the circuit requirements.
void generate(int n, std::ostream& out);
}
