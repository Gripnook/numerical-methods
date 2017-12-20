#include "mesh.h"

#include <iostream>

namespace Circuits {

void generate(int n, std::ostream& out)
{
    out << "1 0 1 1 0" << std::endl;
    for (int i = 0; i < 2 * n + 1; ++i)
    {
        for (int j = 0; j < n + 1; ++j)
        {
            if (j + 1 == n && i == 2 * n)
                out << (1 + j + i * (n + 1)) << " 0 0 1 0" << std::endl;
            else if (j != n)
                out << (1 + j + i * (n + 1)) << " " << (2 + j + i * (n + 1))
                    << " 0 1 0" << std::endl;

            if (j == n && i == 2 * n - 1)
                out << (1 + j + i * (n + 1)) << " 0 0 1 0" << std::endl;
            else if (i != 2 * n)
                out << (1 + j + i * (n + 1)) << " "
                    << (1 + j + (i + 1) * (n + 1)) << " 0 1 0" << std::endl;
        }
    }
}
}
