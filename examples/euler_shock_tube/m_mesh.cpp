#include "m_mesh.h"
#include "m_params.h"
#include <cstdio>

namespace m_mesh {

// Module-level array definitions
FortranArray1D<double> x(m_params::nx);
FortranArray1D<double> xf(m_params::nx + 1);

void init_mesh()
{
    using namespace m_params;

    // Cell faces (uniform)
    for (int i = 1; i <= nx + 1; i++) {
        xf(i) = xmin + (i - 1) * dx;
    }

    // Cell centers
    for (int i = 1; i <= nx; i++) {
        x(i) = 0.5 * (xf(i) + xf(i + 1));
    }

    std::printf(" Mesh: nx=%6d  x=[%8.4f ,%8.4f\n", nx, xmin, xmax);
    std::printf(" dx=%12.4e\n", dx);
}

} // namespace m_mesh
