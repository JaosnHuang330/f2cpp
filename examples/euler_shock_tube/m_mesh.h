#pragma once

#include "runtime/fortran_array.h"

namespace m_mesh {

// Cell centers (nx points)
extern FortranArray1D<double> x;

// Cell faces (nx+1 points)
extern FortranArray1D<double> xf;

void init_mesh();

} // namespace m_mesh
