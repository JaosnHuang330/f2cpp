#pragma once

#include "mod_solver_kfvs.h"
#include "mod_struct_to_array.h"
#include "mod_read_gmsh.h"

namespace solver {

void solver();
void solver_kfvs();
void write_solution_vtk(int iter);

} // namespace solver
