#pragma once

#include <string>
#include "mod_read_gmsh.h"
#include "mod_cell_2d.h"
#include "runtime/fortran_array.h"
#include "runtime/fortran_io.h"

namespace mod_write_vtk {

using mod_read_gmsh::fname;
using mod_cell_2d::obj_cell_2D;

// ---------------------------------------------------------------------------
// Derived type: obj_msh_vis
// ---------------------------------------------------------------------------

struct obj_msh_vis {
    obj_cell_2D* smart_pointer = nullptr;

    // Type-bound procedures
    void write_vtk();      // => write_mesh_vtk
    void write_tecplot();  // => write_mesh_tecplot
};

// ---------------------------------------------------------------------------
// Constructor (interface obj_msh_vis)
// ---------------------------------------------------------------------------
obj_msh_vis obj_msh_vis_constructor(obj_cell_2D& obj_cell_2D_container);

// ---------------------------------------------------------------------------
// Free subroutine declarations
// ---------------------------------------------------------------------------
void write_mesh_vtk(obj_msh_vis& this_);
void write_mesh_tecplot(obj_msh_vis& this_);

} // namespace mod_write_vtk
