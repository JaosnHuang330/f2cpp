#pragma once

#include "mod_read_gmsh.h"
#include "mod_write_vtk.h"
#include "mod_cell_2d.h"
#include "mod_fvm_face_2d.h"

namespace mod_objects {

using mod_cell_2d::obj_cell_2D;
using mod_write_vtk::obj_msh_vis;
using mod_fvm_face_2d::obj_fvm_face_2D;

// ---------------------------------------------------------------------------
// Module-level variables
// ---------------------------------------------------------------------------
extern obj_cell_2D    list_cell;   // type(obj_cell_2D), target
extern obj_msh_vis    visualizer;  // type(obj_msh_vis)
extern obj_fvm_face_2D faces_fvm;  // type(obj_fvm_face_2D)

} // namespace mod_objects
