#pragma once

#include "mod_objects.h"
#include "runtime/fortran_array.h"

namespace mod_struct_to_array {

using mod_objects::list_cell;
using mod_objects::faces_fvm;
using mod_cell_2d::nbfaces;

// ---------------------------------------------------------------------------
// Module-level variables
// ---------------------------------------------------------------------------
extern FortranArray1D<double> vol;
extern FortranArray1D<double> norm_x;
extern FortranArray1D<double> norm_y;
extern FortranArray1D<double> len_norm;
extern FortranArray1D<int>    bc_typ;
extern FortranArray2D<int>    lr_cell;

// ---------------------------------------------------------------------------
// Subroutine declarations
// ---------------------------------------------------------------------------
void struct_to_array();

} // namespace mod_struct_to_array
