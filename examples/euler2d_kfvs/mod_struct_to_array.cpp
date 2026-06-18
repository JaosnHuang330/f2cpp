#include "mod_struct_to_array.h"

namespace mod_struct_to_array {

using mod_objects::list_cell;
using mod_objects::faces_fvm;
using mod_cell_2d::nbfaces;

// ---------------------------------------------------------------------------
// Module-level variable definitions
// ---------------------------------------------------------------------------
FortranArray1D<double> vol;
FortranArray1D<double> norm_x;
FortranArray1D<double> norm_y;
FortranArray1D<double> len_norm;
FortranArray1D<int>    bc_typ;
FortranArray2D<int>    lr_cell;

// ---------------------------------------------------------------------------
// subroutine struct_to_array
// ---------------------------------------------------------------------------
void struct_to_array() {
    int icel, ifac;

    if (!vol.allocated()) {
        vol.allocate(list_cell.nbelm);
    }

    if (!norm_x.allocated()) {
        norm_x.allocate(nbfaces);
    }

    if (!norm_y.allocated()) {
        norm_y.allocate(nbfaces);
    }

    if (!len_norm.allocated()) {
        len_norm.allocate(nbfaces);
    }

    if (!bc_typ.allocated()) {
        bc_typ.allocate(nbfaces);
    }

    if (!lr_cell.allocated()) {
        lr_cell.allocate(nbfaces, 2);
    }

    for (icel = 1; icel <= list_cell.nbelm; icel++) {
        vol(icel) = list_cell.cell(icel).p->vol;
    }

    for (ifac = 1; ifac <= nbfaces; ifac++) {
        norm_x(ifac)   = faces_fvm.face_2D(ifac).f->normal_vector.x;
        norm_y(ifac)   = faces_fvm.face_2D(ifac).f->normal_vector.y;
        len_norm(ifac) = faces_fvm.face_2D(ifac).f->len_nor;
        bc_typ(ifac)   = faces_fvm.face_2D(ifac).f->bc_typ;
    }
}

} // namespace mod_struct_to_array
