#pragma once

#include "mod_cell_2d.h"
#include "mod_vector_algebra.h"
#include "runtime/fortran_array.h"

namespace mod_fvm_face_2d {

using mod_cell_2d::face;
using mod_cell_2d::cell_2D;
using mod_cell_2d::list_cell_2D;
using mod_cell_2d::obj_cell_2D;
using mod_cell_2d::nbfaces;
using mod_vector_algebra::vector_t;

// ---------------------------------------------------------------------------
// Derived types
// ---------------------------------------------------------------------------

// fvm_face_2D extends face
struct fvm_face_2D : public face {
    vector_t  normal_vector;
    vector_t  tangent_vector;
    cell_2D*  left_cell  = nullptr;
    cell_2D*  right_cell = nullptr;
    double    len_nor    = 0.0;
    double    len_tan    = 0.0;
};

struct pt_fvm_face_2D {
    fvm_face_2D* f = nullptr;
};

struct obj_fvm_face_2D {
    FortranArray1D<pt_fvm_face_2D> face_2D;
    obj_cell_2D*                   smart_pointer = nullptr;

    // Type-bound procedures
    void alloc();  // => allocate_face_2D
    void get();    // => assign_face_2D

    // Constructor (from interface obj_fvm_face_2D)
    obj_fvm_face_2D() = default;
    obj_fvm_face_2D(obj_cell_2D& obj_cell_2D_container);
};

// ---------------------------------------------------------------------------
// Free function declarations
// ---------------------------------------------------------------------------
obj_fvm_face_2D obj_fvm_face_2D_constructor(obj_cell_2D& obj_cell_2D_container);
void allocate_face_2D(obj_fvm_face_2D& this_);
void assign_face_2D(obj_fvm_face_2D& this_);

} // namespace mod_fvm_face_2d
