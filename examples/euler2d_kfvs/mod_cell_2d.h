#pragma once

#include "mod_read_gmsh.h"
#include "mod_point_util.h"
#include "mod_vector_algebra.h"
#include "runtime/fortran_array.h"

namespace mod_cell_2d {

using mod_read_gmsh::point;
using mod_read_gmsh::obj_gmsh_reader;
using mod_read_gmsh::node_ident_msh;
using mod_read_gmsh::ptr_node_ident_msh;

// ---------------------------------------------------------------------------
// Derived types
// ---------------------------------------------------------------------------

struct face {
    point   p1;
    point   p2;
    int     bc_typ   = 0;
    int     idface   = 0;
    double  area     = 0.0;
    point   centroid;
};

struct cell_2D {
    int     ident    = 0;
    point   vertex[4];          // vertex(1:4) in Fortran, 0-based in C++
    face    faces[4];           // faces(1:4) in Fortran, 0-based in C++
    cell_2D* neighbor1 = nullptr;
    cell_2D* neighbor2 = nullptr;
    cell_2D* neighbor3 = nullptr;
    cell_2D* neighbor4 = nullptr;
    double  vol      = 0.0;
};

struct list_cell_2D {
    cell_2D* p = nullptr;
};

// obj_cell_2D extends obj_gmsh_reader
struct obj_cell_2D : public obj_gmsh_reader {
    FortranArray1D<list_cell_2D> cell;

    // Type-bound procedures
    void build_cells();       // => construct_cells
    void calcul_vol();        // => calcul_vol_cells
    void calcul_face_info();  // => calcul_area_cent_faces
    void assign_face_id();    // => assign_id_face_v2
    void find_neighbor();     // => detect_neighbor
};

// ---------------------------------------------------------------------------
// Module-level variable
// ---------------------------------------------------------------------------
extern int nbfaces;

// ---------------------------------------------------------------------------
// Free subroutine declarations
// ---------------------------------------------------------------------------
void construct_cells(obj_cell_2D& this_);
void calcul_vol_cells(obj_cell_2D& this_);
void calcul_area_cent_faces(obj_cell_2D& this_);
void sort_vertex(point* vertices, int n);  // private, but declared for internal use
void assign_id_face_v2(obj_cell_2D& this_);
void detect_neighbor(obj_cell_2D& this_);

} // namespace mod_cell_2d
