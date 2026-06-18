#pragma once

#include <string>
#include "runtime/fortran_array.h"
#include "runtime/fortran_io.h"
#include "runtime/fortran_namelist.h"

namespace mod_read_gmsh {

// ---------------------------------------------------------------------------
// Derived types
// ---------------------------------------------------------------------------

struct point {
    int    ident = 0;
    double x     = 0.0;
    double y     = 0.0;
    double z     = 0.0;
};

struct list_point {
    point* p = nullptr;
};

struct node_ident {
    FortranArray1D<int> id_node;
};

// node_ident_msh extends node_ident
struct node_ident_msh : public node_ident {
    int ident    = 0;
    int elem_typ = 0;
    int nb_tags  = 0;
    int tag1     = 0;
    int tag2     = 0;
};

struct ptr_node_ident {
    node_ident* pn = nullptr;
};

struct ptr_node_ident_msh {
    node_ident_msh* pn = nullptr;
};

struct obj_gmsh_reader {
    int nbnode   = 0;
    int nbel_msh = 0;
    int nbelm    = 0;

    FortranArray1D<ptr_node_ident>     id_nodes;
    FortranArray1D<ptr_node_ident_msh> id_nodes_msh;
    FortranArray1D<list_point>         coord_nodes;

    // Type-bound procedures
    void read_msh();          // => read_mesh
    void assign_id_nodes();   // => construct_id_nodes
};

// ---------------------------------------------------------------------------
// Module-level variables (save)
// ---------------------------------------------------------------------------
extern std::string fname;
extern double rho_init;
extern double ux_init;
extern double uy_init;
extern double t_init;

// ---------------------------------------------------------------------------
// Free subroutine declarations (also called as type-bound methods)
// ---------------------------------------------------------------------------
void read_mesh(obj_gmsh_reader& this_);
void construct_id_nodes(obj_gmsh_reader& this_);

} // namespace mod_read_gmsh
