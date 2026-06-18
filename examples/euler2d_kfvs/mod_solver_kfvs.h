#pragma once

#include "mod_cell_2d.h"
#include "mod_fvm_face_2d.h"
#include "mod_struct_to_array.h"
#include "mod_read_gmsh.h"
#include "mod_flux_kfvs.h"
#include "mod_objects.h"
#include "runtime/fortran_array.h"
#include "runtime/fortran_io.h"

namespace mod_solver_kfvs {

using mod_cell_2d::face;
using mod_cell_2d::cell_2D;
using mod_cell_2d::list_cell_2D;
using mod_cell_2d::obj_cell_2D;
using mod_cell_2d::nbfaces;
using mod_fvm_face_2d::fvm_face_2D;
using mod_fvm_face_2d::pt_fvm_face_2D;
using mod_fvm_face_2d::obj_fvm_face_2D;
using mod_struct_to_array::vol;
using mod_struct_to_array::norm_x;
using mod_struct_to_array::norm_y;
using mod_struct_to_array::len_norm;
using mod_struct_to_array::bc_typ;
using mod_struct_to_array::lr_cell;
using mod_objects::list_cell;
using mod_objects::faces_fvm;
using mod_read_gmsh::rho_init;
using mod_read_gmsh::ux_init;
using mod_read_gmsh::uy_init;
using mod_read_gmsh::t_init;
using mod_read_gmsh::fname;

// ---------------------------------------------------------------------------
// Module-level variables
// ---------------------------------------------------------------------------
extern double r_gaz;
extern double invdt;
extern double dt;
extern double tmax;
extern const double cfl;

extern FortranArray1D<double> rho;
extern FortranArray1D<double> ux;
extern FortranArray1D<double> uy;
extern FortranArray1D<double> t;
extern FortranArray1D<double> p;
extern FortranArray1D<double> a;
extern FortranArray1D<double> b;
extern FortranArray1D<double> e;

extern FortranArray2D<double> vect_u;
extern FortranArray2D<double> vect_unew;
extern FortranArray2D<double> flux;
extern FortranArray2D<double> rhs;
extern FortranArray2D<double> rhsdummy_symetrie;
extern FortranArray2D<double> rhsdummy_entree;
extern FortranArray2D<double> rhsdummy_sortie;
extern FortranArray2D<double> rhsdummy_paroi_solid;

extern FortranArray2D<double> vardummy_symetrie;
extern FortranArray2D<double> vardummy_entree;
extern FortranArray2D<double> vardummy_sortie;
extern FortranArray2D<double> vardummy_paroi_solid;

extern int nb_symmetry;
extern int nb_inlet;
extern int nb_outlet;
extern int nb_paroi_solid;
extern int nmax;

// ---------------------------------------------------------------------------
// Forward declaration for write_solution_vtk (defined in namespace solver)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Subroutine declarations
// ---------------------------------------------------------------------------
void donnee_initiale();
void allocate_vardummy();
void allocate_vardummy_multi_elem_airfoil();
void conditions_aux_limites();
void conditions_aux_limites_multi_elem_airfoil();
void calcul_derived_quantities();
void calcul_derived_quantities_multi_elem_airfoil();
void calcul_conservative_vector();
void timestep();
void assign_lr_cell();
void assign_lr_cell_multi_elem_airfoil();
void calcul_flux();
void calcul_flux_multi_elem_airfoil();
void calcul_rhs();
void euler_time_iteration();
void calcul_rho_ux_uy_t();
void write_pressure_coefficient(int iter);
void chk_converge(int iter);

} // namespace mod_solver_kfvs
