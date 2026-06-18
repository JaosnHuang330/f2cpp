#include "mod_solver_kfvs.h"
#include "solver.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <array>

namespace mod_solver_kfvs {

// ---------------------------------------------------------------------------
// Module-level variable definitions
// ---------------------------------------------------------------------------
double r_gaz  = 0.0;
double invdt  = 0.0;
double dt     = 0.0;
double tmax   = 0.0;
const double cfl = 0.9;

FortranArray1D<double> rho;
FortranArray1D<double> ux;
FortranArray1D<double> uy;
FortranArray1D<double> t;
FortranArray1D<double> p;
FortranArray1D<double> a;
FortranArray1D<double> b;
FortranArray1D<double> e;

FortranArray2D<double> vect_u;
FortranArray2D<double> vect_unew;
FortranArray2D<double> flux;
FortranArray2D<double> rhs;
FortranArray2D<double> rhsdummy_symetrie;
FortranArray2D<double> rhsdummy_entree;
FortranArray2D<double> rhsdummy_sortie;
FortranArray2D<double> rhsdummy_paroi_solid;

FortranArray2D<double> vardummy_symetrie;
FortranArray2D<double> vardummy_entree;
FortranArray2D<double> vardummy_sortie;
FortranArray2D<double> vardummy_paroi_solid;

int nb_symmetry    = 0;
int nb_inlet       = 0;
int nb_outlet      = 0;
int nb_paroi_solid = 0;
int nmax           = 0;


// ---------------------------------------------------------------------------
// subroutine donnee_initiale
// ---------------------------------------------------------------------------
void donnee_initiale() {
    if (!rho.allocated()) {
        rho.allocate(list_cell.nbelm);
    }
    if (!ux.allocated()) {
        ux.allocate(list_cell.nbelm);
    }
    if (!uy.allocated()) {
        uy.allocate(list_cell.nbelm);
    }
    if (!t.allocated()) {
        t.allocate(list_cell.nbelm);
    }
    if (!p.allocated()) {
        p.allocate(list_cell.nbelm);
    }
    if (!a.allocated()) {
        a.allocate(list_cell.nbelm);
    }
    if (!b.allocated()) {
        b.allocate(list_cell.nbelm);
    }
    if (!e.allocated()) {
        e.allocate(list_cell.nbelm);
    }
    if (!vect_u.allocated()) {
        vect_u.allocate(list_cell.nbelm, 4);
    }
    if (!vect_unew.allocated()) {
        vect_unew.allocate(list_cell.nbelm, 4);
    }

    for (int i = 1; i <= list_cell.nbelm; ++i) {
        rho(i) = rho_init;
        ux(i)  = ux_init;
        uy(i)  = uy_init;
        t(i)   = t_init;
    }
}

// ---------------------------------------------------------------------------
// subroutine allocate_vardummy
// ---------------------------------------------------------------------------
void allocate_vardummy() {
    if (nb_paroi_solid > 0) return;

    // Allocating vardummy
    for (int i = 1; i <= nbfaces; ++i) {
        fvm_face_2D* pfac = faces_fvm.face_2D(i).f;
        if (pfac->bc_typ == 1) { // Airfoil - paroi solid
            nb_paroi_solid = nb_paroi_solid + 1;
        }
        if (pfac->bc_typ == 2) { // Inflow
            nb_inlet = nb_inlet + 1;
        }
        if (pfac->bc_typ == 3) { // Outflow
            nb_outlet = nb_outlet + 1;
        }
    }

    if (!vardummy_paroi_solid.allocated()) {
        vardummy_paroi_solid.allocate(nb_paroi_solid, 8);
        for (int i = 1; i <= nb_paroi_solid; ++i)
            for (int j = 1; j <= 8; ++j)
                vardummy_paroi_solid(i, j) = 0.0;
    }

    if (!vardummy_entree.allocated()) {
        vardummy_entree.allocate(nb_inlet, 8);
        for (int i = 1; i <= nb_inlet; ++i)
            for (int j = 1; j <= 8; ++j)
                vardummy_entree(i, j) = 0.0;
    }

    if (!vardummy_sortie.allocated()) {
        vardummy_sortie.allocate(nb_outlet, 8);
        for (int i = 1; i <= nb_outlet; ++i)
            for (int j = 1; j <= 8; ++j)
                vardummy_sortie(i, j) = 0.0;
    }
}

// ---------------------------------------------------------------------------
// subroutine allocate_vardummy_multi_elem_airfoil
// ---------------------------------------------------------------------------
void allocate_vardummy_multi_elem_airfoil() {
    //if (nb_symmetry > 0) return;
    if (nb_paroi_solid > 0) return;

    // Allocating vardummy
    for (int i = 1; i <= nbfaces; ++i) {
        fvm_face_2D* pfac = faces_fvm.face_2D(i).f;
        if (pfac->bc_typ == 1) { // inlet
            nb_inlet = nb_inlet + 1;
        }
        if (pfac->bc_typ == 2) { // Airfoil - paroi solid
            nb_paroi_solid = nb_paroi_solid + 1;
        }
        if (pfac->bc_typ == 3) { // Airfoil - paroi solid
            nb_paroi_solid = nb_paroi_solid + 1;
        }
        if (pfac->bc_typ == 4) { // Airfoil - paroi solid
            nb_paroi_solid = nb_paroi_solid + 1;
        }
    }

    //if (!vardummy_symetrie.allocated()) {
    //    vardummy_symetrie.allocate(nb_symmetry, 8);
    //    for (int i = 1; i <= nb_symmetry; ++i)
    //        for (int j = 1; j <= 8; ++j)
    //            vardummy_symetrie(i, j) = 0.0;
    //}

    if (!vardummy_paroi_solid.allocated()) {
        vardummy_paroi_solid.allocate(nb_paroi_solid, 8);
        for (int i = 1; i <= nb_paroi_solid; ++i)
            for (int j = 1; j <= 8; ++j)
                vardummy_paroi_solid(i, j) = 0.0;
    }

    if (!vardummy_entree.allocated()) {
        vardummy_entree.allocate(nb_inlet, 8);
        for (int i = 1; i <= nb_inlet; ++i)
            for (int j = 1; j <= 8; ++j)
                vardummy_entree(i, j) = 0.0;
    }

    //if (!vardummy_sortie.allocated()) {
    //    vardummy_sortie.allocate(nb_outlet, 8);
    //    for (int i = 1; i <= nb_outlet; ++i)
    //        for (int j = 1; j <= 8; ++j)
    //            vardummy_sortie(i, j) = 0.0;
    //}
}


// ---------------------------------------------------------------------------
// subroutine conditions_aux_limites
// ---------------------------------------------------------------------------
void conditions_aux_limites() {
    int cnt_inlet  = 0;
    int cnt_outlet = 0;
    int cnt_wall   = 0;

    for (int icel = 1; icel <= list_cell.nbelm; ++icel) {
        cell_2D* pcel = list_cell.cell(icel).p;
        for (int ifac = 1; ifac <= 4; ++ifac) {
            face* pfac = &pcel->faces[ifac - 1];

            if (pfac->bc_typ == 1) { // Airfoil - paroi solid
                cnt_wall = cnt_wall + 1;
                int idface = pfac->idface;
                double un =  norm_x(idface) * ux(icel) + norm_y(idface) * uy(icel);
                double ut = -norm_y(idface) * ux(icel) + norm_x(idface) * uy(icel);
                vardummy_paroi_solid(cnt_wall, 1) = rho(icel);
                vardummy_paroi_solid(cnt_wall, 2) = norm_x(idface) * (-un) - norm_y(idface) * ut;
                vardummy_paroi_solid(cnt_wall, 3) = norm_y(idface) * (-un) + norm_x(idface) * ut;
                vardummy_paroi_solid(cnt_wall, 4) = t(icel);
            }

            if (pfac->bc_typ == 2) { // Inflow
                cnt_inlet = cnt_inlet + 1;
                vardummy_entree(cnt_inlet, 1) = rho_init;
                vardummy_entree(cnt_inlet, 2) = ux_init;
                vardummy_entree(cnt_inlet, 3) = uy_init;
                vardummy_entree(cnt_inlet, 4) = t_init;
            }

            if (pfac->bc_typ == 3) { // Outflow
                cnt_outlet = cnt_outlet + 1;
                vardummy_sortie(cnt_outlet, 1) = rho(icel);
                vardummy_sortie(cnt_outlet, 2) = ux(icel);
                vardummy_sortie(cnt_outlet, 3) = uy(icel);
                vardummy_sortie(cnt_outlet, 4) = t(icel);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// subroutine conditions_aux_limites_multi_elem_airfoil
// ---------------------------------------------------------------------------
void conditions_aux_limites_multi_elem_airfoil() {
    int cnt_inlet = 0;
    int cnt_wall  = 0;

    for (int icel = 1; icel <= list_cell.nbelm; ++icel) {
        cell_2D* pcel = list_cell.cell(icel).p;
        for (int ifac = 1; ifac <= 4; ++ifac) {
            face* pfac = &pcel->faces[ifac - 1];

            if (pfac->bc_typ == 2 || pfac->bc_typ == 3 || pfac->bc_typ == 4) { // Airfoil - paroi solid
                cnt_wall = cnt_wall + 1;
                int idface = pfac->idface;
                double un =  norm_x(idface) * ux(icel) + norm_y(idface) * uy(icel);
                double ut = -norm_y(idface) * ux(icel) + norm_x(idface) * uy(icel);
                vardummy_paroi_solid(cnt_wall, 1) = rho(icel);
                vardummy_paroi_solid(cnt_wall, 2) = norm_x(idface) * (-un) - norm_y(idface) * ut;
                vardummy_paroi_solid(cnt_wall, 3) = norm_y(idface) * (-un) + norm_x(idface) * ut;
                vardummy_paroi_solid(cnt_wall, 4) = t(icel);
            }

            if (pfac->bc_typ == 1) { // Inflow
                cnt_inlet = cnt_inlet + 1;
                vardummy_entree(cnt_inlet, 1) = 0.2969689477e-4;
                vardummy_entree(cnt_inlet, 2) = 1059.458022;
                vardummy_entree(cnt_inlet, 3) = 0.0;
                vardummy_entree(cnt_inlet, 4) = 1295.646765;
            }
        }
    }
}


// ---------------------------------------------------------------------------
// subroutine calcul_derived_quantities
// ---------------------------------------------------------------------------
void calcul_derived_quantities() {
    r_gaz = 400.0; //1.3806503e-23 / 0.663e-25

    int n = list_cell.nbelm;
    for (int i = 1; i <= n; ++i) {
        p(i) = rho(i) * r_gaz * t(i);
        b(i) = std::sqrt(3.0 * r_gaz * t(i));
        a(i) = rho(i) / (8.0 * b(i) * b(i) * b(i));
        e(i) = 0.5 * rho(i) * (ux(i) * ux(i) + uy(i) * uy(i)) + 3.0 / 2.0 * rho(i) * r_gaz * t(i);
    }

    int nwall = vardummy_paroi_solid.size(1);
    for (int i = 1; i <= nwall; ++i) {
        vardummy_paroi_solid(i, 5) = vardummy_paroi_solid(i, 1) * r_gaz * vardummy_paroi_solid(i, 4);
        vardummy_paroi_solid(i, 7) = std::sqrt(3.0 * r_gaz * vardummy_paroi_solid(i, 4));
        vardummy_paroi_solid(i, 6) = vardummy_paroi_solid(i, 1) / (8.0 * vardummy_paroi_solid(i, 7) * vardummy_paroi_solid(i, 7) * vardummy_paroi_solid(i, 7));
        vardummy_paroi_solid(i, 8) = 0.5 * vardummy_paroi_solid(i, 1) * (vardummy_paroi_solid(i, 2) * vardummy_paroi_solid(i, 2) + vardummy_paroi_solid(i, 3) * vardummy_paroi_solid(i, 3)) +
            3.0 / 2.0 * vardummy_paroi_solid(i, 1) * r_gaz * vardummy_paroi_solid(i, 4);
    }

    int ninlet = vardummy_entree.size(1);
    for (int i = 1; i <= ninlet; ++i) {
        vardummy_entree(i, 5) = vardummy_entree(i, 1) * r_gaz * vardummy_entree(i, 4);
        vardummy_entree(i, 7) = std::sqrt(3.0 * r_gaz * vardummy_entree(i, 4));
        vardummy_entree(i, 6) = vardummy_entree(i, 1) / (8.0 * vardummy_entree(i, 7) * vardummy_entree(i, 7) * vardummy_entree(i, 7));
        vardummy_entree(i, 8) = 0.5 * vardummy_entree(i, 1) * (vardummy_entree(i, 2) * vardummy_entree(i, 2) + vardummy_entree(i, 3) * vardummy_entree(i, 3)) +
            3.0 / 2.0 * vardummy_entree(i, 1) * r_gaz * vardummy_entree(i, 4);
    }

    int noutlet = vardummy_sortie.size(1);
    for (int i = 1; i <= noutlet; ++i) {
        vardummy_sortie(i, 5) = vardummy_sortie(i, 1) * r_gaz * vardummy_sortie(i, 4);
        vardummy_sortie(i, 7) = std::sqrt(3.0 * r_gaz * vardummy_sortie(i, 4));
        vardummy_sortie(i, 6) = vardummy_sortie(i, 1) / (8.0 * vardummy_sortie(i, 7) * vardummy_sortie(i, 7) * vardummy_sortie(i, 7));
        vardummy_sortie(i, 8) = 0.5 * vardummy_sortie(i, 1) * (vardummy_sortie(i, 2) * vardummy_sortie(i, 2) + vardummy_sortie(i, 3) * vardummy_sortie(i, 3)) +
            3.0 / 2.0 * vardummy_sortie(i, 1) * r_gaz * vardummy_sortie(i, 4);
    }
}

// ---------------------------------------------------------------------------
// subroutine calcul_derived_quantities_multi_elem_airfoil
// ---------------------------------------------------------------------------
void calcul_derived_quantities_multi_elem_airfoil() {
    r_gaz = 1.3806503e-23 / 0.663e-25;

    int n = list_cell.nbelm;
    for (int i = 1; i <= n; ++i) {
        p(i) = rho(i) * r_gaz * t(i);
        b(i) = std::sqrt(3.0 * r_gaz * t(i));
        a(i) = rho(i) / (8.0 * b(i) * b(i) * b(i));
        e(i) = 0.5 * rho(i) * (ux(i) * ux(i) + uy(i) * uy(i)) + 3.0 / 2.0 * rho(i) * r_gaz * t(i);
    }

    int nwall = vardummy_paroi_solid.size(1);
    for (int i = 1; i <= nwall; ++i) {
        vardummy_paroi_solid(i, 5) = vardummy_paroi_solid(i, 1) * r_gaz * vardummy_paroi_solid(i, 4);
        vardummy_paroi_solid(i, 7) = std::sqrt(3.0 * r_gaz * vardummy_paroi_solid(i, 4));
        vardummy_paroi_solid(i, 6) = vardummy_paroi_solid(i, 1) / (8.0 * vardummy_paroi_solid(i, 7) * vardummy_paroi_solid(i, 7) * vardummy_paroi_solid(i, 7));
        vardummy_paroi_solid(i, 8) = 0.5 * vardummy_paroi_solid(i, 1) * (vardummy_paroi_solid(i, 2) * vardummy_paroi_solid(i, 2) + vardummy_paroi_solid(i, 3) * vardummy_paroi_solid(i, 3)) +
            3.0 / 2.0 * vardummy_paroi_solid(i, 1) * r_gaz * vardummy_paroi_solid(i, 4);
    }

    int ninlet = vardummy_entree.size(1);
    for (int i = 1; i <= ninlet; ++i) {
        vardummy_entree(i, 5) = vardummy_entree(i, 1) * r_gaz * vardummy_entree(i, 4);
        vardummy_entree(i, 7) = std::sqrt(3.0 * r_gaz * vardummy_entree(i, 4));
        vardummy_entree(i, 6) = vardummy_entree(i, 1) / (8.0 * vardummy_entree(i, 7) * vardummy_entree(i, 7) * vardummy_entree(i, 7));
        vardummy_entree(i, 8) = 0.5 * vardummy_entree(i, 1) * (vardummy_entree(i, 2) * vardummy_entree(i, 2) + vardummy_entree(i, 3) * vardummy_entree(i, 3)) +
            3.0 / 2.0 * vardummy_entree(i, 1) * r_gaz * vardummy_entree(i, 4);
    }
}


// ---------------------------------------------------------------------------
// subroutine calcul_conservative_vector
// ---------------------------------------------------------------------------
void calcul_conservative_vector() {
    int n = list_cell.nbelm;
    for (int i = 1; i <= n; ++i) {
        vect_u(i, 1) = rho(i);
        vect_u(i, 2) = rho(i) * ux(i);
        vect_u(i, 3) = rho(i) * uy(i);
        vect_u(i, 4) = e(i);
    }
    // vect_unew = vect_u
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= 4; ++j) {
            vect_unew(i, j) = vect_u(i, j);
        }
    }
}

// ---------------------------------------------------------------------------
// subroutine timestep
// ---------------------------------------------------------------------------
void timestep() {
    invdt = 0.0;

    for (int i = 1; i <= list_cell.nbelm; ++i) {
        cell_2D* pc = list_cell.cell(i).p;
        double norme_u = std::sqrt(ux(i) * ux(i) + uy(i) * uy(i));

        face* pf1 = &pc->faces[0];
        face* pf2 = &pc->faces[1];
        face* pf3 = &pc->faces[2];
        face* pf4 = &pc->faces[3];

        int face1 = pf1->idface;
        int face2 = pf2->idface;
        int face3 = pf3->idface;
        int face4 = pf4->idface;

        if (face1 == 0) {
            printf(" face1 = 0\n");
            printf(" Please check cell %d\n", i);
        }

        if (face2 == 0) {
            printf(" face2 = 0\n");
        }

        if (face3 == 0) {
            printf(" face3 = 0\n");
        }

        if (face4 == 0) {
            printf(" face4 = 0\n");
            printf(" Please check cell %d\n", i);
        }

        double perimetre = faces_fvm.face_2D(face1).f->len_nor +
                           faces_fvm.face_2D(face2).f->len_nor +
                           faces_fvm.face_2D(face3).f->len_nor +
                           faces_fvm.face_2D(face4).f->len_nor;
        double candidate = norme_u + b(i) * perimetre / pc->vol;
        if (candidate > invdt) invdt = candidate;
    }
    dt = cfl / invdt;
}

// ---------------------------------------------------------------------------
// subroutine assign_lr_cell
// ---------------------------------------------------------------------------
void assign_lr_cell() {
    int cnt_inlet  = 0;
    int cnt_outlet = 0;
    int cnt_wall   = 0;

    // Create left cell - right cell table for boundary faces
    for (int icel = 1; icel <= list_cell.nbelm; ++icel) {
        cell_2D* pcel = list_cell.cell(icel).p;
        for (int ifac = 1; ifac <= 4; ++ifac) {
            face* pfac = &pcel->faces[ifac - 1];

            if (pfac->bc_typ == 1) { // Airfoil - paroi solid
                cnt_wall = cnt_wall + 1;
                int fac = pfac->idface;
                lr_cell(fac, 1) = icel;
                lr_cell(fac, 2) = cnt_wall; // dummy cell for solid wall bc
            }

            if (pfac->bc_typ == 2) { // Inflow
                cnt_inlet = cnt_inlet + 1;
                int fac = pfac->idface;
                lr_cell(fac, 1) = icel;
                lr_cell(fac, 2) = cnt_inlet; // dummy cell for inlet bc
            }

            if (pfac->bc_typ == 3) { // Outflow
                cnt_outlet = cnt_outlet + 1;
                int fac = pfac->idface;
                lr_cell(fac, 1) = icel;
                lr_cell(fac, 2) = cnt_outlet; // dummy cell for outlet bc
            }
        }
    }

    // Create left cell - right cell table for internal faces
    for (int ifac = 1; ifac <= nbfaces; ++ifac) {
        fvm_face_2D* pfac_fvm = faces_fvm.face_2D(ifac).f;
        if (pfac_fvm->left_cell != nullptr && pfac_fvm->right_cell != nullptr) {
            lr_cell(ifac, 1) = pfac_fvm->left_cell->ident;
            lr_cell(ifac, 2) = pfac_fvm->right_cell->ident;
        }
    }
}

// ---------------------------------------------------------------------------
// subroutine assign_lr_cell_multi_elem_airfoil
// ---------------------------------------------------------------------------
void assign_lr_cell_multi_elem_airfoil() {
    int cnt_inlet = 0;
    int cnt_wall  = 0;

    // Create left cell - right cell table for boundary faces
    for (int icel = 1; icel <= list_cell.nbelm; ++icel) {
        cell_2D* pcel = list_cell.cell(icel).p;
        for (int ifac = 1; ifac <= 4; ++ifac) {
            face* pfac = &pcel->faces[ifac - 1];

            if (pfac->bc_typ == 2 || pfac->bc_typ == 3 || pfac->bc_typ == 4) { // Airfoil - paroi solid
                cnt_wall = cnt_wall + 1;
                int fac = pfac->idface;
                lr_cell(fac, 1) = icel;
                lr_cell(fac, 2) = cnt_wall; // dummy cell for solid wall bc
            }

            if (pfac->bc_typ == 1) { // Inflow
                cnt_inlet = cnt_inlet + 1;
                int fac = pfac->idface;
                lr_cell(fac, 1) = icel;
                lr_cell(fac, 2) = cnt_inlet; // dummy cell for inlet bc
            }
        }
    }

    // Create left cell - right cell table for internal faces
    for (int ifac = 1; ifac <= nbfaces; ++ifac) {
        fvm_face_2D* pfac_fvm = faces_fvm.face_2D(ifac).f;
        if (pfac_fvm->left_cell != nullptr && pfac_fvm->right_cell != nullptr) {
            lr_cell(ifac, 1) = pfac_fvm->left_cell->ident;
            lr_cell(ifac, 2) = pfac_fvm->right_cell->ident;
        }
    }
}


// ---------------------------------------------------------------------------
// subroutine calcul_flux
// ---------------------------------------------------------------------------
void calcul_flux() {
    using mod_flux_kfvs::fluxp;
    using mod_flux_kfvs::fluxm;

    if (!flux.allocated()) {
        flux.allocate(nbfaces, 4);
    }

    for (int ifac = 1; ifac <= nbfaces; ++ifac) {
        int left_cell  = lr_cell(ifac, 1);
        int right_cell = lr_cell(ifac, 2);

        if (bc_typ(ifac) == 0) {
            std::array<double,4> flux_plus  = fluxp(rho(left_cell), ux(left_cell), uy(left_cell),
                e(left_cell), p(left_cell), t(left_cell), a(left_cell), b(left_cell),
                norm_x(ifac), norm_y(ifac));
            std::array<double,4> flux_minus = fluxm(rho(right_cell), ux(right_cell), uy(right_cell),
                e(right_cell), p(right_cell), t(right_cell), a(right_cell), b(right_cell),
                norm_x(ifac), norm_y(ifac));
            for (int k = 1; k <= 4; ++k)
                flux(ifac, k) = len_norm(ifac) * (flux_plus[k-1] + flux_minus[k-1]);
        }

        if (bc_typ(ifac) == 1) { // solid wall
            std::array<double,4> flux_plus  = fluxp(rho(left_cell), ux(left_cell), uy(left_cell),
                e(left_cell), p(left_cell), t(left_cell), a(left_cell), b(left_cell),
                norm_x(ifac), norm_y(ifac));
            std::array<double,4> flux_minus = fluxm(
                vardummy_paroi_solid(right_cell, 1), vardummy_paroi_solid(right_cell, 2),
                vardummy_paroi_solid(right_cell, 3), vardummy_paroi_solid(right_cell, 8),
                vardummy_paroi_solid(right_cell, 5), vardummy_paroi_solid(right_cell, 4),
                vardummy_paroi_solid(right_cell, 6), vardummy_paroi_solid(right_cell, 7),
                norm_x(ifac), norm_y(ifac));
            for (int k = 1; k <= 4; ++k)
                flux(ifac, k) = len_norm(ifac) * (flux_plus[k-1] + flux_minus[k-1]);
        }

        if (bc_typ(ifac) == 2) { // Inflow
            std::array<double,4> flux_plus  = fluxp(rho(left_cell), ux(left_cell), uy(left_cell),
                e(left_cell), p(left_cell), t(left_cell), a(left_cell), b(left_cell),
                norm_x(ifac), norm_y(ifac));
            std::array<double,4> flux_minus = fluxm(
                vardummy_entree(right_cell, 1), vardummy_entree(right_cell, 2),
                vardummy_entree(right_cell, 3), vardummy_entree(right_cell, 8),
                vardummy_entree(right_cell, 5), vardummy_entree(right_cell, 4),
                vardummy_entree(right_cell, 6), vardummy_entree(right_cell, 7),
                norm_x(ifac), norm_y(ifac));
            for (int k = 1; k <= 4; ++k)
                flux(ifac, k) = len_norm(ifac) * (flux_plus[k-1] + flux_minus[k-1]);
        }

        if (bc_typ(ifac) == 3) { // Outflow
            std::array<double,4> flux_plus  = fluxp(rho(left_cell), ux(left_cell), uy(left_cell),
                e(left_cell), p(left_cell), t(left_cell), a(left_cell), b(left_cell),
                norm_x(ifac), norm_y(ifac));
            std::array<double,4> flux_minus = fluxm(
                vardummy_sortie(right_cell, 1), vardummy_sortie(right_cell, 2),
                vardummy_sortie(right_cell, 3), vardummy_sortie(right_cell, 8),
                vardummy_sortie(right_cell, 5), vardummy_sortie(right_cell, 4),
                vardummy_sortie(right_cell, 6), vardummy_sortie(right_cell, 7),
                norm_x(ifac), norm_y(ifac));
            for (int k = 1; k <= 4; ++k)
                flux(ifac, k) = len_norm(ifac) * (flux_plus[k-1] + flux_minus[k-1]);
        }
    }
}

// ---------------------------------------------------------------------------
// subroutine calcul_flux_multi_elem_airfoil
// ---------------------------------------------------------------------------
void calcul_flux_multi_elem_airfoil() {
    using mod_flux_kfvs::fluxp;
    using mod_flux_kfvs::fluxm;

    if (!flux.allocated()) {
        flux.allocate(nbfaces, 4);
    }

    for (int ifac = 1; ifac <= nbfaces; ++ifac) {
        int left_cell  = lr_cell(ifac, 1);
        int right_cell = lr_cell(ifac, 2);

        if (bc_typ(ifac) == 0) {
            std::array<double,4> flux_plus  = fluxp(rho(left_cell), ux(left_cell), uy(left_cell),
                e(left_cell), p(left_cell), t(left_cell), a(left_cell), b(left_cell),
                norm_x(ifac), norm_y(ifac));
            std::array<double,4> flux_minus = fluxm(rho(right_cell), ux(right_cell), uy(right_cell),
                e(right_cell), p(right_cell), t(right_cell), a(right_cell), b(right_cell),
                norm_x(ifac), norm_y(ifac));
            for (int k = 1; k <= 4; ++k)
                flux(ifac, k) = len_norm(ifac) * (flux_plus[k-1] + flux_minus[k-1]);
        }

        if (bc_typ(ifac) == 2 || bc_typ(ifac) == 3 || bc_typ(ifac) == 4) { // solid wall
            std::array<double,4> flux_plus  = fluxp(rho(left_cell), ux(left_cell), uy(left_cell),
                e(left_cell), p(left_cell), t(left_cell), a(left_cell), b(left_cell),
                norm_x(ifac), norm_y(ifac));
            std::array<double,4> flux_minus = fluxm(
                vardummy_paroi_solid(right_cell, 1), vardummy_paroi_solid(right_cell, 2),
                vardummy_paroi_solid(right_cell, 3), vardummy_paroi_solid(right_cell, 8),
                vardummy_paroi_solid(right_cell, 5), vardummy_paroi_solid(right_cell, 4),
                vardummy_paroi_solid(right_cell, 6), vardummy_paroi_solid(right_cell, 7),
                norm_x(ifac), norm_y(ifac));
            for (int k = 1; k <= 4; ++k)
                flux(ifac, k) = len_norm(ifac) * (flux_plus[k-1] + flux_minus[k-1]);
        }

        if (bc_typ(ifac) == 1) { // Inflow
            std::array<double,4> flux_plus  = fluxp(rho(left_cell), ux(left_cell), uy(left_cell),
                e(left_cell), p(left_cell), t(left_cell), a(left_cell), b(left_cell),
                norm_x(ifac), norm_y(ifac));
            std::array<double,4> flux_minus = fluxm(
                vardummy_entree(right_cell, 1), vardummy_entree(right_cell, 2),
                vardummy_entree(right_cell, 3), vardummy_entree(right_cell, 8),
                vardummy_entree(right_cell, 5), vardummy_entree(right_cell, 4),
                vardummy_entree(right_cell, 6), vardummy_entree(right_cell, 7),
                norm_x(ifac), norm_y(ifac));
            for (int k = 1; k <= 4; ++k)
                flux(ifac, k) = len_norm(ifac) * (flux_plus[k-1] + flux_minus[k-1]);
        }
    }
}


// ---------------------------------------------------------------------------
// subroutine calcul_rhs
// ---------------------------------------------------------------------------
void calcul_rhs() {
    if (!rhs.allocated()) {
        rhs.allocate(list_cell.nbelm, 4);
    }

    // rhs = 0.0
    for (int i = 1; i <= list_cell.nbelm; ++i)
        for (int j = 1; j <= 4; ++j)
            rhs(i, j) = 0.0;

    for (int ifac = 1; ifac <= nbfaces; ++ifac) {
        int left_cell  = lr_cell(ifac, 1);
        int right_cell = lr_cell(ifac, 2);

        if (bc_typ(ifac) == 0) {
            for (int k = 1; k <= 4; ++k) {
                rhs(left_cell,  k) = rhs(left_cell,  k) - flux(ifac, k);
                rhs(right_cell, k) = rhs(right_cell, k) + flux(ifac, k);
            }
        }

        if (bc_typ(ifac) != 0) {
            for (int k = 1; k <= 4; ++k) {
                rhs(left_cell, k) = rhs(left_cell, k) - flux(ifac, k);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// subroutine euler_time_iteration
// ---------------------------------------------------------------------------
void euler_time_iteration() {
    for (int icel = 1; icel <= list_cell.nbelm; ++icel) {
        for (int k = 1; k <= 4; ++k) {
            vect_unew(icel, k) = vect_u(icel, k) + dt / vol(icel) * rhs(icel, k);
        }
    }
}

// ---------------------------------------------------------------------------
// subroutine calcul_rho_ux_uy_t
// ---------------------------------------------------------------------------
void calcul_rho_ux_uy_t() {
    int nbelm = list_cell.nbelm;

    for (int i = 1; i <= nbelm; ++i) {
        rho(i) = vect_u(i, 1);
        ux(i)  = vect_u(i, 2) / rho(i);
        uy(i)  = vect_u(i, 3) / rho(i);
        t(i)   = 2.0 / (3.0 * r_gaz * rho(i)) * (vect_u(i, 4) - 0.5 * rho(i) * (ux(i) * ux(i) + uy(i) * uy(i)));
    }
}


// ---------------------------------------------------------------------------
// subroutine write_pressure_coefficient(iter)   ! not accurate
// ---------------------------------------------------------------------------
void write_pressure_coefficient(int iter) {
    double rhoref = rho_init;
    double velref = std::sqrt(ux_init * ux_init + uy_init * uy_init);
    double temref = t_init;

    double vinf = velref;
    double pinf = rhoref * r_gaz * temref;

    char citer[8];
    std::snprintf(citer, sizeof(citer), "%07d", iter);

    std::string foutput = std::string("Cp_") + fname + "_" + citer + ".txt";

    fortran_open_unit(20, foutput, "w");
    std::fprintf(fortran_get_unit(20), " X, Y, -Cp\n");

    for (int ifac = 1; ifac <= nbfaces; ++ifac) {
        if (bc_typ(ifac) == 1) { // solid wall
            int left_cell = lr_cell(ifac, 1);
            cell_2D* pcel = list_cell.cell(left_cell).p;

            double ploc = rho(left_cell) * r_gaz * t(left_cell);
            double cp   = 2.0 * (ploc - pinf) / (pinf * velref * velref);

            for (int ilocfac = 1; ilocfac <= 4; ++ilocfac) {
                face* pfac = &pcel->faces[ilocfac - 1];
                if (pfac->bc_typ == 1) {
                    std::fprintf(fortran_get_unit(20), " %g %g %g\n",
                        pfac->centroid.x, pfac->centroid.y, -cp);
                }
            }
        }
    }

    fortran_close_unit(20);
}

// ---------------------------------------------------------------------------
// subroutine chk_converge(iter)
// ---------------------------------------------------------------------------
void chk_converge(int iter) {
    const double tol = 1.0e-9;

    double dr   = 0.0;
    double drho = 0.0;

    for (int i = 1; i <= list_cell.nbelm; ++i) {
        dr   = vect_unew(i, 1) - vect_u(i, 1);
        drho = drho + dr * dr;
    }

    drho = std::sqrt(drho);

    bool lexist = fortran_file_exists("residual.txt");

    if (!lexist) {
        fortran_open_unit(20, "residual.txt", "w");
        std::fprintf(fortran_get_unit(20), "%20s %30s\n", "Iteration", "Residual of density");
        std::fprintf(fortran_get_unit(20), "%20d %30.20f\n", iter, drho);
        fortran_close_unit(20);
    } else if (lexist) {
        fortran_open_unit(20, "residual.txt", "a");
        std::fprintf(fortran_get_unit(20), "%20d %30.20f\n", iter, drho);
        fortran_close_unit(20);
    }

    if (drho <= tol && iter != 1) {
        std::printf("Converged after %9d iterations.\n", iter);
        solver::write_solution_vtk(iter);
        write_pressure_coefficient(iter);
        std::exit(0);
    }
}

} // namespace mod_solver_kfvs
