#include "solver.h"
#include <cmath>
#include <cstdio>
#include <string>

namespace solver {

// ---------------------------------------------------------------------------
// subroutine solver
// ---------------------------------------------------------------------------
void solver() {
    solver_kfvs();
}

// ---------------------------------------------------------------------------
// subroutine solver_kfvs
// ---------------------------------------------------------------------------
void solver_kfvs() {
    using namespace mod_solver_kfvs;
    using namespace mod_struct_to_array;
    using mod_read_gmsh::fname;

    int n;

    struct_to_array();

    donnee_initiale();

    if (fname.substr(0, 3) == "c31") {
        allocate_vardummy_multi_elem_airfoil();
        conditions_aux_limites_multi_elem_airfoil();
        assign_lr_cell_multi_elem_airfoil();
        calcul_derived_quantities_multi_elem_airfoil();
    } else {
        allocate_vardummy();
        conditions_aux_limites();
        assign_lr_cell();
        calcul_derived_quantities();
    }

    calcul_conservative_vector();

    timestep();

    //--- temps de simulation et parametres de sorties
    tmax = 1.0; //0.3828823925d-2
    nmax = 10000; // reduced for testing (was: floor(tmax/dt))

    //--- evolution
    for (n = 1; n <= nmax; ++n) {

        //-- calcul des flux
        if (fname.substr(0, 3) == "c31") {
            calcul_flux_multi_elem_airfoil();
        } else {
            calcul_flux();
        }

        //-- iteration en temps
        calcul_rhs();
        euler_time_iteration();

        //-- verifier si converge
        chk_converge(n);

        //-- mise a jour
        for (int i = 1; i <= vect_u.size(1); ++i)
            for (int k = 1; k <= vect_u.size(2); ++k)
                vect_u(i, k) = vect_unew(i, k);

        //-- calcul de rho,ux,uy,t
        calcul_rho_ux_uy_t();

        //-- mise a jour cl
        if (fname.substr(0, 3) == "c31") {
            conditions_aux_limites_multi_elem_airfoil();
        } else {
            conditions_aux_limites();
        }

        //-- mise a jour des quantites derivees
        if (fname.substr(0, 3) == "c31") {
            calcul_derived_quantities_multi_elem_airfoil();
        } else {
            calcul_derived_quantities();
        }

        //-- sauvegarde resultats (format vtk, lisible par Paraview)
        if ((n % 10000) == 0) {
            printf(" Writing solution file at iteration %d ...\n", n);
            write_solution_vtk(n);
            write_pressure_coefficient(n);
        }
    }
}

// ---------------------------------------------------------------------------
// subroutine write_solution_vtk(iter)
// ---------------------------------------------------------------------------
void write_solution_vtk(int iter) {
    using namespace mod_solver_kfvs;
    using mod_read_gmsh::fname;

    int i;
    char foutput[300];
    char citer[8];

    snprintf(citer, sizeof(citer), "%07d", iter);
    // foutput = trim(fname)//'_'//trim(citer)//'.vtk'
    std::string foutput_str = fname + "_" + std::string(citer) + ".vtk";
    snprintf(foutput, sizeof(foutput), "%s", foutput_str.c_str());

    fortran_open_unit(21, foutput, "w");

    fprintf(fortran_get_unit(21), "# vtk DataFile Version 2.0\n");
    fprintf(fortran_get_unit(21), "VTK format for unstructured mesh\n");
    fprintf(fortran_get_unit(21), "ASCII\n");
    fprintf(fortran_get_unit(21), "DATASET POLYDATA\n");
    // format 1: 'POINTS',i9,' float'
    fprintf(fortran_get_unit(21), "POINTS%9d float\n", list_cell.nbnode);

    for (i = 1; i <= list_cell.nbnode; ++i) {
        fprintf(fortran_get_unit(21), " %g %g %g\n",
                list_cell.coord_nodes(i).p->x,
                list_cell.coord_nodes(i).p->y,
                list_cell.coord_nodes(i).p->z);
    }

    // format 2: 'POLYGONS ',2i9
    fprintf(fortran_get_unit(21), "POLYGONS %9d%9d\n",
            list_cell.nbelm, 5 * list_cell.nbelm);

    for (i = 1; i <= list_cell.nbelm; ++i) {
        // write(21,*) 4, list_cell%id_nodes(i)%pn%id_node(6:9)-1
        fprintf(fortran_get_unit(21), " 4 %d %d %d %d\n",
                list_cell.id_nodes(i).pn->id_node(6) - 1,
                list_cell.id_nodes(i).pn->id_node(7) - 1,
                list_cell.id_nodes(i).pn->id_node(8) - 1,
                list_cell.id_nodes(i).pn->id_node(9) - 1);
    }

    // format 3: 'CELL_DATA',i9
    fprintf(fortran_get_unit(21), "CELL_DATA%9d\n", list_cell.nbelm);
    fprintf(fortran_get_unit(21), "SCALARS CELL_IDENT integer 1\n");
    fprintf(fortran_get_unit(21), "LOOKUP_TABLE default \n");

    for (i = 1; i <= list_cell.nbelm; ++i) {
        fprintf(fortran_get_unit(21), " %d\n",
                list_cell.id_nodes(i).pn->id_node(1));
    }

    fprintf(fortran_get_unit(21), "SCALARS Density float\n");
    fprintf(fortran_get_unit(21), "LOOKUP_TABLE default \n");
    for (i = 1; i <= list_cell.nbelm; ++i) {
        fprintf(fortran_get_unit(21), " %g\n", rho(i));
    }

    fprintf(fortran_get_unit(21), "SCALARS Temperature float\n");
    fprintf(fortran_get_unit(21), "LOOKUP_TABLE default \n");
    for (i = 1; i <= list_cell.nbelm; ++i) {
        fprintf(fortran_get_unit(21), " %g\n", t(i));
    }

    fprintf(fortran_get_unit(21), "SCALARS Pressure float\n");
    fprintf(fortran_get_unit(21), "LOOKUP_TABLE default \n");
    for (i = 1; i <= list_cell.nbelm; ++i) {
        fprintf(fortran_get_unit(21), " %g\n", p(i));
    }

    fprintf(fortran_get_unit(21), "SCALARS Velocity_u float\n");
    fprintf(fortran_get_unit(21), "LOOKUP_TABLE default \n");
    for (i = 1; i <= list_cell.nbelm; ++i) {
        fprintf(fortran_get_unit(21), " %g\n", ux(i));
    }

    fprintf(fortran_get_unit(21), "SCALARS Velocity_v float\n");
    fprintf(fortran_get_unit(21), "LOOKUP_TABLE default \n");
    for (i = 1; i <= list_cell.nbelm; ++i) {
        fprintf(fortran_get_unit(21), " %g\n", uy(i));
    }

    fprintf(fortran_get_unit(21), "SCALARS Velocity_magnitude float\n");
    fprintf(fortran_get_unit(21), "LOOKUP_TABLE default \n");
    for (i = 1; i <= list_cell.nbelm; ++i) {
        fprintf(fortran_get_unit(21), " %g\n",
                std::sqrt(ux(i) * ux(i) + uy(i) * uy(i)));
    }

    fortran_close_unit(21);
}

} // namespace solver
