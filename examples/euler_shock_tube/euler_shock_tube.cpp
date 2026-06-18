//===============================================================================
// euler_shock_tube.cpp - Main driver for 1D Euler shock tube solver
//===============================================================================
#include <cstdio>
#include "m_params.h"
#include "m_mesh.h"
#include "m_solver.h"
#include "m_io.h"

using namespace m_params;
using namespace m_mesh;
using namespace m_solver;
using namespace m_io;

int main(int argc, char** argv)
{
    double rho_norm;

    std::printf("============================================\n");
    std::printf("  1D Euler Shock Tube Solver (HLLC)\n");
    std::printf("============================================\n");

    // Initialize mesh
    init_mesh();

    // Initialize solution (Sod problem)
    init_sod();

    // Print initial state
    rho_norm = density_l2norm();
    print_step(0, 0.0, 0.0, rho_norm);

    // Time integration loop
    while (time_current < t_final && nstep < max_steps) {
        compute_dt();
        advance();

        rho_norm = density_l2norm();

        if (nstep % output_freq == 0) {
            print_step(nstep, time_current, dt, rho_norm);
        }
    }

    // Final output
    print_step(nstep, time_current, dt, rho_norm);
    write_profile(Qc, time_current, output_file);
    write_restart(Qc, time_current, nstep);

    std::printf("============================================\n");
    std::printf(" Done: %6d steps\n", nstep);
    std::printf(" Final time: %10.6f\n", time_current);
    std::printf("============================================\n");

    return 0;
}
