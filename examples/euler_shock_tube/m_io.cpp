#include "m_io.h"
#include "m_params.h"
#include "m_mesh.h"
#include "m_riemann.h"
#include "runtime/fortran_io.h"

#include <cmath>
#include <cstdio>

namespace m_io {

//---------------------------------------------------------------------------
// Write solution profile to text file
//---------------------------------------------------------------------------
void write_profile(const FortranArray2D<double>& Uin, double time,
                   const std::string& filename)
{
    using namespace m_params;
    using namespace m_mesh;
    using namespace m_riemann;

    int iunit = 20;
    fortran_open_unit(iunit, filename, "w");

    fortran_write_unit(iunit, "# Euler 1D solution at t = %10.6f\n", time);
    fortran_write_unit(iunit, "# nx = %6d\n", nx);
    fortran_write_unit(iunit, "# x  rho  u  p  e_int  Mach\n");

    for (int i = 1; i <= nx; ++i) {
        double rho, vel, p;
        cons_to_prim(&Uin(1, i), rho, vel, p);

        double e_int = p / ((gamma_gas - 1.0) * rho);
        double a = std::sqrt(gamma_gas * p / rho);
        double mach = std::abs(vel) / a;

        fortran_write_unit(iunit, "%16.8E %16.8E %16.8E %16.8E %16.8E %16.8E\n",
                           x(i), rho, vel, p, e_int, mach);
    }

    fortran_close_unit(iunit);
    printf(" Solution written to: %s\n", filename.c_str());
}

//---------------------------------------------------------------------------
// Write binary restart file
//---------------------------------------------------------------------------
void write_restart(const FortranArray2D<double>& Uin, double time, int nstep)
{
    using namespace m_params;

    int iunit = 21;
    fortran_open_unit(iunit, restart_file, "wb");

    // Write header: nx, time, nstep (raw binary)
    FILE* fp = fortran_get_unit(iunit);
    std::fwrite(&nx, sizeof(int), 1, fp);
    std::fwrite(&time, sizeof(double), 1, fp);
    std::fwrite(&nstep, sizeof(int), 1, fp);

    // Write the full Uin array (3 x nx doubles, raw binary)
    std::fwrite(&Uin(1, 1), sizeof(double), 3 * nx, fp);

    fortran_close_unit(iunit);
}

//---------------------------------------------------------------------------
// Print step summary to stdout
//---------------------------------------------------------------------------
void print_step(int nstep, double time, double dt, double rho_norm)
{
    printf(" Step %6d  t=%10.6f  dt=%10.3e  ||rho||=%12.5e\n",
           nstep, time, dt, rho_norm);
}

} // namespace m_io
