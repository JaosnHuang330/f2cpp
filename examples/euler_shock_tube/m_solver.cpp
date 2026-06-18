#include "m_solver.h"
#include <cmath>
#include <cstdio>

namespace m_solver {

// Module-level variable definitions
FortranArray2D<double> Qc(3, m_params::nx);
double dt = 0.0;
double time_current = 0.0;
int nstep = 0;

//---------------------------------------------------------------------------
// Initialize with Sod shock tube conditions
//---------------------------------------------------------------------------
void init_sod()
{
    for (int i = 1; i <= m_params::nx; i++) {
        double rho0, u0, p0;
        if (m_mesh::x(i) < 0.5) {
            rho0 = 1.0;
            u0   = 0.0;
            p0   = 1.0;
        } else {
            rho0 = 0.125;
            u0   = 0.0;
            p0   = 0.1;
        }
        double Ucol[3];
        m_riemann::prim_to_cons(rho0, u0, p0, Ucol);
        Qc(1, i) = Ucol[0];
        Qc(2, i) = Ucol[1];
        Qc(3, i) = Ucol[2];
    }

    time_current = 0.0;
    nstep = 0;

    std::printf(" Initialized Sod shock tube problem\n");
    std::printf(" t_final = %8.4f\n", m_params::t_final);
}

//---------------------------------------------------------------------------
// Compute stable time step based on CFL condition
//---------------------------------------------------------------------------
void compute_dt()
{
    double max_speed = 0.0;

    for (int i = 1; i <= m_params::nx; i++) {
        double rho, vel, p;
        double Ucol[3] = { Qc(1, i), Qc(2, i), Qc(3, i) };
        m_riemann::cons_to_prim(Ucol, rho, vel, p);
        double a = std::sqrt(m_params::gamma_gas * p / rho);
        max_speed = std::max(max_speed, std::abs(vel) + a);
    }

    dt = m_params::cfl * m_params::dx / max_speed;

    // Don't overshoot final time
    if (time_current + dt > m_params::t_final) {
        dt = m_params::t_final - time_current;
    }
}

//---------------------------------------------------------------------------
// One time step using forward Euler + HLLC flux
//---------------------------------------------------------------------------
void advance()
{
    FortranArray2D<double> Qnew(3, m_params::nx);

    for (int i = 1; i <= m_params::nx; i++) {
        double rhoL, uL, pL, rhoR, uR, pR;
        double flux_left[3], flux_right[3];

        // Left interface flux
        if (i == 1) {
            // Reflective BC (ghost cell = mirror)
            double UcolL[3] = { Qc(1, 1), Qc(2, 1), Qc(3, 1) };
            m_riemann::cons_to_prim(UcolL, rhoL, uL, pL);
            uL = -uL;  // reflect velocity
            double UcolR[3] = { Qc(1, 1), Qc(2, 1), Qc(3, 1) };
            m_riemann::cons_to_prim(UcolR, rhoR, uR, pR);
        } else {
            double UcolL[3] = { Qc(1, i-1), Qc(2, i-1), Qc(3, i-1) };
            m_riemann::cons_to_prim(UcolL, rhoL, uL, pL);
            double UcolR[3] = { Qc(1, i), Qc(2, i), Qc(3, i) };
            m_riemann::cons_to_prim(UcolR, rhoR, uR, pR);
        }
        m_riemann::hllc_flux(rhoL, uL, pL, rhoR, uR, pR, flux_left);

        // Right interface flux
        if (i == m_params::nx) {
            // Reflective BC
            double UcolL[3] = { Qc(1, m_params::nx), Qc(2, m_params::nx), Qc(3, m_params::nx) };
            m_riemann::cons_to_prim(UcolL, rhoL, uL, pL);
            double UcolR[3] = { Qc(1, m_params::nx), Qc(2, m_params::nx), Qc(3, m_params::nx) };
            m_riemann::cons_to_prim(UcolR, rhoR, uR, pR);
            uR = -uR;
        } else {
            double UcolL[3] = { Qc(1, i), Qc(2, i), Qc(3, i) };
            m_riemann::cons_to_prim(UcolL, rhoL, uL, pL);
            double UcolR[3] = { Qc(1, i+1), Qc(2, i+1), Qc(3, i+1) };
            m_riemann::cons_to_prim(UcolR, rhoR, uR, pR);
        }
        m_riemann::hllc_flux(rhoL, uL, pL, rhoR, uR, pR, flux_right);

        // Update conservative variables
        double factor = dt / m_params::dx;
        Qnew(1, i) = Qc(1, i) - factor * (flux_right[0] - flux_left[0]);
        Qnew(2, i) = Qc(2, i) - factor * (flux_right[1] - flux_left[1]);
        Qnew(3, i) = Qc(3, i) - factor * (flux_right[2] - flux_left[2]);
    }

    // Copy Qnew back to Qc
    for (int i = 1; i <= m_params::nx; i++) {
        Qc(1, i) = Qnew(1, i);
        Qc(2, i) = Qnew(2, i);
        Qc(3, i) = Qnew(3, i);
    }

    time_current = time_current + dt;
    nstep = nstep + 1;
}

//---------------------------------------------------------------------------
// Compute L2 norm of density (for diagnostics)
//---------------------------------------------------------------------------
double density_l2norm()
{
    double norm = 0.0;

    for (int i = 1; i <= m_params::nx; i++) {
        double rho, vel, p;
        double Ucol[3] = { Qc(1, i), Qc(2, i), Qc(3, i) };
        m_riemann::cons_to_prim(Ucol, rho, vel, p);
        norm = norm + rho * rho * m_params::dx;
    }

    norm = std::sqrt(norm);
    return norm;
}

} // namespace m_solver
