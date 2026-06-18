#pragma once

#include "runtime/fortran_array.h"
#include "m_params.h"
#include "m_mesh.h"
#include "m_riemann.h"

namespace m_solver {

// Conservative variables: Qc(1:3, 1:nx)
// Qc(1,:) = rho, Qc(2,:) = rho*u, Qc(3,:) = E
extern FortranArray2D<double> Qc;
extern double dt;
extern double time_current;
extern int nstep;

//---------------------------------------------------------------------------
// Initialize with Sod shock tube conditions
// Left:  rho=1.0, u=0.0, p=1.0    (x < 0.5)
// Right: rho=0.125, u=0.0, p=0.1  (x >= 0.5)
//---------------------------------------------------------------------------
void init_sod();

//---------------------------------------------------------------------------
// Compute stable time step based on CFL condition
//---------------------------------------------------------------------------
void compute_dt();

//---------------------------------------------------------------------------
// One time step using forward Euler + HLLC flux
//---------------------------------------------------------------------------
void advance();

//---------------------------------------------------------------------------
// Compute L2 norm of density (for diagnostics)
//---------------------------------------------------------------------------
double density_l2norm();

} // namespace m_solver
