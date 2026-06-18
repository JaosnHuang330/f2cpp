#pragma once

namespace m_riemann {

//---------------------------------------------------------------------------
// HLL approximate Riemann solver
// Returns the numerical flux at a cell interface given left/right states
//---------------------------------------------------------------------------
void hll_flux(double rhoL, double uL, double pL,
              double rhoR, double uR, double pR,
              double* flux);

//---------------------------------------------------------------------------
// HLLC Riemann solver (more accurate contact resolution)
//---------------------------------------------------------------------------
void hllc_flux(double rhoL, double uL, double pL,
               double rhoR, double uR, double pR,
               double* flux);

//---------------------------------------------------------------------------
// Convert conservative to primitive variables
//---------------------------------------------------------------------------
void cons_to_prim(const double* Uc, double& rho, double& vel, double& p);

//---------------------------------------------------------------------------
// Convert primitive to conservative variables
//---------------------------------------------------------------------------
void prim_to_cons(double rho, double vel, double p, double* Uc);

} // namespace m_riemann
