#include "m_riemann.h"
#include "m_params.h"
#include <cmath>

namespace m_riemann {

using m_params::gamma_gas;

//---------------------------------------------------------------------------
// HLL approximate Riemann solver
//---------------------------------------------------------------------------
void hll_flux(double rhoL, double uL, double pL,
              double rhoR, double uR, double pR,
              double* flux)
{
    double EL, ER, aL, aR, SL, SR;
    double fL[3], fR[3], UUL[3], UUR[3];
    double g1;

    g1 = gamma_gas - 1.0;

    // Sound speeds
    aL = std::sqrt(gamma_gas * pL / rhoL);
    aR = std::sqrt(gamma_gas * pR / rhoR);

    // Total energy
    EL = pL / g1 + 0.5 * rhoL * uL * uL;
    ER = pR / g1 + 0.5 * rhoR * uR * uR;

    // Wave speed estimates (Davis)
    SL = std::min(uL - aL, uR - aR);
    SR = std::max(uL + aL, uR + aR);

    // Conservative variables
    UUL[0] = rhoL;
    UUL[1] = rhoL * uL;
    UUL[2] = EL;

    UUR[0] = rhoR;
    UUR[1] = rhoR * uR;
    UUR[2] = ER;

    // Physical fluxes
    fL[0] = rhoL * uL;
    fL[1] = rhoL * uL * uL + pL;
    fL[2] = (EL + pL) * uL;

    fR[0] = rhoR * uR;
    fR[1] = rhoR * uR * uR + pR;
    fR[2] = (ER + pR) * uR;

    // HLL flux
    if (SL >= 0.0) {
        flux[0] = fL[0];
        flux[1] = fL[1];
        flux[2] = fL[2];
    } else if (SR <= 0.0) {
        flux[0] = fR[0];
        flux[1] = fR[1];
        flux[2] = fR[2];
    } else {
        double denom = 1.0 / (SR - SL);
        flux[0] = (SR * fL[0] - SL * fR[0] + SL * SR * (UUR[0] - UUL[0])) * denom;
        flux[1] = (SR * fL[1] - SL * fR[1] + SL * SR * (UUR[1] - UUL[1])) * denom;
        flux[2] = (SR * fL[2] - SL * fR[2] + SL * SR * (UUR[2] - UUL[2])) * denom;
    }
}

//---------------------------------------------------------------------------
// HLLC Riemann solver (more accurate contact resolution)
//---------------------------------------------------------------------------
void hllc_flux(double rhoL, double uL, double pL,
               double rhoR, double uR, double pR,
               double* flux)
{
    double EL, ER, aL, aR, SL, SR, Sstar;
    double fL[3], fR[3], UUL[3], UUR[3];
    double UstarL[3], UstarR[3];
    double g1, pstar, rhobar, abar;
    double qL, qR;

    g1 = gamma_gas - 1.0;

    // Sound speeds
    aL = std::sqrt(gamma_gas * pL / rhoL);
    aR = std::sqrt(gamma_gas * pR / rhoR);

    // Total energy
    EL = pL / g1 + 0.5 * rhoL * uL * uL;
    ER = pR / g1 + 0.5 * rhoR * uR * uR;

    // Pressure estimate (PVRS)
    rhobar = 0.5 * (rhoL + rhoR);
    abar   = 0.5 * (aL + aR);
    pstar  = 0.5 * (pL + pR) - 0.5 * (uR - uL) * rhobar * abar;
    pstar  = std::max(pstar, 0.0);

    // Wave speed estimates
    if (pstar <= pL) {
        qL = 1.0;
    } else {
        qL = std::sqrt(1.0 + (gamma_gas + 1.0) / (2.0 * gamma_gas) * (pstar / pL - 1.0));
    }
    if (pstar <= pR) {
        qR = 1.0;
    } else {
        qR = std::sqrt(1.0 + (gamma_gas + 1.0) / (2.0 * gamma_gas) * (pstar / pR - 1.0));
    }

    SL = uL - aL * qL;
    SR = uR + aR * qR;
    Sstar = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) /
            (rhoL * (SL - uL) - rhoR * (SR - uR));

    // Conservative variables
    UUL[0] = rhoL;
    UUL[1] = rhoL * uL;
    UUL[2] = EL;

    UUR[0] = rhoR;
    UUR[1] = rhoR * uR;
    UUR[2] = ER;

    // Physical fluxes
    fL[0] = rhoL * uL;
    fL[1] = rhoL * uL * uL + pL;
    fL[2] = (EL + pL) * uL;

    fR[0] = rhoR * uR;
    fR[1] = rhoR * uR * uR + pR;
    fR[2] = (ER + pR) * uR;

    // HLLC flux
    if (SL >= 0.0) {
        flux[0] = fL[0];
        flux[1] = fL[1];
        flux[2] = fL[2];
    } else if (Sstar >= 0.0) {
        // Star left state
        UstarL[0] = rhoL * (SL - uL) / (SL - Sstar);
        UstarL[1] = UstarL[0] * Sstar;
        UstarL[2] = UstarL[0] * (EL / rhoL + (Sstar - uL) *
                    (Sstar + pL / (rhoL * (SL - uL))));
        flux[0] = fL[0] + SL * (UstarL[0] - UUL[0]);
        flux[1] = fL[1] + SL * (UstarL[1] - UUL[1]);
        flux[2] = fL[2] + SL * (UstarL[2] - UUL[2]);
    } else if (SR > 0.0) {
        // Star right state
        UstarR[0] = rhoR * (SR - uR) / (SR - Sstar);
        UstarR[1] = UstarR[0] * Sstar;
        UstarR[2] = UstarR[0] * (ER / rhoR + (Sstar - uR) *
                    (Sstar + pR / (rhoR * (SR - uR))));
        flux[0] = fR[0] + SR * (UstarR[0] - UUR[0]);
        flux[1] = fR[1] + SR * (UstarR[1] - UUR[1]);
        flux[2] = fR[2] + SR * (UstarR[2] - UUR[2]);
    } else {
        flux[0] = fR[0];
        flux[1] = fR[1];
        flux[2] = fR[2];
    }
}

//---------------------------------------------------------------------------
// Convert conservative to primitive variables
//---------------------------------------------------------------------------
void cons_to_prim(const double* Uc, double& rho, double& vel, double& p)
{
    rho = Uc[0];
    vel = Uc[1] / rho;
    p   = (gamma_gas - 1.0) * (Uc[2] - 0.5 * rho * vel * vel);
}

//---------------------------------------------------------------------------
// Convert primitive to conservative variables
//---------------------------------------------------------------------------
void prim_to_cons(double rho, double vel, double p, double* Uc)
{
    Uc[0] = rho;
    Uc[1] = rho * vel;
    Uc[2] = p / (gamma_gas - 1.0) + 0.5 * rho * vel * vel;
}

} // namespace m_riemann
