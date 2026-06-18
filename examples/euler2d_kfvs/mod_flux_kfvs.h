#pragma once

#include <array>

namespace mod_flux_kfvs {

// Returns the positive KFVS flux vector (4 components)
std::array<double,4> fluxp(double rho, double ux, double uy, double e,
                            double p, double t, double a, double b,
                            double nx, double ny);

// Returns the negative KFVS flux vector (4 components)
std::array<double,4> fluxm(double rho, double ux, double uy, double e,
                            double p, double t, double a, double b,
                            double nx, double ny);

} // namespace mod_flux_kfvs
