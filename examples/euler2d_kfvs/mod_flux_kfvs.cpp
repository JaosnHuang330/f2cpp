#include "mod_flux_kfvs.h"

namespace mod_flux_kfvs {

// function fluxp(rho,ux,uy,e,p,t,a,b,nx,ny)
std::array<double,4> fluxp(double rho, double ux, double uy, double e,
                            double p, double t, double a, double b,
                            double nx, double ny)
{
    std::array<double,4> result;

    double un, ut, f1, f2, f3, f4;

    un = nx*ux + ny*uy;
    ut = -ny*ux + nx*uy;

    if (un >= b) {

        f1 = rho*un;
        f2 = rho*un*un + p;
        f3 = ut*f1;
        f4 = (e + p)*un;

    } else if (un <= -b) {

        f1 = 0.0;
        f2 = 0.0;
        f3 = 0.0;
        f4 = 0.0;

    } else {

        f1 = 2.0*a*b*b*(un + b)*(un + b);
        f2 = 4.0/3.0*a*b*b*(un + b)*(un + b)*(un + b);
        f3 = ut*f1;
        f4 = a*b*b/2.0*(un + b)*(un + b)*(un*un + 2.0*un*b + 2*ut*ut + 7.0/3.0*b*b);

    }

    result[0] = f1;
    result[1] = nx*f2 - ny*f3;
    result[2] = ny*f2 + nx*f3;
    result[3] = f4;

    return result;
}

// function fluxm(rho,ux,uy,e,p,t,a,b,nx,ny)
std::array<double,4> fluxm(double rho, double ux, double uy, double e,
                            double p, double t, double a, double b,
                            double nx, double ny)
{
    std::array<double,4> result;

    // fluxm = -fluxp(rho,ux,uy,e,p,t,a,b,-nx,-ny)
    std::array<double,4> fp = fluxp(rho, ux, uy, e, p, t, a, b, -nx, -ny);
    result[0] = -fp[0];
    result[1] = -fp[1];
    result[2] = -fp[2];
    result[3] = -fp[3];

    return result;
}

} // namespace mod_flux_kfvs
