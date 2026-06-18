#pragma once

#include "runtime/fortran_array.h"
#include <string>

namespace m_io {

//---------------------------------------------------------------------------
// Write solution profile to text file
//---------------------------------------------------------------------------
void write_profile(const FortranArray2D<double>& Uin, double time,
                   const std::string& filename);

//---------------------------------------------------------------------------
// Write binary restart file
//---------------------------------------------------------------------------
void write_restart(const FortranArray2D<double>& Uin, double time, int nstep);

//---------------------------------------------------------------------------
// Print step summary to stdout
//---------------------------------------------------------------------------
void print_step(int nstep, double time, double dt, double rho_norm);

} // namespace m_io
