//  euler2d.cpp
//
//  PROGRAM: euler2d
//
//  PURPOSE:  Illustrate the using of gmsh-to-vtk-and-tecplot-Fortran2003.
//

#include "pre_processing.h"
#include "solver.h"

int main(int argc, char** argv) {

    pre_processing::pre_processing(); // refactoring code of gmsh-to-vtk-and-tecplot

    solver::solver();

    return 0;
}
