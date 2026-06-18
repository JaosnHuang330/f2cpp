#include "pre_processing.h"

#include "mod_objects.h"
#include "mod_write_vtk.h"
#include "mod_fvm_face_2d.h"
#include "runtime/fortran_io.h"

#include <chrono>

namespace pre_processing {

using namespace mod_objects;
using mod_write_vtk::obj_msh_vis_constructor;
using mod_fvm_face_2d::obj_fvm_face_2D_constructor;

void pre_processing() {
    float startTime, endTime;

    // call cpu_time(startTime)
    auto t_start = std::chrono::high_resolution_clock::now();
    startTime = static_cast<float>(
        std::chrono::duration<double>(t_start.time_since_epoch()).count());

    list_cell.read_msh();
    list_cell.assign_id_nodes();

    list_cell.build_cells();
    list_cell.find_neighbor();

    visualizer = obj_msh_vis_constructor(list_cell);
    visualizer.write_vtk();

    // call cpu_time(endTime)
    auto t_end = std::chrono::high_resolution_clock::now();
    endTime = static_cast<float>(
        std::chrono::duration<double>(t_end.time_since_epoch()).count());

    fortran_write_unit(6, " Time for treating mesh: %g\n", endTime - startTime);

    visualizer.write_tecplot();

    list_cell.calcul_vol();
    list_cell.calcul_face_info();
    list_cell.assign_face_id();

    faces_fvm = obj_fvm_face_2D_constructor(list_cell);
    faces_fvm.alloc();
    faces_fvm.get();
}

} // namespace pre_processing
