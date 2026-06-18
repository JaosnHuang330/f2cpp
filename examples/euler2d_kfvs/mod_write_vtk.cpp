#include "mod_write_vtk.h"
#include <cstdio>
#include <string>

namespace mod_write_vtk {

using mod_read_gmsh::fname;
using mod_cell_2d::obj_cell_2D;

// ---------------------------------------------------------------------------
// Type-bound procedure dispatchers
// ---------------------------------------------------------------------------

void obj_msh_vis::write_vtk()      { write_mesh_vtk(*this); }
void obj_msh_vis::write_tecplot()  { write_mesh_tecplot(*this); }

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

obj_msh_vis obj_msh_vis_constructor(obj_cell_2D& obj_cell_2D_container) {
    obj_msh_vis obj_msh_vis_constructor_result;
    obj_msh_vis_constructor_result.smart_pointer = &obj_cell_2D_container;
    return obj_msh_vis_constructor_result;
}

// ---------------------------------------------------------------------------
// write_mesh_vtk
// ---------------------------------------------------------------------------

void write_mesh_vtk(obj_msh_vis& this_) {
    int i;
    std::string foutput;

    foutput = std::string(fname) + ".vtk";

    fortran_open_unit(21, foutput, "w");

    fortran_write_unit(21, "%s\n", "# vtk DataFile Version 2.0");
    fortran_write_unit(21, "%s\n", "VTK format for unstructured mesh");
    fortran_write_unit(21, "%s\n", "ASCII");
    fortran_write_unit(21, "%s\n", "DATASET POLYDATA");
    // format 1: 'POINTS' i9 ' float'
    fortran_write_unit(21, "POINTS%9d float\n", this_.smart_pointer->nbnode);

    for (i = 1; i <= this_.smart_pointer->nbnode; i++) {
        fortran_write_unit(21, "%g %g %g\n",
            this_.smart_pointer->coord_nodes(i).p->x,
            this_.smart_pointer->coord_nodes(i).p->y,
            this_.smart_pointer->coord_nodes(i).p->z);
    }

    // format 2: 'POLYGONS ' 2i9
    fortran_write_unit(21, "POLYGONS %9d%9d\n",
        this_.smart_pointer->nbelm,
        5 * this_.smart_pointer->nbelm);

    for (i = 1; i <= this_.smart_pointer->nbelm; i++) {
        // write(21,*) 4, id_node(6:9)-1
        fortran_write_unit(21, "%d %d %d %d %d\n",
            4,
            this_.smart_pointer->id_nodes(i).pn->id_node(6) - 1,
            this_.smart_pointer->id_nodes(i).pn->id_node(7) - 1,
            this_.smart_pointer->id_nodes(i).pn->id_node(8) - 1,
            this_.smart_pointer->id_nodes(i).pn->id_node(9) - 1);
    }

    // format 3: 'CELL_DATA' i9
    fortran_write_unit(21, "CELL_DATA%9d\n", this_.smart_pointer->nbelm);
    fortran_write_unit(21, "%s\n", "SCALARS CELL_IDENT integer 1");
    fortran_write_unit(21, "%s\n", "LOOKUP_TABLE default ");

    for (i = 1; i <= this_.smart_pointer->nbelm; i++) {
        fortran_write_unit(21, "%d\n",
            this_.smart_pointer->id_nodes(i).pn->id_node(1));
    }

    fortran_write_unit(21, "%s\n", "SCALARS NEIGHBOR1 integer 1");
    fortran_write_unit(21, "%s\n", "LOOKUP_TABLE default ");

    for (i = 1; i <= this_.smart_pointer->nbelm; i++) {
        if (this_.smart_pointer->cell(i).p->neighbor1 != nullptr) {
            fortran_write_unit(21, "%d\n",
                this_.smart_pointer->cell(i).p->neighbor1->ident);
        } else {
            fortran_write_unit(21, "%d\n", 0);
        }
    }

    fortran_write_unit(21, "%s\n", "SCALARS NEIGHBOR2 integer 1");
    fortran_write_unit(21, "%s\n", "LOOKUP_TABLE default ");

    for (i = 1; i <= this_.smart_pointer->nbelm; i++) {
        if (this_.smart_pointer->cell(i).p->neighbor2 != nullptr) {
            fortran_write_unit(21, "%d\n",
                this_.smart_pointer->cell(i).p->neighbor2->ident);
        } else {
            fortran_write_unit(21, "%d\n", 0);
        }
    }

    fortran_write_unit(21, "%s\n", "SCALARS NEIGHBOR3 integer 1");
    fortran_write_unit(21, "%s\n", "LOOKUP_TABLE default ");

    for (i = 1; i <= this_.smart_pointer->nbelm; i++) {
        if (this_.smart_pointer->cell(i).p->neighbor3 != nullptr) {
            fortran_write_unit(21, "%d\n",
                this_.smart_pointer->cell(i).p->neighbor3->ident);
        } else {
            fortran_write_unit(21, "%d\n", 0);
        }
    }

    fortran_write_unit(21, "%s\n", "SCALARS NEIGHBOR4 integer 1");
    fortran_write_unit(21, "%s\n", "LOOKUP_TABLE default ");

    for (i = 1; i <= this_.smart_pointer->nbelm; i++) {
        if (this_.smart_pointer->cell(i).p->neighbor4 != nullptr) {
            fortran_write_unit(21, "%d\n",
                this_.smart_pointer->cell(i).p->neighbor4->ident);
        } else {
            fortran_write_unit(21, "%d\n", 0);
        }
    }

    fortran_close_unit(21);
}

// ---------------------------------------------------------------------------
// write_mesh_tecplot
// ---------------------------------------------------------------------------

void write_mesh_tecplot(obj_msh_vis& this_) {
    int i;
    std::string foutput;

    foutput = std::string(fname) + ".dat";

    fortran_open_unit(21, foutput, "w");

    // write(21,'(a)') 'VARIABLES=X,Y'
    fortran_write_unit(21, "%s\n", "VARIABLES=X,Y");
    fortran_write_unit(21, "%s\n", "ZONE T=\"UNSTRUCTURED-COUNTOUR\"");
    fortran_write_unit(21, "%s\n", "ZONETYPE=FEPOLYGON");
    // write(21,*) 'NODES=', nbnode  (list-directed with label)
    fortran_write_unit(21, " NODES= %d\n", this_.smart_pointer->nbnode);
    fortran_write_unit(21, " ELEMENTS= %d\n", this_.smart_pointer->nbelm);
    fortran_write_unit(21, " FACES= %d\n", this_.smart_pointer->nbelm * 4);
    fortran_write_unit(21, " NumConnectedBoundaryFaces= %d\n", 0);
    fortran_write_unit(21, " TotalNumBoundaryConnections= %d\n", 0);

    for (i = 1; i <= this_.smart_pointer->nbnode; i++) {
        fortran_write_unit(21, "%g\n",
            this_.smart_pointer->coord_nodes(i).p->x);
    }

    for (i = 1; i <= this_.smart_pointer->nbnode; i++) {
        fortran_write_unit(21, "%g\n",
            this_.smart_pointer->coord_nodes(i).p->y);
    }

    // Node indexes: 4 faces per element, each face has 2 node indices
    for (i = 1; i <= this_.smart_pointer->nbelm; i++) {
        fortran_write_unit(21, "%d %d\n",
            this_.smart_pointer->id_nodes(i).pn->id_node(6),
            this_.smart_pointer->id_nodes(i).pn->id_node(7));
        fortran_write_unit(21, "%d %d\n",
            this_.smart_pointer->id_nodes(i).pn->id_node(7),
            this_.smart_pointer->id_nodes(i).pn->id_node(8));
        fortran_write_unit(21, "%d %d\n",
            this_.smart_pointer->id_nodes(i).pn->id_node(8),
            this_.smart_pointer->id_nodes(i).pn->id_node(9));
        fortran_write_unit(21, "%d %d\n",
            this_.smart_pointer->id_nodes(i).pn->id_node(9),
            this_.smart_pointer->id_nodes(i).pn->id_node(6));
    }

    for (i = 1; i <= this_.smart_pointer->nbelm; i++) {
        fortran_write_unit(21, "%d %d %d %d\n", i, i, i, i);
    }

    for (i = 1; i <= this_.smart_pointer->nbelm; i++) {
        fortran_write_unit(21, "%d %d %d %d\n", 0, 0, 0, 0);
    }

    fortran_close_unit(21);
}

} // namespace mod_write_vtk
