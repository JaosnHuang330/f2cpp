#include "mod_cell_2d.h"
#include "mod_read_gmsh.h"
#include "mod_point_util.h"
#include "mod_vector_algebra.h"
#include <cmath>
#include <cstdlib>

namespace mod_cell_2d {

using mod_read_gmsh::point;
using mod_read_gmsh::obj_gmsh_reader;
using mod_read_gmsh::node_ident_msh;
using mod_vector_algebra::vector_t;
using mod_vector_algebra::vector_norm;
using mod_point_util::equal;

// ---------------------------------------------------------------------------
// Module-level variable definition
// ---------------------------------------------------------------------------
int nbfaces = 0;

// ---------------------------------------------------------------------------
// Type-bound procedure dispatch methods
// ---------------------------------------------------------------------------
void obj_cell_2D::build_cells()      { construct_cells(*this); }
void obj_cell_2D::calcul_vol()       { calcul_vol_cells(*this); }
void obj_cell_2D::calcul_face_info() { calcul_area_cent_faces(*this); }
void obj_cell_2D::assign_face_id()   { assign_id_face_v2(*this); }
void obj_cell_2D::find_neighbor()    { detect_neighbor(*this); }


// ---------------------------------------------------------------------------
// subroutine construct_cells
// ---------------------------------------------------------------------------
void construct_cells(obj_cell_2D& this_) {
    int i, j;
    int idnode;
    int idnode1, idnode2;
    cell_2D*       pcell = nullptr;
    node_ident_msh* pmsh = nullptr;

    this_.cell.allocate(this_.nbelm);

    for (i = 1; i <= this_.nbelm; i++) {
        this_.cell(i).p = new cell_2D();
    }

    // Assign vertexes of a cell
    for (i = 1; i <= this_.nbelm; i++) {
        this_.cell(i).p->ident = i;

        // vertex 1
        this_.cell(i).p->vertex[0].ident = this_.id_nodes(i).pn->id_node(6);
        idnode = this_.cell(i).p->vertex[0].ident;
        for (j = 1; j <= this_.nbnode; j++) {
            if (idnode == this_.coord_nodes(j).p->ident) {
                this_.cell(i).p->vertex[0].x = this_.coord_nodes(j).p->x;
                this_.cell(i).p->vertex[0].y = this_.coord_nodes(j).p->y;
            }
        }

        // vertex 2
        this_.cell(i).p->vertex[1].ident = this_.id_nodes(i).pn->id_node(7);
        idnode = this_.cell(i).p->vertex[1].ident;
        for (j = 1; j <= this_.nbnode; j++) {
            if (idnode == this_.coord_nodes(j).p->ident) {
                this_.cell(i).p->vertex[1].x = this_.coord_nodes(j).p->x;
                this_.cell(i).p->vertex[1].y = this_.coord_nodes(j).p->y;
            }
        }

        // vertex 3
        this_.cell(i).p->vertex[2].ident = this_.id_nodes(i).pn->id_node(8);
        idnode = this_.cell(i).p->vertex[2].ident;
        for (j = 1; j <= this_.nbnode; j++) {
            if (idnode == this_.coord_nodes(j).p->ident) {
                this_.cell(i).p->vertex[2].x = this_.coord_nodes(j).p->x;
                this_.cell(i).p->vertex[2].y = this_.coord_nodes(j).p->y;
            }
        }

        // vertex 4
        this_.cell(i).p->vertex[3].ident = this_.id_nodes(i).pn->id_node(9);
        idnode = this_.cell(i).p->vertex[3].ident;
        for (j = 1; j <= this_.nbnode; j++) {
            if (idnode == this_.coord_nodes(j).p->ident) {
                this_.cell(i).p->vertex[3].x = this_.coord_nodes(j).p->x;
                this_.cell(i).p->vertex[3].y = this_.coord_nodes(j).p->y;
            }
        }
    }

    // Sort vertexes for each cell
    // sort_vertex failed (commented out in original)
    // do i = 1, nbelm
    //     call sort_vertex(cell(i)%p%vertex)
    // enddo

    // Assign faces for each cell
    for (i = 1; i <= this_.nbelm; i++) {
        this_.cell(i).p->faces[0].p1.x = this_.cell(i).p->vertex[0].x;
        this_.cell(i).p->faces[0].p1.y = this_.cell(i).p->vertex[0].y;
        this_.cell(i).p->faces[0].p2.x = this_.cell(i).p->vertex[1].x;
        this_.cell(i).p->faces[0].p2.y = this_.cell(i).p->vertex[1].y;

        this_.cell(i).p->faces[1].p1.x = this_.cell(i).p->vertex[1].x;
        this_.cell(i).p->faces[1].p1.y = this_.cell(i).p->vertex[1].y;
        this_.cell(i).p->faces[1].p2.x = this_.cell(i).p->vertex[2].x;
        this_.cell(i).p->faces[1].p2.y = this_.cell(i).p->vertex[2].y;

        this_.cell(i).p->faces[2].p1.x = this_.cell(i).p->vertex[2].x;
        this_.cell(i).p->faces[2].p1.y = this_.cell(i).p->vertex[2].y;
        this_.cell(i).p->faces[2].p2.x = this_.cell(i).p->vertex[3].x;
        this_.cell(i).p->faces[2].p2.y = this_.cell(i).p->vertex[3].y;

        this_.cell(i).p->faces[3].p1.x = this_.cell(i).p->vertex[3].x;
        this_.cell(i).p->faces[3].p1.y = this_.cell(i).p->vertex[3].y;
        this_.cell(i).p->faces[3].p2.x = this_.cell(i).p->vertex[0].x;
        this_.cell(i).p->faces[3].p2.y = this_.cell(i).p->vertex[0].y;
    }

    // Assign boundary condition for each face of a cell
    for (i = 1; i <= this_.nbelm; i++) {
        pcell = this_.cell(i).p;

        // face 1: vertex(1) - vertex(2)
        idnode1 = pcell->vertex[0].ident;
        idnode2 = pcell->vertex[1].ident;
        for (j = 1; j <= this_.nbel_msh; j++) {
            pmsh = this_.id_nodes_msh(j).pn;
            if (pmsh->elem_typ == 27) {
                if (idnode1 == pmsh->id_node(6) && idnode2 == pmsh->id_node(7)) {
                    pcell->faces[0].bc_typ = pmsh->tag1;
                    continue;
                } else if (idnode1 == pmsh->id_node(7) && idnode2 == pmsh->id_node(6)) {
                    pcell->faces[0].bc_typ = pmsh->tag1;
                    continue;
                }
            }
        }

        // face 2: vertex(2) - vertex(3)
        idnode1 = pcell->vertex[1].ident;
        idnode2 = pcell->vertex[2].ident;
        for (j = 1; j <= this_.nbel_msh; j++) {
            pmsh = this_.id_nodes_msh(j).pn;
            if (pmsh->elem_typ == 27) {
                if (idnode1 == pmsh->id_node(6) && idnode2 == pmsh->id_node(7)) {
                    pcell->faces[1].bc_typ = pmsh->tag1;
                    continue;
                } else if (idnode1 == pmsh->id_node(7) && idnode2 == pmsh->id_node(6)) {
                    pcell->faces[1].bc_typ = pmsh->tag1;
                    continue;
                }
            }
        }

        // face 3: vertex(3) - vertex(4)
        idnode1 = pcell->vertex[2].ident;
        idnode2 = pcell->vertex[3].ident;
        for (j = 1; j <= this_.nbel_msh; j++) {
            pmsh = this_.id_nodes_msh(j).pn;
            if (pmsh->elem_typ == 27) {
                if (idnode1 == pmsh->id_node(6) && idnode2 == pmsh->id_node(7)) {
                    pcell->faces[2].bc_typ = pmsh->tag1;
                    continue;
                } else if (idnode1 == pmsh->id_node(7) && idnode2 == pmsh->id_node(6)) {
                    pcell->faces[2].bc_typ = pmsh->tag1;
                    continue;
                }
            }
        }

        // face 4: vertex(4) - vertex(1)
        idnode1 = pcell->vertex[3].ident;
        idnode2 = pcell->vertex[0].ident;
        for (j = 1; j <= this_.nbel_msh; j++) {
            pmsh = this_.id_nodes_msh(j).pn;
            if (pmsh->elem_typ == 27) {
                if (idnode1 == pmsh->id_node(6) && idnode2 == pmsh->id_node(7)) {
                    pcell->faces[3].bc_typ = pmsh->tag1;
                    continue;
                } else if (idnode1 == pmsh->id_node(7) && idnode2 == pmsh->id_node(6)) {
                    pcell->faces[3].bc_typ = pmsh->tag1;
                    continue;
                }
            }
        }
    }
}


// ---------------------------------------------------------------------------
// subroutine calcul_vol_cells
// ---------------------------------------------------------------------------
void calcul_vol_cells(obj_cell_2D& this_) {
    int i;
    cell_2D* pc = nullptr;
    double vol, x1, x2, x3, x4, y1, y2, y3, y4;

    for (i = 1; i <= this_.nbelm; i++) {
        pc = this_.cell(i).p;

        x1 = pc->vertex[0].x;
        x2 = pc->vertex[1].x;
        x3 = pc->vertex[2].x;
        x4 = pc->vertex[3].x;

        y1 = pc->vertex[0].y;
        y2 = pc->vertex[1].y;
        y3 = pc->vertex[2].y;
        y4 = pc->vertex[3].y;

        vol = 0.5 * ( (x1 - x3)*(y2 - y4) + (x4 - x2)*(y1 - y3) );
        this_.cell(i).p->vol = vol;
    }
}

// ---------------------------------------------------------------------------
// subroutine calcul_area_cent_faces
// ---------------------------------------------------------------------------
void calcul_area_cent_faces(obj_cell_2D& this_) {
    int i, j;
    cell_2D* pcel = nullptr;
    face*    pfac = nullptr;
    vector_t vec;

    for (i = 1; i <= this_.nbelm; i++) {
        pcel = this_.cell(i).p;
        for (j = 0; j < 4; j++) {
            pfac = &pcel->faces[j];
            vec.x = pfac->p1.x - pfac->p2.x;
            vec.y = pfac->p1.y - pfac->p2.y;
            vec.z = 0.0;
            this_.cell(i).p->faces[j].area       = vector_norm(vec);
            this_.cell(i).p->faces[j].centroid.x = 0.5 * (pfac->p1.x + pfac->p2.x);
            this_.cell(i).p->faces[j].centroid.y = 0.5 * (pfac->p1.y + pfac->p2.y);
        }
    }
}


// ---------------------------------------------------------------------------
// Internal helper: findminimum (internal function of sort_vertex)
// ---------------------------------------------------------------------------
static int findminimum_impl(point* vertices, int start, int endp) {
    int location;
    int i;
    point minimum;

    minimum  = vertices[start - 1];  // Fortran 1-based → C++ 0-based
    location = start;
    for (i = start + 1; i <= endp; i++) {
        if (vertices[i - 1].x < minimum.x) {
            minimum  = vertices[i - 1];
            location = i;
        }
    }
    return location;
}

// ---------------------------------------------------------------------------
// subroutine sort_vertex (private)
// ---------------------------------------------------------------------------
void sort_vertex(point* vertices, int n) {
    // Reference: http://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90
    int i, location;
    point vertices_temp[4];

    for (i = 1; i <= n - 1; i++) {
        location = findminimum_impl(vertices, i, n);
        // swap_point(vertices(i), vertices(location))
        mod_point_util::swap_point(vertices[i - 1], vertices[location - 1]);
    }

    // Copy to temp
    for (i = 0; i < 4; i++) {
        vertices_temp[i] = vertices[i];
    }

    if (vertices_temp[0].y < vertices_temp[1].y) {
        vertices[0] = vertices_temp[0];
        vertices[3] = vertices_temp[1];
    } else {
        // vertices_temp(1)%y >= vertices_temp(2)%y
        vertices[0] = vertices_temp[1];
        vertices[3] = vertices_temp[0];
    }

    if (vertices_temp[2].y < vertices_temp[3].y) {
        vertices[1] = vertices_temp[2];
        vertices[2] = vertices_temp[3];
    } else {
        // vertices_temp(3)%y >= vertices_temp(4)%y
        vertices[1] = vertices_temp[3];
        vertices[2] = vertices_temp[2];
    }
}


// ---------------------------------------------------------------------------
// subroutine assign_id_face_v2
// ---------------------------------------------------------------------------
void assign_id_face_v2(obj_cell_2D& this_) {
    using mod_point_util::equal;
    int i, j, k;
    cell_2D* pc  = nullptr;
    cell_2D* pc2 = nullptr;
    face*    pfac  = nullptr;
    face*    pfac2 = nullptr;

    nbfaces = 0;
    for (i = 1; i <= this_.nbelm; i++) {
        pc = this_.cell(i).p;
        for (j = 0; j < 4; j++) {
            pfac = &pc->faces[j];
            if (pfac->idface == 0) { // a face that has not been indexed yet

                // Check neighbor1
                pc2 = pc->neighbor1;
                if (pc->neighbor1 != nullptr) {
                    for (k = 0; k < 4; k++) {
                        pfac2 = &pc2->faces[k];
                        if (equal(pfac->centroid, pfac2->centroid) &&
                            std::abs(pfac->area - pfac2->area) <= 1.0e-9) {
                            if (pfac2->idface == 0) {
                                nbfaces++;
                                this_.cell(i).p->faces[j].idface           = nbfaces;
                                this_.cell(i).p->neighbor1->faces[k].idface = nbfaces;
                            }
                        }
                    }
                } else {
                    // .not. associated(pc%neighbor1)
                    if (pfac->bc_typ != 0) { // boundary face
                        if (pfac->idface > 0) continue;
                        nbfaces++;
                        this_.cell(i).p->faces[j].idface = nbfaces;
                    }
                }

                // Check neighbor2
                pc2 = pc->neighbor2;
                if (pc->neighbor2 != nullptr) {
                    for (k = 0; k < 4; k++) {
                        pfac2 = &pc2->faces[k];
                        if (equal(pfac->centroid, pfac2->centroid) &&
                            std::abs(pfac->area - pfac2->area) <= 1.0e-9) {
                            if (pfac2->idface == 0) {
                                nbfaces++;
                                this_.cell(i).p->faces[j].idface           = nbfaces;
                                this_.cell(i).p->neighbor2->faces[k].idface = nbfaces;
                            }
                        }
                    }
                } else {
                    // .not. associated(pc%neighbor2)
                    if (pfac->bc_typ != 0) { // boundary face
                        if (pfac->idface > 0) continue;
                        nbfaces++;
                        this_.cell(i).p->faces[j].idface = nbfaces;
                    }
                }

                // Check neighbor3
                pc2 = pc->neighbor3;
                if (pc->neighbor3 != nullptr) {
                    for (k = 0; k < 4; k++) {
                        pfac2 = &pc2->faces[k];
                        if (equal(pfac->centroid, pfac2->centroid) &&
                            std::abs(pfac->area - pfac2->area) <= 1.0e-9) {
                            if (pfac2->idface == 0) {
                                nbfaces++;
                                this_.cell(i).p->faces[j].idface           = nbfaces;
                                this_.cell(i).p->neighbor3->faces[k].idface = nbfaces;
                            }
                        }
                    }
                } else {
                    // .not. associated(pc%neighbor3)
                    if (pfac->bc_typ != 0) { // boundary face
                        if (pfac->idface > 0) continue;
                        nbfaces++;
                        this_.cell(i).p->faces[j].idface = nbfaces;
                    }
                }

                // Check neighbor4
                pc2 = pc->neighbor4;
                if (pc->neighbor4 != nullptr) {
                    for (k = 0; k < 4; k++) {
                        pfac2 = &pc2->faces[k];
                        if (equal(pfac->centroid, pfac2->centroid) &&
                            std::abs(pfac->area - pfac2->area) <= 1.0e-9) {
                            if (pfac2->idface == 0) {
                                nbfaces++;
                                this_.cell(i).p->faces[j].idface           = nbfaces;
                                this_.cell(i).p->neighbor4->faces[k].idface = nbfaces;
                            }
                        }
                    }
                } else {
                    // .not. associated(pc%neighbor4)
                    if (pfac->bc_typ != 0) { // boundary face
                        if (pfac->idface > 0) continue;
                        nbfaces++;
                        this_.cell(i).p->faces[j].idface = nbfaces;
                    }
                }

            } // if (pfac->idface == 0)
        } // for j
    } // for i
}


// ---------------------------------------------------------------------------
// subroutine detect_neighbor
// ---------------------------------------------------------------------------
void detect_neighbor(obj_cell_2D& this_) {
    int i, j, k, t;
    int idnode1, idnode2;
    int cnt;
    cell_2D* pc = nullptr;
    cell_2D* pr = nullptr;

    // Find neighbor1: shares vertex(1) and vertex(2)
    for (i = 1; i <= this_.nbelm; i++) {
        idnode1 = this_.cell(i).p->vertex[0].ident;
        idnode2 = this_.cell(i).p->vertex[1].ident;
        pc = this_.cell(i).p;

        for (j = 1; j <= this_.nbelm; j++) {
            pr = this_.cell(j).p;
            cnt = 0;
            if (pc == pr) continue;
            for (t = 0; t < 4; t++) {
                if (idnode1 == pr->vertex[t].ident) cnt++;
                if (idnode2 == pr->vertex[t].ident) cnt++;
            }
            if (cnt == 2) {
                this_.cell(i).p->neighbor1 = pr;
            }
        }
    }

    // Find neighbor2: shares vertex(2) and vertex(3)
    for (i = 1; i <= this_.nbelm; i++) {
        idnode1 = this_.cell(i).p->vertex[1].ident;
        idnode2 = this_.cell(i).p->vertex[2].ident;
        pc = this_.cell(i).p;

        for (j = 1; j <= this_.nbelm; j++) {
            pr = this_.cell(j).p;
            cnt = 0;
            if (pc == pr) continue;
            for (t = 0; t < 4; t++) {
                if (idnode1 == pr->vertex[t].ident) cnt++;
                if (idnode2 == pr->vertex[t].ident) cnt++;
            }
            if (cnt == 2) {
                this_.cell(i).p->neighbor2 = pr;
            }
        }
    }

    // Find neighbor3: shares vertex(3) and vertex(4)
    for (i = 1; i <= this_.nbelm; i++) {
        idnode1 = this_.cell(i).p->vertex[2].ident;
        idnode2 = this_.cell(i).p->vertex[3].ident;
        pc = this_.cell(i).p;

        for (j = 1; j <= this_.nbelm; j++) {
            pr = this_.cell(j).p;
            cnt = 0;
            if (pc == pr) continue;
            for (t = 0; t < 4; t++) {
                if (idnode1 == pr->vertex[t].ident) cnt++;
                if (idnode2 == pr->vertex[t].ident) cnt++;
            }
            if (cnt == 2) {
                this_.cell(i).p->neighbor3 = pr;
            }
        }
    }

    // Find neighbor4: shares vertex(4) and vertex(1)
    for (i = 1; i <= this_.nbelm; i++) {
        idnode1 = this_.cell(i).p->vertex[3].ident;
        idnode2 = this_.cell(i).p->vertex[0].ident;
        pc = this_.cell(i).p;

        for (j = 1; j <= this_.nbelm; j++) {
            pr = this_.cell(j).p;
            cnt = 0;
            if (pc == pr) continue;
            for (t = 0; t < 4; t++) {
                if (idnode1 == pr->vertex[t].ident) cnt++;
                if (idnode2 == pr->vertex[t].ident) cnt++;
            }
            if (cnt == 2) {
                this_.cell(i).p->neighbor4 = pr;
            }
        }
    }
}

} // namespace mod_cell_2d
