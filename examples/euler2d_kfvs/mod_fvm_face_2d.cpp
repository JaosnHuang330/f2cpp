#include "mod_fvm_face_2d.h"

namespace mod_fvm_face_2d {

using mod_vector_algebra::vector_norm;
using mod_vector_algebra::vector_times_scalar;

// ---------------------------------------------------------------------------
// obj_fvm_face_2D constructor (from interface obj_fvm_face_2D)
// ---------------------------------------------------------------------------
obj_fvm_face_2D::obj_fvm_face_2D(obj_cell_2D& obj_cell_2D_container) {
    smart_pointer = &obj_cell_2D_container;
}

obj_fvm_face_2D obj_fvm_face_2D_constructor(obj_cell_2D& obj_cell_2D_container) {
    obj_fvm_face_2D result;
    result.smart_pointer = &obj_cell_2D_container;
    return result;
}

// ---------------------------------------------------------------------------
// allocate_face_2D (type-bound: alloc)
// ---------------------------------------------------------------------------
void obj_fvm_face_2D::alloc() { allocate_face_2D(*this); }

void allocate_face_2D(obj_fvm_face_2D& this_) {
    int i;

    if (!this_.face_2D.allocated()) {
        this_.face_2D.allocate(nbfaces);
    }

    for (i = 1; i <= nbfaces; i++) {
        this_.face_2D(i).f = new fvm_face_2D();
    }
}

// ---------------------------------------------------------------------------
// assign_face_2D (type-bound: get)
// ---------------------------------------------------------------------------
void obj_fvm_face_2D::get() { assign_face_2D(*this); }

void assign_face_2D(obj_fvm_face_2D& this_) {
    int       icel, ifac, jfac, jfac_nbor;
    cell_2D*  pcel      = nullptr;
    cell_2D*  pcel_nbor = nullptr;
    fvm_face_2D* pfac   = nullptr;
    face*     pfac2     = nullptr;
    face*     pfac2_nbor = nullptr;
    double    norm_vec;

    for (ifac = 1; ifac <= nbfaces; ifac++) {

        pfac = this_.face_2D(ifac).f;
        for (icel = 1; icel <= this_.smart_pointer->nbelm; icel++) {
            pcel = this_.smart_pointer->cell(icel).p;
            for (jfac = 1; jfac <= 4; jfac++) {
                pfac2 = &pcel->faces[jfac - 1];

                if (pfac2->idface == ifac) {

                    if (this_.face_2D(ifac).f->idface > 0) continue; // just to index face_2D one time

                    // Copy data from face to fvm_face_2D
                    this_.face_2D(ifac).f->p1       = this_.smart_pointer->cell(icel).p->faces[jfac - 1].p1;
                    this_.face_2D(ifac).f->p2       = this_.smart_pointer->cell(icel).p->faces[jfac - 1].p2;
                    this_.face_2D(ifac).f->bc_typ   = this_.smart_pointer->cell(icel).p->faces[jfac - 1].bc_typ;
                    this_.face_2D(ifac).f->idface   = this_.smart_pointer->cell(icel).p->faces[jfac - 1].idface;
                    this_.face_2D(ifac).f->area     = this_.smart_pointer->cell(icel).p->faces[jfac - 1].area;
                    this_.face_2D(ifac).f->centroid = this_.smart_pointer->cell(icel).p->faces[jfac - 1].centroid;

                    // Assign left cell and right cell for each face_2D
                    pcel_nbor = pcel->neighbor1;
                    if (pcel_nbor != nullptr) {
                        for (jfac_nbor = 1; jfac_nbor <= 4; jfac_nbor++) {
                            pfac2_nbor = &pcel_nbor->faces[jfac_nbor - 1];
                            if (pfac2_nbor->idface == ifac) {
                                if (pcel->ident <= pcel_nbor->ident) {
                                    this_.face_2D(ifac).f->left_cell  = this_.smart_pointer->cell(icel).p;
                                    this_.face_2D(ifac).f->right_cell = this_.smart_pointer->cell(icel).p->neighbor1;
                                } else if (pcel->ident > pcel_nbor->ident) {
                                    this_.face_2D(ifac).f->left_cell  = this_.smart_pointer->cell(icel).p->neighbor1;
                                    this_.face_2D(ifac).f->right_cell = this_.smart_pointer->cell(icel).p;
                                }
                            }
                        }
                    } else if (pcel_nbor == nullptr) {
                        if (this_.face_2D(ifac).f->left_cell == nullptr) {
                            this_.face_2D(ifac).f->left_cell = this_.smart_pointer->cell(icel).p;
                        }
                    }

                    pcel_nbor = pcel->neighbor2;
                    if (pcel_nbor != nullptr) {
                        for (jfac_nbor = 1; jfac_nbor <= 4; jfac_nbor++) {
                            pfac2_nbor = &pcel_nbor->faces[jfac_nbor - 1];
                            if (pfac2_nbor->idface == ifac) {
                                if (pcel->ident <= pcel_nbor->ident) {
                                    this_.face_2D(ifac).f->left_cell  = this_.smart_pointer->cell(icel).p;
                                    this_.face_2D(ifac).f->right_cell = this_.smart_pointer->cell(icel).p->neighbor2;
                                } else if (pcel->ident > pcel_nbor->ident) {
                                    this_.face_2D(ifac).f->left_cell  = this_.smart_pointer->cell(icel).p->neighbor2;
                                    this_.face_2D(ifac).f->right_cell = this_.smart_pointer->cell(icel).p;
                                }
                            }
                        }
                    } else if (pcel_nbor == nullptr) {
                        if (this_.face_2D(ifac).f->left_cell == nullptr) {
                            this_.face_2D(ifac).f->left_cell = this_.smart_pointer->cell(icel).p;
                        }
                    }

                    pcel_nbor = pcel->neighbor3;
                    if (pcel_nbor != nullptr) {
                        for (jfac_nbor = 1; jfac_nbor <= 4; jfac_nbor++) {
                            pfac2_nbor = &pcel_nbor->faces[jfac_nbor - 1];
                            if (pfac2_nbor->idface == ifac) {
                                if (pcel->ident <= pcel_nbor->ident) {
                                    this_.face_2D(ifac).f->left_cell  = this_.smart_pointer->cell(icel).p;
                                    this_.face_2D(ifac).f->right_cell = this_.smart_pointer->cell(icel).p->neighbor3;
                                } else if (pcel->ident > pcel_nbor->ident) {
                                    this_.face_2D(ifac).f->left_cell  = this_.smart_pointer->cell(icel).p->neighbor3;
                                    this_.face_2D(ifac).f->right_cell = this_.smart_pointer->cell(icel).p;
                                }
                            }
                        }
                    } else if (pcel_nbor == nullptr) {
                        if (this_.face_2D(ifac).f->left_cell == nullptr) {
                            this_.face_2D(ifac).f->left_cell = this_.smart_pointer->cell(icel).p;
                        }
                    }

                    pcel_nbor = pcel->neighbor4;
                    if (pcel_nbor != nullptr) {
                        for (jfac_nbor = 1; jfac_nbor <= 4; jfac_nbor++) {
                            pfac2_nbor = &pcel_nbor->faces[jfac_nbor - 1];
                            if (pfac2_nbor->idface == ifac) {
                                if (pcel->ident <= pcel_nbor->ident) {
                                    this_.face_2D(ifac).f->left_cell  = this_.smart_pointer->cell(icel).p;
                                    this_.face_2D(ifac).f->right_cell = this_.smart_pointer->cell(icel).p->neighbor4;
                                } else if (pcel->ident > pcel_nbor->ident) {
                                    this_.face_2D(ifac).f->left_cell  = this_.smart_pointer->cell(icel).p->neighbor4;
                                    this_.face_2D(ifac).f->right_cell = this_.smart_pointer->cell(icel).p;
                                }
                            }
                        }
                    } else if (pcel_nbor == nullptr) {
                        if (this_.face_2D(ifac).f->left_cell == nullptr) {
                            this_.face_2D(ifac).f->left_cell = this_.smart_pointer->cell(icel).p;
                        }
                    }

                    // Assign tangent and normal vectors
                    this_.face_2D(ifac).f->tangent_vector = vector_t{
                        this_.face_2D(ifac).f->p2.x - this_.face_2D(ifac).f->p1.x,
                        this_.face_2D(ifac).f->p2.y - this_.face_2D(ifac).f->p1.y,
                        0.0
                    };
                    norm_vec = vector_norm(this_.face_2D(ifac).f->tangent_vector);
                    this_.face_2D(ifac).f->len_tan        = norm_vec;
                    this_.face_2D(ifac).f->tangent_vector = vector_times_scalar(this_.face_2D(ifac).f->tangent_vector, 1.0 / norm_vec);

                    if (this_.face_2D(ifac).f->right_cell == this_.smart_pointer->cell(icel).p) {
                        this_.face_2D(ifac).f->normal_vector = vector_t{
                            this_.face_2D(ifac).f->p1.y - this_.face_2D(ifac).f->p2.y,
                            this_.face_2D(ifac).f->p2.x - this_.face_2D(ifac).f->p1.x,
                            0.0
                        };
                    } else if (this_.face_2D(ifac).f->left_cell == this_.smart_pointer->cell(icel).p) {
                        this_.face_2D(ifac).f->normal_vector = vector_t{
                            -this_.face_2D(ifac).f->p1.y + this_.face_2D(ifac).f->p2.y,
                            -this_.face_2D(ifac).f->p2.x + this_.face_2D(ifac).f->p1.x,
                            0.0
                        };
                    }
                    norm_vec = vector_norm(this_.face_2D(ifac).f->normal_vector);
                    this_.face_2D(ifac).f->len_nor        = norm_vec;
                    this_.face_2D(ifac).f->normal_vector  = vector_times_scalar(this_.face_2D(ifac).f->normal_vector, 1.0 / norm_vec);
                }
            }
        }
    }
}

} // namespace mod_fvm_face_2d
