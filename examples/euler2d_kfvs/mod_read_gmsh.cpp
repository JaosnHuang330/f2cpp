#include "mod_read_gmsh.h"

#include <cstdio>
#include <cstring>
#include <string>

namespace mod_read_gmsh {

// ---------------------------------------------------------------------------
// Module-level variable definitions
// ---------------------------------------------------------------------------
std::string fname;
double rho_init = 0.0;
double ux_init  = 0.0;
double uy_init  = 0.0;
double t_init   = 0.0;

// ---------------------------------------------------------------------------
// Type-bound method dispatch
// ---------------------------------------------------------------------------
void obj_gmsh_reader::read_msh()        { read_mesh(*this); }
void obj_gmsh_reader::assign_id_nodes() { construct_id_nodes(*this); }


// ---------------------------------------------------------------------------
// subroutine read_mesh(this)
// ---------------------------------------------------------------------------
void read_mesh(obj_gmsh_reader& this_) {
    // Namelist /mesh_file/ fname
    // Namelist /init_valu/ rho_init, ux_init, uy_init, t_init

    std::string finput;
    std::string cline;
    int i;
    int elem_typ;

    // open(unit=100, file='input.dat')
    fortran_open_unit(100, "input.dat", "r");
    {
        NamelistGroup nml_mesh_file("mesh_file");
        nml_mesh_file.add("fname", fname);
        nml_mesh_file.read(100);
    }
    {
        NamelistGroup nml_init_valu("init_valu");
        nml_init_valu.add("rho_init", rho_init);
        nml_init_valu.add("ux_init",  ux_init);
        nml_init_valu.add("uy_init",  uy_init);
        nml_init_valu.add("t_init",   t_init);
        nml_init_valu.read(100);
    }
    fortran_close_unit(100);

    // finput = trim(fname)//'.msh'
    finput = fname + ".msh";

    // write(6,'(a)') 'Input mesh file: ',finput
    fprintf(fortran_get_unit(6), "Input mesh file: %s\n", finput.c_str());

    // open(unit=20, file=finput, status='old')
    fortran_open_unit(20, finput, "r");

    // read(20,*) x4 — skip 4 lines
    fortran_read_list(20);
    fortran_read_list(20);
    fortran_read_list(20);
    fortran_read_list(20);

    // read(20,*) this%nbnode
    fortran_read_list(20, &this_.nbnode);

    // allocate(this%coord_nodes(1:this%nbnode))
    this_.coord_nodes.allocate(this_.nbnode);

    // do i=1,nbnode; allocate(this%coord_nodes(i)%p); end do
    for (i = 1; i <= this_.nbnode; i++) {
        this_.coord_nodes(i).p = new point();
    }

    // do i=1,nbnode; read(20,*) ident, x, y, z; end do
    for (i = 1; i <= this_.nbnode; i++) {
        fortran_read_list(20,
            &this_.coord_nodes(i).p->ident,
            &this_.coord_nodes(i).p->x,
            &this_.coord_nodes(i).p->y,
            &this_.coord_nodes(i).p->z);
    }

    // read(20,*) — skip 2 lines
    fortran_read_list(20);
    fortran_read_list(20);

    // read(20,*) this%nbel_msh
    fortran_read_list(20, &this_.nbel_msh);

    // allocate(this%id_nodes_msh(1:this%nbel_msh))
    this_.id_nodes_msh.allocate(this_.nbel_msh);

    // do i=1,nbel_msh; allocate(this%id_nodes_msh(i)%pn); end do
    for (i = 1; i <= this_.nbel_msh; i++) {
        this_.id_nodes_msh(i).pn = new node_ident_msh();
    }

    this_.nbelm = 0;

    // do i=1,nbel_msh; read header fields + select case on elem_typ
    for (i = 1; i <= this_.nbel_msh; i++) {
        fortran_read_list(20,
            &this_.id_nodes_msh(i).pn->ident,
            &this_.id_nodes_msh(i).pn->elem_typ,
            &this_.id_nodes_msh(i).pn->nb_tags,
            &this_.id_nodes_msh(i).pn->tag1,
            &this_.id_nodes_msh(i).pn->tag2);

        elem_typ = this_.id_nodes_msh(i).pn->elem_typ;

        switch (elem_typ) {
            case 1:  // 2-node line
                this_.id_nodes_msh(i).pn->id_node.allocate(7);
                break;
            case 3:  // 4-node quadrangle
                this_.id_nodes_msh(i).pn->id_node.allocate(9);
                this_.nbelm = this_.nbelm + 1;
                break;
            case 15: // 1-node point
                this_.id_nodes_msh(i).pn->id_node.allocate(6);
                break;
            case 37: // 5-node edge quadrangle
                this_.id_nodes_msh(i).pn->id_node.allocate(30);
                this_.nbelm = this_.nbelm + 1;
                break;
            case 27: // boundary 5-node edge
                this_.id_nodes_msh(i).pn->id_node.allocate(10);
                break;
            default:
                printf("comming soon !\n");
                std::exit(1);
                break;
        }
    }

    // close(unit=20)
    fortran_close_unit(20);

    // Reopen and scan for "$Elements"
    // open(unit=20, file=finput, status='old')
    fortran_open_unit(20, finput, "r");

    // do; read(20,*) cline; cline=trim(cline); if (cline=='$Elements') exit; enddo
    do {
        fortran_read_list(20, &cline);
        // trim is implicit in std::string from fortran_read_list
        // strip trailing whitespace just in case
        while (!cline.empty() && (cline.back() == ' ' || cline.back() == '\r' || cline.back() == '\n')) {
            cline.pop_back();
        }
    } while (cline != "$Elements");

    // read(20,*) — skip the element count line
    fortran_read_list(20);

    // do i=1,nbel_msh; read(20,*) id_node(:); end do
    for (i = 1; i <= this_.nbel_msh; i++) {
        int sz = this_.id_nodes_msh(i).pn->id_node.size();
        fortran_read_list_array(20, this_.id_nodes_msh(i).pn->id_node, 1, sz);
    }

    // close(unit=20)
    fortran_close_unit(20);
}



// ---------------------------------------------------------------------------
// subroutine construct_id_nodes(this)
// ---------------------------------------------------------------------------
void construct_id_nodes(obj_gmsh_reader& this_) {
    int i, j;

    // allocate(this%id_nodes(1:this%nbelm))
    this_.id_nodes.allocate(this_.nbelm);

    // do i=1,nbelm; allocate pn; allocate id_node based on elem_typ
    for (i = 1; i <= this_.nbelm; i++) {
        this_.id_nodes(i).pn = new node_ident();

        if (this_.id_nodes_msh(i).pn->elem_typ == 3) {
            this_.id_nodes(i).pn->id_node.allocate(9);
        }

        if (this_.id_nodes_msh(i).pn->elem_typ == 37) {
            this_.id_nodes(i).pn->id_node.allocate(30);
        }
    }

    // j=0; loop over nbel_msh, copy id_node for elem_typ 3 and 37
    j = 0;
    for (i = 1; i <= this_.nbel_msh; i++) {
        if (this_.id_nodes_msh(i).pn->elem_typ == 3) {
            j = j + 1;
            // this%id_nodes(j)%pn%id_node(:) = this%id_nodes_msh(i)%pn%id_node(:)
            int sz = this_.id_nodes_msh(i).pn->id_node.size();
            for (int k = 1; k <= sz; k++) {
                this_.id_nodes(j).pn->id_node(k) = this_.id_nodes_msh(i).pn->id_node(k);
            }
        }

        if (this_.id_nodes_msh(i).pn->elem_typ == 37) {
            j = j + 1;
            // this%id_nodes(j)%pn%id_node(:) = this%id_nodes_msh(i)%pn%id_node(:)
            int sz = this_.id_nodes_msh(i).pn->id_node.size();
            for (int k = 1; k <= sz; k++) {
                this_.id_nodes(j).pn->id_node(k) = this_.id_nodes_msh(i).pn->id_node(k);
            }
        }
    }
}


} // namespace mod_read_gmsh
