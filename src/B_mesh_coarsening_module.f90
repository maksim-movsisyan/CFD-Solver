module mesh_coarsening_module
use mesh_module
use mg_coarse_mesh_building_module
use mg_aglomeration_module
use output_fields_module
implicit none
!=======================================================================
!=============== AGGLOMERATION INFO DATA STRUCTURE =====================
!=======================================================================
type :: agglomeration_info_t
	integer, allocatable :: fine_to_coarse(:)							!=== MAP: FINE_IDX -> COARSE_IDX, DIM = (fine_ncells) ===
	integer, allocatable :: fine_cells_ptr(:), fine_cells(:)			!=== FINE CELLS OF COARSE CELL, DIM = (coarse_ncells+1) / DIM = (sum(nfine_per_coarse)) ===
	integer :: coarse_ncells

	contains
	procedure :: agglomerate => perform_mesh_agglomeration
	procedure :: get_coarse => get_coarse_mesh
end type



contains
!=======================================================================
!===================== AGGLOMERATION INFO METHODS ======================
!=======================================================================
	subroutine perform_mesh_agglomeration(this, fine_mesh, target_fine_per_coarse)
		implicit none
		class(agglomeration_info_t), intent(inout) :: this
		type(mesh_t), intent(in) :: fine_mesh
		integer, intent(in) :: target_fine_per_coarse
		
		call mg_perform_agglomeration(fine_mesh%dim, fine_mesh%ncells,&
									  fine_mesh%nfaces, fine_mesh%nbfaces,&
									  fine_mesh%face_left_cell, fine_mesh%face_right_cell,&
									  fine_mesh%cell_faces_ptr, fine_mesh%cell_faces,&
									  fine_mesh%face_area, fine_mesh%cell_center,&
									  target_fine_per_coarse,&
									  this%coarse_ncells, this%fine_to_coarse,&
									  this%fine_cells_ptr, this%fine_cells) 
	end subroutine

	subroutine get_coarse_mesh(this, fine_mesh, coarse_mesh, output_mesh, filename)
		implicit none
		class(agglomeration_info_t), intent(in) :: this
		logical, intent(in) :: output_mesh
		type(mesh_t), intent(in) :: fine_mesh
		type(mesh_t), intent(inout) :: coarse_mesh
		character(len=*), optional, intent(in) :: filename
		
		character(len=256) :: filename2
		
		coarse_mesh%dim = fine_mesh%dim
		coarse_mesh%ncells = this%coarse_ncells
		
		call mg_build_coarse_mesh(fine_mesh%dim, fine_mesh%ncells,&
								  fine_mesh%nfaces, fine_mesh%nbfaces, fine_mesh%nnodes,&
								  fine_mesh%face_left_cell, fine_mesh%face_right_cell,&
								  fine_mesh%face_zone, fine_mesh%face_type,&
								  fine_mesh%face_nodes_ptr, fine_mesh%face_nodes,&
								  fine_mesh%face_center, fine_mesh%face_normal, fine_mesh%face_area,&
								  fine_mesh%node_coords, fine_mesh%cell_volume, fine_mesh%cell_center,&
								  this%fine_to_coarse, this%fine_cells_ptr, this%fine_cells,&
								  coarse_mesh%ncells, coarse_mesh%nfaces,&
								  coarse_mesh%nbfaces, coarse_mesh%nnodes,&
								  coarse_mesh%face_center, coarse_mesh%face_normal, coarse_mesh%face_area, coarse_mesh%face_weight,&
								  coarse_mesh%cell_volume, coarse_mesh%cell_center, coarse_mesh%node_coords,&
								  coarse_mesh%face_zone, coarse_mesh%face_type, coarse_mesh%face_bidx,&
								  coarse_mesh%face_left_cell, coarse_mesh%face_right_cell,&
								  coarse_mesh%face_nodes_ptr, coarse_mesh%face_nodes,&
								  coarse_mesh%cell_nodes_ptr, coarse_mesh%cell_nodes,&
								  coarse_mesh%cell_faces_ptr, coarse_mesh%cell_faces, coarse_mesh%cell_type)
								  
		coarse_mesh%zone_names = fine_mesh%zone_names
		coarse_mesh%zone_subnames = fine_mesh%zone_subnames
		
		filename2 = 'agglomerated_mesh.vtk'
		if (present(filename)) filename2 = trim(filename)
		
		if (output_mesh) then
			if (coarse_mesh%dim == 2) then
				call write_mesh_faces_vtk(trim(adjustl(filename2)), coarse_mesh%dim,&
										  coarse_mesh%nnodes, coarse_mesh%node_coords,&
										  coarse_mesh%nfaces, coarse_mesh%face_nodes_ptr,&
										  coarse_mesh%face_nodes, coarse_mesh%face_zone)
			end if
		end if
	end subroutine








end module
