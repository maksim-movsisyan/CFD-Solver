program mesh_generator
	use mesh_module
	use output_fields_module
	use structured_mesh_module
	implicit none
	type(mesh_t) :: mesh
	type(structured_mesh_t) strd_mesh

	call strd_mesh%read_map('mesh.txt')
	mesh%dim = 2
	call strd_mesh%assemble(mesh%ncells, mesh%nfaces, mesh%nnodes, mesh%nbfaces,&
							mesh%node_coords, mesh%face_nodes, mesh%face_nodes_ptr,&
							mesh%face_left_cell, mesh%face_right_cell, mesh%face_zone,&	
							mesh%face_type, mesh%cell_type, mesh%zone_names, mesh%zone_subnames)
							
	call mesh%write_fluent('..\mesh\strcd_mesh.cas')

	
end program
