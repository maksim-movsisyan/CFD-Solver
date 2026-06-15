program mg_initialization
	use mesh_module
	use mesh_coarsening_module
	use solver_config_module
	use density_based_solver_module
	use mg_operators_module
	use output_fields_module
	implicit none
	character(len=128) :: meshfile
	type(mesh_t) :: mesh, coarse_mesh
	type(agglomeration_info_t) :: agglom_info
	type(multigrid_config_t) :: mg_settings
	type(density_based_solver_t) :: SOLVER
	real(8) :: time_start, time_end
	integer :: iunit, idx
	real(8), allocatable :: Q_fine(:, :)
	
	call cpu_time(time_start)
	
	open(newunit=iunit, file='..\input\B_mesh_filename.txt', status='old', action='read')
	read(iunit, '(a)') meshfile
	idx = index(meshfile, 'data')
	write(meshfile, '(a, a)') '..' // trim(adjustl(meshfile(idx+4:)))
	close(iunit)
		
	call mg_settings%read_cfg('..\input\multigrid_settings.txt')
	
	!=== MESH READING: ===
	mesh%verbosity=100
	call mesh%read_fluent(trim(adjustl(meshfile)))
	call mesh%compute_metric()
	
	!=== AGGLOMERATION: ===
	call agglom_info%agglomerate(mesh, mg_settings%target_size)
	call agglom_info%get_coarse(mesh, coarse_mesh, mg_settings%output_mesh)
	
	!=== SOLUTION ON COARSE MESH: ===
	call SOLVER%initialize(coarse_mesh, 'MGINIT_input.txt',&
										'..\input\B_boundary_conditions.txt',&
										'..\input\B_physical_propertis.txt',&
										'..\input\linear_solver_settings.txt')
	SOLVER%opm%output_vtk = .false.
	SOLVER%opm%output_dat = .false.
	call SOLVER%run()	
	
	!=== EXTRAPOLATION BACK TO FINE MESH: ===
	allocate(Q_fine(mesh%dim+2, mesh%ncells))
	Q_fine = 0.d0
	
	if (mg_settings%p_order == 1) then
		call mg_prolongation1_vector_cell(mesh%ncells, mesh%dim+2,&
										  agglom_info%fine_to_coarse,&
										  SOLVER%fldsm%registry(SOLVER%QIDX)%values, Q_fine)
	else
		call mg_prolongation2_vector_cell(mesh%dim, mesh%ncells, mesh%dim+2,&
										  agglom_info%fine_to_coarse,&
										  mesh%cell_center, coarse_mesh%cell_center,&
										  SOLVER%fldsm%registry(SOLVER%QIDX)%values,&
										  SOLVER%fldsm%registry(SOLVER%GRADQIDX)%values,&
										  Q_fine)
	end if
	
	!=== OUTPUT INITIALIZATION FILE: ===
	call output_fields_vtk(Q_fine, SOLVER%opm%field_names,&
						   'init_file\SOLUTION_0000.vtk', 'INITIALIZATION_FILE',&
						   mesh%dim, mesh%ncells, mesh%nfaces, mesh%nnodes,&
						   mesh%cell_nodes, mesh%cell_nodes_ptr, mesh%node_coords)
	
	call cpu_time(time_end)
	print*, 'Total executional time = ', time_end - time_start, ' sec'
end program
