program main
	use mesh_module
	use density_based_solver_module
	implicit none
	character(len=128) :: meshfile
	type(density_based_solver_t) :: SOLVER
	type(mesh_t) :: mesh
	real(8) :: time_start, time_end
	integer :: iunit
	
	call cpu_time(time_start)
	
	open(newunit=iunit, file='data\input\B_mesh_filename.txt', status='old', action='read')
	read(iunit, '(a)') meshfile
	close(iunit)
	
	mesh%verbosity=100
	call mesh%read_fluent(meshfile)
	call mesh%compute_metric()
	call SOLVER%initialize(mesh)
	call SOLVER%run()
	call mesh%write_fluent('run\output\mesh.cas')
	
	
	
	call cpu_time(time_end)
	print*, 'Total executional time = ', time_end - time_start, ' sec'
end program
