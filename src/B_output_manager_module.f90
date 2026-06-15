module output_manager_module
use output_fields_module
use mesh_module, only: mesh_t
use bc_manager_module, only: bc_manager_t
use solver_config_module, only: db_solver_config_t
use bc_wall_module, only: bc_wall_t
implicit none
!=======================================================================
!==================== OUTPUT MANAGER DATA STRUCTURE ====================
!=======================================================================
type :: output_manager_t
	type(mesh_t), pointer :: mesh => null()
	type(bc_manager_t), pointer :: bcm => null()
	type(db_solver_config_t), pointer :: settings => null()
	
	!fields:
	character(len=128) :: output_dir = 'run\output\'					!=== OUTPUT DIRECTORY ===
	character(len=128) :: name = 'SOLUTION'								!=== OUTPUT NAME ==
	character(len=128) :: path											!=== FULL OUTPUT PATH = DIR+NAME ===
	
	character(len=32), allocatable :: field_names(:)
	integer, allocatable :: field_idx(:)
	
	logical :: output_vtk = .true.
	logical :: output_dat = .true.	
	
	integer :: current_step = 0
	integer :: SAVE_NITER
	
	!log file:
	integer :: log_iunit
	character(len=128) :: log_file = 'run\log\log_file.txt'				!=== LOG FILE ===
	logical :: IS_PRINT, IS_LOG
	integer :: PRINT_NITER, LOG_NITER
	logical :: is_log_header_write = .false.					
	
	!monitor points file:
	integer :: mp_iunit
	character(len=128) :: mp_file = 'run\monitors\mp_file.txt' 			!=== MONITOR POINTS FILE ===
	
	!integral characteritics file:
	integer :: ic_iunit
	character(len=128) :: ic_file = 'run\integrals\ic_file.txt' 		!=== INTERGRAL CHARACTERISTICS FILE ===
	
	contains
	procedure :: initialize => opm_initialize
	!fields:
	procedure :: output_fields => opm_output_fields
	procedure :: input_fields => opm_input_fields
	!log:
	procedure :: print_log => opm_print_log
	procedure :: write_log => opm_write_log
	!ic:
	procedure :: write_forces => opm_write_forces
	procedure :: write_wfluxes => opm_write_wfluxes
end type


contains
!=======================================================================
!======================= OUTPUT MANAGER METHODS ========================
!=======================================================================
	subroutine opm_initialize(this, mesh, bcm, settings, field_names)
        implicit none
        class(output_manager_t), intent(inout) :: this
        type(mesh_t), intent(in), target :: mesh
        type(bc_manager_t), intent(in), target :: bcm
        type(db_solver_config_t), intent(in), target :: settings
        character(len=*), intent(in) :: field_names(:)
        
        integer :: stat
        character(len=512) :: cmd
        
        this%mesh => mesh
        this%bcm => bcm
        this%settings => settings
        
        this%field_names = field_names
            
        this%path = trim(this%output_dir) // trim(this%name)
        
        cmd = 'if not exist "' // trim(this%path) // '" mkdir "' // trim(this%path) // '"'
        
        call execute_command_line(trim(cmd), exitstat=stat)

        if (stat /= 0) then
            print *, "Error executing command: ", trim(cmd)
        else
            print *, "Output directory is ready: ", trim(this%path)
        end if
    end subroutine
	
	subroutine opm_output_fields(this, iter, fields)
		implicit none
		class(output_manager_t), intent(inout) :: this
		integer, intent(in) :: iter
		real(8), intent(in) :: fields(:, :)								!=== DIM = (nvars, ncells); reshape([p, u, v, w, t], [5, size(p)])) ===
		
		character(len=256) :: filename
		character(len=128) :: step_title
		
		if (mod(iter, this%settings%SAVE_NITER) .ne. 0) return
		
		this%current_step = this%current_step + 1
		
		if (this%output_vtk) then
			!=== DIR NAME: ===
			write(filename, '(A, "\", A, "_", I4.4, ".vtk")') &
					trim(this%path), trim(this%name), this%current_step
			
			!=== FILE TITLE: ===
			write(step_title, '(A, " Step: ", I0)') trim(this%name), this%current_step
			
			!=== OUTPUT VTK: ===
			call output_fields_vtk(fields, this%field_names, trim(filename), trim(step_title),&
								   this%mesh%dim, this%mesh%ncells, this%mesh%nfaces, this%mesh%nnodes,&
								   this%mesh%cell_nodes, this%mesh%cell_nodes_ptr, this%mesh%node_coords)
		end if
		
		
		if (this%output_dat) then
			!=== DIR NAME: ===
			write(filename, '(A, "\", A, "_", I4.4, ".dat")') &
					trim(this%path), trim(this%name), this%current_step
			
			!=== FILE TITLE: ===
			write(step_title, '(A, " Step: ", I0)') trim(this%name), this%current_step
			
			!=== OUTPUT DAT: ===
			call output_fields_dat(fields, this%field_names, trim(filename), trim(step_title),&
								   this%mesh%ncells, this%mesh%nfaces, this%mesh%nbfaces,&
								   this%mesh%face_zone, this%mesh%face_type,&
								   this%mesh%face_left_cell, this%mesh%face_right_cell,&
								   this%settings%USE_GHOST_CELLS)	
		end if
		
		
		
	end subroutine
	
	subroutine opm_input_fields(this, fields, is_init)
		implicit none
		class(output_manager_t), intent(inout) :: this
		real(8), intent(inout) :: fields(:, :)								!=== DIM = (nvars, ncells); reshape([p, u, v, w, t], [5, size(p)])) ===
		logical, intent(inout) :: is_init
		
		character(len=256) :: filename
		
		is_init = .false.
		
		!=== .VTK INITIALIZATION: ===
		if (this%settings%INITIALIZATION_TYPE == 2) then
			deallocate(this%field_names)
			this%current_step = 1
			write(filename, '(A, "\", A, "_", I4.4, ".vtk")') &
					trim(this%path), trim(this%name), this%current_step
					
			call input_fields_vtk(fields, this%field_names, filename, this%mesh%dim, this%mesh%ncells)
			is_init = .true.
		end if 
		
		!=== .DAT INITIALIZATION: ===
		if (this%settings%INITIALIZATION_TYPE == 3) then
			deallocate(this%field_names)
			this%current_step = 1
			write(filename, '(A, "\", A, "_", I4.4, ".dat")') &
					trim(this%path), trim(this%name), this%current_step
					
			call input_fields_dat(fields, this%field_names, filename, this%mesh%ncells)
			is_init = .true.
		end if 
		
		!=== MULTIGRID INITIALIZATION: ===
		if (this%settings%INITIALIZATION_TYPE == 4) then
			deallocate(this%field_names)
			this%current_step = 0
			write(filename, '(A, "\", A, "_", I4.4, ".vtk")') &
					'data\mg_initialization\init_file', trim(this%name), this%current_step
					
			call input_fields_vtk(fields, this%field_names, filename, this%mesh%dim, this%mesh%ncells)
			is_init = .true.
		end if 
		
		
				
	end subroutine

	subroutine opm_print_log(this, iter, header_line, data_line)
		implicit none
		class(output_manager_t), intent(inout) :: this
		integer, intent(in) :: iter
		character(len=*), intent(in) :: header_line, data_line
		
		if (.not. this%settings%IS_PRINT) return
		
		if (mod(iter, this%settings%PRINT_NITER) == 0 .or. iter == 1) then
			if (mod(iter, 500) == 0 .or. iter == 1) then
				write(*, *) ''
				write(*, '(120("-"))')
				write(*, '(A)') trim(header_line)
				write(*, '(120("-"))')
			end if
			write(*, '(A)') trim(data_line)
		end if
	end subroutine 
	
	subroutine opm_write_log(this, iter, header_line, data_line)
		implicit none
		class(output_manager_t), intent(inout) :: this
		integer, intent(in) :: iter
		character(len=*), intent(in) :: header_line, data_line
		
		if (.not. this%settings%IS_LOG) return
		
		if (.not. this%is_log_header_write) then
			write(this%log_iunit, '(120("-"))')
			write(this%log_iunit, '(A)') trim(header_line)
			write(this%log_iunit, '(120("-"))')
			this%is_log_header_write = .true.
		end if
		
		if (mod(iter, this%settings%LOG_NITER) == 0 .or. iter == 1) then
			write(this%log_iunit, '(A)') trim(data_line)
		end if
	end subroutine
	
	subroutine opm_write_forces(this, iter, values_ptr, values_grad_ptr, USE_GHOST_CELLS)
		implicit none
		class(output_manager_t), intent(inout) :: this
		integer, intent(in) :: iter
		real(8), intent(in), contiguous :: values_ptr(:, :)
		real(8), intent(in), contiguous :: values_grad_ptr(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		
		real(8) :: force(this%mesh%dim)
		integer :: bc_id, i
		
		if (mod(iter, this%settings%LOG_NITER) .ne. 0 .and. iter .ne. 1) return
		
		write(this%ic_iunit, '(i0)', advance='no') iter
		
		do bc_id = 1, this%bcm%num_bc
			select type(bc => this%bcm%bc_list(bc_id)%bc)
			type is (bc_wall_t)
				call bc%compute_force(values_ptr, values_grad_ptr, USE_GHOST_CELLS, force)
				
				do i = 1, this%mesh%dim
					write(this%ic_iunit, '(ES15.7)', advance='no') force(i)
				end do
			end select
		end do
		write(this%ic_iunit, *) ''
	end subroutine
	
	subroutine opm_write_wfluxes(this, iter, values_ptr, values_grad_ptr, Pr, cp, USE_GHOST_CELLS)
		implicit none
		class(output_manager_t), intent(inout) :: this
		integer, intent(in) :: iter
		real(8), intent(in) :: Pr, cp
		real(8), intent(in), contiguous :: values_ptr(:, :)
		real(8), intent(in), contiguous :: values_grad_ptr(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		
		integer :: bc_id
		character(len=256) :: filename
		
		if (mod(iter, this%settings%SAVE_NITER) .ne. 0) return
						
		do bc_id = 1, this%bcm%num_bc
			select type(bc => this%bcm%bc_list(bc_id)%bc)
			type is (bc_wall_t)
				write(filename, '(A, "\", A, i0, A, "_", I4.4, ".txt")') &
					trim('run\integrals'), 'wall', bc_id, 'Fluxes', this%current_step
				call bc%compute_fluxes(Pr, cp, values_ptr, values_grad_ptr, USE_GHOST_CELLS, filename)
				
			end select
		end do
	end subroutine

end module
