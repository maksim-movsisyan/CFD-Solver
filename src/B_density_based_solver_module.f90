module density_based_solver_module
use mesh_module, only: mesh_t
use field_manager_module, only: field_manager_t, LOC_CELL
use bc_manager_module, only: bc_manager_t
use fluxes_manager_module, only: flux_manager_t
use physical_properties_manager_module, only: phys_prop_manager_t
use output_manager_module, only: output_manager_t
use solver_config_module, only: db_solver_config_t
use diff_operator_module, only: gradient_t
use diff_operator_factory_module, only: gradient_factory
use density_based_iterators_module
use density_based_iterator_factory_module
use time_step_calculation_module
implicit none

!=======================================================================
!================= DENSITY-BASED SOLVER DATA STRUCTURE =================
!=======================================================================
type :: density_based_solver_t
	character(len=128) :: cfg_file =&
								'data\input\COUPLED_solver_settings.txt'!=== INPUT FILE ===
	character(len=128) :: log_file = 'run\log\COUPLED_log_file.txt'		!=== LOG FILE ===
	character(len=128) :: bc_filename =&
								'data\input\B_boundary_conditions.txt'	!=== BC FILE ===
	character(len=128) :: pp_filename =&
								'data\input\B_physical_propertis.txt'	!=== PHYSICAL PROPERTIES FILE ===
	
	
	!== SOLVER PARAMETERS: ===
	type(db_solver_config_t) :: settings
	
	!=== SOLVER DATA: ===
	real(8), allocatable :: ptime_step(:)								!=== DIM = (ncells) ===
	real(8), allocatable :: cell_lenght(:)								!=== CELL CHARACTERISTIC LENGTH: DIM = (ncells) ===
	real(8), allocatable :: Q_FOS(:, :)									!=== IF IS_FOS => DIM = (nvars, ncells) ===
	
	real(8), allocatable :: eq_res(:)									!=== DIM = (nvars) ===
	real(8), allocatable :: init_eq_res(:)								!=== DIM = (nvars) ===
	integer :: num_iter = 0
	
	integer :: QIDX, GRADQIDX, FLUXESQIDX
	integer :: UIDX = -1
	
	type(field_manager_t) :: fldsm
	type(bc_manager_t) :: bcm
	type(flux_manager_t) :: flxsm
	type(phys_prop_manager_t) :: ppm
	type(output_manager_t) :: opm
	
	type(mesh_t), pointer :: mesh => null()
	class(gradient_t), pointer :: gradient => null()
	class(density_based_iterator_t), pointer :: iterator => null()
	contains
	procedure :: initialize => dbs_initialize
	procedure :: ptime_step_calc => dbs_calculate_ptime_step
	procedure :: set_initial_guess => dbs_set_initial_guess
	procedure :: run => dbs_run_calculation
end type


contains
!=======================================================================
!===================== DENSITY-BASED SOLVER METHODS ====================
!=======================================================================
	subroutine dbs_initialize(this, mesh, filename, bc_filename, pp_filename, ls_filename)
		implicit none
		class(density_based_solver_t), intent(inout) :: this
		type(mesh_t), target, intent(in) :: mesh
		character(len=*), optional, intent(in)  :: filename, bc_filename, pp_filename, ls_filename
		
		
		!=== READING SOLVER AND PHYSICAL PARAMETERS: ===
		if (present(filename)) this%cfg_file = filename
		if (present(bc_filename)) this%bc_filename = bc_filename
		if (present(pp_filename)) this%pp_filename = pp_filename
		call this%settings%read_cfg(this%cfg_file)
		call this%ppm%initialize(this%pp_filename)
			
		!=== MESH POINTER: ===
		this%mesh => mesh
		
		!=== GRADIENT POINTER: ===
		this%gradient => gradient_factory%create(mesh, this%settings%GRADIENT_TYPE)
		
		!=== FIELDS MANAGER: ===
		call this%fldsm%initialize(mesh)
			!primitive variables [P U V (W) T]:
		call this%fldsm%add_field('Q',  .true., LOC_CELL, mesh%dim+2, mesh%ncells, mesh%nbfaces)
		this%QIDX = this%fldsm%get_idx('Q')
			!primitive variables gradients [gradP_x, gradP_y, ... , gradT_x, gradT_y, (gradT_z)]:
		call this%fldsm%add_field('gradQ',  .false., LOC_CELL, mesh%dim*(mesh%dim+2), mesh%ncells, mesh%nbfaces)
		this%GRADQIDX = this%fldsm%get_idx('gradQ')
			!connecting variables with gradeints:
		call this%fldsm%registry(this%QIDX)%set_gradient(this%fldsm%registry(this%GRADQIDX)%values)
			!fluxes:
		call this%fldsm%add_field('fluxesQ',  .false., LOC_CELL, mesh%dim+2, mesh%ncells, mesh%nbfaces)
		this%FLUXESQIDX = this%fldsm%get_idx('fluxesQ')
		
			!consrvaive variables [Ro RoU RoV (RoW) RoE]:
		if (this%settings%USE_CONSERVATIVE_VARS) then
			call this%fldsm%add_field('U',  .false., LOC_CELL, mesh%dim+2, mesh%ncells, mesh%nbfaces)
			this%UIDX = this%fldsm%get_idx('U')
		end if
		
		!=== BOUNDARY CONDITION MANAGER: ===
		call this%bcm%initialize(mesh, this%fldsm, this%ppm, this%bc_filename)
		
		!=== FLUXES MANAGER: ===
		call this%flxsm%initialize(mesh, this%fldsm, this%ppm, this%settings)
		
		!=== ARRAYS ALLOCATION: ===
		allocate(this%ptime_step(mesh%ncells),&
				 this%cell_lenght(mesh%ncells),&
				 this%eq_res(mesh%dim+2), this%init_eq_res(mesh%dim+2))
				 
		!=== CELL CHARACTERISTIC LENGTH COMPUTATION: ===
		call compute_cell_lenght()
		
		!=== ITERATOR INITIALIZATION: ===
		this%iterator => db_iterator_factory%create(this%settings%SCHEME)
		if (present(ls_filename)) this%iterator%ls_filename = ls_filename	
		call this%iterator%initialize(this%mesh, this%fldsm, this%flxsm, this%bcm,&
									  this%ppm, this%gradient, this%settings, this%ptime_step,&
									  this%QIDX, this%UIDX, this%GRADQIDX, this%FLUXESQIDX)					  
		
		
		!=== OUTPUT MANAGER INITIALIZATION: ===
		if (mesh%dim == 2) then
			call this%opm%initialize(this%mesh, this%bcm, this%settings,&
							     [character(len=32) :: "Pressure", "Velocity_X", "Velocity_Y", "Temperature"])
		else
			call this%opm%initialize(this%mesh, this%bcm, this%settings,&
							     [character(len=32) :: "Pressure", "Velocity_X", "Velocity_Y", "Velocity_Z", "Temperature"])
		end if
		
		contains
		subroutine compute_cell_lenght()
			integer :: cell_idx, face_idx
			integer :: i, pos1, pos2
			real(8) :: dx, area
			
			do cell_idx = 1, this%mesh%ncells		
				!=== CHARACTERISTIC LENGTH: ===
				dx = 0.d0
				pos1 = this%mesh%cell_faces_ptr(cell_idx)
				pos2 = this%mesh%cell_faces_ptr(cell_idx + 1) - 1
				do i = pos1, pos2
					face_idx = this%mesh%cell_faces(i)
					area = this%mesh%face_area(face_idx)
					
					dx = dx + area
				end do			
				dx = this%mesh%cell_volume(cell_idx)/dx
				this%cell_lenght(cell_idx) = dx
			end do
		end subroutine
	end subroutine
	
	subroutine dbs_calculate_ptime_step(this)
		implicit none
		class(density_based_solver_t), intent(inout), target :: this
		
		real(8), pointer, contiguous :: vars_ptr(:, :)
		integer :: cell_idx, ncells, dim
		integer :: i, d
		real(8) :: P, V(this%mesh%dim), T, R_gas, k
		real(8) :: dx, dtau_inv, dtau_visc
		real(8) :: CFL
		integer :: MODEL
		
		vars_ptr => this%fldsm%registry(this%QIDX)%values
		ncells = this%mesh%ncells
		dim = this%mesh%dim
		R_gas = this%ppm%R_gas; k = this%ppm%k
		CFL = this%settings%CFL; MODEL = this%settings%MODEL
		
		do cell_idx = 1, ncells
			!=== CELL VARIABLES: ===
			P = vars_ptr(1, cell_idx)
			do d = 1, dim
				V(d) = vars_ptr(1+d, cell_idx)
			end do
			T = vars_ptr(dim+2, cell_idx)
			
			!=== CHARACTERISTIC LENGTH: ===
			dx = this%cell_lenght(cell_idx)
			
			!=== PTIME-STEP CALCULATION: ===
			call compute_ptime_step(dim, CFL, P, V, T, k, R_gas,&
								    dx, dtau_inv, dtau_visc)
			
			if (MODEL == 2) then
				this%ptime_step(cell_idx) = 1.d0/(1.d0/dtau_inv + 1.d0/dtau_visc)
			else
				this%ptime_step(cell_idx) = dtau_inv
			end if
		end do
	
	end subroutine

	subroutine dbs_set_initial_guess(this)
		implicit none
		class(density_based_solver_t), intent(inout), target :: this
		
		real(8), pointer, contiguous :: vars_ptr(:, :)
		real(8), pointer, contiguous :: cvars_ptr(:, :)
		real(8) :: P0, V0(3), T0, Ro0, E0
		integer :: i, d
		logical :: is_init
		
		vars_ptr => this%fldsm%registry(this%QIDX)%values
		call this%opm%input_fields(vars_ptr, is_init)
			
		if (.not. is_init) then	
			call this%bcm%get_ref_values(V0, P0, T0)

			do i = 1, this%mesh%ncells
				vars_ptr(1, i) = P0
				do d = 1, this%mesh%dim
					vars_ptr(1+d, i) = V0(d)
				end do
				vars_ptr(this%mesh%dim+2, i) = T0
			end do
		end if
		
		if (this%settings%USE_CONSERVATIVE_VARS) then
			cvars_ptr => this%fldsm%registry(this%UIDX)%values
			
			do i = 1, this%mesh%ncells
				Ro0 = P0/(this%ppm%R_gas*T0)
				cvars_ptr(1, i) = Ro0
				do d = 1, this%mesh%dim
					cvars_ptr(1+d, i) = Ro0*V0(d)
				end do
				cvars_ptr(this%mesh%dim+2, i) = Ro0*(this%ppm%cv*T0 + 0.5d0*dot_product(V0(1:this%mesh%dim), V0(1:this%mesh%dim)))
			end do
		end if
		
	end subroutine
	
	subroutine dbs_run_calculation(this)
		implicit none
		class(density_based_solver_t), intent(inout), target :: this
		
		integer :: iter, i, d, ncells, dim, TRUE_ORDER
		real(8) :: time_start, TRUE_CFL, weight
		real(8), pointer, contiguous :: vars_ptr(:, :), grad_vars_prt(:, :),&
										flxs_ptr(:, :), cvars_ptr(:, :)
		
		ncells = this%mesh%ncells
		dim = this%mesh%dim
		
		!=== SETTING INIIAL GUESS: ===
		call this%set_initial_guess()
		
		!=== SETTING FIELDS POINTERS: ===
		vars_ptr => this%fldsm%registry(this%QIDX)%values
		grad_vars_prt => this%fldsm%registry(this%GRADQIDX)%values
		flxs_ptr => this%fldsm%registry(this%FLUXESQIDX)%values
		if (this%UIDX > 0) cvars_ptr => this%fldsm%registry(this%UIDX)%values
		
		call cpu_time(time_start)
		open(newunit=this%opm%log_iunit, file=this%opm%log_file, status='replace')
		open(newunit=this%opm%ic_iunit, file=this%opm%ic_file, status='replace')
		!=== FIRST ORDER SMOOTHING LOOP: ===
		TRUE_ORDER = this%settings%ORDER
		TRUE_CFL = this%settings%CFL
		allocate(this%Q_FOS(dim+2, ncells))
		
		!=== FOS_NITER/2 ITERATIONS WITH 1ST ORDER AND LOW CFL: ===
		this%settings%CFL = this%settings%FOS_CFL
		this%settings%ORDER = 1
		do iter = 1, this%settings%FOS_NITER/2
			call itern()
		end do
		
		!=== FOS_NITER/2 ITERATIONS WITH WEIGHTNING 1ST ORDER	 ===
		!=== AND HIGH ORDER AND WITH LOW CFL: ===
		do iter = this%settings%FOS_NITER/2 + 1, this%settings%FOS_NITER
			weight = real(iter - this%settings%FOS_NITER/2, 8)/real(this%settings%FOS_NITER/2, 8)
			
			!=== FIRST ORDER STEP: ===
			this%settings%ORDER = 1
			call itern()
			do i = 1, ncells
				do d = 1, dim+2
					this%Q_FOS(d, i) = vars_ptr(d, i)
				end do
			end do
			
			!=== HIGH ORDER STEP: ===
			this%settings%ORDER = TRUE_ORDER
			call itern()
			
			do i = 1, ncells
				vars_ptr(:, i) = weight*vars_ptr(:, i) + (1.d0 - weight)*this%Q_FOS(:, i)
			end do
			
			!=== UPDATING CONSERVATIVE VARIABLES: ===
			if (this%settings%USE_CONSERVATIVE_VARS) then
				do i = 1, ncells
					!density:
					cvars_ptr(1, i) = vars_ptr(1, i)/(this%ppm%R_gas*vars_ptr(dim+2, i))
					
					!momentum:
					do d = 1, dim
						cvars_ptr(1+d, i) = cvars_ptr(1, i)*vars_ptr(1+d, i)
					end do
					
					!energy:
					cvars_ptr(dim+2, i) = cvars_ptr(1, i)*(this%ppm%cv*vars_ptr(dim+2, i)&
										  + 0.5d0*dot_product(vars_ptr(2:dim+1, i), vars_ptr(2:dim+1, i)))
				end do	
			end if 				
		end do
		
		this%settings%CFL = TRUE_CFL
		this%settings%ORDER = TRUE_ORDER
		deallocate(this%Q_FOS)
		
		!=== NON LINEAR ITERATIONS LOOP: ==
		do iter = this%settings%FOS_NITER + 1, this%settings%NITER
			call itern()		
		end do
		close(this%opm%log_iunit)
		close(this%opm%ic_iunit)
			
		!=== OUTPUT FIELDS: ===
		call this%opm%output_fields(this%settings%SAVE_NITER, vars_ptr)
		call this%opm%write_wfluxes(this%settings%SAVE_NITER, vars_ptr, grad_vars_prt,&
									this%ppm%Pr, this%ppm%cp,&
									this%settings%USE_GHOST_CELLS)
		contains
		subroutine itern()
			!=== UPDATING BOUNDARY FACES/GHOST CELLS VALUES: ===
			call this%bcm%update_boundary_values('Q', .true., this%settings%USE_GHOST_CELLS, this%QIDX)
			
			!=== COMPUTING GRADIENTS: ===
			call this%gradient%apply_vector_cellfield(vars_ptr, grad_vars_prt, this%mesh%dim + 2, this%settings%USE_GHOST_CELLS)
			
			!=== UPDATING BOUNDARY FACES/GHOST CELLS GRADIENTS: ===
			call this%bcm%update_boundary_values('Q', .true., this%settings%USE_GHOST_CELLS, this%QIDX)
			
			!=== COMPUTING PSEUDO-TIME STEP: ===
			call this%ptime_step_calc()
			
			!=== COMPUTING FLUXES: ===
			call this%flxsm%compute_fluxes_cmprs(this%QIDX, this%FLUXESQIDX)
			
			!=== UPDATING SOLUTION: ===
			call this%iterator%iterate(iter, time_start)
			
			!=== RECOMPUTING PRIMITIVE VARIABLES (IF USE_CV): ===
			if (this%settings%USE_CONSERVATIVE_VARS) then
				do i = 1, ncells
					!velocity:
					do d = 2, dim + 1
						vars_ptr(d, i) = cvars_ptr(d, i)/cvars_ptr(1, i)
					end do
					
					!temperature:
					vars_ptr(dim+2, i) = 1.d0/(this%ppm%cv)*(cvars_ptr(dim+2, i)/cvars_ptr(1, i) - 0.5d0*dot_product(vars_ptr(2:dim+1, i), vars_ptr(2:dim+1, i)))
					
					!pressure:
					vars_ptr(1, i) = cvars_ptr(1, i)*vars_ptr(dim+2, i)*this%ppm%R_gas
					
				end do
			end if
			
			!=== PRINT INFORMATION: ===
			call this%opm%print_log(iter, this%iterator%log_header, this%iterator%log_data)
			
			!=== LOG INFORMATION: ===
			call this%opm%write_log(iter, this%iterator%log_header, this%iterator%log_data)
			
			!=== LOG IC: ===
			call this%opm%write_forces(iter, vars_ptr, grad_vars_prt, this%settings%USE_GHOST_CELLS)
			
			!=== OUTPUT FIELDS: ===
			call this%opm%output_fields(iter, vars_ptr)
			call this%opm%write_wfluxes(iter, vars_ptr, grad_vars_prt,&
										this%ppm%Pr, this%ppm%cp,&
										this%settings%USE_GHOST_CELLS)
			
			this%num_iter = this%num_iter + 1	
		end subroutine
	end subroutine










end module
