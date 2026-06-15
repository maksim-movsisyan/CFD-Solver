module mg_mesh_hierarchy_module
use mesh_module
use solver_config_module
use mesh_coarsening_module
use linear_operator_module
use mg_operators_module

use preconditioner_bcsr_ilu0_module
use matrix_initialization_module
use db_mtrx_assembling_module

!=== FOR DB SOLVERS: ===
use field_manager_module, only: field_manager_t, LOC_CELL
use bc_manager_module, only: bc_manager_t
use fluxes_manager_module, only: flux_manager_t
use physical_properties_manager_module, only: phys_prop_manager_t
use diff_operator_module
use solver_config_module, only: db_solver_config_t
use time_step_calculation_module

implicit none
!=======================================================================
!============= MULTI GRID HIERARCHY LEVEL DATA STRUCTURE ===============
!=======================================================================
type :: level_data_t
	type(mesh_t), pointer :: mesh => null()
	type(agglomeration_info_t) :: agglom_info
	
	integer :: level_number                   
	logical :: is_finest                       
	logical :: is_coarsest   
	
	real(8), allocatable :: phi(:), rhs(:), tmp(:)
	class(linear_operator_t), pointer :: linear_operator => null()
	
	real(8), allocatable :: res(:), cor(:)
	
	type(bcsr_operator_t) :: LU
	integer, allocatable :: iw(:)
	
	!=== ADDITIONAL DATA FOR DB SOLVERS: ===
	type(field_manager_t), pointer :: fldsm => null()
	type(bc_manager_t), pointer :: bcm => null()
	type(flux_manager_t), pointer :: flxsm => null()
	class(gradient_t), pointer :: gradient => null()
	type(phys_prop_manager_t), pointer :: ppm => null()	
	type(db_solver_config_t), pointer :: settings => null()
	
	integer :: QIDX, GRADQIDX, FLUXESQIDX
	integer :: UIDX = -1
	
	real(8), pointer, contiguous :: ptime_step(:) => null()					 
	real(8), pointer, contiguous :: cell_lenght(:) => null()
	
	real(8), pointer, contiguous :: vars_ptr(:, :) => null(),&
									grad_vars_ptr(:, :) => null(),&
									flxs_ptr(:, :) => null(),&
									cvars_ptr(:, :) => null()
									
									
	integer, allocatable :: map_LL(:), map_LR(:), map_RL(:),&
							map_RR(:), map_LB(:)			
									
	contains
	procedure :: initialize_db => level_data_init_dens_bsd		 
	procedure :: compute_ptime => level_data_compute_ptime_step
	procedure :: smooth_db => level_data_smooth_db
end type


!=======================================================================
!============ MULTI GRID MULTILEVEL HIERARCHY DATA STRUCTURE ===========
!=======================================================================
type :: multilevel_hierarchy_t
	type(level_data_t), allocatable :: levels(:)
	type(multigrid_config_t) :: settings     
	integer :: current_max_level       
	
	contains
	procedure :: initialize => mg_hierarchy_initialize
	procedure :: add_level => mg_hierarchy_add_level
	procedure :: create => mg_hierarchy_create
	procedure :: v_cycle => mg_perform_vcycle
	
	!=== FOR DB SOLVER: ===
	procedure :: initialize_db => mg_hierarchy_initialize_db
end type

contains
!=======================================================================
!================= MULTI GRID HIERARCHY LEVEL METHODS ==================
!=======================================================================
subroutine level_data_init_dens_bsd(this, ppm, settings, bc_filename)
	implicit none
	class(level_data_t), intent(inout) :: this
	type(phys_prop_manager_t), intent(in), target :: ppm
	type(db_solver_config_t), intent(in), target :: settings
	character(len=*), intent(in) :: bc_filename
	
	!=== SETTINGS POINTER: ===
	this%settings => settings
	
	!=== PHYSICAL PROPERTIES MANAGER POINTER: ===
	this%ppm => ppm
	
	!=== GRADIENT: ===
	allocate(LSQ_gradient_t :: this%gradient)
	call this%gradient%set_mesh(this%mesh)
	
	!=== FIELDS MANAGER: ===
	allocate(field_manager_t::this%fldsm)
	call this%fldsm%initialize(this%mesh)
		!primitive variables [P U V (W) T]:
	call this%fldsm%add_field('Q',  .true., LOC_CELL, this%mesh%dim+2, this%mesh%ncells, this%mesh%nbfaces)
	this%QIDX = this%fldsm%get_idx('Q')
	this%vars_ptr => this%fldsm%registry(this%QIDX)%values
		!primitive variables gradients [gradP_x, gradP_y, ... , gradT_x, gradT_y, (gradT_z)]:
	call this%fldsm%add_field('gradQ',  .false., LOC_CELL, this%mesh%dim*(this%mesh%dim+2), this%mesh%ncells, this%mesh%nbfaces)
	this%GRADQIDX = this%fldsm%get_idx('gradQ')
	this%grad_vars_ptr => this%fldsm%registry(this%GRADQIDX)%values
		!connecting variables with gradeints:
	call this%fldsm%registry(this%QIDX)%set_gradient(this%fldsm%registry(this%GRADQIDX)%values)
		!fluxes:
	call this%fldsm%add_field('fluxesQ',  .false., LOC_CELL, this%mesh%dim+2, this%mesh%ncells, this%mesh%nbfaces)
	this%FLUXESQIDX = this%fldsm%get_idx('fluxesQ')
	this%flxs_ptr => this%fldsm%registry(this%FLUXESQIDX)%values
	
		!consrvaive variables [Ro RoU RoV (RoW) RoE]:
	if (this%settings%USE_CONSERVATIVE_VARS) then
		call this%fldsm%add_field('U',  .false., LOC_CELL, this%mesh%dim+2, this%mesh%ncells, this%mesh%nbfaces)
		this%UIDX = this%fldsm%get_idx('U')
		this%cvars_ptr => this%fldsm%registry(this%UIDX)%values
	end if
	
	!=== BOUNDARY CONDITION MANAGER: ===
	allocate(bc_manager_t::this%bcm)
	call this%bcm%initialize(this%mesh, this%fldsm, this%ppm, bc_filename)
	
	!=== FLUXES MANAGER: ===
	allocate(flux_manager_t::this%flxsm)
	call this%flxsm%initialize(this%mesh, this%fldsm, this%ppm, this%settings)
	
	!=== ARRAYS ALLOCATION: ===
	allocate(this%ptime_step(this%mesh%ncells),&
			 this%cell_lenght(this%mesh%ncells))
			 
	!=== CELL CHARACTERISTIC LENGTH COMPUTATION: ===
	call compute_cell_lenght()
		
	
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

subroutine level_data_compute_ptime_step(this)
	implicit none
	class(level_data_t), intent(inout), target :: this
		
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

subroutine level_data_smooth_db(this, niter, omega)
	implicit none
	class(level_data_t), intent(inout) :: this
	integer, intent(in) :: niter
	real(8), intent(in) :: omega 
	
	integer :: k
	
	do k = 1, niter
		if (norm2(this%phi) > 1e-10) then
			call this%linear_operator%apply(this%phi, this%res)
		else
			this%res = 0.d0
		end if
		this%res = this%rhs - this%res
		
		call apply_bcsr_ilu0_preconditioner(this%LU%n,&
											this%LU%bs,&
											this%LU%values,&
											this%LU%col_indices,&
											this%LU%diag_indices,&
											this%LU%row_ptr, this%res, this%tmp)
											
		this%phi = this%phi + omega*this%tmp
	end do
	
end subroutine

!=======================================================================
!============= MULTI GRID MULTILEVEL HIERARCHY METHODS =================
!=======================================================================
subroutine mg_hierarchy_initialize(this, mesh, filename)
	implicit none
	class(multilevel_hierarchy_t), intent(inout) :: this
	type(mesh_t), target, intent(in) :: mesh
	character(len=*), optional, intent(in) :: filename
	
	integer :: dim, ncells, nvars
	
	if (present(filename)) then
		call this%settings%read_cfg(filename)
	else 
		call this%settings%read_cfg('data\input\multigrid_settings.txt')
	end if
		
	allocate(this%levels(this%settings%nlevels))
	
	this%levels(1)%mesh => mesh
	this%levels(1)%level_number = 1
	this%levels(1)%is_finest = .true.
	this%levels(1)%is_coarsest = .false.
	
	
	dim = mesh%dim
	ncells = mesh%ncells
	nvars = dim + 2
	
	allocate(this%levels(1)%phi(nvars*ncells))
	allocate(this%levels(1)%rhs(nvars*ncells))
	allocate(this%levels(1)%res(nvars*ncells))
	allocate(this%levels(1)%cor(nvars*ncells))
	allocate(this%levels(1)%tmp(nvars*ncells))
	
	this%current_max_level = 1
end subroutine

subroutine mg_hierarchy_add_level(this)
	implicit none
	class(multilevel_hierarchy_t), intent(inout) :: this
	
	integer :: new_level
	integer :: dim, ncells, nvars
	character(len=256) :: filename
	
	
	if (this%current_max_level >= this%settings%nlevels) return
	
	new_level = this%current_max_level + 1
	
	allocate(mesh_t :: this%levels(new_level)%mesh)
	
	write(filename, '(a, i0, ".vtk")') 'run\output\MG_HIERARCHY\mesh', new_level
	
	call this%levels(new_level)%agglom_info%agglomerate(this%levels(new_level-1)%mesh, this%settings%target_size)
	call this%levels(new_level)%agglom_info%get_coarse(this%levels(new_level-1)%mesh, this%levels(new_level)%mesh,&
													   this%settings%output_mesh, filename)
	
	this%levels(new_level)%level_number = new_level
	this%levels(new_level)%is_finest = .false.
	this%levels(new_level)%is_coarsest = .false.
	
	this%levels(new_level-1)%is_coarsest = .false.
	
	this%current_max_level = new_level
	
	
	dim = this%levels(new_level)%mesh%dim
	ncells = this%levels(new_level)%mesh%ncells
	nvars = dim + 2
	
	allocate(this%levels(new_level)%phi(nvars*ncells))
	allocate(this%levels(new_level)%rhs(nvars*ncells))
	allocate(this%levels(new_level)%res(nvars*ncells))
	allocate(this%levels(new_level)%cor(nvars*ncells))
	allocate(this%levels(new_level)%tmp(nvars*ncells))

	this%levels(new_level)%is_coarsest = .true.
	
end subroutine

subroutine mg_hierarchy_create(this, mesh, filename)
	implicit none
	class(multilevel_hierarchy_t), intent(inout) :: this
	type(mesh_t), target, intent(in) :: mesh
	character(len=*), optional, intent(in) :: filename
	
	integer :: level
	
	call this%initialize(mesh)       
	
	do level = 1, this%settings%nlevels - 1
		call this%add_level()
	end do
	
end subroutine

subroutine mg_hierarchy_initialize_db(this, linear_operator, flxsm, fldsm, bcm, ppm, settings,&
								      gradient, ptime_step,&
								      QIDX, GRADQIDX, FLUXESQIDX, UIDX, bc_filename)
	implicit none
	class(multilevel_hierarchy_t), intent(inout) :: this
	class(linear_operator_t), intent(in), target :: linear_operator
	type(flux_manager_t), intent(in), target :: flxsm
	type(field_manager_t), intent(in), target :: fldsm
	type(bc_manager_t), intent(in), target :: bcm
	type(phys_prop_manager_t), intent(in), target :: ppm
	type(db_solver_config_t), intent(in), target :: settings
	class(gradient_t), intent(in), target :: gradient
	real(8), intent(inout), target :: ptime_step(:)
	integer, intent(in) :: QIDX, GRADQIDX, FLUXESQIDX, UIDX
	character(len=*), intent(in) :: bc_filename
	
	integer :: i
	
	!=== FINES LEVEL DATA: ===
	this%levels(1)%flxsm => flxsm
	this%levels(1)%fldsm => fldsm
	this%levels(1)%bcm => bcm
	this%levels(1)%ppm => ppm
	this%levels(1)%settings => settings
	this%levels(1)%gradient => gradient
	
	this%levels(1)%ptime_step => ptime_step
	
	this%levels(1)%QIDX = QIDX
	this%levels(1)%UIDX = UIDX
	this%levels(1)%GRADQIDX = GRADQIDX
	this%levels(1)%FLUXESQIDX = FLUXESQIDX
	
	this%levels(1)%vars_ptr => this%levels(1)%fldsm%registry(QIDX)%values
	this%levels(1)%grad_vars_ptr => this%levels(1)%fldsm%registry(GRADQIDX)%values
	this%levels(1)%flxs_ptr => this%levels(1)%fldsm%registry(FLUXESQIDX)%values
	
	
	if (UIDX > 0) then
		this%levels(1)%cvars_ptr => this%levels(1)%fldsm%registry(UIDX)%values
	end if
	
	this%levels(1)%linear_operator => linear_operator
	
	!=== OTHER LEVELS DATA: ===
	do i = 2, this%current_max_level
		call this%levels(i)%initialize_db(ppm, settings, bc_filename)
	end do
	
	do i = 2, this%current_max_level
		select type(lo => linear_operator)
			type is (bcsr_operator_t)
				allocate(bcsr_operator_t::this%levels(i)%linear_operator)
				select type(lo2 => this%levels(i)%linear_operator)
				type is (bcsr_operator_t)
					lo2%bs = this%levels(i)%mesh%dim+2		 
					call bcsr_mtrx_initialize1(this%levels(i)%mesh%ncells, this%levels(i)%mesh%nfaces,&
											   this%levels(i)%mesh%nbfaces,lo2%bs,&
											   this%levels(i)%mesh%cell_faces_ptr, this%levels(i)%mesh%cell_faces,&
											   this%levels(i)%mesh%face_left_cell, this%levels(i)%mesh%face_right_cell,&
											   lo2%n, lo2%nnz,&
											   lo2%values,&
											   lo2%col_indices,&
											   lo2%row_ptr,&
											   lo2%diag_indices)
				end select
				
			type is (jacfree_operator_t)
				allocate(jacfree_operator_t::this%levels(i)%linear_operator)
				select type(lo2 => this%levels(i)%linear_operator)
				type is (jacfree_operator_t)
					call lo2%initialize(this%levels(i)%mesh,&
									    this%levels(i)%fldsm,&
									    this%levels(i)%flxsm,&
									    this%levels(i)%bcm,&
									    this%levels(i)%ppm,&
									    this%levels(i)%gradient,&
									    this%levels(i)%ptime_step,&
										this%levels(i)%QIDX,&
										this%levels(i)%UIDX,&
										this%levels(i)%FLUXESQIDX, this%levels(i)%settings%USE_GHOST_CELLS,&
										this%levels(i)%settings%USE_CONSERVATIVE_VARS)
				end select
				
				
			type is (algdiff_operator_t)
				allocate(algdiff_operator_t::this%levels(i)%linear_operator)
				select type(lo2 => this%levels(i)%linear_operator)
				type is (algdiff_operator_t)							 
					call lo2%initialize(this%levels(i)%mesh,&
									    this%levels(i)%fldsm,&
									    this%levels(i)%flxsm,&
									    this%levels(i)%bcm,&
									    this%levels(i)%ppm,&
									    this%levels(i)%gradient,&
									    this%levels(i)%ptime_step,&
										this%levels(i)%QIDX,&
										this%levels(i)%UIDX,&
										this%levels(i)%FLUXESQIDX, this%levels(i)%settings%USE_GHOST_CELLS,&
										this%levels(i)%settings%USE_CONSERVATIVE_VARS)
				end select
				
			end select
	end do
	
	!=== SMOOTHING MATRIX INITIALIZATION: ===
	do i = 1, this%current_max_level
		this%levels(i)%LU%bs = this%levels(i)%mesh%dim+2
		call bcsr_mtrx_initialize1(this%levels(i)%mesh%ncells, this%levels(i)%mesh%nfaces,&
								   this%levels(i)%mesh%nbfaces, this%levels(i)%LU%bs,&
								   this%levels(i)%mesh%cell_faces_ptr, this%levels(i)%mesh%cell_faces,&
								   this%levels(i)%mesh%face_left_cell, this%levels(i)%mesh%face_right_cell,&
								   this%levels(i)%LU%n, this%levels(i)%LU%nnz,&
								   this%levels(i)%LU%values,&
								   this%levels(i)%LU%col_indices,&
								   this%levels(i)%LU%row_ptr,&
								   this%levels(i)%LU%diag_indices)
		
		call db_create_mtrx_fill_map(this%levels(i)%mesh%nfaces, this%levels(i)%mesh%nbfaces,&
									 this%levels(i)%mesh%face_left_cell, this%levels(i)%mesh%face_right_cell,&
									 this%levels(i)%mesh%face_bidx, this%levels(i)%LU%col_indices,&
									 this%levels(i)%LU%row_ptr,&
									 this%levels(i)%map_LL, this%levels(i)%map_LR, this%levels(i)%map_RL,&
									 this%levels(i)%map_RR, this%levels(i)%map_LB)
		
		allocate(this%levels(i)%iw(this%levels(i)%mesh%ncells))		
	end do
	
end subroutine

recursive subroutine mg_perform_vcycle(this, level_fine)
	implicit none
	class(multilevel_hierarchy_t), intent(inout), target :: this
	integer, intent(in) :: level_fine

	integer :: i, level_coarse
	real(8), pointer, contiguous :: tmp_ptr_fine(:, :), tmp_ptr_coarse(:, :)
		
	
	
	!=== PRE SMOOTHING: ===
	call this%levels(level_fine)%smooth_db(this%settings%nu1, this%settings%omega)
	 
	!=== FINE MESH RESIDUAL: ===
	if (norm2(this%levels(level_fine)%phi) > 1e-8) then
		call this%levels(level_fine)%linear_operator%apply(this%levels(level_fine)%phi, this%levels(level_fine)%res)
	else
		this%levels(level_fine)%res = 0.0d0 
	end if 
	this%levels(level_fine)%res = this%levels(level_fine)%rhs - this%levels(level_fine)%res

	!=== TRANSITION ON THE NEXT LEVEL: ===
	level_coarse = level_fine + 1
	
	!=== RESIDUALS RESTRICTION: ===
	tmp_ptr_fine(1:this%levels(level_fine)%mesh%dim+2, 1:this%levels(level_fine)%mesh%ncells) => this%levels(level_fine)%res
	tmp_ptr_coarse(1:this%levels(level_coarse)%mesh%dim+2, 1:this%levels(level_coarse)%mesh%ncells) => this%levels(level_coarse)%res
	
	call mg_restriction_vector_cell(this%levels(level_coarse)%mesh%ncells,&
									this%levels(level_coarse)%mesh%dim+2,&
									this%levels(level_coarse)%mesh%cell_volume,&
									this%levels(level_fine)%mesh%cell_volume,&
									this%levels(level_coarse)%agglom_info%fine_cells_ptr,&
									this%levels(level_coarse)%agglom_info%fine_cells,&
									tmp_ptr_fine, tmp_ptr_coarse)
	
	this%levels(level_coarse)%rhs = this%levels(level_coarse)%res						  
	
	!=== RECURSION: ===
	if (level_coarse == this%current_max_level) then
		this%levels(level_coarse)%phi = 0.d0
		call this%levels(level_coarse)%smooth_db(this%settings%max_iter, this%settings%omega)

	else
		this%levels(level_coarse)%phi = 0.d0
		call this%v_cycle(level_coarse)
	end if
	
	
	!=== CORRECTION PROLONGATION: ===
	tmp_ptr_fine(1:this%levels(level_fine)%mesh%dim+2, 1:this%levels(level_fine)%mesh%ncells) => this%levels(level_fine)%cor
	tmp_ptr_coarse(1:this%levels(level_coarse)%mesh%dim+2, 1:this%levels(level_coarse)%mesh%ncells) => this%levels(level_coarse)%phi
	
	call mg_prolongation1_vector_cell(this%levels(level_fine)%mesh%ncells,&
									  this%levels(level_fine)%mesh%dim+2,&
									  this%levels(level_coarse)%agglom_info%fine_to_coarse,&
									  tmp_ptr_coarse, tmp_ptr_fine)
		
	
	this%levels(level_fine)%phi = this%levels(level_fine)%phi + this%levels(level_fine)%cor
	
	
	!=== POST SMOOTHING: ===
	call this%levels(level_fine)%smooth_db(this%settings%nu2, this%settings%omega)
end subroutine
 



end module
