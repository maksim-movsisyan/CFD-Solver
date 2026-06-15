module density_based_iterators_module
use mesh_module, only: mesh_t
use field_manager_module, only: field_manager_t
use fluxes_manager_module, only: flux_manager_t
use bc_manager_module, only: bc_manager_t
use physical_properties_manager_module, only: phys_prop_manager_t
use matrix_manager_module, only: matrix_manager_t
use solver_config_module, only: db_solver_config_t
use linear_solver_module, only: linear_solver_t
use linear_solver_factory_module, only: linear_solver_factory
use linear_operator_module, only: bcsr_operator_t, jacfree_operator_t, algdiff_operator_t
use diff_operator_module, only: gradient_t
use preconditioner_module, only: bcsr_ilu0_preconditioner_t, blusgs_preconditioner_t, gmg_preconditioner_t
use mg_mesh_hierarchy_module
use BLUSGS_module
implicit none
!=======================================================================
!=============== DENSITY-BASED ITERATOR DATA STRUCTURE =================
!=======================================================================
type, abstract :: density_based_iterator_t
	type(mesh_t), pointer :: mesh => null()
	type(field_manager_t), pointer :: fldsm => null()
	type(flux_manager_t), pointer :: flxsm => null()
	type(phys_prop_manager_t), pointer :: ppm => null()
	type(bc_manager_t), pointer :: bcm => null()
	
	type(db_solver_config_t), pointer :: settings => null()
	
	real(8), pointer, contiguous :: ptime_step(:) => null()
	
	real(8), pointer, contiguous :: vars_ptr(:, :) => null()			!=== PRIMITIVE VARIALBLES ARRAY ===
	real(8), pointer, contiguous :: cvars_ptr(:, :) => null()			!=== CONSERVATIVE VARIALBLES ARRAY ===
	real(8), pointer, contiguous :: flxs_ptr(:, :) => null()			!=== FLUXES ARRAY ===
	
	integer :: QIDX, GRADQIDX, FLUXESQIDX
	integer :: UIDX = -1
	
	character(len=200) :: dynamic_fmt_header, dynamic_fmt_data,&
						  log_header, log_data							!=== LOG UNITS ===
						  
						  
	character(len=128) :: ls_filename = 'data\input\linear_solver_settings.txt'
	
    contains
    procedure(dbi_initialize_interface), deferred :: initialize
    procedure(dbi_iterate_interface), deferred :: iterate
end type

abstract interface
    subroutine dbi_initialize_interface(this, mesh, fldsm, flxsm, bcm, ppm, gradient, settings, ptime_step,&
										QIDX, UIDX, GRADQIDX, FLUXESQIDX)
        import :: density_based_iterator_t, mesh_t, field_manager_t, flux_manager_t,&
				  phys_prop_manager_t, gradient_t, db_solver_config_t, bc_manager_t
        class(density_based_iterator_t), intent(inout) :: this
        type(mesh_t), intent(in), target :: mesh
        type(field_manager_t), intent(in), target :: fldsm
        type(flux_manager_t), intent(in), target :: flxsm
        type(bc_manager_t), intent(in), target :: bcm
        type(phys_prop_manager_t), intent(in), target :: ppm
        class(gradient_t), intent(in), target :: gradient
        type(db_solver_config_t), intent(in), target :: settings
        real(8), intent(in), target :: ptime_step(:)
        integer, intent(in) :: QIDX, UIDX, GRADQIDX, FLUXESQIDX
    end subroutine
    
    subroutine dbi_iterate_interface(this, iter, time_start)
        import :: density_based_iterator_t
        class(density_based_iterator_t), intent(inout), target :: this
        integer, intent(in) :: iter
        real(8), intent(in) :: time_start
    end subroutine
    
end interface

!=======================================================================
!=========== DENSITY-BASED EXPLICIT ITERATOR DATA STRUCTURE ============
!=======================================================================
type, extends(density_based_iterator_t) :: db_explicit_iterator_t
	contains
	procedure :: initialize => db_explt_initialize
    procedure :: iterate => db_explt_iterate
end type

!=======================================================================
!======= DENSITY-BASED CLASSICAL IMPLICIT ITERATOR DATA STRUCTURE ======
!=======================================================================
type, extends(density_based_iterator_t) :: db_implicit_iterator_t
	type(bcsr_operator_t) :: matrix										!=== SYSTEM MATRIX ===
	type(matrix_manager_t) :: mtrxm										!=== MATRIX MANAGER ===
	real(8), allocatable :: rhs(:), phi(:)								!=== DIM = (nvars*ncells) ===
	class(linear_solver_t), pointer :: linear_solver => null()			!=== LINEAR SOLVER ===
	
	contains
	procedure :: initialize => db_implt_initialize
    procedure :: iterate => db_implt_iterate
end type

!=======================================================================
!==== DENSITY-BASED JACOBIAN FREE IMPLICIT ITERATOR DATA STRUCTURE =====
!=======================================================================
type, extends(density_based_iterator_t) :: db_jacfree_iterator_t
	type(jacfree_operator_t) :: jacobianFree							!=== SYSTEM OPERATOR ===
	real(8), allocatable :: rhs(:), phi(:)								!=== DIM = (nvars*ncells) ===
	class(linear_solver_t), pointer :: linear_solver => null()			!=== LINEAR SOLVER ===

	type(bcsr_operator_t), pointer :: precond => null()					!=== PRECONDITIONER MATRIX ===
	type(matrix_manager_t) :: mtrxm										!=== MATRIX MANAGER FOR PRECONDITIONER MATRIX ===
							
	contains
	procedure :: initialize => db_jacfree_initialize
    procedure :: iterate => db_jacfree_iterate
end type

!=======================================================================
!====== DENSITY-BASED ALGH DIFF IMPLICIT ITERATOR DATA STRUCTURE =======
!=======================================================================
type, extends(density_based_iterator_t) :: db_algdiff_iterator_t
	type(algdiff_operator_t) :: algDiff									!=== SYSTEM OPERATOR ===
	real(8), allocatable :: rhs(:), phi(:)								!=== DIM = (nvars*ncells) ===
	class(linear_solver_t), pointer :: linear_solver => null()			!=== LINEAR SOLVER ===

	type(bcsr_operator_t), pointer :: precond => null()					!=== PRECONDITIONER MATRIX ===
	type(matrix_manager_t) :: mtrxm										!=== MATRIX MANAGER FOR PRECONDITIONER MATRIX ===
							
	contains
	procedure :: initialize => db_algdiff_initialize
    procedure :: iterate => db_algdiff_iterate
end type

!=======================================================================
!======== DENSITY-BASED BLUSGS IMPLICIT ITERATOR DATA STRUCTURE ========
!=======================================================================
type, extends(density_based_iterator_t) :: db_blusgs_iterator_t
	real(8), allocatable :: sprf(:)										!=== SPECTRAL RADIUS ON FACES, DIM = (nfaces) ===
	real(8), allocatable :: J_f(:, :, :)								!=== JACOBIANS ON FACES, DIM = (dim+2, dim+2, nfaces) ===
	real(8), allocatable :: J_d(:, :, :)								!=== DIAGONAL JACOBIANS, DIM = (dim+2, dim+2, ncells) ===
	
	real(8), allocatable :: phi(:)										!=== DIM = (nvars*ncells) ===
	real(8) :: w = 0.5d0, fluxes_sign = 1.d0			
	contains
	procedure :: initialize => db_blusgs_initialize
    procedure :: iterate => db_blusgs_iterate
end type

contains
!=======================================================================
!============== DENSITY-BASED EXPLICIT ITERATOR METHODS ================
!=======================================================================
	subroutine db_explt_initialize(this, mesh, fldsm, flxsm, bcm, ppm, gradient, settings, ptime_step,&
								   QIDX, UIDX, GRADQIDX, FLUXESQIDX)
        implicit none
        class(db_explicit_iterator_t), intent(inout) :: this
        type(mesh_t), intent(in), target :: mesh
        type(field_manager_t), intent(in), target :: fldsm
        type(flux_manager_t), intent(in), target :: flxsm
        type(bc_manager_t), intent(in), target :: bcm
        type(phys_prop_manager_t), intent(in), target :: ppm
        class(gradient_t), intent(in), target :: gradient
        type(db_solver_config_t), intent(in), target :: settings
        real(8), intent(in), target :: ptime_step(:)
        integer, intent(in) :: QIDX, UIDX, GRADQIDX, FLUXESQIDX
		
		character(len=20)  :: dim_str
		
		!=== POINTERS: ===
		this%mesh => mesh
		this%fldsm => fldsm
		this%flxsm => flxsm
		this%bcm => bcm
		this%ppm => ppm
	
		this%settings => settings
	
		this%ptime_step => ptime_step
		
		!=== VARIABLES INDICES: ===
		this%QIDX = QIDX
		this%UIDX = UIDX
		this%GRADQIDX = GRADQIDX
		this%FLUXESQIDX = FLUXESQIDX
		
		!=== VARIABLES POINTERS: ===
		this%vars_ptr => fldsm%registry(this%QIDX)%values
		this%flxs_ptr => fldsm%registry(this%FLUXESQIDX)%values
		if (UIDX > 0) this%cvars_ptr => fldsm%registry(this%UIDX)%values
		
		!=== LOG DATA SETTINGS: ===
		write(dim_str, '(I0)') this%mesh%dim + 2 
		this%dynamic_fmt_header = '(''  Iter'', ' // trim(dim_str) // 'A15, ''   Time(s)'')'
		this%dynamic_fmt_data   = '(I8, ' // trim(dim_str) // 'ES15.7, F10.3)'
		if (this%mesh%dim == 2) then
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Energy'
		else
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Mom-Z', 'Energy'
		end if
									
	end subroutine
	
	subroutine db_explt_iterate(this, iter, time_start)
		implicit none
		class(db_explicit_iterator_t), intent(inout), target :: this
		integer, intent(in) :: iter
		real(8), intent(in) :: time_start
		
		real(8), pointer, contiguous :: vars_ptr(:, :), cvars_ptr(:, :), flxs_ptr(:, :)	
		real(8) :: time_end								
		integer :: i, d
		integer :: ncells, dim
		
		ncells = this%mesh%ncells
		dim = this%mesh%dim
		
		vars_ptr => this%vars_ptr
		flxs_ptr => this%flxs_ptr
		cvars_ptr => this%cvars_ptr
		
		!=== UPDATING SOLUTION: ===
		if (.not. this%settings%USE_CONSERVATIVE_VARS) then
			!primitive variables:
			do i = 1, ncells
				do d = 1, dim + 2
					vars_ptr(d, i) = vars_ptr(d, i) - this%ptime_step(i)/this%mesh%cell_volume(i)*flxs_ptr(d, i)
				end do
			end do
		else
			!conservarive variables:
			do i = 1, ncells
				do d = 1, dim + 2
					cvars_ptr(d, i) = cvars_ptr(d, i) - this%ptime_step(i)/this%mesh%cell_volume(i)*flxs_ptr(d, i)
				end do
			end do
		end if
		
		!=== CPU TIME: ===
		call cpu_time(time_end)
	
		!=== LOG OUTPUT: ===
		if (this%settings%IS_PRINT .or. this%settings%IS_LOG) then
			if (mod(iter, this%settings%PRINT_NITER) == 0 .or. iter == 1) then
				write(this%log_data, this%dynamic_fmt_data) iter,&
					(maxval(dabs(flxs_ptr(i, 1:ncells))), i = 1, dim + 2),&
					time_end - time_start
			end if
		end if
    end subroutine


!=======================================================================
!============== DENSITY-BASED IMPLICIT ITERATOR METHODS ================
!=======================================================================
	subroutine db_implt_initialize(this, mesh, fldsm, flxsm, bcm, ppm, gradient, settings, ptime_step,&
								   QIDX, UIDX, GRADQIDX, FLUXESQIDX)
        implicit none
        class(db_implicit_iterator_t), intent(inout) :: this
        type(mesh_t), intent(in), target :: mesh
        type(field_manager_t), intent(in), target :: fldsm
        type(flux_manager_t), intent(in), target :: flxsm
        type(bc_manager_t), intent(in), target :: bcm
        type(phys_prop_manager_t), intent(in), target :: ppm
        class(gradient_t), intent(in), target :: gradient
        type(db_solver_config_t), intent(in), target :: settings
        real(8), intent(in), target :: ptime_step(:)
        integer, intent(in) :: QIDX, UIDX, GRADQIDX, FLUXESQIDX
		
		character(len=20)  :: dim_str
		
		!=== POINTERS: ===
		this%mesh => mesh
		this%fldsm => fldsm
		this%flxsm => flxsm
		this%bcm => bcm
		this%ppm => ppm
	
		this%settings => settings
	
		this%ptime_step => ptime_step
		
		!=== VARIABLES INDICES: ===
		this%QIDX = QIDX
		this%UIDX = UIDX
		this%GRADQIDX = GRADQIDX
		this%FLUXESQIDX = FLUXESQIDX
		
		!=== VARIABLES POINTERS: ===
		this%vars_ptr => fldsm%registry(this%QIDX)%values
		this%flxs_ptr => fldsm%registry(this%FLUXESQIDX)%values
		if (UIDX > 0) this%cvars_ptr => fldsm%registry(this%UIDX)%values
		
		!=== MATRIX MANAGER INITIALIZATION: ===
		call this%mtrxm%initialize(mesh, fldsm, bcm, ppm, settings, this%matrix,&
								  ptime_step, QIDX, GRADQIDX)
		
		
		!=== SYSTEM MATRIX INITIALIZATION: ===
		allocate(this%rhs(mesh%ncells*(mesh%dim + 2)),&
				 this%phi(mesh%ncells*(mesh%dim + 2)))
		call this%mtrxm%initialize_bcsr_mtrx()
		
		!=== LINEAR SOLVER INITIALIZATION: ===		 
		this%linear_solver => linear_solver_factory%create(this%ls_filename)
		call this%linear_solver%set_operator(this%matrix)		 
		
		!=== PRECONDITIONER SETUP ===
		if (associated(this%linear_solver%M)) then
			call this%linear_solver%M%setup(this%matrix, this%mesh)
			
			select type(prec => this%linear_solver%M)
			type is (blusgs_preconditioner_t)
				prec%ppm => ppm
				prec%settings => settings
				prec%Q => fldsm%registry(this%QIDX)%values
				prec%ptime_step => ptime_step
				
			type is (gmg_preconditioner_t)
				call prec%gmg_hierarchy%initialize_db(this%matrix, this%flxsm, this%fldsm, this%bcm, this%ppm, this%settings,&
													  gradient, this%ptime_step,&
													  this%QIDX, this%GRADQIDX, this%FLUXESQIDX, this%UIDX, 'data\input\B_boundary_conditions.txt')
				
				
			end select
		end if
		
		!=== LOG DATA SETTINGS: ===
		write(dim_str, '(I0)') this%mesh%dim + 2 
		this%dynamic_fmt_header = '(''  Iter'', ' // trim(dim_str) // 'A15, ''    L-Iter'', ''     L-Resid'', ''   Time(s)'')'
		this%dynamic_fmt_data   = '(I8, ' // trim(dim_str) // 'ES15.7, I10, ES15.7, F10.3)'
		if (this%mesh%dim == 2) then
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Energy'
		else
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Mom-Z', 'Energy'
		end if
	end subroutine
	
	subroutine db_implt_iterate(this, iter, time_start)
		implicit none
		class(db_implicit_iterator_t), intent(inout), target :: this
		integer, intent(in) :: iter
		real(8), intent(in) :: time_start
		
		real(8), pointer, contiguous :: vars_ptr(:, :), flxs_ptr(:, :), cvars_ptr(:, :)
		real(8) :: time_end					
		integer :: i, d, off
		integer :: dim, ncells
		
		dim = this%mesh%dim
		ncells = this%mesh%ncells
		
		vars_ptr => this%vars_ptr
		flxs_ptr => this%flxs_ptr
		cvars_ptr => this%cvars_ptr
		
		!=== RIGHT HAND SIDE INITIALIZATION: ===
		do i = 1, this%mesh%ncells
			off = (i-1)*(dim+2)
			do d = 1, dim+2
				this%rhs(off + d) = -flxs_ptr(d, i)
			end do
		end do
		
		!=== SYSTEM MATRIX ASSEMBLING: ===
		if (iter < 20 .or. mod(iter, this%settings%MAS_NITER) == 0) then
			call this%mtrxm%db_assmeble_bcsr_mtrx()
		end if
											 
		!=== PRECONDITIONER UPDATE: ===
		if (associated(this%linear_solver%M)) then
			if (iter < 20 .or. mod(iter, this%settings%PAS_NITER) == 0) then
				call this%linear_solver%M%update()
			end if
		end if
			
		!=== SLAE SOLUTION: ===
		this%phi = 0.d0
		call this%linear_solver%solve(this%rhs, this%phi)
		
		!=== UPDATING SOLUTION: ===
		if (.not. this%settings%USE_CONSERVATIVE_VARS) then
			!primitive variables:
			do i = 1, ncells
				off = (i-1)*(dim+2)
				vars_ptr(:, i) = vars_ptr(:, i) + this%phi(1 + off:dim+2 + off)
			end do
		else
			!conservative variables:
			do i = 1, ncells
				off = (i-1)*(dim+2)
				cvars_ptr(:, i) = cvars_ptr(:, i) + this%phi(1 + off:dim+2 + off)
			end do
		end if	
			
		
		!=== CPU TIME: ===
		call cpu_time(time_end)
		
		!=== LOG OUTPUT: ===
		if (this%settings%IS_PRINT .or. this%settings%IS_LOG) then
			if (mod(iter, this%settings%PRINT_NITER) == 0 .or. iter == 1) then					
				write(this%log_data, this%dynamic_fmt_data) iter,&
					(maxval(dabs(flxs_ptr(i, 1:ncells))), i = 1, dim + 2),&
					this%linear_solver%iterations,&
					this%linear_solver%residual_drop,&
					time_end - time_start
			end if
		end if
    end subroutine


!=======================================================================
!============= DENSITY-BASED JACOBIAN FREE ITERATOR METHODS ============
!=======================================================================
	subroutine db_jacfree_initialize(this, mesh, fldsm, flxsm, bcm, ppm, gradient, settings, ptime_step,&
								     QIDX, UIDX, GRADQIDX, FLUXESQIDX)
        implicit none
        class(db_jacfree_iterator_t), intent(inout) :: this
        type(mesh_t), intent(in), target :: mesh
        type(field_manager_t), intent(in), target :: fldsm
        type(flux_manager_t), intent(in), target :: flxsm
        type(bc_manager_t), intent(in), target :: bcm
        type(phys_prop_manager_t), intent(in), target :: ppm
        class(gradient_t), intent(in), target :: gradient
        type(db_solver_config_t), intent(in), target :: settings
        real(8), intent(in), target :: ptime_step(:)
        integer, intent(in) :: QIDX, UIDX, GRADQIDX, FLUXESQIDX
		
		character(len=20)  :: dim_str
		
		!=== POINTERS: ===
		this%mesh => mesh
		this%fldsm => fldsm
		this%flxsm => flxsm
		this%bcm => bcm
		this%ppm => ppm

		this%settings => settings
	
		this%ptime_step => ptime_step
		
		!=== VARIABLES INDICES: ===
		this%QIDX = QIDX
		this%UIDX = UIDX
		this%GRADQIDX = GRADQIDX
		this%FLUXESQIDX = FLUXESQIDX
		
		!=== VARIABLES POINTERS: ===
		this%vars_ptr => fldsm%registry(this%QIDX)%values
		this%flxs_ptr => fldsm%registry(this%FLUXESQIDX)%values	
		if (UIDX > 0) this%cvars_ptr => fldsm%registry(this%UIDX)%values
		
		!=== LINEAR SYSTEM INITIALIZATION: ===
		allocate(this%rhs(mesh%ncells*(mesh%dim + 2)),&
				 this%phi(mesh%ncells*(mesh%dim + 2)))
		call this%jacobianFree%initialize(mesh, fldsm, flxsm, bcm, ppm, gradient, ptime_step,&
										  QIDX, UIDX, FLUXESQIDX, this%settings%USE_GHOST_CELLS,&
										  this%settings%USE_CONSERVATIVE_VARS)
		
		!=== LINEAR SOLVER INITIALIZATION: ===		 
		this%linear_solver => linear_solver_factory%create(this%ls_filename)
		call this%linear_solver%set_operator(this%jacobianFree)		 
		
		!=== PRECONDITIONER SETUP ===
		if (associated(this%linear_solver%M)) then
			!call this%linear_solver%M%setup(mesh=this%mesh)			!=== EXPLICIT SETUP ===
			
			!=== IMPLICIT MATRIX ASSOCIATING: ===
			select type(prec => this%linear_solver%M)
			type is (bcsr_ilu0_preconditioner_t)
				this%precond => prec%LU
				allocate(prec%iw(mesh%ncells))
				!=== PRECOND MATRIX MANAGER INITIALIZATION: ===
				call this%mtrxm%initialize(mesh, fldsm, bcm, ppm, settings, this%precond,&
										  ptime_step, QIDX, GRADQIDX)
				call this%mtrxm%initialize_bcsr_mtrx()
			type is (blusgs_preconditioner_t)
				prec%ppm => ppm
				prec%settings => settings
				prec%Q => fldsm%registry(this%QIDX)%values
				prec%ptime_step => ptime_step
				call this%linear_solver%M%setup(mesh=this%mesh)	
			
			type is (gmg_preconditioner_t)
				call this%linear_solver%M%setup(mesh=this%mesh)
				call prec%gmg_hierarchy%initialize_db(this%jacobianFree, this%flxsm, this%fldsm, this%bcm, this%ppm, this%settings,&
													  gradient, this%ptime_step,&
													  this%QIDX, this%GRADQIDX, this%FLUXESQIDX, this%UIDX, 'data\input\B_boundary_conditions.txt')
				
				
			end select
		end if
		
		!=== LOG DATA SETTINGS: ===
		write(dim_str, '(I0)') this%mesh%dim + 2 
		this%dynamic_fmt_header = '(''  Iter'', ' // trim(dim_str) // 'A15, ''    L-Iter'', ''     L-Resid'', ''   Time(s)'')'
		this%dynamic_fmt_data   = '(I8, ' // trim(dim_str) // 'ES15.7, I10, ES15.7, F10.3)'
		if (this%mesh%dim == 2) then
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Energy'
		else
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Mom-Z', 'Energy'
		end if
	end subroutine
	
	subroutine db_jacfree_iterate(this, iter, time_start)
		implicit none
		class(db_jacfree_iterator_t), intent(inout), target :: this
		integer, intent(in) :: iter
		real(8), intent(in) :: time_start
		
		real(8), pointer, contiguous :: vars_ptr(:, :), flxs_ptr(:, :), cvars_ptr(:, :)
		real(8) :: time_end					
		integer :: i, d, off
		integer :: dim, ncells
		
		dim = this%mesh%dim
		ncells = this%mesh%ncells
		
		vars_ptr => this%vars_ptr
		flxs_ptr => this%flxs_ptr
		cvars_ptr => this%cvars_ptr
		
		!=== RIGHT HAND SIDE INITIALIZATION: ===
		do i = 1, this%mesh%ncells
			off = (i-1)*(dim+2)
			do d = 1, dim+2
				this%rhs(off + d) = -flxs_ptr(d, i)
			end do
		end do
											 
		!=== PRECONDITIONER UPDATE: ===
		if (associated(this%linear_solver%M)) then
			if (iter < 20 .or. mod(iter, this%settings%PAS_NITER) == 0) then
				call this%mtrxm%db_assmeble_bcsr_mtrx()
				call this%linear_solver%M%update()
			end if
		end if
			
		!=== SLAE SOLUTION: ===
		this%phi = 0.d0
		call this%linear_solver%solve(this%rhs, this%phi)
		
		!=== UPDATING SOLUTION: ===
		if (.not. this%settings%USE_CONSERVATIVE_VARS) then
			!primitive variables:
			do i = 1, ncells
				off = (i-1)*(dim+2)
				vars_ptr(:, i) = vars_ptr(:, i) + this%phi(1 + off:dim+2 + off)
			end do
		else
			!conservative variables:
			do i = 1, ncells
				off = (i-1)*(dim+2)
				cvars_ptr(:, i) = cvars_ptr(:, i) + this%phi(1 + off:dim+2 + off)
			end do
		end if	
		
		!=== CPU TIME: ===
		call cpu_time(time_end)
	
		!=== LOG OUTPUT: ===
		if (this%settings%IS_PRINT .or. this%settings%IS_LOG) then
			if (mod(iter, this%settings%PRINT_NITER) == 0 .or. iter == 1) then
				write(this%log_data, this%dynamic_fmt_data) iter,&
					(maxval(dabs(flxs_ptr(i, 1:ncells))), i = 1, dim + 2),&
					this%linear_solver%iterations,&
					this%linear_solver%residual_drop,&
					time_end - time_start
			end if
		end if
		
    end subroutine


!=======================================================================
!=============== DENSITY-BASED ALGH DIFF ITERATOR METHODS ==============
!=======================================================================
	subroutine db_algdiff_initialize(this, mesh, fldsm, flxsm, bcm, ppm, gradient, settings, ptime_step,&
								     QIDX, UIDX, GRADQIDX, FLUXESQIDX)
        implicit none
        class(db_algdiff_iterator_t), intent(inout) :: this
        type(mesh_t), intent(in), target :: mesh
        type(field_manager_t), intent(in), target :: fldsm
        type(flux_manager_t), intent(in), target :: flxsm
        type(bc_manager_t), intent(in), target :: bcm
        type(phys_prop_manager_t), intent(in), target :: ppm
        class(gradient_t), intent(in), target :: gradient
        type(db_solver_config_t), intent(in), target :: settings
        real(8), intent(in), target :: ptime_step(:)
        integer, intent(in) :: QIDX, UIDX, GRADQIDX, FLUXESQIDX
		
		character(len=20)  :: dim_str
		
		!=== POINTERS: ===
		this%mesh => mesh
		this%fldsm => fldsm
		this%flxsm => flxsm
		this%bcm => bcm
		this%ppm => ppm

		this%settings => settings
	
		this%ptime_step => ptime_step
		
		!=== VARIABLES INDICES: ===
		this%QIDX = QIDX
		this%UIDX = UIDX
		this%GRADQIDX = GRADQIDX
		this%FLUXESQIDX = FLUXESQIDX
		
		!=== VARIABLES POINTERS: ===
		this%vars_ptr => fldsm%registry(this%QIDX)%values
		this%flxs_ptr => fldsm%registry(this%FLUXESQIDX)%values	
		if (UIDX > 0) this%cvars_ptr => fldsm%registry(this%UIDX)%values
		
		!=== LINEAR SYSTEM INITIALIZATION: ===
		allocate(this%rhs(mesh%ncells*(mesh%dim + 2)),&
				 this%phi(mesh%ncells*(mesh%dim + 2)))
		call this%algDiff%initialize(mesh, fldsm, flxsm, bcm, ppm, gradient, ptime_step,&
									 QIDX, UIDX, FLUXESQIDX, this%settings%USE_GHOST_CELLS,&
									 this%settings%USE_CONSERVATIVE_VARS)
		
		!=== LINEAR SOLVER INITIALIZATION: ===		 
		this%linear_solver => linear_solver_factory%create(this%ls_filename)
		call this%linear_solver%set_operator(this%algDiff)		 
		
		!=== PRECONDITIONER SETUP ===
		if (associated(this%linear_solver%M)) then
			!call this%linear_solver%M%setup(mesh=this%mesh)			!=== EXPLICIT SETUP ===
			
			!=== IMPLICIT MATRIX ASSOCIATING: ===
			select type(prec => this%linear_solver%M)
			type is (bcsr_ilu0_preconditioner_t)
				this%precond => prec%LU
				allocate(prec%iw(mesh%ncells))
				!=== PRECOND MATRIX MANAGER INITIALIZATION: ===
				call this%mtrxm%initialize(mesh, fldsm, bcm, ppm, settings, this%precond,&
										  ptime_step, QIDX, GRADQIDX)
				call this%mtrxm%initialize_bcsr_mtrx()
			
			type is (blusgs_preconditioner_t)
				prec%ppm => ppm
				prec%settings => settings
				prec%Q => fldsm%registry(this%QIDX)%values
				prec%ptime_step => ptime_step
				call this%linear_solver%M%setup(mesh=this%mesh)
				
			type is (gmg_preconditioner_t)
				call this%linear_solver%M%setup(mesh=this%mesh)
				call prec%gmg_hierarchy%initialize_db(this%algDiff, this%flxsm, this%fldsm, this%bcm, this%ppm, this%settings,&
													  gradient, this%ptime_step,&
													  this%QIDX, this%GRADQIDX, this%FLUXESQIDX, this%UIDX, 'data\input\B_boundary_conditions.txt')
				
				
			end select
		end if
		
		!=== LOG DATA SETTINGS: ===
		write(dim_str, '(I0)') this%mesh%dim + 2 
		this%dynamic_fmt_header = '(''  Iter'', ' // trim(dim_str) // 'A15, ''    L-Iter'', ''     L-Resid'', ''   Time(s)'')'
		this%dynamic_fmt_data   = '(I8, ' // trim(dim_str) // 'ES15.7, I10, ES15.7, F10.3)'
		if (this%mesh%dim == 2) then
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Energy'
		else
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Mom-Z', 'Energy'
		end if
	end subroutine
	
	subroutine db_algdiff_iterate(this, iter, time_start)
		implicit none
		class(db_algdiff_iterator_t), intent(inout), target :: this
		integer, intent(in) :: iter
		real(8), intent(in) :: time_start
		
		real(8), pointer, contiguous :: vars_ptr(:, :), flxs_ptr(:, :), cvars_ptr(:, :)
		real(8) :: time_end					
		integer :: i, d, off
		integer :: dim, ncells
		
		dim = this%mesh%dim
		ncells = this%mesh%ncells
		
		vars_ptr => this%vars_ptr
		flxs_ptr => this%flxs_ptr
		cvars_ptr => this%cvars_ptr
		
		!=== RIGHT HAND SIDE INITIALIZATION: ===
		do i = 1, this%mesh%ncells
			off = (i-1)*(dim+2)
			do d = 1, dim+2
				this%rhs(off + d) = -flxs_ptr(d, i)
			end do
		end do
											 
		!=== PRECONDITIONER UPDATE: ===
		if (associated(this%linear_solver%M)) then
			if (iter < 20 .or. mod(iter, this%settings%PAS_NITER) == 0) then
				call this%mtrxm%db_assmeble_bcsr_mtrx()
				call this%linear_solver%M%update()
			end if
		end if
			
		!=== SLAE SOLUTION: ===
		this%phi = 0.d0
		call this%linear_solver%solve(this%rhs, this%phi)
		
		!=== UPDATING SOLUTION: ===
		if (.not. this%settings%USE_CONSERVATIVE_VARS) then
			!primitive variables:
			do i = 1, ncells
				off = (i-1)*(dim+2)
				vars_ptr(:, i) = vars_ptr(:, i) + this%phi(1 + off:dim+2 + off)
			end do
		else
			!conservative variables:
			do i = 1, ncells
				off = (i-1)*(dim+2)
				cvars_ptr(:, i) = cvars_ptr(:, i) + this%phi(1 + off:dim+2 + off)
			end do
		end if	
		
		!=== CPU TIME: ===
		call cpu_time(time_end)
	
		!=== LOG OUTPUT: ===
		if (this%settings%IS_PRINT .or. this%settings%IS_LOG) then
			if (mod(iter, this%settings%PRINT_NITER) == 0 .or. iter == 1) then
				write(this%log_data, this%dynamic_fmt_data) iter,&
					(maxval(dabs(flxs_ptr(i, 1:ncells))), i = 1, dim + 2),&
					this%linear_solver%iterations,&
					this%linear_solver%residual_drop,&
					time_end - time_start
			end if
		end if
		
    end subroutine


!=======================================================================
!=============== DENSITY-BASED ALGH DIFF ITERATOR METHODS ==============
!=======================================================================
	subroutine db_blusgs_initialize(this, mesh, fldsm, flxsm, bcm, ppm, gradient, settings, ptime_step,&
								     QIDX, UIDX, GRADQIDX, FLUXESQIDX)
        implicit none
        class(db_blusgs_iterator_t), intent(inout) :: this
        type(mesh_t), intent(in), target :: mesh
        type(field_manager_t), intent(in), target :: fldsm
        type(flux_manager_t), intent(in), target :: flxsm
        type(bc_manager_t), intent(in), target :: bcm
        type(phys_prop_manager_t), intent(in), target :: ppm
        class(gradient_t), intent(in), target :: gradient
        type(db_solver_config_t), intent(in), target :: settings
        real(8), intent(in), target :: ptime_step(:)
        integer, intent(in) :: QIDX, UIDX, GRADQIDX, FLUXESQIDX
		
		character(len=20)  :: dim_str
		
		!=== POINTERS: ===
		this%mesh => mesh
		this%fldsm => fldsm
		this%flxsm => flxsm
		this%bcm => bcm
		this%ppm => ppm

		this%settings => settings
	
		this%ptime_step => ptime_step
		
		!=== VARIABLES INDICES: ===
		this%QIDX = QIDX
		this%UIDX = UIDX
		this%GRADQIDX = GRADQIDX
		this%FLUXESQIDX = FLUXESQIDX
		
		!=== VARIABLES POINTERS: ===
		this%vars_ptr => fldsm%registry(this%QIDX)%values
		this%flxs_ptr => fldsm%registry(this%FLUXESQIDX)%values	
		if (UIDX > 0) this%cvars_ptr => fldsm%registry(this%UIDX)%values
		
		!=== LINEAR SYSTEM INITIALIZATION: ===
		allocate(this%phi((mesh%dim+2)*mesh%ncells),&
				 this%sprf(mesh%nfaces),&
				 this%J_f(mesh%dim+2, mesh%dim+2, mesh%nfaces),&
				 this%J_d(mesh%dim+2, mesh%dim+2, mesh%ncells))
		
		
		!=== LOG DATA SETTINGS: ===
		write(dim_str, '(I0)') this%mesh%dim + 2 
		this%dynamic_fmt_header = '(''  Iter'', ' // trim(dim_str) // 'A15, ''   Time(s)'')'
		this%dynamic_fmt_data   = '(I8, ' // trim(dim_str) // 'ES15.7, F10.3)'
		if (this%mesh%dim == 2) then
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Energy'
		else
			write(this%log_header, this%dynamic_fmt_header) 'Mass', 'Mom-X', 'Mom-Y', 'Mom-Z', 'Energy'
		end if
	end subroutine
	
	subroutine db_blusgs_iterate(this, iter, time_start)
		implicit none
		class(db_blusgs_iterator_t), intent(inout), target :: this
		integer, intent(in) :: iter
		real(8), intent(in) :: time_start
		
		real(8), pointer, contiguous :: vars_ptr(:, :), flxs_ptr(:, :), cvars_ptr(:, :)
		real(8) :: time_end					
		integer :: i, d, off
		integer :: dim, ncells
		
		dim = this%mesh%dim
		ncells = this%mesh%ncells
		
		vars_ptr => this%vars_ptr
		flxs_ptr => this%flxs_ptr
		cvars_ptr => this%cvars_ptr
		
											 
		!=== SYSTEM MATRIX ASSEMBLING: ===
		if (iter < 20 .or. mod(iter, this%settings%MAS_NITER) == 0) then
			call blusgs_setup(this%mesh%dim, this%mesh%ncells,&
							  this%mesh%nfaces, this%mesh%nbfaces,&
							  this%mesh%face_left_cell, this%mesh%face_right_cell,&
							  this%mesh%cell_faces, this%mesh%cell_faces_ptr,&
							  this%mesh%face_area, this%mesh%face_normal,&
							  this%mesh%face_center, this%mesh%cell_center,&
							  this%mesh%cell_volume, this%ptime_step,&
							  this%ppm%k, this%ppm%R_gas, this%ppm%cv,&
							  this%ppm%cp, this%ppm%Pr, this%w,&
							  vars_ptr, this%sprf, this%J_f, this%J_d,&
							  this%settings%MODEL,&
							  this%settings%USE_GHOST_CELLS,&
							  this%settings%USE_CONSERVATIVE_VARS)
		end if
			
		!=== SLAE SOLUTION: ===
		this%phi = 0.d0
		call blusgs_apply(this%mesh%dim, this%mesh%ncells, this%w, this%fluxes_sign,&
						  this%mesh%face_left_cell, this%mesh%face_right_cell,&
						  this%mesh%cell_faces, this%mesh%cell_faces_ptr,&
						  this%phi, flxs_ptr, this%sprf, this%J_f, this%J_d)
		
		!=== UPDATING SOLUTION: ===
		if (.not. this%settings%USE_CONSERVATIVE_VARS) then
			!primitive variables:
			do i = 1, ncells
				off = (i-1)*(dim+2)
				vars_ptr(:, i) = vars_ptr(:, i) + this%phi(1 + off:dim+2 + off)
			end do
		else
			!conservative variables:
			do i = 1, ncells
				off = (i-1)*(dim+2)
				cvars_ptr(:, i) = cvars_ptr(:, i) + this%phi(1 + off:dim+2 + off)
			end do
		end if	
		
		!=== CPU TIME: ===
		call cpu_time(time_end)
	
		!=== LOG OUTPUT: ===
		if (this%settings%IS_PRINT .or. this%settings%IS_LOG) then
			if (mod(iter, this%settings%PRINT_NITER) == 0 .or. iter == 1) then
				write(this%log_data, this%dynamic_fmt_data) iter,&
					(maxval(dabs(flxs_ptr(i, 1:ncells))), i = 1, dim + 2),&
					time_end - time_start
			end if
		end if
		
    end subroutine

end module
