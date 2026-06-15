module linear_operator_module
use linear_operator_dense_module
use linear_operator_csr_module
use linear_operator_bcsr_module

!for matrix free:
use field_manager_module, only: field_manager_t, LOC_CELL
use fluxes_manager_module, only: flux_manager_t, CONS_TO_PRIM_D
use bc_manager_module, only: bc_manager_t
use physical_properties_manager_module, only: phys_prop_manager_t
use diff_operator_module, only: gradient_t
use mesh_module, only: mesh_t
implicit none
!=======================================================================
!================= LINEAR OPERATOR BASIC DATA STRUCTURE ================
!=======================================================================
type, abstract :: linear_operator_t
	integer :: nrows, ncols
	character(len=64) :: name = 'LinearOperator'
	logical :: is_matrix_free = .false.
contains
	procedure(apply_interface), deferred :: apply
	procedure(get_diagonal_interface), deferred :: get_diagonal
	procedure(update_interface), deferred :: update
end type
    
abstract interface
	subroutine apply_interface(this, x, y)
		import :: linear_operator_t
		class(linear_operator_t), intent(in) :: this
		real(8), intent(in), contiguous :: x(:)
		real(8), intent(inout), contiguous :: y(:)
	end subroutine
	
	subroutine get_diagonal_interface(this, diag)
		import :: linear_operator_t
		class(linear_operator_t), intent(in) :: this
		real(8), intent(inout), contiguous :: diag(:)
	end subroutine
	
	subroutine update_interface(this, params)
		import :: linear_operator_t
		class(linear_operator_t), intent(inout) :: this
		real(8), intent(in), optional :: params(:)
	end subroutine
end interface
    
 

!=======================================================================
!==================== DEFAULT MATRIX DATA STRUCTURE ====================   
!=======================================================================
type, extends(linear_operator_t) :: dense_operator_t
	real(8), allocatable :: matrix(:,:)
contains
	procedure :: apply => dense_apply
	procedure :: get_diagonal => dense_get_diagonal
	procedure :: update => dense_update
end type


!=======================================================================
!======================= CSR MATRIX DATA STRUCTURE =====================   
!=======================================================================
type, extends(linear_operator_t) :: csr_operator_t
	integer :: n														!=== NUMBER OF ROWS ===
	integer :: nnz														!=== NUMBER OF NONZERO ELEMENTS ===
	integer :: bs														!=== BLOCK SIZE (if necesery) ===	
	real(8), allocatable :: values(:)									!=== DIM = (1:nnz) ===
	integer, allocatable :: col_indices(:)								!=== DIM = (1:nnz) ===
	integer, allocatable :: diag_indices(:)								!=== DIM = (1:n) ===
	integer, allocatable :: row_ptr(:)									!=== DIM = (1:n + 1) ===
contains
	procedure :: apply => csr_apply
	procedure :: get_diagonal => csr_get_diagonal
	procedure :: update => csr_update
end type
    

!=======================================================================
!==================== BlockCSR MATRIX DATA STRUCTURE ===================
!=======================================================================
type, extends(linear_operator_t) :: bcsr_operator_t
	integer :: n														!=== NUMBER OF ROWS (<-> NUMBER OF CELLS) ===
	integer :: nnz														!=== NUMBER OF NONZERO BLOCKS ===
	integer :: bs														!=== BLOCK SIZE ===	
	real(8), allocatable :: values(:,:,:)								!=== DIM = (bs, bs, nnz) ===
	integer, allocatable :: col_indices(:)								!=== DIM = (1:nnz) ===
	integer, allocatable :: diag_indices(:)								!=== DIM = (1:n) ===
	integer, allocatable :: row_ptr(:)									!=== DIM = (1:n + 1) ===
contains
	procedure :: apply => bcsr_apply
	procedure :: get_diagonal => bcsr_get_diagonal
	procedure :: update => bcsr_update
end type
  
      
!=======================================================================
!================= MATRIX FREE1 OPERATOR DATA STRUCTURE ================
!=======================================================================
type, extends(linear_operator_t) :: jacfree_operator_t
	type(mesh_t), pointer :: mesh => null()
	type(field_manager_t), pointer :: fldsm => null()
	type(flux_manager_t), pointer :: flxsm => null()
	type(bc_manager_t), pointer :: bcm => null()
	type(phys_prop_manager_t), pointer :: ppm => null()
	class(gradient_t), pointer :: gradient => null()
	real(8), pointer, contiguous :: ptime_step(:) => null()
	
	real(8), contiguous, pointer :: Q(:, :), Q_eps(:, :),&
									gradQ_eps(:, :),&
									Fluxes(:, :), Fluxes_eps(:, :),&
									U(:, :), U_eps(:, :)
	integer :: QIDX, FLUXESQIDX, QIDX_EPS, GRADQIDX_EPS, FLUXESQIDX_EPS
	integer :: UIDX = -1, UIDX_EPS = -1
	
	logical :: USE_GHOST_CELLS = .false.
	logical :: USE_CONSERVATIVE_VARS = .false.
	
	contains
	procedure :: initialize => jacfree_operator_initialize
	procedure :: apply => jacfree_apply
	procedure :: get_diagonal => jacfree_get_diagonal
	procedure :: update => jacfree_update
end type

!=======================================================================
!================= MATRIX FREE2 OPERATOR DATA STRUCTURE ================
!=======================================================================
type, extends(linear_operator_t) :: algdiff_operator_t
	type(mesh_t), pointer :: mesh => null()
	type(field_manager_t), pointer :: fldsm => null()
	type(flux_manager_t), pointer :: flxsm => null()
	type(bc_manager_t), pointer :: bcm => null()
	type(phys_prop_manager_t), pointer :: ppm => null()
	class(gradient_t), pointer :: gradient => null()
	real(8), pointer, contiguous :: ptime_step(:) => null()
	
	real(8), contiguous, pointer :: Q(:, :), Q_eps(:, :),&
									gradQ_eps(:, :),&
									Fluxes(:, :), Fluxes_eps(:, :),&
									U(:, :), U_eps(:, :)
	integer :: QIDX, FLUXESQIDX, QIDX_EPS, GRADQIDX_EPS, FLUXESQIDX_EPS
	integer :: UIDX = -1, UIDX_EPS = -1
	
	logical :: USE_GHOST_CELLS = .false.
	logical :: USE_CONSERVATIVE_VARS = .false.
	
	contains
	procedure :: initialize => algdiff_operator_initialize
	procedure :: apply => algdiff_apply
	procedure :: get_diagonal => algdiff_get_diagonal
	procedure :: update => algdiff_update
end type
        
contains
!=======================================================================
!======================== DEFAULT MATRIX METHODS =======================   
!=======================================================================
	subroutine dense_apply(this, x, y)
		implicit none
		class(dense_operator_t), intent(in) :: this
		real(8), intent(in), contiguous :: x(:)
		real(8), intent(inout), contiguous :: y(:)
		
		call dense_mtrx_matvec(size(this%matrix, dim=1), this%matrix, x, y)
	end subroutine

	subroutine dense_get_diagonal(this, diag)
		implicit none
		class(dense_operator_t), intent(in) :: this
		real(8), intent(inout), contiguous :: diag(:)
		
		call dense_mtrx_get_diagonal(size(this%matrix, dim=1), this%matrix, diag)
	end subroutine

	subroutine dense_update(this, params)
		implicit none
		class(dense_operator_t), intent(inout) :: this
		real(8), intent(in), optional :: params(:)
		
		print*, 'Not done yet :('
	end subroutine
	

!=======================================================================
!========================== CSR MATRIX METHODS =========================   
!=======================================================================
	subroutine csr_apply(this, x, y)
		implicit none
		class(csr_operator_t), intent(in) :: this
		real(8), intent(in), contiguous :: x(:)
		real(8), intent(inout), contiguous :: y(:)
		
		call csr_mtrx_matvec(this%n,&
							 this%values,&
							 this%col_indices,&
							 this%row_ptr, x, y)
	end subroutine

	subroutine csr_get_diagonal(this, diag)
		implicit none
		class(csr_operator_t), intent(in) :: this
		real(8), intent(inout), contiguous :: diag(:)
		
		call csr_mtrx_get_diagonal(this%n,&
								   this%values,&
								   this%diag_indices,&
								   diag)
	end subroutine

	subroutine csr_update(this, params)
		implicit none
		class(csr_operator_t), intent(inout) :: this
		real(8), intent(in), optional :: params(:)
		
		print*, 'Not done yet :('
	end subroutine
	

!=======================================================================
!======================= BlockCSR MATRIX METHODS =======================
!=======================================================================
	subroutine bcsr_apply(this, x, y)
		implicit none
		class(bcsr_operator_t), intent(in) :: this
		real(8), intent(in), contiguous :: x(:)
		real(8), intent(inout), contiguous :: y(:)
		
		call bcsr_mtrx_matvec(this%n, this%bs,&
							  this%values,&
							  this%col_indices,&
							  this%row_ptr, x, y)
	end subroutine

	subroutine bcsr_get_diagonal(this, diag)
		implicit none
		class(bcsr_operator_t), intent(in) :: this
		real(8), intent(inout), contiguous :: diag(:)
		
		call bcsr_mtrx_get_main_diagonal(this%n, this%bs,&
										 this%values,&
										 this%diag_indices, diag)
	end subroutine

	subroutine bcsr_update(this, params)
		implicit none
		class(bcsr_operator_t), intent(inout) :: this
		real(8), intent(in), optional :: params(:)
		
		print*, 'Not done yet :('
	end subroutine
	

!=======================================================================
!===================== MATRIX FREE1 OPERATOR METHODS ===================
!=======================================================================
	subroutine jacfree_operator_initialize(this, mesh, fldsm, flxsm, bcm, ppm, gradient, ptime_step,&
										   QIDX, UIDX, FLUXESQIDX, USE_GHOST_CELLS, USE_CONSERVATIVE_VARS)
		implicit none
		class(jacfree_operator_t), intent(inout) :: this
		type(mesh_t), intent(in), target :: mesh
		type(field_manager_t), intent(in), target :: fldsm
		type(flux_manager_t), intent(in), target :: flxsm
		type(bc_manager_t), intent(in), target :: bcm
		type(phys_prop_manager_t), intent(in), target :: ppm
		class(gradient_t), intent(in), target :: gradient
		real(8), intent(in), target :: ptime_step(:)
		integer, intent(in) :: QIDX, UIDX, FLUXESQIDX
		logical, intent(in) :: USE_GHOST_CELLS, USE_CONSERVATIVE_VARS
		
		this%USE_GHOST_CELLS = USE_GHOST_CELLS
		this%USE_CONSERVATIVE_VARS = USE_CONSERVATIVE_VARS
		this%QIDX = QIDX
		this%UIDX = UIDX
		this%FLUXESQIDX = FLUXESQIDX
		
		this%mesh => mesh
		this%fldsm => fldsm
		this%flxsm => flxsm
		this%bcm => bcm
		this%ppm => ppm
		this%gradient => gradient
		this%ptime_step => ptime_step
		
		!=== EPS PRIMITIVE VARIABLES [P U V (W) T]: ===
		call this%fldsm%add_field('Q_eps',  .false., LOC_CELL, mesh%dim+2, mesh%ncells, mesh%nbfaces)
		this%QIDX_EPS = this%fldsm%get_idx('Q_eps')
		!=== EPS PRIMITIVE VARIABLES GRADIENTS [gradP_x, gradP_y, ... , gradT_x, gradT_y, (gradT_z)]: ===
		call this%fldsm%add_field('gradQ_eps',  .false., LOC_CELL, mesh%dim*(mesh%dim+2), mesh%ncells, mesh%nbfaces)
		this%GRADQIDX_EPS = this%fldsm%get_idx('gradQ_eps')
		!=== CONNECTING VARIABLES WITH GRADIENTS: ===
		call this%fldsm%registry(this%QIDX_EPS)%set_gradient(this%fldsm%registry(this%GRADQIDX_EPS)%values)
		!=== EPS FLUXES: ===
		call this%fldsm%add_field('fluxesQ_eps',  .false., LOC_CELL, mesh%dim+2, mesh%ncells, mesh%nbfaces)
		this%FLUXESQIDX_EPS = this%fldsm%get_idx('fluxesQ_eps')
		
		!=== EPS CONSERVATIVE VARIABLES: ===
		if (USE_CONSERVATIVE_VARS) then
			call this%fldsm%add_field('U_eps',  .false., LOC_CELL, mesh%dim+2, mesh%ncells, mesh%nbfaces)
			this%UIDX_EPS = this%fldsm%get_idx('U_eps')
			
			this%U => fldsm%registry(UIDX)%values
			this%U_eps => fldsm%registry(this%UIDX_EPS)%values			
		end if
				
		this%Q => fldsm%registry(QIDX)%values
		this%Fluxes => fldsm%registry(FLUXESQIDX)%values
		this%Q_eps => fldsm%registry(this%QIDX_EPS)%values
		this%gradQ_eps => fldsm%registry(this%QIDX_EPS)%grad
		this%Fluxes_eps => fldsm%registry(this%FLUXESQIDX_EPS)%values
	end subroutine

	subroutine jacfree_apply(this, x, y)
		implicit none
		class(jacfree_operator_t), intent(in) :: this
		real(8), intent(in), contiguous :: x(:)
		real(8), intent(inout), contiguous :: y(:)
		
		real(8) :: eps, eps_inv, diag, x_norm
		integer :: cell_idx, offset, ncells, dim, nvars
		integer :: d
		
		ncells = this%mesh%ncells
		dim = this%mesh%dim
		nvars = dim + 2
		
		x_norm = norm2(x)
		if (x_norm < 1e-12) then
			y = 0.d0
			return
		end if
		
		eps = 1e-8*(1.0d0 + norm2(this%Q(:, 1:ncells)))/x_norm
		eps_inv = 1.d0/eps
		
		!=== VARIABLES VARIATION: ===
		if (.not. this%USE_CONSERVATIVE_VARS) then
			!primitive variables:
			do cell_idx = 1, ncells
				offset = nvars*(cell_idx - 1)
				this%Q_eps(:, cell_idx) = this%Q(:, cell_idx) + eps*x(1+offset:nvars+offset)
			end do
		else
			!conservative variables: 
			do cell_idx = 1, ncells
				offset = nvars*(cell_idx - 1)
				this%U_eps(:, cell_idx) = this%U(:, cell_idx) + eps*x(1+offset:nvars+offset)
			end do
			
			!recomputing primitive variables:
			do cell_idx = 1, ncells
				!velocity:
				do d = 2, dim + 1
					this%Q_eps(d, cell_idx) = this%U_eps(d, cell_idx)/this%U_eps(1, cell_idx)
				end do
				
				!temperature:
				this%Q_eps(dim+2, cell_idx) = 1.d0/(this%ppm%cv)*(this%U_eps(dim+2, cell_idx)/this%U_eps(1, cell_idx)&
											  - 0.5d0*dot_product(this%Q_eps(2:dim+1, cell_idx), this%Q_eps(2:dim+1, cell_idx)))
				
				!pressure:
				this%Q_eps(1, cell_idx) = this%U_eps(1, cell_idx)*this%Q_eps(dim+2, cell_idx)*this%ppm%R_gas
			end do
		end if
		
				
		!=== UPDATING BOUNDARY FACES/GHOST CELLS VALUES: ===
		call this%bcm%update_boundary_values('Q', .true., this%USE_GHOST_CELLS, this%QIDX_EPS)
		
		!=== COMPUTING GRADIENTS: ===
		call this%gradient%apply_vector_cellfield(this%Q_eps, this%gradQ_eps,&
												  nvars, this%USE_GHOST_CELLS)
												  
		!=== UPDATING BOUNDARY FACES/GHOST CELLS GRADIENTS: ===
		call this%bcm%update_boundary_values('Q', .true., this%USE_GHOST_CELLS, this%QIDX_EPS)
		
		!=== COMPUTING FLUXES: ===
		call this%flxsm%compute_fluxes_cmprs(this%QIDX_EPS, this%FLUXESQIDX_EPS)
		
		!=== MATRIX-VECTOR PORDUCT: ===
		do cell_idx = 1, ncells
			offset = nvars*(cell_idx - 1)
			diag = this%mesh%cell_volume(cell_idx)/this%ptime_step(cell_idx)
			
			y(1+offset:nvars+offset) = diag*x(1+offset:nvars+offset) +&
									   (this%Fluxes_eps(:, cell_idx) - this%Fluxes(:, cell_idx))*eps_inv
		end do
		
	end subroutine

	subroutine jacfree_get_diagonal(this, diag)
		implicit none
		class(jacfree_operator_t), intent(in) :: this
		real(8), intent(inout), contiguous :: diag(:)
		
		print*, 'Not done yet :('
	end subroutine

	subroutine jacfree_update(this, params)
		implicit none
		class(jacfree_operator_t), intent(inout) :: this
		real(8), intent(in), optional :: params(:)
		
		print*, 'Not done yet :('
	end subroutine


!=======================================================================
!===================== MATRIX FREE2 OPERATOR METHODS ===================
!=======================================================================
	subroutine algdiff_operator_initialize(this, mesh, fldsm, flxsm, bcm, ppm, gradient, ptime_step,&
										   QIDX, UIDX, FLUXESQIDX, USE_GHOST_CELLS, USE_CONSERVATIVE_VARS)
		implicit none
		class(algdiff_operator_t), intent(inout) :: this
		type(mesh_t), intent(in), target :: mesh
		type(field_manager_t), intent(in), target :: fldsm
		type(flux_manager_t), intent(in), target :: flxsm
		type(bc_manager_t), intent(in), target :: bcm
		type(phys_prop_manager_t), intent(in), target :: ppm
		class(gradient_t), intent(in), target :: gradient
		real(8), intent(in), target :: ptime_step(:)
		integer, intent(in) :: QIDX, UIDX, FLUXESQIDX
		logical, intent(in) :: USE_GHOST_CELLS, USE_CONSERVATIVE_VARS
		
		this%USE_GHOST_CELLS = USE_GHOST_CELLS
		this%USE_CONSERVATIVE_VARS = USE_CONSERVATIVE_VARS
		this%QIDX = QIDX
		this%UIDX = UIDX
		this%FLUXESQIDX = FLUXESQIDX
		
		this%mesh => mesh
		this%fldsm => fldsm
		this%flxsm => flxsm
		this%bcm => bcm
		this%ppm => ppm
		this%gradient => gradient
		this%ptime_step => ptime_step
		
		!=== EPS PRIMITIVE VARIABLES [P U V (W) T]: ===
		call this%fldsm%add_field('Q_eps',  .false., LOC_CELL, mesh%dim+2, mesh%ncells, mesh%nbfaces)
		this%QIDX_EPS = this%fldsm%get_idx('Q_eps')
		!=== EPS PRIMITIVE VARIABLES GRADIENTS [gradP_x, gradP_y, ... , gradT_x, gradT_y, (gradT_z)]: ===
		call this%fldsm%add_field('gradQ_eps',  .false., LOC_CELL, mesh%dim*(mesh%dim+2), mesh%ncells, mesh%nbfaces)
		this%GRADQIDX_EPS = this%fldsm%get_idx('gradQ_eps')
		!=== CONNECTING VARIABLES WITH GRADIENTS: ===
		call this%fldsm%registry(this%QIDX_EPS)%set_gradient(this%fldsm%registry(this%GRADQIDX_EPS)%values)
		!=== EPS FLUXES: ===
		call this%fldsm%add_field('fluxesQ_eps',  .false., LOC_CELL, mesh%dim+2, mesh%ncells, mesh%nbfaces)
		this%FLUXESQIDX_EPS = this%fldsm%get_idx('fluxesQ_eps')
		
		!=== EPS CONSERVATIVE VARIABLES: ===
		if (USE_CONSERVATIVE_VARS) then
			call this%fldsm%add_field('U_eps',  .false., LOC_CELL, mesh%dim+2, mesh%ncells, mesh%nbfaces)
			this%UIDX_EPS = this%fldsm%get_idx('U_eps')
			
			this%U => fldsm%registry(UIDX)%values
			this%U_eps => fldsm%registry(this%UIDX_EPS)%values			
		end if
				
		this%Q => fldsm%registry(QIDX)%values
		this%Fluxes => fldsm%registry(FLUXESQIDX)%values
		this%Q_eps => fldsm%registry(this%QIDX_EPS)%values
		this%gradQ_eps => fldsm%registry(this%QIDX_EPS)%grad
		this%Fluxes_eps => fldsm%registry(this%FLUXESQIDX_EPS)%values
	end subroutine

	subroutine algdiff_apply(this, x, y)
		implicit none
		class(algdiff_operator_t), intent(in) :: this
		real(8), intent(in), contiguous :: x(:)
		real(8), intent(inout), contiguous :: y(:)
		
		real(8) :: diag, x_norm
		integer :: cell_idx, offset, ncells, dim, nvars
		integer :: d
		
		ncells = this%mesh%ncells
		dim = this%mesh%dim
		nvars = dim + 2
		
		x_norm = norm2(x)
		if (x_norm < 1e-12) then
			y = 0.d0
			return
		end if
		
		!=== VARIABLES VARIATION: ===
		if (.not. this%USE_CONSERVATIVE_VARS) then
			!primitive variables:
			do cell_idx = 1, ncells
				offset = nvars*(cell_idx - 1)
				this%Q_eps(:, cell_idx) = x(1+offset:nvars+offset)
			end do
		else
			!conservative variables: 
			do cell_idx = 1, ncells
				offset = nvars*(cell_idx - 1)
				this%U_eps(:, cell_idx) = x(1+offset:nvars+offset)
			end do
			
			!recomputing primitive variables:
			call CONS_TO_PRIM_D(dim, ncells, this%ppm%cv, this%ppm%R_gas, this%U, this%U_eps, this%Q, this%Q_eps)
		end if
		
				
		!=== UPDATING BOUNDARY FACES/GHOST CELLS VALUES: ===
		call this%bcm%ADupdate_boundary_values('Q', .true., this%USE_GHOST_CELLS, this%QIDX, this%QIDX_EPS)
		
		!=== COMPUTING GRADIENTS: ===
		call this%gradient%ADapply_vector_cellfield(this%Q_eps, this%gradQ_eps,&
												    nvars, this%USE_GHOST_CELLS)
												  
		!=== UPDATING BOUNDARY FACES/GHOST CELLS GRADIENTS: ===
		call this%bcm%ADupdate_boundary_values('Q', .true., this%USE_GHOST_CELLS, this%QIDX, this%QIDX_EPS)
		
		!=== COMPUTING FLUXES: ===
		call this%flxsm%compute_ADdiff_fluxes_cmprs(this%QIDX, this%QIDX_EPS, this%FLUXESQIDX_EPS)
		
		!=== MATRIX-VECTOR PORDUCT: ===
		do cell_idx = 1, ncells
			offset = nvars*(cell_idx - 1)
			diag = this%mesh%cell_volume(cell_idx)/this%ptime_step(cell_idx)
			
			y(1+offset:nvars+offset) = diag*x(1+offset:nvars+offset) +&
									   this%Fluxes_eps(:, cell_idx)
		end do
		
	end subroutine

	subroutine algdiff_get_diagonal(this, diag)
		implicit none
		class(algdiff_operator_t), intent(in) :: this
		real(8), intent(inout), contiguous :: diag(:)
		
		print*, 'Not done yet :('
	end subroutine

	subroutine algdiff_update(this, params)
		implicit none
		class(algdiff_operator_t), intent(inout) :: this
		real(8), intent(in), optional :: params(:)
		
		print*, 'Not done yet :('
	end subroutine
		
		
end module
