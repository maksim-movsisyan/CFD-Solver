module preconditioner_module
use mesh_module, only: mesh_t
use linear_operator_module
use preconditioner_block_jacobi_module
use preconditioner_csr_ilu0_module
use preconditioner_bcsr_ilu0_module
use string_utils_module

!for matrix free:
use matrix_initialization_module

!for blusgs:
use BLUSGS_module
use physical_properties_manager_module, only: phys_prop_manager_t
use solver_config_module, only: db_solver_config_t

!for gmg:
use mg_mesh_hierarchy_module
use mg_operators_module
use db_mtrx_assembling_module

implicit none

!=======================================================================
!================ PRECONDITIONER BASIC DATA STRUCTURE ==================
!=======================================================================
type, abstract :: preconditioner_t
	class(linear_operator_t), pointer :: A => null()
	type(mesh_t), pointer :: mesh => null()
	character(len=64) :: name = 'Undefined Preconditioner'
	
	logical :: is_setup = .false.
	logical :: is_updated = .false.
	logical :: needs_matrix = .true.   
	logical :: needs_mesh = .false.   
	 
	integer :: verbosity = -1
contains
	procedure(precon_setup_interface), deferred :: setup
	procedure(precon_apply_interface), deferred :: apply
	procedure(precon_update_interface), deferred :: update
	procedure(precon_destroy_interface), deferred :: destroy
end type
    
abstract interface
	subroutine precon_setup_interface(this, A, mesh, params)
		import :: preconditioner_t, mesh_t, linear_operator_t
		class(preconditioner_t), intent(inout) :: this
		class(linear_operator_t), target, optional, intent(in) :: A
		type(mesh_t), target, optional, intent(in) :: mesh
		character(len=*), optional, intent(in) :: params				!pramas = "key1=value1;key2=value2..."
	end subroutine
	
	subroutine precon_apply_interface(this, r, z)
		import :: preconditioner_t
		class(preconditioner_t), intent(inout) :: this
		real(8), intent(in), contiguous, target :: r(:)
		real(8), intent(inout), contiguous :: z(:)
	end subroutine
	
	subroutine precon_update_interface(this, A)
		import :: preconditioner_t, linear_operator_t
		class(preconditioner_t), intent(inout) :: this
		class(linear_operator_t), optional, intent(in) :: A
		
	end subroutine
	
	subroutine precon_destroy_interface(this)
		import :: preconditioner_t
		class(preconditioner_t), intent(inout) :: this
	end subroutine
end interface

!=======================================================================
!============== BLOCK JACOBI PRECONDITIONER DATA STRUCTURE =============
!=======================================================================
type, extends(preconditioner_t) :: block_jacobi_preconditioner_t
	integer :: n, bs													!=== MATRIX DIMENSION ===
	real(8), allocatable :: invDIAG(:, :, :)							!=== DIM = (bs, bs, n) ===
	real(8) :: omega = 1.0d0											!=== RELAXATION FACTOR ===

contains
	procedure :: setup => setup_block_jacobi
	procedure :: apply => apply_block_jacobi
	procedure :: update => update_block_jacobi
	procedure :: destroy => destroy_block_jacobi
	
end type

!=======================================================================
!============== CSR ILU(0) PRECONDITIONER DATA STRUCTURE ===============
!=======================================================================
type, extends(preconditioner_t) :: csr_ilu0_preconditioner_t
	type(csr_operator_t) :: LU											!=== LU MATRIX ===
	integer, allocatable :: iw(:)										!=== INTEGER WORK ARRAY ===
contains
	procedure :: setup => setup_csr_ilu0
	procedure :: apply => apply_csr_ilu0
	procedure :: update => update_csr_ilu0
	procedure :: destroy => destroy_csr_ilu0
	
end type

!=======================================================================
!============= BCSR ILU(0) PRECONDITIONER DATA STRUCTURE ===============
!=======================================================================
type, extends(preconditioner_t) :: bcsr_ilu0_preconditioner_t
	type(bcsr_operator_t) :: LU											!=== LU MATRIX ===
	integer, allocatable :: iw(:)										!=== INTEGER WORK ARRAY ===

contains
	procedure :: setup => setup_bcsr_ilu0
	procedure :: apply => apply_bcsr_ilu0
	procedure :: update => update_bcsr_ilu0
	procedure :: destroy => destroy_bcsr_ilu0
	
end type

!=======================================================================
!================ BLUSGS PRECONDITIONER DATA STRUCTURE =================
!=======================================================================
type, extends(preconditioner_t) :: blusgs_preconditioner_t
	type(phys_prop_manager_t), pointer :: ppm => null()
	type(db_solver_config_t), pointer :: settings => null()
	
	real(8), allocatable :: sprf(:)										!=== SPECTRAL RADIUS ON FACES, DIM = (nfaces) ===
	real(8), allocatable :: J_f(:, :, :)								!=== JACOBIANS ON FACES, DIM = (dim+2, dim+2, nfaces) ===
	real(8), allocatable :: J_d(:, :, :)								!=== DIAGONAL JACOBIANS, DIM = (dim+2, dim+2, ncells) ===
	
	real(8), pointer, contiguous :: Q(:, :) => null()					!=== DIM = (nvars, ncells) ===
	real(8), pointer, contiguous :: ptime_step(:) => null()				!=== DIM = (nvars, ncells) ===
	real(8) :: w = 0.5d0	
contains
	procedure :: setup => setup_blusgs_preconditioner
	procedure :: apply => apply_blusgs_preconditioner
	procedure :: update => update_blusgs_preconditioner
	procedure :: destroy => destroy_blusgs_preconditioner
	
end type


!=======================================================================
!======= GEOMETRIC MULTIGRID PRECONDITIONER DATA STRUCTURE =============
!=======================================================================
type, extends(preconditioner_t) :: gmg_preconditioner_t
	type(multilevel_hierarchy_t) :: gmg_hierarchy
contains
	procedure :: setup => setup_gmg_preconditioner
	procedure :: apply => apply_gmg_preconditioner
	procedure :: update => update_gmg_preconditioner
	procedure :: destroy => destroy_gmg_preconditioner
	
end type
        
contains
!=======================================================================
!================== BLOCK JACOBI PRECONDITIONER METHODS ================
!=======================================================================
	subroutine setup_block_jacobi(this, A, mesh, params)
		implicit none
        class(block_jacobi_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), target, optional, intent(in) :: A
        type(mesh_t), target, optional, intent(in):: mesh
        character(len=*), optional, intent(in) :: params				!=== params = [n=n_val;bs=bs_val] ===
        
        class(bcsr_operator_t), pointer :: bcsr_ptr
        logical :: found
		real(8) :: val_r
        
        call this%destroy()
        this%name = 'Block-Jacobi preconditioner'
        if (present(mesh)) this%mesh => mesh
        
        if (present(A)) then
			!=== MATRIX INITIALIZATION: ===
			select type(A)
			type is (bcsr_operator_t)
				bcsr_ptr => A
				this%A => bcsr_ptr
				
				this%n = bcsr_ptr%n
				this%bs = bcsr_ptr%bs
				
				allocate(this%invDIAG(this%bs, this%bs, this%n))
				this%is_setup = .true.	
			end select
        else if (present(params)) then
			!=== NON-MATRIX INITIALIZATION: ===
			call get_key_value_real(params, 'bs', val_r, found)
			if (found) this%bs = nint(val_r)

			call get_key_value_real(params, 'n', val_r, found)
			if (found) this%n = nint(val_r)

			call get_key_value_real(params, 'omega', val_r, found)
			if (found) this%omega = val_r
			
			allocate(this%invDIAG(this%bs, this%bs, this%n))
			this%is_setup = .true.
        end if
        
        if (.not. this%is_setup) then
			error stop 'Block-Jacobi precondtioner stting up error'
        end if
    end subroutine
    
    subroutine apply_block_jacobi(this, r, z)
		implicit none
        class(block_jacobi_preconditioner_t), intent(inout) :: this
        real(8), intent(in), contiguous, target :: r(:)
        real(8), intent(inout), contiguous :: z(:)
        
       call apply_block_jacobi_preconditioner(this%n, this%bs, this%invDIAG, r, z)
    end subroutine
    
    subroutine update_block_jacobi(this, A)
		implicit none
        class(block_jacobi_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), optional, intent(in) :: A
		
		class(bcsr_operator_t), pointer :: bcsr_ptr
		
        if (present(A)) then
            call this%destroy()
            call this%setup(A)
        end if
        
        if (associated(this%A)) then
			!=== MATRIX UPDATING, invDIAG CAN BE NON FILLED: ===
			select type (mtrx => this%A)
			type is (bcsr_operator_t)
				bcsr_ptr => mtrx
			end select
			call update_block_jacobi_preconditioner2(this%n, this%bs,&
													 bcsr_ptr%values,&
													 bcsr_ptr%diag_indices,&
													 this%invDIAG)
			this%is_updated = .true.
        else
			!=== MATRIX UPDATING, invDIAG MUST CONTAIN DIAGONAL MATRIX BLOCKS: ===
			call update_block_jacobi_preconditioner1(this%n, this%bs, this%invDIAG)
			this%is_updated = .true.
        end if
    end subroutine
    
    subroutine destroy_block_jacobi(this)
		implicit none
        class(block_jacobi_preconditioner_t), intent(inout) :: this
        
        nullify(this%A)
        nullify(this%mesh)
        this%name = 'Undefined Preconditioner'
        
        if (allocated(this%invDIAG)) deallocate(this%invDIAG)
    end subroutine


!=======================================================================
!================== CSR ILU(0) PRECONDITIONER METHODS ==================
!=======================================================================    
    subroutine setup_csr_ilu0(this, A, mesh, params)
		implicit none
        class(csr_ilu0_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), target, optional, intent(in) :: A
        type(mesh_t), target, optional, intent(in):: mesh
        character(len=*), optional, intent(in) :: params
        
        class(csr_operator_t), pointer :: csr_ptr
        		
        call this%destroy()
        this%name = 'CSR-ILU(0) preconditioner'
        if (present(mesh)) this%mesh => mesh
        
        if (present(A)) then
			!=== MATRIX INITIALIZATION: ===
			select type(A)
			type is (csr_operator_t)
				csr_ptr => A
				this%A => csr_ptr
				
				this%LU%n = csr_ptr%n
				this%LU%nnz = csr_ptr%nnz
				this%LU%bs = csr_ptr%bs
				
				allocate(this%LU%values(this%LU%nnz), this%LU%col_indices(this%LU%nnz),&
						 this%LU%diag_indices(this%LU%n), this%LU%row_ptr(this%LU%n + 1),&
						 this%iw(this%LU%n))
				
				this%LU%col_indices = csr_ptr%col_indices
				this%LU%diag_indices = csr_ptr%diag_indices
				this%LU%row_ptr = csr_ptr%row_ptr
			
				this%is_setup = .true.
			end select
        else if (present(params)) then
			!=== NON-MATRIX INITIALIZATION: ===
			!=== NOT DONE YET, MESH BASED INITIALIZATION TO BE DONE! ===
        end if
        
        if (.not. this%is_setup) then
			error stop 'CSR-ILU(0) precondtioner stting up error'
        end if
    end subroutine
    
    subroutine apply_csr_ilu0(this, r, z)
		implicit none
        class(csr_ilu0_preconditioner_t), intent(inout) :: this
        real(8), intent(in), contiguous, target :: r(:)
        real(8), intent(inout), contiguous :: z(:)
        
        call apply_csr_ilu0_preconditioner(this%LU%n,&
										   this%LU%values,&
										   this%LU%col_indices,&
										   this%LU%diag_indices,&
										   this%LU%row_ptr, r, z)
    end subroutine
    
    subroutine update_csr_ilu0(this, A)
		implicit none
        class(csr_ilu0_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), optional, intent(in) :: A
		
		class(csr_operator_t), pointer :: csr_ptr
		
        if (present(A)) then
            call this%destroy()
            call this%setup(A)
        end if
        
        if (associated(this%A)) then
			!=== MATRIX UPDATING, LU%VALUES CAN BE NON FILLED: ===
			select type (mtrx => this%A)
			type is (csr_operator_t)
				csr_ptr => mtrx
			end select
			this%LU%values = csr_ptr%values
			call update_csr_ilu0_preconditioner(this%LU%n,&
											    this%LU%values,&
											    this%LU%col_indices,&
												this%LU%diag_indices,&
												this%LU%row_ptr, this%iw)
			this%is_updated = .true.
        else
			!=== MATRIX UPDATING, LU%VALUES MUST CONTAIN MATRIX VALUES (LU == A): ===
			call update_csr_ilu0_preconditioner(this%LU%n,&
											    this%LU%values,&
											    this%LU%col_indices,&
												this%LU%diag_indices,&
												this%LU%row_ptr, this%iw)
			this%is_updated = .true.
        end if
    end subroutine
    
    subroutine destroy_csr_ilu0(this)
		implicit none
        class(csr_ilu0_preconditioner_t), intent(inout) :: this
        
        nullify(this%A)
        nullify(this%mesh)
        this%name = 'Undefined Preconditioner'
        
        if (allocated(this%LU%values)) deallocate(this%LU%values)
        if (allocated(this%LU%col_indices)) deallocate(this%LU%col_indices)
        if (allocated(this%LU%diag_indices)) deallocate(this%LU%diag_indices)
        if (allocated(this%LU%row_ptr)) deallocate(this%LU%row_ptr)
        if (allocated(this%iw)) deallocate(this%iw)
    end subroutine
    
    
!=======================================================================
!================== BCSR ILU(0) PRECONDITIONER METHODS =================
!=======================================================================    
    subroutine setup_bcsr_ilu0(this, A, mesh, params)
		implicit none
        class(bcsr_ilu0_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), target, optional, intent(in) :: A
        type(mesh_t), target, optional, intent(in):: mesh
        character(len=*), optional, intent(in) :: params
        
        class(bcsr_operator_t), pointer :: bcsr_ptr
                		
        call this%destroy()
        this%name = 'BCSR-ILU(0) preconditioner'
        if (present(mesh)) this%mesh => mesh
        
        if (present(A)) then
			!=== MATRIX INITIALIZATION: ===
			select type(A)
			type is (bcsr_operator_t)
				bcsr_ptr => A
				this%A => bcsr_ptr
				
				this%LU%n = bcsr_ptr%n
				this%LU%nnz = bcsr_ptr%nnz
				this%LU%bs = bcsr_ptr%bs
				
				allocate(this%LU%values(this%LU%bs, this%LU%bs, this%LU%nnz),&
				 this%LU%col_indices(this%LU%nnz),&
				 this%LU%diag_indices(this%LU%n),&
				 this%LU%row_ptr(this%LU%n + 1),&
				 this%iw(this%LU%n))
        
				this%LU%col_indices = bcsr_ptr%col_indices
				this%LU%diag_indices = bcsr_ptr%diag_indices
				this%LU%row_ptr = bcsr_ptr%row_ptr
        	
				this%is_setup = .true.
			end select
        else if (present(mesh)) then
			!=== NON-MATRIX MESH BASED INITIALIZATION ===
			this%LU%bs = this%mesh%dim + 2		 
			call bcsr_mtrx_initialize1(this%mesh%ncells, this%mesh%nfaces,&
									   this%mesh%nbfaces, this%LU%bs,&
									   this%mesh%cell_faces_ptr, this%mesh%cell_faces,&
									   this%mesh%face_left_cell, this%mesh%face_right_cell,&
									   this%LU%n, this%LU%nnz,&
									   this%LU%values,&
									   this%LU%col_indices,&
									   this%LU%row_ptr,&
									   this%LU%diag_indices)
			allocate(this%iw(this%LU%n))
			this%is_setup = .true.				   
        end if
        
        if (.not. this%is_setup) then
			error stop 'BCSR-ILU(0) precondtioner stting up error'
        end if
    end subroutine
    
    subroutine apply_bcsr_ilu0(this, r, z)
		implicit none
        class(bcsr_ilu0_preconditioner_t), intent(inout) :: this
        real(8), intent(in), contiguous, target :: r(:)
        real(8), intent(inout), contiguous :: z(:)
       
        call apply_bcsr_ilu0_preconditioner(this%LU%n,&
											this%LU%bs,&
											this%LU%values,&
											this%LU%col_indices,&
											this%LU%diag_indices,&
											this%LU%row_ptr, r, z)
    end subroutine
    
    subroutine update_bcsr_ilu0(this, A)
		implicit none
        class(bcsr_ilu0_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), optional, intent(in) :: A
		
		class(bcsr_operator_t), pointer :: bcsr_ptr
				
        if (present(A)) then
            call this%destroy()
            call this%setup(A)
        end if
        
        if (associated(this%A)) then
			!=== MATRIX UPDATING, LU%VALUES CAN BE NON FILLED: ===
			select type (mtrx => this%A)
			type is (bcsr_operator_t)
				bcsr_ptr => mtrx
			end select
			this%LU%values = bcsr_ptr%values
			
			call update_bcsr_ilu0_preconditioner(this%LU%n,&
												 this%LU%bs,&
												 this%LU%values,&
												 this%LU%col_indices,&
												 this%LU%diag_indices,&
												 this%LU%row_ptr, this%iw)
			this%is_updated = .true.
        else
			!=== MATRIX UPDATING, LU%VALUES MUST CONTAIN MATRIX VALUES (LU == A): ===
			call update_bcsr_ilu0_preconditioner(this%LU%n,&
												 this%LU%bs,&
												 this%LU%values,&
												 this%LU%col_indices,&
												 this%LU%diag_indices,&
												 this%LU%row_ptr, this%iw)
			this%is_updated = .true.
        end if
    end subroutine
    
    subroutine destroy_bcsr_ilu0(this)
		implicit none
        class(bcsr_ilu0_preconditioner_t), intent(inout) :: this
        
        nullify(this%A)
        nullify(this%mesh)
        this%name = 'Undefined Preconditioner'
        
        if (allocated(this%LU%values)) deallocate(this%LU%values)
        if (allocated(this%LU%col_indices)) deallocate(this%LU%col_indices)
        if (allocated(this%LU%diag_indices)) deallocate(this%LU%diag_indices)
        if (allocated(this%LU%row_ptr)) deallocate(this%LU%row_ptr)
        if (allocated(this%iw)) deallocate(this%iw)
    end subroutine


!=======================================================================
!==================== BLUSGS PRECONDITIONER METHODS ====================
!=======================================================================    
    subroutine setup_blusgs_preconditioner(this, A, mesh, params)
		implicit none
        class(blusgs_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), target, optional, intent(in) :: A
        type(mesh_t), target, optional, intent(in):: mesh
        character(len=*), optional, intent(in) :: params
        
                		
        call this%destroy()
        this%name = 'BLUSGS preconditioner'
        if (present(mesh)) this%mesh => mesh
        
        allocate(this%sprf(mesh%nfaces),&
				 this%J_f(mesh%dim+2, mesh%dim+2, mesh%nfaces),&
				 this%J_d(mesh%dim+2, mesh%dim+2, mesh%ncells))
				 
        this%is_setup = .true.	
        
        if (.not. this%is_setup) then
			error stop 'BLUSGS precondtioner stting up error'
        end if
    end subroutine
    
    subroutine apply_blusgs_preconditioner(this, r, z)
		implicit none
        class(blusgs_preconditioner_t), intent(inout) :: this
        real(8), intent(in), contiguous, target :: r(:)
        real(8), intent(inout), contiguous :: z(:)
		
		real(8), pointer, contiguous :: RHS(:, :)
		
		RHS(1:this%mesh%dim+2, 1:this%mesh%ncells) => r(:)
		z = 0.d0
        call blusgs_apply(this%mesh%dim, this%mesh%ncells, this%w, -1.d0,&
						  this%mesh%face_left_cell,&
						  this%mesh%face_right_cell,&
						  this%mesh%cell_faces,&
						  this%mesh%cell_faces_ptr,&
						  z, RHS, this%sprf, this%J_f, this%J_d)
    end subroutine
    
    subroutine update_blusgs_preconditioner(this, A)
		implicit none
        class(blusgs_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), optional, intent(in) :: A
		
		call blusgs_setup(this%mesh%dim, this%mesh%ncells,&
						  this%mesh%nfaces, this%mesh%nbfaces,&
						  this%mesh%face_left_cell, this%mesh%face_right_cell,&
						  this%mesh%cell_faces, this%mesh%cell_faces_ptr,&
						  this%mesh%face_area, this%mesh%face_normal,&
						  this%mesh%face_center, this%mesh%cell_center,&
						  this%mesh%cell_volume, this%ptime_step,&
						  this%ppm%k, this%ppm%R_gas, this%ppm%cv,&
						  this%ppm%cp, this%ppm%Pr, this%w,&
						  this%Q, this%sprf, this%J_f, this%J_d,&
						  this%settings%MODEL,&
						  this%settings%USE_GHOST_CELLS,&
						  this%settings%USE_CONSERVATIVE_VARS)
		
        this%is_updated = .true.
    end subroutine
    
    subroutine destroy_blusgs_preconditioner(this)
		implicit none
        class(blusgs_preconditioner_t), intent(inout) :: this
        
        nullify(this%mesh)
        this%name = 'Undefined Preconditioner'
        
        if (allocated(this%sprf)) deallocate(this%sprf)
        if (allocated(this%J_f)) deallocate(this%J_f)
        if (allocated(this%J_d)) deallocate(this%J_d)
    end subroutine



!=======================================================================
!===================== GMG PRECONDITIONER METHODS ======================
!=======================================================================    
    subroutine setup_gmg_preconditioner(this, A, mesh, params)
		implicit none
        class(gmg_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), target, optional, intent(in) :: A
        type(mesh_t), target, optional, intent(in):: mesh
        character(len=*), optional, intent(in) :: params
        
                		
        call this%destroy()
        this%name = 'GMG preconditioner'
        if (present(mesh)) this%mesh => mesh
        
        call this%gmg_hierarchy%create(mesh)
				 
        this%is_setup = .true.	
        
        if (.not. this%is_setup) then
			error stop 'GMG precondtioner stting up error'
        end if
    end subroutine
    
    subroutine apply_gmg_preconditioner(this, r, z)
		implicit none
        class(gmg_preconditioner_t), intent(inout) :: this
        real(8), intent(in), contiguous, target :: r(:)
        real(8), intent(inout), contiguous :: z(:)
		
		this%gmg_hierarchy%levels(1)%phi = 0.d0
		this%gmg_hierarchy%levels(1)%rhs = r
		call this%gmg_hierarchy%v_cycle(1)
		z = this%gmg_hierarchy%levels(1)%phi
		
    end subroutine
    
    subroutine update_gmg_preconditioner(this, A)
		implicit none
        class(gmg_preconditioner_t), intent(inout) :: this
        class(linear_operator_t), optional, intent(in) :: A
		
		logical :: tmp_flag = .false.
		integer :: i, j, d, dim
		
		dim = this%gmg_hierarchy%levels(1)%mesh%dim
		
		do i = 2, this%gmg_hierarchy%current_max_level
			!=== PRIMITIVE VARIABLES RESTRICTION: ===
			call mg_restriction_vector_cell(this%gmg_hierarchy%levels(i)%mesh%ncells,&
											this%gmg_hierarchy%levels(i)%mesh%dim+2,&
											this%gmg_hierarchy%levels(i)%mesh%cell_volume,&
											this%gmg_hierarchy%levels(i-1)%mesh%cell_volume,&
											this%gmg_hierarchy%levels(i)%agglom_info%fine_cells_ptr,&
											this%gmg_hierarchy%levels(i)%agglom_info%fine_cells,&
											this%gmg_hierarchy%levels(i-1)%vars_ptr,&
											this%gmg_hierarchy%levels(i)%vars_ptr)
											
			!=== UPDATING BOUNDARY FACES/GHOST CELLS VALUES: ===
			call this%gmg_hierarchy%levels(i)%bcm%update_boundary_values('Q', .true., this%gmg_hierarchy%levels(i)%settings%USE_GHOST_CELLS,&
																		 this%gmg_hierarchy%levels(i)%QIDX)
			
			!=== COMPUTING GRADIENTS: ===
			call this%gmg_hierarchy%levels(i)%gradient%apply_vector_cellfield(this%gmg_hierarchy%levels(i)%vars_ptr,&
																			  this%gmg_hierarchy%levels(i)%grad_vars_ptr,&
																			  this%gmg_hierarchy%levels(i)%mesh%dim + 2,&
																			  this%gmg_hierarchy%levels(i)%settings%USE_GHOST_CELLS)
			
			!=== UPDATING BOUNDARY FACES/GHOST CELLS GRADIENTS: ===
			call this%gmg_hierarchy%levels(i)%bcm%update_boundary_values('Q', .true., this%gmg_hierarchy%levels(i)%settings%USE_GHOST_CELLS,&
											     this%gmg_hierarchy%levels(i)%QIDX)
			
			!=== COMPUTING PSEUDO-TIME STEP: ===
			call this%gmg_hierarchy%levels(i)%compute_ptime()
			
			!=== COMPUTING FLUXES: ===
			call this%gmg_hierarchy%levels(i)%flxsm%compute_fluxes_cmprs(this%gmg_hierarchy%levels(i)%QIDX,&
																		 this%gmg_hierarchy%levels(i)%FLUXESQIDX)
																		 
			select type(lo => this%gmg_hierarchy%levels(i)%linear_operator)
			type is (bcsr_operator_t)
				call bcsr_mtrx_cmpt(lo)
			end select 
																		 															
		end do
		
		!=== SMOOTHING MATRIX COPYING (IF NOT MATRIX-FREE): ===
		do i = 1, this%gmg_hierarchy%current_max_level
			select type(lo => this%gmg_hierarchy%levels(i)%linear_operator)
				type is (bcsr_operator_t)
				this%gmg_hierarchy%levels(i)%LU%values = lo%values
				tmp_flag = .true.
			end select
		end do
	
		!=== SMOOTHING MATRIX COMPUTATION (IF MATRIX-FREE): ===
		if (.not. tmp_flag) then
			do i = 1, this%gmg_hierarchy%current_max_level
				call bcsr_mtrx_cmpt(this%gmg_hierarchy%levels(i)%LU)
			end do
		end if
		
		
		!=== BILU(0) DECOMPOSING: ===
		do i = 1, this%gmg_hierarchy%current_max_level
			call update_bcsr_ilu0_preconditioner(this%gmg_hierarchy%levels(i)%LU%n,&
												 this%gmg_hierarchy%levels(i)%LU%bs,&
												 this%gmg_hierarchy%levels(i)%LU%values,&
												 this%gmg_hierarchy%levels(i)%LU%col_indices,&
												 this%gmg_hierarchy%levels(i)%LU%diag_indices,&
												 this%gmg_hierarchy%levels(i)%LU%row_ptr, this%gmg_hierarchy%levels(i)%iw)
		end do
		
        this%is_updated = .true.
        
        contains
        subroutine bcsr_mtrx_cmpt(lo)
			type(bcsr_operator_t) :: lo
			if (.not. this%gmg_hierarchy%levels(i)%settings%USE_CONSERVATIVE_VARS) then
				!=== PRIMITIVE VARIABLES: ===
				call db_assemble_bcsr_matrix_inv(this%gmg_hierarchy%levels(i)%mesh%dim,&
												 this%gmg_hierarchy%levels(i)%mesh%ncells,&
												 this%gmg_hierarchy%levels(i)%mesh%nfaces,&
												 this%gmg_hierarchy%levels(i)%mesh%nbfaces,&
												 this%gmg_hierarchy%levels(i)%mesh%face_left_cell,&
												 this%gmg_hierarchy%levels(i)%mesh%face_right_cell,&
												 this%gmg_hierarchy%levels(i)%mesh%face_bidx,&
												 this%gmg_hierarchy%levels(i)%mesh%face_normal,&
												 this%gmg_hierarchy%levels(i)%mesh%face_area,&
												 this%gmg_hierarchy%levels(i)%mesh%face_weight,&
												 this%gmg_hierarchy%levels(i)%mesh%cell_volume,&
												 this%gmg_hierarchy%levels(i)%ptime_step,&
												 this%gmg_hierarchy%levels(i)%vars_ptr,&
												 lo%values,&
												 lo%diag_indices,&
												 this%gmg_hierarchy%levels(i)%map_LL,&
												 this%gmg_hierarchy%levels(i)%map_LR,&
												 this%gmg_hierarchy%levels(i)%map_RL,&
												 this%gmg_hierarchy%levels(i)%map_RR,&
												 this%gmg_hierarchy%levels(i)%map_LB,&
												 this%gmg_hierarchy%levels(i)%ppm%k,&
												 this%gmg_hierarchy%levels(i)%ppm%R_gas,&
												 this%gmg_hierarchy%levels(i)%ppm%cp,&
												 this%gmg_hierarchy%levels(i)%ppm%cv,&
												 this%gmg_hierarchy%levels(i)%settings%STABILITY_OPERATOR_TYPE,&
												 this%gmg_hierarchy%levels(i)%settings%USE_GHOST_CELLS, .true.)
				if (this%gmg_hierarchy%levels(i)%settings%MODEL == 2) then
					call db_assemble_bcsr_matrix_visc(this%gmg_hierarchy%levels(i)%mesh%dim,&
													  this%gmg_hierarchy%levels(i)%mesh%ncells,&
													  this%gmg_hierarchy%levels(i)%mesh%nfaces,&
													  this%gmg_hierarchy%levels(i)%mesh%nbfaces,&
													  this%gmg_hierarchy%levels(i)%mesh%face_left_cell,&
													  this%gmg_hierarchy%levels(i)%mesh%face_right_cell,&
													  this%gmg_hierarchy%levels(i)%mesh%face_bidx,&
													  this%gmg_hierarchy%levels(i)%mesh%face_normal,&
													  this%gmg_hierarchy%levels(i)%mesh%face_area,&
													  this%gmg_hierarchy%levels(i)%mesh%face_weight,&   
													  this%gmg_hierarchy%levels(i)%mesh%cell_center,&  
													  this%gmg_hierarchy%levels(i)%vars_ptr,&
													  lo%values,&
													  this%gmg_hierarchy%levels(i)%map_LL,&
													  this%gmg_hierarchy%levels(i)%map_LR,&
													  this%gmg_hierarchy%levels(i)%map_RL,&
													  this%gmg_hierarchy%levels(i)%map_RR,&
													  this%gmg_hierarchy%levels(i)%map_LB,&
													  this%gmg_hierarchy%levels(i)%ppm%k,&
													  this%gmg_hierarchy%levels(i)%ppm%R_gas,&
													  this%gmg_hierarchy%levels(i)%ppm%cp,&
													  this%gmg_hierarchy%levels(i)%ppm%cv,&
													  this%gmg_hierarchy%levels(i)%ppm%Pr,&   
													  this%gmg_hierarchy%levels(i)%settings%USE_GHOST_CELLS, .false.)
				end if
			else
				!=== CONSERVATIVE VARIABLES: ===
				call db_assemble_bcsr_matrix_inv_cv(this%gmg_hierarchy%levels(i)%mesh%dim,&
													this%gmg_hierarchy%levels(i)%mesh%ncells,&
													this%gmg_hierarchy%levels(i)%mesh%nfaces,&
													this%gmg_hierarchy%levels(i)%mesh%nbfaces,&
													this%gmg_hierarchy%levels(i)%mesh%face_left_cell,&
													this%gmg_hierarchy%levels(i)%mesh%face_right_cell,&
													this%gmg_hierarchy%levels(i)%mesh%face_bidx,&
													this%gmg_hierarchy%levels(i)%mesh%face_normal,&
													this%gmg_hierarchy%levels(i)%mesh%face_area,&
													this%gmg_hierarchy%levels(i)%mesh%face_weight,&
													this%gmg_hierarchy%levels(i)%mesh%cell_volume,&
													this%gmg_hierarchy%levels(i)%ptime_step,&
													this%gmg_hierarchy%levels(i)%vars_ptr,&
													lo%values,&
													lo%diag_indices,&
													this%gmg_hierarchy%levels(i)%map_LL,&
													this%gmg_hierarchy%levels(i)%map_LR,&
													this%gmg_hierarchy%levels(i)%map_RL,&
													this%gmg_hierarchy%levels(i)%map_RR,&
													this%gmg_hierarchy%levels(i)%map_LB,&
													this%gmg_hierarchy%levels(i)%ppm%k,&
													this%gmg_hierarchy%levels(i)%ppm%R_gas,&
													this%gmg_hierarchy%levels(i)%ppm%cp,&
													this%gmg_hierarchy%levels(i)%ppm%cv,&
													this%gmg_hierarchy%levels(i)%settings%STABILITY_OPERATOR_TYPE,&
													this%gmg_hierarchy%levels(i)%settings%USE_GHOST_CELLS, .true.)
				if (this%gmg_hierarchy%levels(i)%settings%MODEL == 2) then
					call db_assemble_bcsr_matrix_visc_cv(this%gmg_hierarchy%levels(i)%mesh%dim,&
														 this%gmg_hierarchy%levels(i)%mesh%ncells,&
														 this%gmg_hierarchy%levels(i)%mesh%nfaces,&
														 this%gmg_hierarchy%levels(i)%mesh%nbfaces,&
														 this%gmg_hierarchy%levels(i)%mesh%face_left_cell,&
														 this%gmg_hierarchy%levels(i)%mesh%face_right_cell,&
														 this%gmg_hierarchy%levels(i)%mesh%face_bidx,&
														 this%gmg_hierarchy%levels(i)%mesh%face_normal,&
														 this%gmg_hierarchy%levels(i)%mesh%face_area,&
														 this%gmg_hierarchy%levels(i)%mesh%face_weight,&   
														 this%gmg_hierarchy%levels(i)%mesh%cell_center,&  
														 this%gmg_hierarchy%levels(i)%vars_ptr,&
														 lo%values,&
														 this%gmg_hierarchy%levels(i)%map_LL,&
														 this%gmg_hierarchy%levels(i)%map_LR,&
														 this%gmg_hierarchy%levels(i)%map_RL,&
														 this%gmg_hierarchy%levels(i)%map_RR,&
														 this%gmg_hierarchy%levels(i)%map_LB,&
														 this%gmg_hierarchy%levels(i)%ppm%k,&
														 this%gmg_hierarchy%levels(i)%ppm%R_gas,&
														 this%gmg_hierarchy%levels(i)%ppm%cp,&
														 this%gmg_hierarchy%levels(i)%ppm%cv,&
														 this%gmg_hierarchy%levels(i)%ppm%Pr,&    
														 this%gmg_hierarchy%levels(i)%settings%USE_GHOST_CELLS, .false.)
				end if
			end if
			
			!=== BOUNDARY CONDITIONS: ===
			call this%gmg_hierarchy%levels(i)%bcm%add_boundary_jacobians(lo%values,&
																		 this%gmg_hierarchy%levels(i)%map_LB,&
																		 'Q', this%gmg_hierarchy%levels(i)%settings%MODEL,&
																		 this%gmg_hierarchy%levels(i)%settings%USE_GHOST_CELLS,&
																		 this%gmg_hierarchy%levels(i)%settings%USE_CONSERVATIVE_VARS,&
																		 this%gmg_hierarchy%levels(i)%QIDX)
        end subroutine
    end subroutine
    
    subroutine destroy_gmg_preconditioner(this)
		implicit none
        class(gmg_preconditioner_t), intent(inout) :: this
        
        nullify(this%mesh)
        this%name = 'Undefined Preconditioner'
        
        print*, 'GMG PREC DSTR TBD'
    end subroutine
                
end module
