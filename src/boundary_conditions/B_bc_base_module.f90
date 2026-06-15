module bc_base_module
use mesh_module, only: mesh_t
implicit none
private
!=======================================================================
!==================== BC LOCAL TYPES NUMERATION ========================
!=======================================================================
integer, parameter, public :: BC_TYPE_UNDEFINED = 0
integer, parameter, public :: BC_TYPE_WALL = 1
integer, parameter, public :: BC_TYPE_INLET = 2
integer, parameter, public :: BC_TYPE_OUTLET = 3
integer, parameter, public :: BC_TYPE_SYMMETRY = 4
integer, parameter, public :: BC_TYPE_FARFIELD = 5
integer, parameter, public :: BC_TYPE_PRESSURE_INLET = 6
integer, parameter, public :: BC_TYPE_PRESSURE_OUTLET = 7
integer, parameter, public :: BC_TYPE_MASS_FLOW_INLET = 8
integer, parameter, public :: BC_TYPE_PERIODIC = 9

!=== WALL SUBTYPES ===
integer, parameter, public :: BC_WALL_SLIP_REIMAN = 100
integer, parameter, public :: BC_WALL_NO_SLIP = 101
integer, parameter, public :: BC_WALL_SLIP = 102
integer, parameter, public :: BC_WALL_ADIABATIC = 103
integer, parameter, public :: BC_WALL_ISOTHERMAL = 104
integer, parameter, public :: BC_WALL_HEAT_FLUX = 105

!=== INLET SUBTYPES ===
integer, parameter, public :: BC_INLET_SUBSONIC = 201
integer, parameter, public :: BC_INLET_SUPERSONIC = 202
integer, parameter, public :: BC_INLET_TOTAL = 203

!=== OUTLET SUBTYPES ===
integer, parameter, public :: BC_OUTLET_SUBSONIC = 301
integer, parameter, public :: BC_OUTLET_SUPERSONIC = 302


!=======================================================================
!=========== BOUNDARY CONDITION ABSTRACT TYPE DATA STRUCTURE ===========
!=======================================================================
type, abstract, public :: bc_base_t
	character(len=32) :: bc_name											
	integer :: zone_id													!=== BOUNDARY PATCH FACES ZONE (mesh%face_zone) ===	
	
	integer :: dim														!=== PROBLEM DIMENSION ===
	type(mesh_t), pointer :: mesh => null()								!=== MESH POINTER ===
	
	integer :: nfaces													!=== NUMBER OF FACES IN PATCH ===
	integer :: bc_type = BC_TYPE_UNDEFINED								!=== BC TYPE ===
	integer :: bc_subtype = BC_TYPE_UNDEFINED							!=== BC SUBTYPE ===
	
	
	integer, allocatable :: face_indices(:)								!=== GLOBAL BOUNDARY FACE INDICES, DIM = (nfaces) ===
	integer, allocatable :: cell_indices(:)								!=== GLOBAL BOUNDARY CELL INDICES, DIM = (nfaces) ===
    integer, allocatable :: b_indices(:)								!=== LOCAL BOUNDARY FACE INDICES, DIM = (nfaces)  ===
    
	integer :: verbosity = 0											!=== OUTPUT INFORMATION LVL ===
	
contains
	procedure :: initialize => bc_initialize							!=== BC PATCH INITIALIZATION SUBROUITNE ===
	procedure(apply_patch_interface), deferred :: apply_patch			!=== MAIN SUBROUTINE FOR APPLYING BOUNDARY CONDITION ON PATCH face_indices ===
	procedure(ADdiff_apply_patch_interface), deferred ::&
											ADdiff_apply_patch			!=== ALG DIFF MAIN SUBROUTINE FOR APPLYING BOUNDARY CONDITION ON PATCH face_indices ===
	procedure(apply_ghosts_patch_interface), deferred ::&
												apply_ghosts_patch		!=== SUBROUTINE FOR UPDATING GHOST CELLS VALUES ===
	procedure(ADdiff_apply_ghosts_patch_interface), deferred ::&
											ADdiff_apply_ghosts_patch	!=== ALG DIFF SUBROUTINE FOR UPDATING GHOST CELLS VALUES ===
	procedure(apply_patch_jmatrix_inv_interface), deferred ::&
												apply_jmatrix_inv_patch	!=== FOR DENSITY-BASED SOLVER: SUBROUITNE FOR ADDING INVISCID CONTRIBUTION OF BOUNDARY IN JACOBIAN MATRIX ===
	procedure(apply_patch_jmatrix_visc_interface), deferred ::&
											   apply_jmatrix_visc_patch !=== FOR DENSITY-BASED SOLVER: SUBROUITNE FOR ADDING VISCOUS CONTRIBUTION OF BOUNDARY IN JACOBIAN MATRIX ===
end type

abstract interface
	subroutine apply_patch_interface(this, R_gas, k, values_ptr, values_grad_ptr, name, update_grad)
		import :: bc_base_t
		implicit none
		class(bc_base_t), intent(inout) :: this
		real(8), intent(in) :: R_gas, k
		real(8), intent(inout), contiguous :: values_ptr(:, :)
		real(8), intent(inout), contiguous :: values_grad_ptr(:, :)
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
	end subroutine
	
	subroutine ADdiff_apply_patch_interface(this, R_gas, k, values_ptr, values_ptrd,&
											values_grad_ptr, values_grad_ptrd, name, update_grad)
		import :: bc_base_t
		implicit none
		class(bc_base_t), intent(inout) :: this
		real(8), intent(in) :: R_gas, k
		real(8), intent(in), contiguous :: values_ptr(:, :)
		real(8), intent(in), contiguous :: values_grad_ptr(:, :)
		real(8), intent(inout), contiguous :: values_ptrd(:, :)
		real(8), intent(inout), contiguous :: values_grad_ptrd(:, :)
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
	end subroutine
	
	subroutine apply_ghosts_patch_interface(this, values_ptr, values_grad_ptr, name, update_grad)
		import :: bc_base_t
		implicit none
		class(bc_base_t), intent(inout) :: this
		real(8), intent(inout), contiguous :: values_ptr(:, :)
		real(8), intent(inout), contiguous :: values_grad_ptr(:, :)
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
	end subroutine
	
	subroutine ADdiff_apply_ghosts_patch_interface(this, values_ptr, values_ptrd,&
												   values_grad_ptr, values_grad_ptrd, name, update_grad)
		import :: bc_base_t
		implicit none
		class(bc_base_t), intent(inout) :: this
		real(8), intent(in), contiguous :: values_ptr(:, :)
		real(8), intent(in), contiguous :: values_grad_ptr(:, :)
		real(8), intent(inout), contiguous :: values_ptrd(:, :)
		real(8), intent(inout), contiguous :: values_grad_ptrd(:, :)
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
	end subroutine
	
	subroutine apply_patch_jmatrix_inv_interface(this, R_gas, k, values_ptr, values, map_LB, name)
		import :: bc_base_t
		implicit none
		class(bc_base_t), intent(inout) :: this
		real(8), intent(in) :: R_gas, k
		real(8), intent(in), contiguous :: values_ptr(:, :)
		integer, intent(in), contiguous :: map_LB(:)
		real(8), intent(inout), contiguous :: values(:, :, :)
		character(len=*), intent(in) :: name
	end subroutine
	
	subroutine apply_patch_jmatrix_visc_interface(this, R_gas, k, cp, Pr, values_ptr, values, map_LB, name)
		import :: bc_base_t
		implicit none
		class(bc_base_t), intent(inout) :: this
		real(8), intent(in) :: R_gas, k, cp, Pr
		real(8), intent(in), contiguous :: values_ptr(:, :)
		integer, intent(in), contiguous :: map_LB(:)
		real(8), intent(inout), contiguous :: values(:, :, :)
		character(len=*), intent(in) :: name
	end subroutine
	
	
end interface

contains
!=======================================================================
!=============== BOUNDARY CONDITION ABSTRACT TYPE METHODS ==============
!=======================================================================
	subroutine bc_initialize(this, mesh, args)
		implicit none
		class(bc_base_t), intent(inout) :: this
		type(mesh_t), intent(in), target :: mesh
		real(8), optional, intent(in) :: args(:)
		
		integer :: i, counter
		integer :: left_cell, right_cell
		
		this%mesh => mesh
		this%dim = mesh%dim
		
		counter = 0
		do i = mesh%nfaces - mesh%nbfaces + 1, mesh%nfaces
			left_cell = mesh%face_left_cell(i)
			right_cell = mesh%face_right_cell(i)
			
			if (mesh%face_zone(i) == this%zone_id) then
				if (left_cell > this%mesh%ncells .or. right_cell > this%mesh%ncells) then
					counter = counter + 1
				end if
			end if
		end do
		this%nfaces = counter
		
		if (counter > 0) then
			allocate(this%face_indices(counter))
			allocate(this%cell_indices(counter))
			allocate(this%b_indices(counter))
		end if 
		
		counter = 0
		do i = mesh%nfaces - mesh%nbfaces + 1, mesh%nfaces
			left_cell = mesh%face_left_cell(i)
			right_cell = mesh%face_right_cell(i)
			
			if (mesh%face_zone(i) == this%zone_id) then
				if (left_cell > this%mesh%ncells .or. right_cell > this%mesh%ncells)  then
					counter = counter + 1
					this%face_indices(counter) = i
					this%b_indices(counter) = mesh%face_bidx(i)
				  
					
					if (left_cell > this%mesh%ncells) then
						this%cell_indices(counter) = right_cell
					else
						this%cell_indices(counter) = left_cell
					end if
				end if
			end if
		end do
		
		if (this%verbosity > 0) then
			print*, 'BOUNDARY CONDITION <', this%bc_name, '> INITIALIZED SUCCESSFULLY, NUM FACES: ', this%nfaces
		end if
		
	end subroutine
end module
