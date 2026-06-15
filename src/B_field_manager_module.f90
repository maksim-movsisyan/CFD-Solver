module field_manager_module
use mesh_module, only: mesh_t
implicit none
private
!=======================================================================
!===================== DATA LOCATION NUMERATION ========================
!=======================================================================
integer, parameter, public :: LOC_CELL = 1 								!=== CELL CENTERS ===
integer, parameter, public :: LOC_FACE = 2								!=== FACE CENTERS ===
integer, parameter, public :: LOC_NODE = 3								!=== NODE CENTERS ===
integer, parameter, public :: LOC_B_FACE = 4 							!=== BOUNDARY FACE CENTERS ===


!=======================================================================
!=================== ABSTRACT FIELD DATA STRUCTURE =====================
!=======================================================================
type, public :: field_t
	character(len=32) :: name		     								!=== FIELD NAME ===
	integer :: location = 0												!=== FIELD LOCATION (CELL/FACE/NODE/BOUNDARY FACE) ===
	integer :: n_comp = 0												!=== NUMBER OF COMPONETNS: SCALAR = 1; VECTORS = dim; TENSOR = dim^2; ... ===
	integer :: n_elem = 0												!=== NUMBER OF ELEMENTS (CELLS/FACES/NODES/BOUNDARY FACES) ===
	
	integer :: n_internal = 0											!=== NUMBER OF INTERNAL VALUES ===
	integer :: n_external = 0											!=== NUMBER OF EXTERNAL VALUES ===
	
    logical :: is_output = .false.  									!=== OUTPUT FLAG ===

	real(8), allocatable :: values(:,:)									!=== FIELD VALUES, DIM = (n_comp, n_elem) <-> default: AoS ===		
	logical :: has_bound = .false.										!=== INCLUDE BOUNDARY FLAG ===
	
	real(8), pointer, contiguous :: grad(:, :) => null()				!=== GRADIENT POINTER ===
	logical :: has_gradient = .false.									!=== INCLUDE GRADIENT FLAG ===
	
contains
	procedure :: get_ptr => field_get_ptr
	procedure :: set_gradient => field_set_gradinet_ptr
end type

!=======================================================================
!===================== FIELD MANAGER DATA STRUCTURE ====================
!=======================================================================
type, public :: field_manager_t
	type(field_t) :: registry(100)										!=== FIELDS REGISTRY ===
	integer :: count = 0												!=== NUMBER OF FIELDS ===
	
	type(mesh_t), pointer :: mesh => null()								!=== MESH POINTER ===
	integer :: verbosity = 0											!=== OUTPUT INFORMATION LVL ===
	
contains
	procedure :: add_field => fm_add_field								!=== SUBROUTINE FOR ADDING FIELD IN REGISTRY ===
	procedure :: get_ptr => fm_get_ptr									!=== FUNCTION FOR GETTING FIELD POINTER (BASIC FUNCTION) ===
	procedure :: get_idx => fm_get_idx									!=== FUNCTION FOR GETTING FIELD INDEX IN REGISTRY ===
	
	
	procedure :: get_scalar_ptr => fm_get_scalar_ptr					!=== FUNCTION FOR GETTING FIELD POINTER (SCALAR FIELD: PTR(:)) ===
	procedure :: get_vector_ptr => fm_get_vector_ptr					!=== FUNCTION FOR GETTING FIELD POINTER (VECTOR FIELD: PTR(:, :)) ===	
	procedure :: get_component_ptr => fm_get_component_ptr				!=== FUNCTION FOR GETTING FIELD POINTER (VECTOR COMPONENT SLICE!!!: PTR(:, :)) ===
		
    
    procedure :: initialize => fm_initialize							!=== FIELD MANAGER T CONSTRUCTOR ===
	procedure :: finalize => fm_finalize								!=== FIELD MANAGER T DESTRUCTOR ===
end type


contains
!=======================================================================
!======================= ABSTRACT FIELD METHODS ========================
!=======================================================================
	function field_get_ptr(this) result(ptr)
		implicit none
		class(field_t), intent(in), target :: this
		real(8), pointer, contiguous :: ptr(:,:)
		ptr => this%values
	end function

	subroutine field_set_gradinet_ptr(this, grad_ptr)
		implicit none
        class(field_t), intent(inout), target :: this
        real(8), intent(in), contiguous, target :: grad_ptr(:, :)
               
       this%grad => grad_ptr
       this%has_gradient = .true. 
	end subroutine

!=======================================================================
!========================= FIELD MANAGER METHODS =======================
!=======================================================================
	!=== ADD FIELD TO REGISTRY: ===
	subroutine fm_add_field(this, name, is_output, location, n_comp, n_internal, n_external)
		implicit none
        class(field_manager_t), intent(inout) :: this
        character(len=*), intent(in) :: name
        logical, intent(in) :: is_output
        integer, intent(in) :: location, n_comp, n_internal, n_external
        
        integer :: i
		
        this%count = this%count + 1
        
        associate(f => this%registry(this%count))
            f%name = name
            f%is_output	= is_output
            f%location = location
            f%n_comp = n_comp
            f%n_internal = n_internal
            f%n_external = n_external
            f%n_elem = n_internal + n_external
            allocate(f%values(f%n_comp, f%n_elem), source=0.0d0)
            
            if (n_external > 0) f%has_bound = .true.
        end associate
        
        if (this%verbosity > 1) then
			write(*, '(a, a, a)') 'Field ', trim(name), ' added to field manager successfully!'
		end if
			
    end subroutine
    
	!=== GET FIELD IDX (IN REGISTRY): ===
    function fm_get_idx(this, name) result(idx)
		implicit none
        class(field_manager_t), intent(in) :: this
        character(len=*), intent(in) :: name
        character(len=32) :: search_name
        
        integer :: idx, i
        
        search_name = name
        idx = -1
        
        do i = 1, this%count
            if (this%registry(i)%name == search_name) then
                idx = i
                return
            end if
        end do
    end function

	!=== GET FIELD POINTER, BASIC: ===
    function fm_get_ptr(this, name) result(ptr)
		implicit none
        class(field_manager_t), intent(in), target :: this
        character(len=*), intent(in) :: name
        real(8), pointer, contiguous :: ptr(:,:)
        integer :: idx
        
        idx = this%get_idx(name)
        if (idx > 0) then
            ptr => this%registry(idx)%values
        else
            ptr => null()
        end if
    end function


	!=== GET FIELD POINTER, SCALAR FIELD: ===
	function fm_get_scalar_ptr(this, name) result(ptr)
		implicit none
        class(field_manager_t), intent(in), target :: this
        character(len=*), intent(in) :: name
        real(8), pointer, contiguous :: ptr(:)
        
        integer :: idx
        
        idx = this%get_idx(name)
        if (idx > 0) then
            ptr => this%registry(idx)%values(1, :)
        else
            ptr => null()
        end if
    end function
    
    
    !=== GET FIELD POINTER, VECTOR FIELD: ===
    function fm_get_vector_ptr(this, name) result(ptr)
		implicit none
        class(field_manager_t), intent(in), target :: this
        character(len=*), intent(in) :: name
        real(8), pointer, contiguous :: ptr(:,:)
        ptr => this%get_ptr(name)
    end function
    
  
    !=== GET FIELD POINTER, VECTOR COMPONETNT FIELD: ===
    function fm_get_component_ptr(this, name, comp) result(ptr)
		implicit none
        class(field_manager_t), intent(in), target :: this
        character(len=*), intent(in) :: name
        integer, intent(in) :: comp
        real(8), pointer :: ptr(:)
        integer :: idx
        
        idx = this%get_idx(name)
        if (idx > 0 .and. comp <= this%registry(idx)%n_comp) then
            ptr => this%registry(idx)%values(comp, :)
        else
            ptr => null()
        end if
    end function
    
  
	!=== DESTRUCTOR ===
	subroutine fm_finalize(this)
		implicit none
		class(field_manager_t), intent(inout) :: this
		
		integer :: i
		character(len=32) :: name
	
		do i = 1, this%count
			name = this%registry(i)%name
			
			this%registry(i)%name = 'None'
			this%registry(i)%is_output = .false.
			this%registry(i)%location = 0
			this%registry(i)%n_comp = 0
			this%registry(i)%n_internal = 0
			this%registry(i)%n_external = 0
			this%registry(i)%n_elem = 0
			deallocate(this%registry(i)%values)
			
			if (this%verbosity > 1) then
				write(*, '(a, a, a)') 'Field ', trim(name), ' removed from field manager successfully!'
			end if
		end do
			
		if (this%verbosity > 0) then
			write(*, '(a, i0)') 'Field Manager successfully destroyed! Number of destroyed fields = ', this%count
		end if
		this%count = 0	
		
	end subroutine
	
	!=== SET MESH ===
	subroutine fm_initialize(this, mesh)
		implicit none
		class(field_manager_t), intent(inout) :: this
		type(mesh_t), intent(in), target :: mesh
		
		this%mesh => mesh
	end subroutine
	
end module


