module diff_operator_module
use mesh_module, only: mesh_t
use gradient_computation_module
use ADdiff_gradient_computation_module
implicit none
!=======================================================================
!=============== GRADIENT OPERATOR BASIC DATA STRUCTURE ================
!=======================================================================
type, abstract :: gradient_t
	integer :: dim = 2
	type(mesh_t), pointer :: mesh => null()
	character(len=64) :: name = 'GradientDefault'
	
	contains
	procedure(grad_apply_scalar_cellfield_interface), deferred :: apply_scalar_cellfield
	procedure(grad_apply_vector_cellfield_interface), deferred :: apply_vector_cellfield
	procedure(grad_apply_scalar_cell_interface), deferred :: apply_scalar_cell
	procedure(grad_apply_vector_cell_interface), deferred :: apply_vector_cell
	
	procedure(ADgrad_apply_scalar_cellfield_interface), deferred :: ADapply_scalar_cellfield
	procedure(ADgrad_apply_vector_cellfield_interface), deferred :: ADapply_vector_cellfield
	procedure :: set_mesh => gradient_set_mesh
end type

abstract interface
	subroutine grad_apply_scalar_cellfield_interface(this, Field, Field_grad, USE_GHOST_CELLS)
		import :: gradient_t
		class(gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:)
		logical, intent(in) :: USE_GHOST_CELLS
		real(8), intent(inout), contiguous :: Field_grad(:, :)
	end subroutine
	
	subroutine grad_apply_vector_cellfield_interface(this, Field, Field_grad, nvars, USE_GHOST_CELLS)
		import :: gradient_t
		class(gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: nvars
		real(8), intent(inout), contiguous :: Field_grad(:, :)
	end subroutine
	
	subroutine grad_apply_scalar_cell_interface(this, cell_idx, Field, Cell_grad, USE_GHOST_CELLS)
		import :: gradient_t
		class(gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: cell_idx
		real(8), intent(inout), contiguous :: Cell_grad(:)
	end subroutine
	
	subroutine grad_apply_vector_cell_interface(this, cell_idx, Field, Cell_grad, nvars, USE_GHOST_CELLS)
		import :: gradient_t
		class(gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: cell_idx
		integer, intent(in) :: nvars
		real(8), intent(inout), contiguous :: Cell_grad(:)
	end subroutine
	
	
	
	subroutine ADgrad_apply_scalar_cellfield_interface(this, dField, dField_grad, USE_GHOST_CELLS)
		import :: gradient_t
		class(gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: dField(:)
		logical, intent(in) :: USE_GHOST_CELLS
		real(8), intent(inout), contiguous :: dField_grad(:, :)
	end subroutine
	
	subroutine ADgrad_apply_vector_cellfield_interface(this, dField, dField_grad, nvars, USE_GHOST_CELLS)
		import :: gradient_t
		class(gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: dField(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: nvars
		real(8), intent(inout), contiguous :: dField_grad(:, :)
	end subroutine
end interface
    
    
!=======================================================================
!=========== GREEN GAUSS GRADIENT OPERATOR DATA STRUCTURE ==============
!=======================================================================
type, extends(gradient_t) :: GG_gradient_t

	contains
	procedure :: apply_scalar_cellfield => GG_apply_scalar_cellfield
	procedure :: apply_vector_cellfield => GG_apply_vector_cellfield
	procedure :: apply_scalar_cell => GG_apply_scalar_cell
	procedure :: apply_vector_cell => GG_apply_vector_cell
	
	procedure :: ADapply_scalar_cellfield => GG_ADapply_scalar_cellfield
	procedure :: ADapply_vector_cellfield => GG_ADapply_vector_cellfield
end type


!=======================================================================
!=========== LEAST SQUARES GRADIENT OPERATOR DATA STRUCTURE ============
!=======================================================================
type, extends(gradient_t) :: LSQ_gradient_t
	real(8), allocatable :: lsq_w(:, :)
	logical :: is_updated = .false.
	contains
	procedure :: apply_scalar_cellfield => LSQ_apply_scalar_cellfield
	procedure :: apply_vector_cellfield => LSQ_apply_vector_cellfield
	procedure :: apply_scalar_cell => LSQ_apply_scalar_cell
	procedure :: apply_vector_cell => LSQ_apply_vector_cell
	
	procedure :: ADapply_scalar_cellfield => LSQ_ADapply_scalar_cellfield
	procedure :: ADapply_vector_cellfield => LSQ_ADapply_vector_cellfield
end type

contains
!=======================================================================
!=================== GRADIENT OPERATOR BASIC METHIDS ===================
!=======================================================================
	subroutine gradient_set_mesh(this, mesh)
		implicit none
		class(gradient_t), intent(inout) :: this
		type(mesh_t), intent(in), target :: mesh
		
		this%mesh => mesh
		this%dim = mesh%dim
	end subroutine	


!=======================================================================
!================ GREEN GAUSS GRADIENT OPERATOR METHODS ================
!=======================================================================
	subroutine GG_apply_scalar_cellfield(this, Field, Field_grad, USE_GHOST_CELLS)
		implicit none
		class(GG_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:)
		logical, intent(in) :: USE_GHOST_CELLS
		real(8), intent(inout), contiguous :: Field_grad(:, :)
		
		if (.not. associated(this%mesh)) return
		call compute_grad_scalarfield_gg(Field, Field_grad,&
										 this%mesh%face_left_cell,&
										 this%mesh%face_right_cell,&
										 this%mesh%face_normal,&
										 this%mesh%face_area,&
										 this%mesh%cell_volume,&
										 this%mesh%face_weight,&
										 this%mesh%ncells,&
										 this%mesh%nfaces,&
										 this%mesh%nbfaces,&
										 this%mesh%dim,&
										 USE_GHOST_CELLS)			        
	end subroutine

	subroutine GG_apply_vector_cellfield(this, Field, Field_grad, nvars, USE_GHOST_CELLS)
		implicit none
		class(GG_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: nvars
		real(8), intent(inout), contiguous :: Field_grad(:, :)
		
		if (.not. associated(this%mesh)) return
		
		call compute_grad_vectorfield_gg(Field, Field_grad,&
										 this%mesh%face_left_cell,&
										 this%mesh%face_right_cell,&
										 this%mesh%face_normal,&
										 this%mesh%face_area,&
										 this%mesh%cell_volume,&
										 this%mesh%face_weight,&
										 this%mesh%ncells,&
										 this%mesh%nfaces,&
										 this%mesh%nbfaces,&
										 this%mesh%dim, nvars,&
										 USE_GHOST_CELLS)						 
	end subroutine
	
	subroutine GG_apply_scalar_cell(this, cell_idx, Field, Cell_grad, USE_GHOST_CELLS)
		implicit none
		class(GG_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: cell_idx
		real(8), intent(inout), contiguous :: Cell_grad(:)
		
		if (.not. associated(this%mesh)) return
		
		call compute_grad_scalarcell_gg(Field, Cell_grad, cell_idx,&
										this%mesh%ncells,&
									    this%mesh%dim,&
										this%mesh%face_left_cell,&
										this%mesh%face_right_cell,&
										this%mesh%cell_faces,&
										this%mesh%cell_faces_ptr,&
										this%mesh%face_normal,&
										this%mesh%face_area,&
										this%mesh%cell_volume,&
										this%mesh%face_weight,&
										USE_GHOST_CELLS)
	end subroutine

	subroutine GG_apply_vector_cell(this, cell_idx, Field, Cell_grad, nvars, USE_GHOST_CELLS)
		implicit none
		class(GG_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: cell_idx
		integer, intent(in) :: nvars
		real(8), intent(inout), contiguous :: Cell_grad(:)
		
		if (.not. associated(this%mesh)) return
		call compute_grad_vectorcell_gg(Field, Cell_grad, cell_idx,&
										this%mesh%ncells,&
										this%mesh%dim, nvars,&
										this%mesh%face_left_cell,&
										this%mesh%face_right_cell,&
										this%mesh%cell_faces,&
										this%mesh%cell_faces_ptr,&
										this%mesh%face_normal,&
										this%mesh%face_area,&
										this%mesh%cell_volume,&
										this%mesh%face_weight,&
										USE_GHOST_CELLS)

	end subroutine



	subroutine GG_ADapply_scalar_cellfield(this, dField, dField_grad, USE_GHOST_CELLS)
		implicit none
		class(GG_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: dField(:)
		logical, intent(in) :: USE_GHOST_CELLS
		real(8), intent(inout), contiguous :: dField_grad(:, :)
		
		if (.not. associated(this%mesh)) return
		
		call COMPUTE_GRAD_SCALARFIELD_GG_D(dField, dField_grad,&
										   this%mesh%face_left_cell,&
										   this%mesh%face_right_cell,&
										   this%mesh%face_normal,&
										   this%mesh%face_area,&
										   this%mesh%cell_volume,&
										   this%mesh%face_weight,&
										   this%mesh%ncells,&
										   this%mesh%nfaces,&
										   this%mesh%nbfaces,&
										   this%mesh%dim,&
										   USE_GHOST_CELLS)				        
	end subroutine

	subroutine GG_ADapply_vector_cellfield(this, dField, dField_grad, nvars, USE_GHOST_CELLS)
		implicit none
		class(GG_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: dField(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: nvars
		real(8), intent(inout), contiguous :: dField_grad(:, :)
		
		if (.not. associated(this%mesh)) return
		
		call COMPUTE_GRAD_VECTORFIELD_GG_D(dField, dField_grad,&
										   this%mesh%face_left_cell,&
										   this%mesh%face_right_cell,&
										   this%mesh%face_normal,&
										   this%mesh%face_area,&
										   this%mesh%cell_volume,&
										   this%mesh%face_weight,&
										   this%mesh%ncells,&
										   this%mesh%nfaces,&
										   this%mesh%nbfaces,&
										   this%mesh%dim, nvars,&
										   USE_GHOST_CELLS)						 
	end subroutine
	
	
!=======================================================================
!============== LEAST SQUARES GRADIENT OPERATOR METHODS ================
!=======================================================================
	subroutine LSQ_apply_scalar_cellfield(this, Field, Field_grad, USE_GHOST_CELLS)
		implicit none
		class(LSQ_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:)
		logical, intent(in) :: USE_GHOST_CELLS
		real(8), intent(inout), contiguous :: Field_grad(:, :)
		
		if (.not. associated(this%mesh)) return
		
		if (.not. this%is_updated) then
			if (.not. allocated(this%lsq_w)) allocate(this%lsq_w(this%mesh%dim, size(this%mesh%cell_faces)))
			
			call precompute_lsq_weights(this%mesh%dim,&
										this%mesh%ncells,&
										this%mesh%cell_faces_ptr,&
										this%mesh%cell_faces,&
										this%mesh%face_left_cell,&
										this%mesh%face_right_cell,&
										this%mesh%cell_center,&
										this%mesh%face_center, this%lsq_w)
            this%is_updated = .true.
		end if
		
		call compute_grad_scalarfield_lsq(Field, Field_grad,&
										  this%mesh%ncells,&
										  this%mesh%dim,&
										  this%mesh%cell_faces_ptr,&
										  this%mesh%cell_faces,&
                                          this%mesh%face_left_cell,&
                                          this%mesh%face_right_cell,&
                                          this%lsq_w)
	end subroutine

	subroutine LSQ_apply_vector_cellfield(this, Field, Field_grad, nvars, USE_GHOST_CELLS)
		implicit none
		class(LSQ_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: nvars
		real(8), intent(inout), contiguous :: Field_grad(:, :)
		
		if (.not. associated(this%mesh)) return
		
		if (.not. this%is_updated) then
			if (.not. allocated(this%lsq_w)) allocate(this%lsq_w(this%mesh%dim, size(this%mesh%cell_faces)))
			
			call precompute_lsq_weights(this%mesh%dim,&
										this%mesh%ncells,&
										this%mesh%cell_faces_ptr,&
										this%mesh%cell_faces,&
										this%mesh%face_left_cell,&
										this%mesh%face_right_cell,&
										this%mesh%cell_center,&
										this%mesh%face_center, this%lsq_w)
            this%is_updated = .true.
		end if
                                    
		call compute_grad_vectorfield_lsq(Field, Field_grad,&
										  this%mesh%ncells, nvars,&
										  this%mesh%dim,&
										  this%mesh%cell_faces_ptr,&
										  this%mesh%cell_faces,&
										  this%mesh%face_left_cell,&
										  this%mesh%face_right_cell,&
										  this%lsq_w)
	end subroutine
	
	subroutine LSQ_apply_scalar_cell(this, cell_idx, Field, Cell_grad, USE_GHOST_CELLS)
		implicit none
		class(LSQ_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: cell_idx
		real(8), intent(inout), contiguous :: Cell_grad(:)
		
		if (.not. associated(this%mesh)) return
													
		call compute_grad_scalarcell_lsq(Field, Cell_grad,&
										 cell_idx, this%mesh%ncells,&
									     this%mesh%dim,&
									     this%mesh%cell_faces_ptr,&
									     this%mesh%cell_faces,&
									     this%mesh%face_left_cell,&
									     this%mesh%face_right_cell,&
									     this%mesh%cell_center,&
									     this%mesh%face_center)
	end subroutine

	subroutine LSQ_apply_vector_cell(this, cell_idx, Field, Cell_grad, nvars, USE_GHOST_CELLS)
		implicit none
		class(LSQ_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: Field(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: cell_idx
		integer, intent(in) :: nvars
		real(8), intent(inout), contiguous :: Cell_grad(:)
		
		if (.not. associated(this%mesh)) return
		call compute_grad_vectorcell_lsq(Field, Cell_grad,&
										 cell_idx, this%mesh%ncells,&
										 this%mesh%dim, nvars,&
										 this%mesh%cell_faces_ptr,&
										 this%mesh%cell_faces,&
										 this%mesh%face_left_cell,&
										 this%mesh%face_right_cell,&
										 this%mesh%cell_center,&
										 this%mesh%face_center)
	end subroutine

	
	
	subroutine LSQ_ADapply_scalar_cellfield(this, dField, dField_grad, USE_GHOST_CELLS)
		implicit none
		class(LSQ_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: dField(:)
		logical, intent(in) :: USE_GHOST_CELLS
		real(8), intent(inout), contiguous :: dField_grad(:, :)
		
		if (.not. associated(this%mesh)) return
		
		if (.not. this%is_updated) then
			if (.not. allocated(this%lsq_w)) allocate(this%lsq_w(this%mesh%dim, size(this%mesh%cell_faces)))
			
			call precompute_lsq_weights(this%mesh%dim,&
										this%mesh%ncells,&
										this%mesh%cell_faces_ptr,&
										this%mesh%cell_faces,&
										this%mesh%face_left_cell,&
										this%mesh%face_right_cell,&
										this%mesh%cell_center,&
										this%mesh%face_center, this%lsq_w)
            this%is_updated = .true.
		end if
		
		call COMPUTE_GRAD_SCALARFIELD_LSQ_D(dField, dField_grad,&
										  this%mesh%ncells,&
										  this%mesh%dim,&
										  this%mesh%cell_faces_ptr,&
										  this%mesh%cell_faces,&
                                          this%mesh%face_left_cell,&
                                          this%mesh%face_right_cell,&
                                          this%lsq_w)
	end subroutine

	subroutine LSQ_ADapply_vector_cellfield(this, dField, dField_grad, nvars, USE_GHOST_CELLS)
		implicit none
		class(LSQ_gradient_t), intent(inout) :: this
		real(8), intent(in), contiguous :: dField(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in) :: nvars
		real(8), intent(inout), contiguous :: dField_grad(:, :)
		
		if (.not. associated(this%mesh)) return
		
		if (.not. this%is_updated) then
			if (.not. allocated(this%lsq_w)) allocate(this%lsq_w(this%mesh%dim, size(this%mesh%cell_faces)))
			
			call precompute_lsq_weights(this%mesh%dim,&
										this%mesh%ncells,&
										this%mesh%cell_faces_ptr,&
										this%mesh%cell_faces,&
										this%mesh%face_left_cell,&
										this%mesh%face_right_cell,&
										this%mesh%cell_center,&
										this%mesh%face_center, this%lsq_w)
            this%is_updated = .true.
		end if
                                    
		call COMPUTE_GRAD_VECTORFIELD_LSQ_D(dField, dField_grad,&
										    this%mesh%ncells, nvars,&
										    this%mesh%dim,&
										    this%mesh%cell_faces_ptr,&
										    this%mesh%cell_faces,&
										    this%mesh%face_left_cell,&
										    this%mesh%face_right_cell,&
										    this%lsq_w)
	end subroutine
	
end module
