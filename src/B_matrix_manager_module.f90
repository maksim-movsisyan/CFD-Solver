module matrix_manager_module
use mesh_module, only: mesh_t
use field_manager_module, only: field_manager_t
use linear_operator_module, only: bcsr_operator_t
use bc_manager_module, only: bc_manager_t
use physical_properties_manager_module, only: phys_prop_manager_t
use solver_config_module, only: db_solver_config_t
use matrix_initialization_module
use db_mtrx_assembling_module

implicit none
!=======================================================================
!==================== FLUXES MANAGER DATA STRUCTURE ====================
!=======================================================================
type, public :: matrix_manager_t
    type(mesh_t), pointer :: mesh => null()
    type(field_manager_t), pointer :: fldsm => null()
    type(bc_manager_t), pointer :: bcm => null()
    type(phys_prop_manager_t), pointer :: ppm => null()
    type(db_solver_config_t), pointer :: settings => null()
    real(8), pointer, contiguous :: ptime_step(:) => null()
    
    type(bcsr_operator_t), pointer :: mtrx => null()
    real(8), pointer, contiguous :: vars_ptr(:, :) => null()
    
    integer, allocatable :: map_LL(:), map_LR(:), map_RL(:),&			!=== MAPPING FOR MATRIX ASSEMBLING ===
							map_RR(:), map_LB(:)						!=== DIM = (nfaces - nbfaces), DIMB = (nbfaces) ===
							
    
    integer :: QIDX, GRADQIDX
   
contains
    procedure :: initialize => mtrxm_initialize
    procedure :: initialize_bcsr_mtrx => mtrxm_initialize_bcsr_mtrx
    procedure :: db_assmeble_bcsr_mtrx => mtrxm_db_assmeble_bcsr_mtrx
end type

contains
!=======================================================================
!======================= FLUXES MANAGER METHODS ========================
!=======================================================================
	subroutine mtrxm_initialize(this, mesh, fldsm, bcm, ppm, settings, mtrx, ptime_step, QIDX, GRADQIDX)
		implicit none
		class(matrix_manager_t), intent(inout) :: this
		type(mesh_t), intent(in), target:: mesh
		type(field_manager_t), intent(in), target :: fldsm
		type(bc_manager_t), intent(in), target :: bcm
		type(phys_prop_manager_t), intent(in), target :: ppm
		type(db_solver_config_t), intent(in), target :: settings
		type(bcsr_operator_t), intent(in), target :: mtrx
		real(8), intent(in), target :: ptime_step(:)
		integer, intent(in) :: QIDX, GRADQIDX
		
		this%mesh => mesh
		this%fldsm => fldsm
		this%bcm => bcm
		this%ppm => ppm
		this%settings => settings
		this%mtrx => mtrx
		this%ptime_step => ptime_step
		
		
		this%QIDX = QIDX
		this%GRADQIDX = GRADQIDX
		
		this%vars_ptr => fldsm%registry(this%QIDX)%values
	end subroutine

	subroutine mtrxm_initialize_bcsr_mtrx(this)
		implicit none
		class(matrix_manager_t), intent(inout) :: this
		
		this%mtrx%bs = this%mesh%dim+2		 
		call bcsr_mtrx_initialize1(this%mesh%ncells, this%mesh%nfaces,&
								   this%mesh%nbfaces, this%mtrx%bs,&
								   this%mesh%cell_faces_ptr, this%mesh%cell_faces,&
								   this%mesh%face_left_cell, this%mesh%face_right_cell,&
								   this%mtrx%n, this%mtrx%nnz,&
								   this%mtrx%values,&
								   this%mtrx%col_indices,&
								   this%mtrx%row_ptr,&
								   this%mtrx%diag_indices)
								   
		call db_create_mtrx_fill_map(this%mesh%nfaces, this%mesh%nbfaces,&
									 this%mesh%face_left_cell, this%mesh%face_right_cell,&
									 this%mesh%face_bidx, this%mtrx%col_indices,&
									 this%mtrx%row_ptr,&
									 this%map_LL, this%map_LR, this%map_RL,&
									 this%map_RR, this%map_LB)
	end subroutine
	
	subroutine mtrxm_db_assmeble_bcsr_mtrx(this)
		implicit none
		class(matrix_manager_t), intent(inout) :: this
		
		if (.not. associated(this%mtrx)) return
		
		if (.not. this%settings%USE_CONSERVATIVE_VARS) then
			!=== PRIMITIVE VARIABLES: ===
			call db_assemble_bcsr_matrix_inv(this%mesh%dim,&
											 this%mesh%ncells,&
											 this%mesh%nfaces,&
											 this%mesh%nbfaces,&
											 this%mesh%face_left_cell,&
											 this%mesh%face_right_cell,&
											 this%mesh%face_bidx,&
											 this%mesh%face_normal,&
											 this%mesh%face_area,&
											 this%mesh%face_weight,&
											 this%mesh%cell_volume,&
											 this%ptime_step,&
											 this%vars_ptr,&
											 this%mtrx%values,&
											 this%mtrx%diag_indices,&
											 this%map_LL, this%map_LR,&
											 this%map_RL, this%map_RR, this%map_LB,&
											 this%ppm%k, this%ppm%R_gas,&
											 this%ppm%cp, this%ppm%cv,&
											 this%settings%STABILITY_OPERATOR_TYPE,&
											 this%settings%USE_GHOST_CELLS, .true.)
			if (this%settings%MODEL == 2) then
				call db_assemble_bcsr_matrix_visc(this%mesh%dim,&
												  this%mesh%ncells,&
												  this%mesh%nfaces,&
												  this%mesh%nbfaces,&
												  this%mesh%face_left_cell,&
												  this%mesh%face_right_cell,&
												  this%mesh%face_bidx,&
												  this%mesh%face_normal,&
												  this%mesh%face_area,&
												  this%mesh%face_weight,&   
												  this%mesh%cell_center,&  
												  this%vars_ptr,&
												  this%mtrx%values,&
												  this%map_LL, this%map_LR,&
												  this%map_RL, this%map_RR, this%map_LB,&
												  this%ppm%k, this%ppm%R_gas,&
												  this%ppm%cp, this%ppm%cv, this%ppm%Pr,&   
												  this%settings%USE_GHOST_CELLS, .false.)
			end if
		else
			!=== CONSERVATIVE VARIABLES: ===
			call db_assemble_bcsr_matrix_inv_cv(this%mesh%dim,&
											    this%mesh%ncells,&
											    this%mesh%nfaces,&
											    this%mesh%nbfaces,&
											    this%mesh%face_left_cell,&
											    this%mesh%face_right_cell,&
											    this%mesh%face_bidx,&
											    this%mesh%face_normal,&
											    this%mesh%face_area,&
											    this%mesh%face_weight,&
											    this%mesh%cell_volume,&
											    this%ptime_step,&
											    this%vars_ptr,&
											    this%mtrx%values,&
											    this%mtrx%diag_indices,&
											    this%map_LL, this%map_LR,&
											    this%map_RL, this%map_RR, this%map_LB,&
											    this%ppm%k, this%ppm%R_gas,&
											    this%ppm%cp, this%ppm%cv,&
											    this%settings%STABILITY_OPERATOR_TYPE,&
											    this%settings%USE_GHOST_CELLS, .true.)
			if (this%settings%MODEL == 2) then
				call db_assemble_bcsr_matrix_visc_cv(this%mesh%dim,&
												     this%mesh%ncells,&
												     this%mesh%nfaces,&
												     this%mesh%nbfaces,&
												     this%mesh%face_left_cell,&
												     this%mesh%face_right_cell,&
												     this%mesh%face_bidx,&
												     this%mesh%face_normal,&
												     this%mesh%face_area,&
												     this%mesh%face_weight,&   
												     this%mesh%cell_center,&  
												     this%vars_ptr,&
												     this%mtrx%values,&
												     this%map_LL, this%map_LR,&
												     this%map_RL, this%map_RR, this%map_LB,&
												     this%ppm%k, this%ppm%R_gas,&
												     this%ppm%cp, this%ppm%cv, this%ppm%Pr,&   
												     this%settings%USE_GHOST_CELLS, .false.)
			end if
		end if
		
		!=== BOUNDARY CONDITIONS: ===
		call this%bcm%add_boundary_jacobians(this%mtrx%values, this%map_LB,&
											 'Q', this%settings%MODEL,&
											 this%settings%USE_GHOST_CELLS,&
											 this%settings%USE_CONSERVATIVE_VARS,&
											 this%QIDX)
											 
	end subroutine
	
end module
