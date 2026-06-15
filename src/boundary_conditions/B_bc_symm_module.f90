module bc_symm_module
use bc_base_module
use bc_apply_faces_module
use bc_apply_ghost_cells_module

use ADdiff_bc_kernels_module
use ADdiff_bc_ghost_kernels_module
use physical_properties_module
use stability_operator_module, only: get_Aq_matrix
implicit none
!=======================================================================
!=============== SYMM BOUNDARY CONDITION DATA STRUCTURE ================
!=======================================================================
type, extends(bc_base_t), public :: bc_symm_t

contains
	procedure :: apply_patch => apply_symm
	procedure :: ADdiff_apply_patch => ADdiff_apply_symm
	procedure :: apply_ghosts_patch => apply_symm_ghosts
	procedure :: ADdiff_apply_ghosts_patch => ADdiff_apply_symm_ghosts
	procedure :: apply_jmatrix_inv_patch => apply_symm_jmatrix_inv
	procedure :: apply_jmatrix_visc_patch => apply_symm_jmatrix_visc
end type

contains
!=======================================================================
!=================== SYMM BOUNDARY CONDITION METHODS ===================
!=======================================================================
    subroutine apply_symm(this, R_gas, k, values_ptr, values_grad_ptr, name, update_grad)
		implicit none
		class(bc_symm_t), intent(inout) :: this
		real(8), intent(in) :: R_gas, k
		real(8), intent(inout), contiguous :: values_ptr(:, :)
		real(8), intent(inout), contiguous :: values_grad_ptr(:, :)
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
        
        integer :: dim, ncells
        
		dim = this%dim
		ncells = this%mesh%ncells
		
        select case(trim(name))
        case('Vel')														!=== VELOCITY FIELD ===
			call apply_slip_vector_loop(dim, ncells, this%nfaces, 1, dim, update_grad,&
									    values_ptr, values_grad_ptr,&
									    this%b_indices, this%cell_indices, this%face_indices,&
									    this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
			
        case('Prs')														!=== PRESSURE FIELD ===
			call apply_extrapolation0_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
										   values_ptr, values_grad_ptr,&
										   this%b_indices, this%cell_indices, this%face_indices)
		    
        case('Temp')													!=== TEMPERATURE FIELD ===
            call apply_fixed_gradient_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
										   [0.d0], values_ptr, values_grad_ptr,&
									       this%b_indices, this%cell_indices, this%face_indices,&
									       this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
									       
        case('Q')														!=== PRIMITIVE VARIABLES [P U V W T] FIELD ===
			call apply_symm_Q_loop(dim, ncells, this%nfaces, update_grad,&
								   values_ptr, values_grad_ptr,&
								   this%b_indices, this%cell_indices, this%face_indices,&
								   this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
											 
        case default
			print*, 'UNSUPPORTABLE FIELD FOR BC UPDATING...'
			stop
        end select
    end subroutine
    
    subroutine apply_symm_ghosts(this, values_ptr, values_grad_ptr, name, update_grad)
		implicit none
		class(bc_symm_t), intent(inout) :: this
		real(8), intent(inout), contiguous :: values_ptr(:, :)
		real(8), intent(inout), contiguous :: values_grad_ptr(:, :)
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
        
        integer :: dim, ncells
        
		dim = this%dim
		ncells = this%mesh%ncells
		
        select case(trim(name))
        case('Vel')														!=== VELOCITY FIELD ===
			call apply_slip_vector_ghost_loop(dim, ncells, this%nfaces, 1, dim, update_grad,&
										  	  values_ptr, values_grad_ptr,&
											  this%b_indices, this%cell_indices, this%face_indices,&
											  this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal) 
			
        case('Prs')														!=== PRESSURE FIELD ===
			call apply_extrapolation0_ghost_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
											     values_ptr, values_grad_ptr,&
											     this%b_indices, this%cell_indices, this%face_indices)
		    
        case('Temp')													!=== TEMPERATURE FIELD ===
            call apply_fixed_gradient_ghost_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
											     [0.d0], values_ptr, values_grad_ptr,&
											     this%b_indices, this%cell_indices, this%face_indices,&
											     this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
        case('Q')														!=== PRIMITIVE VARIABLES [P U V W T] FIELD ===
			call apply_symm_Q_ghost_loop(dim, ncells, this%nfaces, update_grad,&
									     values_ptr, values_grad_ptr,&
									     this%b_indices, this%cell_indices, this%face_indices,&
									     this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
        case default
			print*, 'UNSUPPORTABLE FIELD FOR BC UPDATING...'
			stop
        end select
    end subroutine
    
    subroutine apply_symm_jmatrix_inv(this, R_gas, k, values_ptr, values, map_LB, name)
		implicit none
		class(bc_symm_t), intent(inout) :: this
		real(8), intent(in) :: R_gas, k
		real(8), intent(in), contiguous :: values_ptr(:, :)
		integer, intent(in), contiguous :: map_LB(:)
		real(8), intent(inout), contiguous :: values(:, :, :)
		character(len=*), intent(in) :: name
		
		integer :: dim, ncells, off
        integer :: i, d, face_idx, b_idx
        real(8) :: n(1:this%dim), area
		real(8) :: Aq(this%mesh%dim+2, this%mesh%dim+2), mask(this%mesh%dim+2)
		real(8) :: P_b, V_b(this%mesh%dim), T_b, Ro_b
        		
		dim = this%dim
		ncells = this%mesh%ncells
		
		select case(trim(name))
		case('Q')														!=== PRIMITIVE VARIABLES [P U V W T] FIELD ===
			!=== CREATING MASK Q_b = Q_b(Q_L): ===
			mask = 0.d0
			mask(1) = 1.d0
			mask(dim+2) = 1.d0
			
			!=== UPDATING VALUES: ===
			do i = 1, this%nfaces					
				b_idx = this%b_indices(i)
				face_idx = this%face_indices(i)
				
				n = this%mesh%face_normal(:, face_idx)
				area = this%mesh%face_area(face_idx)
				
				P_b = values_ptr(1, ncells+b_idx)
				V_b = values_ptr(2:dim+1, ncells+b_idx)
				T_b = values_ptr(dim+2, ncells+b_idx)
				Ro_b = P_b/(R_gas*T_b)
				
				call get_Aq_matrix(dim, R_gas, k, V_b, n, Ro_b, T_b, Aq)
				
				do d = 1, dim+2
					Aq(:, d) = Aq(:, d)*mask(d)
				end do
				
				values(:, :, map_LB(b_idx)) = values(:, :, map_LB(b_idx)) + Aq*area				  		
			end do
				
		
		case default
			print*, 'UNSUPPORTABLE FIELD FOR BC UPDATING...'
			stop
		end select
	end subroutine
	
	subroutine apply_symm_jmatrix_visc(this, R_gas, k, cp, Pr, values_ptr, values, map_LB, name)
		implicit none
		class(bc_symm_t), intent(inout) :: this
		real(8), intent(in) :: R_gas, k, cp, Pr
		real(8), intent(in), contiguous :: values_ptr(:, :)
		integer, intent(in), contiguous :: map_LB(:)
		real(8), intent(inout), contiguous :: values(:, :, :)
		character(len=*), intent(in) :: name
		
		integer :: dim, ncells, off
        integer :: i, d, face_idx, cell_idx, b_idx
        real(8) :: n(1:this%dim), ksi(1:this%dim), area
		real(8) :: Mtrx_L(this%mesh%dim+2, this%mesh%dim+2), mask(this%mesh%dim+2)
		real(8) :: P_b, V_b(this%mesh%dim), T_b, Ro_b, A_b
		real(8) :: teta, temp, RL
		real(8) :: mu, lambda
        		
		dim = this%dim
		ncells = this%mesh%ncells
		
		select case(trim(name))
		case('Q')														!=== PRIMITIVE VARIABLES [P U V W T] FIELD ===
			!=== CREATING MASK Q_b = Q_b(Q_L): ===
			mask = 0.d0			
		
		case default
			print*, 'UNSUPPORTABLE FIELD FOR BC UPDATING...'
			stop
		end select
	end subroutine
	
	SUBROUTINE ADdiff_apply_symm(this, R_gas, k, values_ptr, values_ptrd,&
								 values_grad_ptr, values_grad_ptrd, name, update_grad)
		IMPLICIT NONE
		CLASS(bc_symm_t), INTENT(INOUT) :: this
		REAL(8), INTENT(IN) :: R_gas, k
		REAL(8), INTENT(IN), CONTIGUOUS :: values_ptr(:, :)
		REAL(8), INTENT(IN), CONTIGUOUS :: values_grad_ptr(:, :)
		REAL(8), INTENT(INOUT), CONTIGUOUS :: values_ptrd(:, :)
		REAL(8), INTENT(INOUT), CONTIGUOUS :: values_grad_ptrd(:, :)
		CHARACTER(LEN=*), INTENT(IN) :: name
		LOGICAL, INTENT(IN) :: update_grad
		        
        INTEGER :: dim, ncells, off
        INTEGER :: i, d, cell_idx, face_idx, b_idx
        REAL(8) :: n(1:this%dim), ksi(1:this%dim), dist
        		
		dim = this%dim
		ncells = this%mesh%ncells
		DO i=1,this%nfaces
			b_idx = this%b_indices(i)
			cell_idx = this%cell_indices(i)
			face_idx = this%face_indices(i)
			n = this%mesh%face_normal(:, face_idx)
			ksi = this%mesh%face_center(:, face_idx) - this%mesh%cell_center(:, cell_idx)
			dist = 0.d0
			DO d=1,dim
				dist = dist + ksi(d)*n(d)
			END DO
			!=== PRESSURE: ===
			CALL APPLY_EXTRAPOLATION0_BC_D(dim, values_ptrd(1, cell_idx),&
										   values_grad_ptrd(1:dim, cell_idx),&
			&                              values_ptrd(1, ncells+b_idx), &
			&                              values_grad_ptrd(1:dim, ncells+b_idx), &
			&                              update_grad)
			!=== VELOCITY: ===
			off = dim*dim
			CALL APPLY_SLIP_VECTOR_BC_D(dim, n, values_ptrd(2:dim+1, cell_idx),&
			&                           values_grad_ptrd(1+dim:dim+off, cell_idx), &
			&                           values_ptrd(2:dim+1, ncells+b_idx), &
			&                           values_grad_ptrd(1+dim:dim+off, ncells+b_idx), dist, update_grad)
			!=== TEMPERATURE: ===
			off = dim*(dim+1)
			CALL APPLY_EXTRAPOLATION0_BC_D(dim, values_ptrd(dim+2, cell_idx), &
			&                              values_grad_ptrd(1+off:dim+off, cell_idx),&
			&                              values_ptrd(dim+2, ncells+b_idx), &
			&                              values_grad_ptrd(1+off:dim+off, ncells+b_idx), .false.)
			
			IF (update_grad) CALL COPY_GRAD_T_FROM_CELL_D(dim, n, &
			&                                             values_grad_ptrd(1+off:dim+off, cell_idx), &
			&                                             values_grad_ptrd(1+off:dim+off, ncells+b_idx))
		END DO
        
    END SUBROUTINE
	
	SUBROUTINE ADdiff_apply_symm_ghosts(this, values_ptr, values_ptrd,&
								 values_grad_ptr, values_grad_ptrd, name, update_grad)
		IMPLICIT NONE
		CLASS(bc_symm_t), INTENT(INOUT) :: this
		REAL(8), INTENT(IN), CONTIGUOUS :: values_ptr(:, :)
		REAL(8), INTENT(IN), CONTIGUOUS :: values_grad_ptr(:, :)
		REAL(8), INTENT(INOUT), CONTIGUOUS :: values_ptrd(:, :)
		REAL(8), INTENT(INOUT), CONTIGUOUS :: values_grad_ptrd(:, :)
		CHARACTER(LEN=*), INTENT(IN) :: name
		LOGICAL, INTENT(IN) :: update_grad
		        
        INTEGER :: dim, ncells, off
        INTEGER :: i, d, cell_idx, face_idx, b_idx
        REAL(8) :: n(1:this%dim), ksi(1:this%dim), dist
        		
		dim = this%dim
		ncells = this%mesh%ncells
		DO i=1,this%nfaces
			b_idx = this%b_indices(i)
			cell_idx = this%cell_indices(i)
			face_idx = this%face_indices(i)
			n = this%mesh%face_normal(:, face_idx)
			ksi = this%mesh%face_center(:, face_idx) - this%mesh%cell_center(:, cell_idx)
			dist = 0.d0
			DO d=1,dim
				dist = dist + ksi(d)*n(d)
			END DO
			!=== PRESSURE: ===
			CALL APPLY_EXTRAPOLATION0_GHOST_BC_D(dim, values_ptrd(1, cell_idx),&
										   values_grad_ptrd(1:dim, cell_idx),&
			&                              values_ptrd(1, ncells+b_idx), &
			&                              values_grad_ptrd(1:dim, ncells+b_idx), &
			&                              update_grad)
			!=== VELOCITY: ===
			off = dim*dim
			CALL APPLY_SLIP_VECTOR_GHOST_BC_D(dim, n, values_ptrd(2:dim+1, cell_idx),&
			&                           values_grad_ptrd(1+dim:dim+off, cell_idx), &
			&                           values_ptrd(2:dim+1, ncells+b_idx), &
			&                           values_grad_ptrd(1+dim:dim+off, ncells+b_idx), update_grad)
			!=== TEMPERATURE: ===
			off = dim*(dim+1)
			CALL APPLY_FIXED_GRADIENT_GHOST_BC_D(dim, n, values_ptrd(dim+2, cell_idx), &
												 values_grad_ptrd(1+off:dim+off, cell_idx),&
												 values_ptrd(dim+2, ncells+b_idx),&
												 values_grad_ptrd(1+off:dim+off, ncells+b_idx), update_grad)			
		END DO
        
    END SUBROUTINE
	
end module
