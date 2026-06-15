module bc_wall_module
use bc_base_module
use bc_apply_faces_module
use bc_apply_ghost_cells_module

use ADdiff_bc_kernels_module
use ADdiff_bc_ghost_kernels_module
use physical_properties_module
use integral_characteristics_module
use stability_operator_module, only: get_Aq_matrix
implicit none
!=======================================================================
!=============== WALL BOUNDARY CONDITION DATA STRUCTURE ================
!=======================================================================
type, extends(bc_base_t), public :: bc_wall_t
	real(8) :: wall_velocity(3) = 0.d0
	real(8) :: wall_temperature = 300.0
	real(8) :: wall_heat_flux = 0.d0
contains
	procedure :: apply_patch => apply_wall
	procedure :: ADdiff_apply_patch => ADdiff_apply_wall
	procedure :: apply_ghosts_patch => apply_wall_ghosts
	procedure :: ADdiff_apply_ghosts_patch => ADdiff_apply_wall_ghosts
	procedure :: apply_jmatrix_inv_patch => apply_wall_jmatrix_inv
	procedure :: apply_jmatrix_visc_patch => apply_wall_jmatrix_visc
	
	procedure :: compute_force => compute_force_wall
	procedure :: compute_fluxes => compute_fluxes_wall
end type

contains
!=======================================================================
!=================== WALL BOUNDARY CONDITION METHODS ===================
!=======================================================================
    subroutine apply_wall(this, R_gas, k, values_ptr, values_grad_ptr, name, update_grad)
		implicit none
		class(bc_wall_t), intent(inout) :: this
		real(8), intent(in) :: R_gas, k
		real(8), intent(inout), contiguous :: values_ptr(:, :)
		real(8), intent(inout), contiguous :: values_grad_ptr(:, :)
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
        
        integer :: dim, ncells
        real(8) :: fxd_values(5)
		
		dim = this%dim
		ncells = this%mesh%ncells
		
        select case(trim(name))
        case('Vel')														!=== VELOCITY FIELD ===
			select case (this%bc_subtype)
			case(BC_WALL_SLIP, BC_WALL_SLIP_REIMAN)
				call apply_slip_vector_loop(dim, ncells, this%nfaces, 1, dim, update_grad,&
											values_ptr, values_grad_ptr,&
											this%b_indices, this%cell_indices, this%face_indices,&
											this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
				
			case default
				call apply_fixed_value_loop(dim, ncells, this%nfaces, 1, dim, update_grad,&
											this%wall_velocity(1:dim), values_ptr, values_grad_ptr,&
											this%b_indices, this%cell_indices, this%face_indices,&
											this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
				
			end select

        case('Prs')														!=== PRESSURE FIELD ===
			call apply_extrapolation0_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
										   values_ptr, values_grad_ptr,&
										   this%b_indices, this%cell_indices, this%face_indices)				
		    
        case('Temp')													!=== TEMPERATURE FIELD ===
            select case (this%bc_subtype)
			case(BC_WALL_ISOTHERMAL)
				fxd_values(1) = this%wall_temperature
				call apply_fixed_value_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
											fxd_values, values_ptr, values_grad_ptr,&
											this%b_indices, this%cell_indices, this%face_indices,&
											this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
			
			case(BC_WALL_HEAT_FLUX, BC_WALL_ADIABATIC)
				fxd_values(1) = this%wall_heat_flux
				call apply_fixed_gradient_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
											   fxd_values, values_ptr, values_grad_ptr,&
											   this%b_indices, this%cell_indices, this%face_indices,&
											   this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
				
			case default
				call apply_extrapolation0_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
											   values_ptr, values_grad_ptr,&
											   this%b_indices, this%cell_indices, this%face_indices)
				
			end select
            
        case('Q')														!=== PRIMITIVE VARIABLES [P U V W T] FIELD ===
			select case (this%bc_subtype)
			case(BC_WALL_ISOTHERMAL)
				fxd_values(1:dim) = this%wall_velocity(1:dim)
				fxd_values(dim+1) = this%wall_temperature
				call apply_isothermal_noslip_wall_Q_loop(dim, ncells, this%nfaces, update_grad,&
														 fxd_values, values_ptr, values_grad_ptr,&
														 this%b_indices, this%cell_indices, this%face_indices,&
														 this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
				
			case(BC_WALL_HEAT_FLUX, BC_WALL_ADIABATIC)
				fxd_values(1:dim) = this%wall_velocity(1:dim)
				fxd_values(dim+1) = this%wall_heat_flux
				call apply_fixed_flux_noslip_wall_Q_loop(dim, ncells, this%nfaces, update_grad,&
														 fxd_values, values_ptr, values_grad_ptr,&
														 this%b_indices, this%cell_indices, this%face_indices,&
														 this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
			
			case(BC_WALL_SLIP)
				call apply_symm_Q_loop(dim, ncells, this%nfaces, update_grad,&
									   values_ptr, values_grad_ptr,&
									   this%b_indices, this%cell_indices, this%face_indices,&
									   this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
			
			case(BC_WALL_SLIP_REIMAN)
				call apply_reiman_slip_wall_Q_loop(dim, ncells, this%nfaces, update_grad, k, R_gas,&
												   values_ptr, values_grad_ptr,&
												   this%b_indices, this%cell_indices, this%face_indices,&
												   this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
				
			case default
				fxd_values(1:dim) = this%wall_velocity(1:dim)
				fxd_values(dim+1) = 0.d0
				call apply_fixed_flux_noslip_wall_Q_loop(dim, ncells, this%nfaces, update_grad,&
														 fxd_values, values_ptr, values_grad_ptr,&
														 this%b_indices, this%cell_indices, this%face_indices,&
														 this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
			end select
        case default
			print*, 'UNSUPPORTABLE FIELD FOR BC UPDATING...'
			stop
        end select
    end subroutine
	
	subroutine apply_wall_ghosts(this, values_ptr, values_grad_ptr, name, update_grad)
		implicit none
		class(bc_wall_t), intent(inout) :: this
		real(8), intent(inout), contiguous :: values_ptr(:, :)
		real(8), intent(inout), contiguous :: values_grad_ptr(:, :)
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
        
        integer :: dim, ncells
        real(8) :: fxd_values(5)
		
		dim = this%dim
		ncells = this%mesh%ncells
		
        select case(trim(name))
        case('Vel')														!=== VELOCITY FIELD ===
			select case (this%bc_subtype)
			case(BC_WALL_SLIP, BC_WALL_SLIP_REIMAN)
				call apply_slip_vector_ghost_loop(dim, ncells, this%nfaces, 1, dim, update_grad,&
												  values_ptr, values_grad_ptr,&
												  this%b_indices, this%cell_indices, this%face_indices,&
												  this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
				
			case default
				call apply_fixed_value_ghost_loop(dim, ncells, this%nfaces, 1, dim, update_grad,&
												  this%wall_velocity(1:dim), values_ptr, values_grad_ptr,&
												  this%b_indices, this%cell_indices, this%face_indices,&
												  this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
				
			end select

        case('Prs')														!=== PRESSURE FIELD ===
			call apply_extrapolation0_ghost_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
											     values_ptr, values_grad_ptr,&
											     this%b_indices, this%cell_indices, this%face_indices)	
		    
        case('Temp')													!=== TEMPERATURE FIELD ===
            select case (this%bc_subtype)
			case(BC_WALL_ISOTHERMAL)
				fxd_values(1) = this%wall_temperature
				call apply_fixed_value_ghost_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
											  	  fxd_values, values_ptr, values_grad_ptr,&
												  this%b_indices, this%cell_indices, this%face_indices,&
												  this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
			
			case(BC_WALL_HEAT_FLUX, BC_WALL_ADIABATIC)
				fxd_values(1) = this%wall_heat_flux
				call apply_fixed_gradient_ghost_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
												     fxd_values, values_ptr, values_grad_ptr,&
												     this%b_indices, this%cell_indices, this%face_indices,&
												     this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
				
			case default
				call apply_extrapolation0_ghost_loop(dim, ncells, this%nfaces, 1, 1, update_grad,&
												     values_ptr, values_grad_ptr,&
												     this%b_indices, this%cell_indices, this%face_indices)
				
			end select
            
        case('Q')														!=== PRIMITIVE VARIABLES [P U V W T] FIELD ===
			select case (this%bc_subtype)
			case(BC_WALL_ISOTHERMAL)
				fxd_values(1:dim) = this%wall_velocity(1:dim)
				fxd_values(dim+1) = this%wall_temperature
				call apply_isothermal_noslip_wall_Q_ghost_loop(dim, ncells, this%nfaces, update_grad,&
															   fxd_values, values_ptr, values_grad_ptr,&
															   this%b_indices, this%cell_indices, this%face_indices,&
															   this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
				
			case(BC_WALL_HEAT_FLUX, BC_WALL_ADIABATIC)
				fxd_values(1:dim) = this%wall_velocity(1:dim)
				fxd_values(dim+1) = this%wall_heat_flux
				call apply_fixed_flux_noslip_wall_Q_ghost_loop(dim, ncells, this%nfaces, update_grad,&
															   fxd_values, values_ptr, values_grad_ptr,&
															   this%b_indices, this%cell_indices, this%face_indices,&
															   this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
			
			case(BC_WALL_SLIP, BC_WALL_SLIP_REIMAN)
				call apply_symm_Q_ghost_loop(dim, ncells, this%nfaces, update_grad,&
										     values_ptr, values_grad_ptr,&
										     this%b_indices, this%cell_indices, this%face_indices,&
										     this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
		
			case default
				fxd_values(1:dim) = this%wall_velocity(1:dim)
				fxd_values(dim+1) = 0.d0
				call apply_fixed_flux_noslip_wall_Q_ghost_loop(dim, ncells, this%nfaces, update_grad,&
															   fxd_values, values_ptr, values_grad_ptr,&
															   this%b_indices, this%cell_indices, this%face_indices,&
															   this%mesh%cell_center, this%mesh%face_center, this%mesh%face_normal)
			end select
        case default
			print*, 'UNSUPPORTABLE FIELD FOR BC UPDATING...'
			stop
        end select
    end subroutine

	subroutine apply_wall_jmatrix_inv(this, R_gas, k, values_ptr, values, map_LB, name)
		implicit none
		class(bc_wall_t), intent(inout) :: this
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
			mask(1) = 1.0d0

			select case (this%bc_subtype)
			case(BC_WALL_ISOTHERMAL)
				mask(2:dim+2) = 0.0d0
			case(BC_WALL_HEAT_FLUX, BC_WALL_ADIABATIC, BC_WALL_SLIP, BC_WALL_SLIP_REIMAN)
				mask(2:dim+1) = 0.0d0
				mask(dim+2) = 1.0d0
			case default
				mask(2:dim+1) = 0.0d0
				mask(dim+2)   = 1.0d0
			end select
			
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
	
	subroutine apply_wall_jmatrix_visc(this, R_gas, k, cp, Pr, values_ptr, values, map_LB, name)
		implicit none
		class(bc_wall_t), intent(inout) :: this
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
			select case (this%bc_subtype)
			case(BC_WALL_SLIP_REIMAN, BC_WALL_SLIP)
				mask = 0.d0
			case(BC_WALL_NO_SLIP, BC_WALL_ISOTHERMAL)
				mask = 1.d0
				mask(1) = 0.d0
			case(BC_WALL_ADIABATIC, BC_WALL_HEAT_FLUX)
				mask = 1.d0
				mask(1) = 0.d0
				mask(dim+2) = 0.d0
			end select
			
			!=== UPDATING VALUES: ===
			do i = 1, this%nfaces					
				b_idx = this%b_indices(i)
				cell_idx = this%cell_indices(i)
				face_idx = this%face_indices(i)
				
				n = this%mesh%face_normal(:, face_idx)
				ksi = this%mesh%face_center(:, face_idx) - this%mesh%cell_center(:, cell_idx)
				area = this%mesh%face_area(face_idx)
				
				!=== BOUNDARY VALUES: ===
				P_b = values_ptr(1, ncells+b_idx)
				V_b = values_ptr(2:dim+1, ncells+b_idx)
				T_b = values_ptr(dim+2, ncells+b_idx)
				Ro_b = P_b/(R_gas*T_b)
				A_b = dsqrt(k*R_gas*T_b)
				
				!=== PHYSICAL TRANSPORT PROPERTIES: ===
				mu = MU_air(T_b)
				lambda = lambda_air(T_b, Pr, cp)
				
				RL = norm2(ksi)
				teta = 1.d0/A_b**2
				temp = 1.d0/(RL*(T_b*cp*teta - 1.d0))
				
				!=== BOUNDARY JACOBIAN MATRIX: ===
				Mtrx_L = 0.d0
				do d = 1, dim
					Mtrx_L(1, d+1) = -mu*V_b(d)*temp
					Mtrx_L(d+1, d+1) = mu/(RL*Ro_b)
					Mtrx_L(dim+2, d+1) = -mu*T_b*V_b(d)*teta*temp/Ro_b
				end do
				Mtrx_L(1, dim+2) = lambda*temp
				Mtrx_L(dim+2, dim+2) = lambda*T_b*teta*temp/Ro_b
					
				do d = 1, dim+2
					Mtrx_L(:, d) = Mtrx_L(:, d)*mask(d)
				end do
				
				values(:, :, map_LB(b_idx)) = values(:, :, map_LB(b_idx)) + Mtrx_L*area				  		
			end do
				
		
		case default
			print*, 'UNSUPPORTABLE FIELD FOR BC UPDATING...'
			stop
		end select
	end subroutine

	SUBROUTINE ADdiff_apply_wall(this, R_gas, k, values_ptr, values_ptrd,&
								 values_grad_ptr, values_grad_ptrd, name, update_grad)
		IMPLICIT NONE
		CLASS(bc_wall_t), INTENT(INOUT) :: this
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
		SELECT CASE (this%bc_subtype)
			CASE (BC_WALL_ISOTHERMAL)
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
					CALL EXTPRLT_PRS_D()
					CALL FXD_VEL_D()
					CALL FXD_TEMP_D()
				END DO
			CASE (BC_WALL_HEAT_FLUX, BC_WALL_ADIABATIC)
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
					CALL EXTPRLT_PRS_D()
					CALL FXD_VEL_D()
					CALL FXD_TFLUX_D()
				END DO
			CASE (BC_WALL_SLIP)
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
					CALL EXTPRLT_PRS_D()
					CALL SLP_VEL_D()
					CALL EXTPRLT_TEMP_D()
				END DO
			CASE (BC_WALL_SLIP_REIMAN)
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
					CALL RMN_PRS_D()
					CALL SLP_VEL_D()
					CALL EXTPRLT_TEMP_D()
				END DO
			CASE DEFAULT
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
					CALL EXTPRLT_PRS_D()
					CALL FXD_VEL_D()
					CALL EXTPRLT_TEMP_D()
				END DO
		END SELECT

		CONTAINS
		!  Differentiation of extprlt_prs in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		!=== PRIMITIVE VARIABLES AUX SUBROUTINES: ===
		SUBROUTINE EXTPRLT_PRS_D()
			CALL APPLY_EXTRAPOLATION0_BC_D(dim, values_ptrd(1, cell_idx),&
			&                              values_grad_ptrd(1:dim , cell_idx),&
			&                              values_ptrd(1, ncells+b_idx), &
			&                              values_grad_ptrd(1:dim, ncells+b_idx), &
			&                              update_grad)
		END SUBROUTINE

		!  Differentiation of fxd_vel in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE FXD_VEL_D()
			DO d=1,dim
				off = dim*d
				CALL APPLY_FIXED_VALUE_BC_D(dim, n, values_ptrd(1+d, cell_idx), &
				&                             values_grad_ptrd(1+off:dim+off, cell_idx),&
				&                             values_grad_ptrd(1+off:dim+off, ncells+b_idx), dist, update_grad)
			END DO
		END SUBROUTINE

		!  Differentiation of slp_vel in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE SLP_VEL_D()
			off = dim*dim
			CALL APPLY_SLIP_VECTOR_BC_D(dim, n, values_ptrd(2:dim+1, cell_idx), &
			&                           values_grad_ptrd(1+dim:dim+off, cell_idx), &
			&                           values_ptrd(2:dim+1, ncells+b_idx), &
			&                           values_grad_ptrd(1+dim:dim+off, ncells+b_idx), dist, update_grad)
		END SUBROUTINE

		!  Differentiation of fxd_temp in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE FXD_TEMP_D()
			off = dim*(dim+1)
			CALL APPLY_FIXED_VALUE_BC_D(dim, n, values_ptrd(dim+2, cell_idx), &
			&                           values_grad_ptrd(1+off:dim+off, cell_idx), &
			&                           values_grad_ptrd(1+off:dim+off, ncells+b_idx), dist, update_grad)
		END SUBROUTINE

		!  Differentiation of fxd_tflux in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE FXD_TFLUX_D()
			off = dim*(dim+1)
			CALL APPLY_FIXED_GRADIENT_BC_D(dim, n, values_ptrd(dim+2, cell_idx), &
			&                              values_grad_ptrd(1+off:dim+off, cell_idx),&
			&                              values_ptrd(dim+2, ncells+b_idx)&
			&                              , values_grad_ptrd(1+off:dim+off, ncells+ b_idx), update_grad)
		END SUBROUTINE

		!  Differentiation of extprlt_temp in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE EXTPRLT_TEMP_D()
			off = dim*(dim+1)
			CALL APPLY_EXTRAPOLATION0_BC_D(dim, values_ptrd(dim+2, cell_idx), &
			&                              values_grad_ptrd(1+off:dim+off,cell_idx),&
			&                              values_ptrd(dim+2, ncells+b_idx), &
			&                              values_grad_ptrd(1+off:dim+off, ncells+b_idx), update_grad)
		END SUBROUTINE 

		!  Differentiation of rmn_prs in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE RMN_PRS_D()
			CALL APPLY_REIMAN_PRESSURE_BC_D(dim, values_ptr(1, cell_idx), &
			&                               values_ptrd(1, cell_idx), &
			&                               values_grad_ptr(1:dim, cell_idx), &
			&                               values_grad_ptrd(1:dim, cell_idx), &
			&                               values_ptr(2:dim+1, cell_idx), &
			&                               values_ptrd(2:dim+1, cell_idx), &
			&                               values_ptr(dim+2, cell_idx), values_ptrd&
			&                               (dim+2, cell_idx), n, k, r_gas, &
			&                               values_ptrd(1, ncells+b_idx),&
			&                               values_grad_ptrd(1:dim, ncells+b_idx), update_grad)
		END SUBROUTINE
        
    END SUBROUTINE

	SUBROUTINE ADdiff_apply_wall_ghosts(this, values_ptr, values_ptrd,&
								 values_grad_ptr, values_grad_ptrd, name, update_grad)
		IMPLICIT NONE
		CLASS(bc_wall_t), INTENT(INOUT) :: this
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
		SELECT CASE (this%bc_subtype)
			CASE (BC_WALL_ISOTHERMAL)
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
					CALL EXTPRLT_PRS_D()
					CALL FXD_VEL_D()
					CALL FXD_TEMP_D()
				END DO
			CASE (BC_WALL_HEAT_FLUX, BC_WALL_ADIABATIC)
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
					CALL EXTPRLT_PRS_D()
					CALL FXD_VEL_D()
					CALL FXD_TFLUX_D()
				END DO
			CASE (BC_WALL_SLIP, BC_WALL_SLIP_REIMAN)
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
					CALL EXTPRLT_PRS_D()
					CALL SLP_VEL_D()
					CALL EXTPRLT_TEMP_D()
				END DO
			CASE DEFAULT
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
					CALL EXTPRLT_PRS_D()
					CALL FXD_VEL_D()
					CALL EXTPRLT_TEMP_D()
				END DO
		END SELECT

		CONTAINS
		!  Differentiation of extprlt_prs in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		!=== PRIMITIVE VARIABLES AUX SUBROUTINES: ===
		SUBROUTINE EXTPRLT_PRS_D()
			CALL APPLY_EXTRAPOLATION0_GHOST_BC_D(dim, values_ptrd(1, cell_idx),&
			&                              values_grad_ptrd(1:dim , cell_idx),&
			&                              values_ptrd(1, ncells+b_idx), &
			&                              values_grad_ptrd(1:dim, ncells+b_idx), &
			&                              update_grad)
		END SUBROUTINE

		!  Differentiation of fxd_vel in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE FXD_VEL_D()
			DO d=1,dim
				off = dim*d
				CALL APPLY_FIXED_VALUE_GHOST_BC_D(dim, n, values_ptrd(1+d, cell_idx), &
				&                             values_grad_ptrd(1+off:dim+off, cell_idx),&
											  values_ptrd(1+d, ncells+b_idx), &
				&                             values_grad_ptrd(1+off:dim+off, ncells+b_idx), dist, update_grad)
			END DO
		END SUBROUTINE

		!  Differentiation of slp_vel in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE SLP_VEL_D()
			off = dim*dim
			CALL APPLY_SLIP_VECTOR_GHOST_BC_D(dim, n, values_ptrd(2:dim+1, cell_idx), &
			&                           values_grad_ptrd(1+dim:dim+off, cell_idx), &
			&                           values_ptrd(2:dim+1, ncells+b_idx), &
			&                           values_grad_ptrd(1+dim:dim+off, ncells+b_idx), update_grad)
		END SUBROUTINE

		!  Differentiation of fxd_temp in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE FXD_TEMP_D()
			off = dim*(dim+1)
			CALL APPLY_FIXED_VALUE_GHOST_BC_D(dim, n, values_ptrd(dim+2, cell_idx), &
			&                           values_grad_ptrd(1+off:dim+off, cell_idx), &
										values_ptrd(dim+2, ncells+b_idx), &
			&                           values_grad_ptrd(1+off:dim+off, ncells+b_idx), dist, update_grad)
		END SUBROUTINE

		!  Differentiation of fxd_tflux in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE FXD_TFLUX_D()
			off = dim*(dim+1)
			CALL APPLY_FIXED_GRADIENT_GHOST_BC_D(dim, n, values_ptrd(dim+2, cell_idx), &
			&                              values_grad_ptrd(1+off:dim+off, cell_idx),&
			&                              values_ptrd(dim+2, ncells+b_idx)&
			&                              , values_grad_ptrd(1+off:dim+off, ncells+ b_idx), update_grad)
		END SUBROUTINE

		!  Differentiation of extprlt_temp in forward (tangent) mode:
		!   variations   of useful results: values_grad_ptr values_ptr
		!   with respect to varying inputs: values_grad_ptr values_ptr
		SUBROUTINE EXTPRLT_TEMP_D()
			off = dim*(dim+1)
			CALL APPLY_EXTRAPOLATION0_GHOST_BC_D(dim, values_ptrd(dim+2, cell_idx), &
			&                              values_grad_ptrd(1+off:dim+off,cell_idx),&
			&                              values_ptrd(dim+2, ncells+b_idx), &
			&                              values_grad_ptrd(1+off:dim+off, ncells+b_idx), update_grad)
		END SUBROUTINE         
    END SUBROUTINE



!=======================================================================
!============ WALL BOUNDARY CONDITION POSTPROCESSING METHODS ===========
!=======================================================================
	subroutine compute_force_wall(this, values_ptr, values_grad_ptr, USE_GHOST_CELLS, force)
		implicit none
		class(bc_wall_t), intent(inout) :: this
		real(8), intent(in), contiguous :: values_ptr(:, :)
		real(8), intent(in), contiguous :: values_grad_ptr(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		real(8), intent(inout) :: force(:)
		
		call compute_force_patch(this%face_indices, this%cell_indices, this%b_indices,&
								 this%mesh%dim, this%mesh%ncells,&
								 this%mesh%face_area, this%mesh%face_normal,&
								 values_ptr, values_grad_ptr, force, USE_GHOST_CELLS)
	end subroutine

	subroutine compute_fluxes_wall(this, Pr, cp, values_ptr, values_grad_ptr, USE_GHOST_CELLS, filename)
		implicit none
		class(bc_wall_t), intent(inout) :: this
		real(8), intent(in) :: Pr, cp
		real(8), intent(in), contiguous :: values_ptr(:, :)
		real(8), intent(in), contiguous :: values_grad_ptr(:, :)
		logical, intent(in) :: USE_GHOST_CELLS
		character(len=*), intent(in) :: filename
		
		call compute_wfluxes_patch(filename, this%mesh%dim, this%mesh%ncells, Pr, cp,&
								   this%face_indices, this%cell_indices, this%b_indices,&
								   this%mesh%face_center, this%mesh%face_normal,&
								   values_ptr, values_grad_ptr, USE_GHOST_CELLS)
	end subroutine

end module
