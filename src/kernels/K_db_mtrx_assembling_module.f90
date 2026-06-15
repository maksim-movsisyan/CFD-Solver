module db_mtrx_assembling_module
use array_sorting_module
use stability_operator_module
use physical_properties_module
implicit none
contains
!=======================================================================
!=== SUBROUTINE FOR ASSEMBLING BCSR MATRIX FOR DENSITY-BASED SOLVER  ===
!=== IN PRIMITIVE VARIABLES											 ===
!=======================================================================	
	pure subroutine db_assemble_bcsr_matrix_inv(dim, ncells, nfaces, nbfaces,&
										        face_left_cell, face_right_cell, face_bidx,&
										        face_normal, face_area, face_weight,&
										        cell_volume, pseudo_time, Q,&
										        values, diag_indices,&
										        map_LL, map_LR, map_RL, map_RR, map_LB,&
										        k, R_gas, cp, cv,&
										        STABILITY_OPERATOR_TYPE, USE_GHOST_CELLS, IS_VANISH)
		implicit none
		integer, intent(in) :: dim, ncells, nfaces, nbfaces
		integer, intent(in) :: STABILITY_OPERATOR_TYPE
		logical, intent(in) :: USE_GHOST_CELLS, IS_VANISH
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
										   face_bidx(:), diag_indices(:)
		integer, intent(in), contiguous :: map_LL(:), map_LR(:), map_RL(:), map_RR(:), map_LB(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:),&
										   face_weight(:), cell_volume(:), pseudo_time(:)
		real(8), intent(in), contiguous :: Q(:, :)
		real(8), intent(in) :: k, R_gas, cp, cv
		real(8), intent(inout), contiguous :: values(:, :, :)
									
		integer :: cell_idx, face_idx, b_idx, diag_idx, d
		integer :: L, R
		real(8) :: n(dim), area, g, val
		real(8) :: Mtrx_L(dim+2,dim+2), Mtrx_R(dim+2,dim+2)
		
		if (IS_VANISH) values = 0.d0
		
		!=== INTERNAL FACES: ===
		select case(STABILITY_OPERATOR_TYPE)
			case(1)
#				define STABILITY_OPERATOR_NAME HLL_stability_universal
#				include "inviscid_jacobian_kernel.inc"
#				undef STABILITY_OPERATOR_NAME
			case(2)
#				define STABILITY_OPERATOR_NAME HLLC_stability_universal
#				include "inviscid_jacobian_kernel.inc"
#				undef STABILITY_OPERATOR_NAME
			case(3)
#				define STABILITY_OPERATOR_NAME ROE_stability_universal
#				include "inviscid_jacobian_kernel.inc"
#				undef STABILITY_OPERATOR_NAME
			case(4)
#				define STABILITY_OPERATOR_NAME AUSM_stability_universal
#				include "inviscid_jacobian_kernel.inc"
#				undef STABILITY_OPERATOR_NAME
			case(5)
#				define STABILITY_OPERATOR_NAME AUSM_plus_stability_universal
#				include "inviscid_jacobian_kernel.inc"
#				undef STABILITY_OPERATOR_NAME
		end select
		
		
		!=== EXTERNAL FACES (IS GHOST CELLS USED): ===
		if (USE_GHOST_CELLS) then
			select case(STABILITY_OPERATOR_TYPE)
				case(1)
#					define STABILITY_OPERATOR_NAME HLL_stability_universal
#					include "inviscid_jacobian_ghosts_kernel.inc"
#					undef STABILITY_OPERATOR_NAME
				case(2)
#					define STABILITY_OPERATOR_NAME HLLC_stability_universal
#					include "inviscid_jacobian_ghosts_kernel.inc"
#					undef STABILITY_OPERATOR_NAME
				case(3)
#					define STABILITY_OPERATOR_NAME ROE_stability_universal
#					include "inviscid_jacobian_ghosts_kernel.inc"
#					undef STABILITY_OPERATOR_NAME
				case(4)
#					define STABILITY_OPERATOR_NAME AUSM_stability_universal
#					include "inviscid_jacobian_ghosts_kernel.inc"
#					undef STABILITY_OPERATOR_NAME
				case(5)
#					define STABILITY_OPERATOR_NAME AUSM_plus_stability_universal
#					include "inviscid_jacobian_ghosts_kernel.inc"
#					undef STABILITY_OPERATOR_NAME
			end select
		end if
		
		
		!=== PSEUDO-TIME ADD: ===
		do cell_idx = 1, ncells
			val = cell_volume(cell_idx)/pseudo_time(cell_idx)
			diag_idx = diag_indices(cell_idx)
			
			do d = 1, dim+2
				values(d, d, diag_idx) = values(d, d, diag_idx) + val
			end do			
		end do	
	end subroutine
	
	pure subroutine db_assemble_bcsr_matrix_visc(dim, ncells, nfaces, nbfaces,&
										         face_left_cell, face_right_cell, face_bidx,&
										         face_normal, face_area, face_weight,&
										         cell_center, Q,&
										         values, map_LL, map_LR, map_RL, map_RR, map_LB,&
										         k, R_gas, cp, cv, Pr,&
										         USE_GHOST_CELLS, IS_VANISH)
		implicit none
		integer, intent(in) :: dim, ncells, nfaces, nbfaces
		logical, intent(in) :: USE_GHOST_CELLS, IS_VANISH
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
										   face_bidx(:)
		integer, intent(in), contiguous :: map_LL(:), map_LR(:), map_RL(:), map_RR(:), map_LB(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:),&
										   face_weight(:), cell_center(:, :)
		real(8), intent(in), contiguous :: Q(:, :)
		real(8), intent(in) :: k, R_gas, cp, cv, Pr
		real(8), intent(inout), contiguous :: values(:, :, :)
									
		integer :: cell_idx, face_idx, b_idx, diag_idx, d
		integer :: L, R
		real(8) :: n(dim), ksi(dim), area, g, val
		real(8) :: mu, lambda
		real(8) :: Mtrx_L(dim+2,dim+2), Mtrx_R(dim+2,dim+2)
		
		if (IS_VANISH) values = 0.d0
		
		!=== INTERNAL FACES: ===
		do face_idx = 1, nfaces - nbfaces
			L = face_left_cell(face_idx)
			R = face_right_cell(face_idx)
			
			n = face_normal(:, face_idx)
			ksi = cell_center(:, R) - cell_center(:, L)
			area = face_area(face_idx)
			g = face_weight(face_idx)
			
			!=== PHYSICAL TRANSPORT PROPERTIES: ===
			mu = MU_air(g*Q(dim+2, R) + (1.d0 - g)*Q(dim+2, L))
			lambda = lambda_air(g*Q(dim+2, R) + (1.d0 - g)*Q(dim+2, L), Pr, cp)
			
			!=== STABILITY OPERATOR: ===																
			call VISCOUS_stability_universal(dim, k, R_gas, mu, lambda, cp,&
											 Q(1, L), Q(2, L), Q(3, L), Q(4, L), Q(dim+2, L),&
											 Q(1, R), Q(2, R), Q(3, R), Q(4, R), Q(dim+2, R), g, n, ksi, Mtrx_L, Mtrx_R)
			
			!=== JACOBIAN UPDATING: ===
			values(:, :, map_LL(face_idx)) = values(:, :, map_LL(face_idx)) + Mtrx_L*area
			values(:, :, map_LR(face_idx)) = values(:, :, map_LR(face_idx)) + Mtrx_R*area
			values(:, :, map_RL(face_idx)) = values(:, :, map_RL(face_idx)) - Mtrx_L*area
			values(:, :, map_RR(face_idx)) = values(:, :, map_RR(face_idx)) - Mtrx_R*area
			
		end do
		
		
		!=== EXTERNAL FACES (IS GHOST CELLS USED): ===
		if (USE_GHOST_CELLS) then
			do face_idx = nfaces - nbfaces + 1, nfaces
				L = face_left_cell(face_idx)
				R = face_right_cell(face_idx)
				b_idx = face_bidx(face_idx)
				
				n = face_normal(:, face_idx)
				ksi = cell_center(:, R) - cell_center(:, L)
				area = face_area(face_idx)
				g = face_weight(face_idx)
				
				!=== PHYSICAL TRANSPORT PROPERTIES: ===
				mu = MU_air(g*Q(dim+2, R) + (1.d0 - g)*Q(dim+2, L))
				lambda = lambda_air(g*Q(dim+2, R) + (1.d0 - g)*Q(dim+2, L), Pr, cp)
			
				!=== STABILITY OPERATOR: ===
				call VISCOUS_stability_universal(dim, k, R_gas, mu, lambda, cp,&
											 Q(1, L), Q(2, L), Q(3, L), Q(4, L), Q(dim+2, L),&
											 Q(1, R), Q(2, R), Q(3, R), Q(4, R), Q(dim+2, R), g, n, ksi, Mtrx_L, Mtrx_R)

				!=== JACOBIAN UPDATING: ===
				values(:, :, map_LB(b_idx)) = values(:, :, map_LB(b_idx)) + Mtrx_L*area
			end do
		end if
	
	end subroutine
	
	
	
!=======================================================================
!=== SUBROUTINE FOR ASSEMBLING BCSR MATRIX FOR DENSITY-BASED SOLVER  ===
!=== IN CONSERVARIVE VARIABLES (CV) 								 ===
!=======================================================================		
	pure subroutine db_assemble_bcsr_matrix_inv_cv(dim, ncells, nfaces, nbfaces,&
												   face_left_cell, face_right_cell, face_bidx,&
										           face_normal, face_area, face_weight,&
										           cell_volume, pseudo_time, Q,&
										           values, diag_indices,&
										           map_LL, map_LR, map_RL, map_RR, map_LB,&
										           k, R_gas, cp, cv,&
										           STABILITY_OPERATOR_TYPE, USE_GHOST_CELLS, IS_VANISH)
		implicit none
		integer, intent(in) :: dim, ncells, nfaces, nbfaces
		integer, intent(in) :: STABILITY_OPERATOR_TYPE
		logical, intent(in) :: USE_GHOST_CELLS, IS_VANISH
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
										   face_bidx(:), diag_indices(:)
		integer, intent(in), contiguous :: map_LL(:), map_LR(:), map_RL(:), map_RR(:), map_LB(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:),&
										   face_weight(:), cell_volume(:), pseudo_time(:)
		real(8), intent(in), contiguous :: Q(:, :)
		real(8), intent(in) :: k, R_gas, cp, cv
		real(8), intent(inout), contiguous :: values(:, :, :)
									
		integer :: cell_idx, face_idx, b_idx, diag_idx, d
		integer :: L, R
		real(8) :: n(dim), area, g, val
		real(8) :: Mtrx_L(dim+2,dim+2), Mtrx_R(dim+2,dim+2)
		
		if (IS_VANISH) values = 0.d0
		
		!=== INTERNAL FACES: ===
		select case(STABILITY_OPERATOR_TYPE)
			case(1)
#				define STABILITY_OPERATOR_NAME HLL_stability_conservative
#				include "inviscid_jacobian_kernel.inc"
#				undef STABILITY_OPERATOR_NAME
			case(2)
#				define STABILITY_OPERATOR_NAME HLLC_stability_conservative
#				include "inviscid_jacobian_kernel.inc"
#				undef STABILITY_OPERATOR_NAME
			case(3)
#				define STABILITY_OPERATOR_NAME ROE_stability_conservative
#				include "inviscid_jacobian_kernel.inc"
#				undef STABILITY_OPERATOR_NAME
		end select
		
		
		!=== EXTERNAL FACES (IF GHOST CELLS USED): ===
		if (USE_GHOST_CELLS) then
			select case(STABILITY_OPERATOR_TYPE)
				case(1)
#					define STABILITY_OPERATOR_NAME HLL_stability_conservative
#					include "inviscid_jacobian_ghosts_kernel.inc"
#					undef STABILITY_OPERATOR_NAME
				case(2)
#					define STABILITY_OPERATOR_NAME HLLC_stability_conservative
#					include "inviscid_jacobian_ghosts_kernel.inc"
#					undef STABILITY_OPERATOR_NAME
				case(3)
#					define STABILITY_OPERATOR_NAME ROE_stability_conservative
#					include "inviscid_jacobian_ghosts_kernel.inc"
#					undef STABILITY_OPERATOR_NAME
			end select
		end if
		
		
		!=== PSEUDO-TIME ADD: ===
		do cell_idx = 1, ncells
			val = cell_volume(cell_idx)/pseudo_time(cell_idx)
			diag_idx = diag_indices(cell_idx)
			
			do d = 1, dim+2
				values(d, d, diag_idx) = values(d, d, diag_idx) + val
			end do			
		end do	
	end subroutine

	pure subroutine db_assemble_bcsr_matrix_visc_cv(dim, ncells, nfaces, nbfaces,&
										            face_left_cell, face_right_cell, face_bidx,&
										            face_normal, face_area, face_weight,&
										            cell_center, Q,&
										            values, map_LL, map_LR, map_RL, map_RR, map_LB,&
										            k, R_gas, cp, cv, Pr,&
										            USE_GHOST_CELLS, IS_VANISH)
		implicit none
		integer, intent(in) :: dim, ncells, nfaces, nbfaces
		logical, intent(in) :: USE_GHOST_CELLS, IS_VANISH
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
										   face_bidx(:)
		integer, intent(in), contiguous :: map_LL(:), map_LR(:), map_RL(:), map_RR(:), map_LB(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:),&
										   face_weight(:), cell_center(:, :)
		real(8), intent(in), contiguous :: Q(:, :)
		real(8), intent(in) :: k, R_gas, cp, cv, Pr
		real(8), intent(inout), contiguous :: values(:, :, :)
									
		integer :: cell_idx, face_idx, b_idx, diag_idx, d
		integer :: L, R
		real(8) :: n(dim), ksi(dim), area, g, val
		real(8) :: mu, lambda
		real(8) :: Mtrx_L(dim+2,dim+2), Mtrx_R(dim+2,dim+2)
		
		if (IS_VANISH) values = 0.d0
		
		!=== INTERNAL FACES: ===
		do face_idx = 1, nfaces - nbfaces
			L = face_left_cell(face_idx)
			R = face_right_cell(face_idx)
			
			n = face_normal(:, face_idx)
			ksi = cell_center(:, R) - cell_center(:, L)
			area = face_area(face_idx)
			g = face_weight(face_idx)
			
			!=== PHYSICAL TRANSPORT PROPERTIES: ===
			mu = MU_air(g*Q(dim+2, R) + (1.d0 - g)*Q(dim+2, L))
			lambda = lambda_air(g*Q(dim+2, R) + (1.d0 - g)*Q(dim+2, L), Pr, cp)
			
			!=== STABILITY OPERATOR: ===																
			call VISCOUS_stability_conservative(dim, k, R_gas, mu, lambda, cp,&
												Q(1, L), Q(2, L), Q(3, L), Q(4, L), Q(dim+2, L),&
												Q(1, R), Q(2, R), Q(3, R), Q(4, R), Q(dim+2, R), g, n, ksi, Mtrx_L, Mtrx_R)
			
			!=== JACOBIAN UPDATING: ===
			values(:, :, map_LL(face_idx)) = values(:, :, map_LL(face_idx)) + Mtrx_L*area
			values(:, :, map_LR(face_idx)) = values(:, :, map_LR(face_idx)) + Mtrx_R*area
			values(:, :, map_RL(face_idx)) = values(:, :, map_RL(face_idx)) - Mtrx_L*area
			values(:, :, map_RR(face_idx)) = values(:, :, map_RR(face_idx)) - Mtrx_R*area
			
		end do
		
		
		!=== EXTERNAL FACES (IS GHOST CELLS USED): ===
		if (USE_GHOST_CELLS) then
			do face_idx = nfaces - nbfaces + 1, nfaces
				L = face_left_cell(face_idx)
				R = face_right_cell(face_idx)
				b_idx = face_bidx(face_idx)
				
				n = face_normal(:, face_idx)
				ksi = cell_center(:, R) - cell_center(:, L)
				area = face_area(face_idx)
				g = face_weight(face_idx)
				
				!=== PHYSICAL TRANSPORT PROPERTIES: ===
				mu = MU_air(g*Q(dim+2, R) + (1.d0 - g)*Q(dim+2, L))
				lambda = lambda_air(g*Q(dim+2, R) + (1.d0 - g)*Q(dim+2, L), Pr, cp)
			
				!=== STABILITY OPERATOR: ===
				call VISCOUS_stability_conservative(dim, k, R_gas, mu, lambda, cp,&
												    Q(1, L), Q(2, L), Q(3, L), Q(4, L), Q(dim+2, L),&
												    Q(1, R), Q(2, R), Q(3, R), Q(4, R), Q(dim+2, R), g, n, ksi, Mtrx_L, Mtrx_R)

				!=== JACOBIAN UPDATING: ===
				values(:, :, map_LB(b_idx)) = values(:, :, map_LB(b_idx)) + Mtrx_L*area
			end do
		end if
	
	end subroutine
	





!=======================================================================
!================= MATRIX FILL-MAP GENERATION SUBROUTINE ===============	
!=======================================================================
	pure subroutine db_create_mtrx_fill_map(nfaces, nbfaces,&
											face_left_cell, face_right_cell,&
											face_bidx, col_indices, row_ptr,&
											map_LL, map_LR, map_RL, map_RR, map_LB)
		implicit none
		integer, intent(in) :: nfaces, nbfaces
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
										   face_bidx(:), col_indices(:), row_ptr(:)
		integer, allocatable, intent(inout) :: map_LL(:), map_LR(:), map_RL(:), map_RR(:), map_LB(:)									
		
		integer :: face_idx, b_idx
		integer :: left_cell, right_cell
		integer :: pos1l, pos1r, pos2l, pos2r
		
		
		allocate(map_LL(nfaces - nbfaces),&
				 map_LR(nfaces - nbfaces),&
				 map_RL(nfaces - nbfaces),&
				 map_RR(nfaces - nbfaces),&
				 map_LB(nbfaces))
		
		!=== INTERNAL FACES: ===
		do face_idx = 1, nfaces - nbfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			pos1l = row_ptr(left_cell)
			pos2l = row_ptr(left_cell + 1) - 1
			map_LL(face_idx) = binary_search(col_indices, pos1l, pos2l, left_cell)
			map_LR(face_idx) = binary_search(col_indices, pos1l, pos2l, right_cell)
			
			pos1r = row_ptr(right_cell)
			pos2r = row_ptr(right_cell + 1) - 1
			map_RL(face_idx) = binary_search(col_indices, pos1r, pos2r, left_cell)
			map_RR(face_idx) = binary_search(col_indices, pos1r, pos2r, right_cell)
		end do
		
		!=== EXTERNAL FACES: ===
		do face_idx = nfaces - nbfaces + 1, nfaces
			left_cell = face_left_cell(face_idx)
			b_idx = face_bidx(face_idx)
			
			pos1l = row_ptr(left_cell)
			pos2l = row_ptr(left_cell + 1) - 1
			map_LB(b_idx) = binary_search(col_indices, pos1l, pos2l, left_cell)
		end do
	end subroutine
	
end module
