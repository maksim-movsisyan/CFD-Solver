module fluxes_calculation_module
use reiman_solver_module
use reconstruction_module
use physical_properties_module
use stability_operator_module, only: get_J_mtrx
implicit none
contains
!=======================================================================
!=========== CMPRS INVISCID FLUXES CALCLUCATION SUBROUTINE: ============
!=======================================================================
	pure subroutine compute_inviscid_fluxes_cmprs(dim, nfaces, nbfaces,&
												  face_left_cell, face_right_cell,&
												  face_normal, face_area,&
												  face_center, cell_center,&
												  Q, gradQ, Fluxes,&
												  k, R_gas, cp, cv,&
												  ORDER, LIMITER,&
												  REIMAN_SOLVER_TYPE, USE_GHOST_CELLS)
		implicit none
		integer, intent(in) :: dim, nfaces, nbfaces
		integer, intent(in) :: ORDER, LIMITER, REIMAN_SOLVER_TYPE
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:),&
										   face_center(:, :), cell_center(:, :)
		real(8), intent(in), contiguous :: Q(:, :), gradQ(:, :)
		real(8), intent(in) :: k, R_gas, cp, cv
		real(8), intent(inout), contiguous :: Fluxes(:, :)
		
		integer :: n_loop_faces		
		integer :: face_idx
		integer :: left_cell, right_cell
		real(8) :: n(dim), area 
		real(8) :: Var_L(dim + 2), Var_R(dim + 2)						!=== [P U V (W) T] ===
		real(8) :: face_fluxes(dim + 2)
		
		real(8) :: P, V(dim), T, Ro, H, V_n, V2
		integer :: d
		
		n_loop_faces = nfaces - nbfaces
		if (USE_GHOST_CELLS) n_loop_faces = nfaces
						
		!=== INTERNAL FACES: ===
		select case(REIMAN_SOLVER_TYPE)
			case(1)
#				define REIMAN_SOLVER_NAME HLL_reiman_solver
#				include "inviscid_flux_kernel.inc"
#				undef REIMAN_SOLVER_NAME
			case(2)
#				define REIMAN_SOLVER_NAME HLLC_reiman_solver
#				include "inviscid_flux_kernel.inc"
#				undef REIMAN_SOLVER_NAME
			case(3)
#				define REIMAN_SOLVER_NAME ROE_reiman_solver
#				include "inviscid_flux_kernel.inc"
#				undef REIMAN_SOLVER_NAME
			case(4)
#				define REIMAN_SOLVER_NAME AUSM_reiman_solver
#				include "inviscid_flux_kernel.inc"
#				undef REIMAN_SOLVER_NAME
			case(5)
#				define REIMAN_SOLVER_NAME AUSMplus_reiman_solver
#				include "inviscid_flux_kernel.inc"
#				undef REIMAN_SOLVER_NAME
		end select
		
									
		!=== EXTERNAL FACES: ===
		if (.not. USE_GHOST_CELLS) then
			!=== UPDATING FLUXES USING BOUNDARY FACES VALUES: ===
			do face_idx = nfaces - nbfaces + 1, nfaces
				left_cell = face_left_cell(face_idx)
				right_cell = face_right_cell(face_idx)			

				n = face_normal(:, face_idx)
				area = face_area(face_idx)
							
				V2 = 0.d0
				P = Q(1, right_cell)
				do d = 1, dim
					V(d) = Q(1+d, right_cell)
					V2 = V2 + V(d)**2
				end do
				T = Q(dim+2, right_cell)
				
				V_n = DOT_PRODUCT(V, n)
				Ro = P/(R_gas*T); H = cp*T + 0.5d0*V2
				
				face_fluxes(1) = Ro*V_n
				do d = 1, dim
					face_fluxes(1+d) = Ro*V(d)*V_n + P*n(d)
				end do
				face_fluxes(dim+2) = Ro*H*V_n
							
				!=== FLUXES UPDATING: ===
				Fluxes(:, left_cell) = Fluxes(:, left_cell) + face_fluxes*area
			end do
		end if
							
	end subroutine

!=======================================================================
!=========== CMPRS VISCOUS FLUXES CALCLUCATION SUBROUTINE: =============
!=======================================================================
	pure subroutine compute_viscous_fluxes_cmprs(dim, nfaces, nbfaces,&
											 	 face_left_cell, face_right_cell,&
												 face_normal, face_area, face_weight,&
												 face_center, cell_center,&
												 Q, gradQ, Fluxes,&
												 k, R_gas, cp, cv, Pr, USE_GHOST_CELLS) 
		implicit none
		integer, intent(in) :: dim, nfaces, nbfaces
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:), face_weight(:),&
										   face_center(:, :), cell_center(:, :)
		real(8), intent(in), contiguous :: Q(:, :), gradQ(:, :)
		real(8), intent(in) :: k, R_gas, cp, cv, Pr
		real(8), intent(inout), contiguous :: Fluxes(:, :)
		
		integer :: n_loop_faces
		integer :: face_idx, d, v
		integer :: left_cell, right_cell
		integer :: offset, v_offset
		real(8) :: w, n(dim), ksi(dim), area, dist
		real(8) :: Q_f(dim+2), gradQ_f(dim*(dim+2))
		real(8) :: divV, lambda, mu
		real(8) :: qn_f, taun_f(dim)
		real(8) :: face_fluxes(dim+2)
		
		n_loop_faces = nfaces - nbfaces
		if (USE_GHOST_CELLS) n_loop_faces = nfaces
		
		!=== INTERNAL FACES: ===
		do face_idx = 1, n_loop_faces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			w = face_weight(face_idx)
			
			n = face_normal(:, face_idx)
			area = face_area(face_idx)
			
			ksi = cell_center(:, right_cell) - cell_center(:, left_cell)
			dist = norm2(ksi)
			ksi = ksi/dist
						
			!=== LINEAR INTERPOLATION ON FACE: ===
			Q_f = w*Q(:, right_cell) + (1.d0 - w)*Q(:, left_cell)
			gradQ_f = w*gradQ(:, right_cell) + (1.d0 - w)*gradQ(:, left_cell)
						
			!=== PHYSICAL TRANSPORT PROPERTIES: ===
			mu = MU_air(Q_f(dim+2))
			lambda = lambda_air(Q_f(dim+2), Pr, cp)
				
			!=== ENERGY DIFFUSION: ===
			offset = dim*(dim + 1)	
			qn_f = -lambda*((Q(dim+2, right_cell) - Q(dim+2, left_cell))/dist +&
							dot_product(gradQ_f(1+offset:dim+offset), n - ksi))
			
			!=== MOMENTUM DIFFUSION: ===
			divV = 0.d0
			do d = 1, dim
				offset = d*dim
				taun_f(d) = (Q(d+1, right_cell) - Q(d+1, left_cell))/dist +&
								dot_product(gradQ_f(1+offset:dim+offset), n - ksi)							
				do v = 1, dim
					v_offset = v*dim
					taun_f(d) = taun_f(d) + gradQ_f(d + v_offset)*n(v)
				end do
				divV = divV + gradQ_f(d + offset)
			end do
			
			do d = 1, dim
				taun_f(d) = mu*(taun_f(d) - 2.0d0/3.0d0*divV*n(d))
			end do
			
			!=== DIFFUSION FLUXES: ===
			face_fluxes = 0.d0
			do d = 1, dim
				face_fluxes(d+1) = -taun_f(d)
			end do
			face_fluxes(dim+2) = -dot_product(taun_f, Q_f(2:dim+1)) + qn_f
			
			!=== FLUXES UPDATING: ===
			Fluxes(:, left_cell) = Fluxes(:, left_cell) + face_fluxes*area
			Fluxes(:, right_cell) = Fluxes(:, right_cell) - face_fluxes*area
		end do
		
		!=== EXTERNAL FACES: ===
		if (.not. USE_GHOST_CELLS) then
			!=== UPDATING FLUXES USING BOUNDARY FACES VALUES: ===
			do face_idx = nfaces - nbfaces + 1, nfaces
				left_cell = face_left_cell(face_idx)
				right_cell = face_right_cell(face_idx)
							
				n = face_normal(:, face_idx)
				area = face_area(face_idx)
							
				!=== FACE VALUES: ===
				Q_f = Q(:, right_cell)
				gradQ_f = gradQ(:, right_cell)
						
				!=== PHYSICAL TRANSPORT PROPERTIES: ===
				mu = MU_air(Q_f(dim+2))
				lambda = lambda_air(Q_f(dim+2), Pr, cp)
				
				!=== ENERGY DIFFUSION: ===
				offset = dim*(dim + 1)	
				qn_f = -lambda*dot_product(gradQ_f(1+offset:dim+offset), n)
				
				!=== MOMENTUM DIFFUSION: ===
				taun_f = 0.d0
				divV = 0.d0
				do d = 1, dim
					offset = d*dim
					do v = 1, dim
						v_offset = v*dim
						taun_f(d) = taun_f(d) + (gradQ_f(v + offset) + gradQ_f(d + v_offset))*n(v)
					end do
					divV = divV + gradQ_f(d + offset)
				end do
				
				do d = 1, dim
					taun_f(d) = mu*(taun_f(d) - 2.d0/3.d0*divV*n(d))
				end do
				
				!=== DIFFUSION FLUXES: ===
				face_fluxes = 0.d0
				do d = 1, dim
					face_fluxes(d+1) = -taun_f(d)
				end do
				face_fluxes(dim+2) = -dot_product(taun_f, Q_f(2:dim+1)) + qn_f
				
				!=== FLUXES UPDATING: ===
				Fluxes(:, left_cell) = Fluxes(:, left_cell) + face_fluxes*area
			end do
		end if					
	end subroutine	
		

!=======================================================================
!===== CONSERVATIVE VARIABLES TO PRIMITIVE VARIABLES SUBROUTINE: =======
!=======================================================================
	pure subroutine cons_fluxes_to_prim_fluxes(dim, ncells, Q, Fluxes, k, R_gas, cv, cp)
		implicit none
		integer, intent(in) :: dim, ncells
		real(8), intent(in) :: k, R_gas, cv, cp
		real(8), intent(in), contiguous :: Q(:, :)
		real(8), intent(inout), contiguous :: Fluxes(:, :)
		
		real(8) :: P, V(3), V2, T, Ro, H, A
		integer :: i, d, j, m
		real(8) :: tmp_res(dim + 2), Jac(dim + 2, dim + 2)
		

		V = 0.d0						
		do i = 1, ncells
			P = Q(1, i)
			do d = 1, dim
				V(d) = Q(1+d, i)
			end do
			T = Q(dim+2, i)
			Ro = P/(R_gas*T)
			H = cp*T + 0.5d0*(V(1)**2 + V(2)**2 + (V(3)**2*(dim-2)))
			A = dsqrt(k*R_gas*T)
			
			call get_J_mtrx(dim, cp, k, A, V, H, Ro, T, Jac)
			do j = 1, dim + 2
				tmp_res(j) = 0.0d0
				do m = 1, dim + 2
					tmp_res(j) = tmp_res(j) + Jac(j, m)*Fluxes(m, i)
				end do
			end do
			Fluxes(:, i) = tmp_res
		end do
				
	end subroutine
end module
