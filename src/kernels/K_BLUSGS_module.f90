module BLUSGS_module
use physical_properties_module
use stability_operator_module
implicit none
contains
!=======================================================================
!=== SUBROUTINE FOR COMPUTING SPECTARL RADIUS AND JACOBIANS ON FACES ===
!=======================================================================
pure subroutine blusgs_compute_face_values(dim, nfaces, nbfaces,&
										   face_left_cell, face_right_cell,&
										   face_area, face_normal, face_center,&
										   cell_center, k, R_gas, cv, cp, Pr,&
										   Q, lambda_face, Jacobian_face,&
										   MODEL, USE_GHOST_CELLS, USE_CONSERVATIVE_VARS)
	implicit none
	integer, intent(in) :: MODEL
    integer, intent(in) :: dim, nfaces, nbfaces
    real(8), intent(in) :: k, R_gas, cv, cp, Pr
    logical, intent(in) :: USE_GHOST_CELLS, USE_CONSERVATIVE_VARS
    real(8), intent(in), contiguous :: Q(:, :)
    real(8), intent(in), contiguous :: face_area(:), face_normal(:, :),&
									   face_center(:, :), cell_center(:, :)
    integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:)
    real(8), intent(inout), contiguous :: lambda_face(:), Jacobian_face(:, :, :)
    
    integer :: n_loop_faces, face_idx, d
    integer :: left_cell, right_cell
    real(8) :: P_L, V_L(3), T_L, V2_L, A_L, H_L, Ro_L
    real(8) :: P_R, V_R(3), T_R, V2_R, A_R, H_R, E_R, Ro_R, Vn_R
    real(8) :: Ro_tilda, V_tilda(3), V2_tilda, Vn_tilda, H_tilda, A_tilda, T_tilda, P_tilda, E_tilda
    real(8) :: sqrtRoL, sqrtRoR, invSumSqrt
    real(8) :: n(dim), ksi(dim), area, mu, lambda, RL
    real(8) :: Jac(dim+2, dim+2)
 
    V_tilda = 0.d0; V_R = 0.d0; V_L = 0.d0
    n_loop_faces = nfaces - nbfaces
    if (USE_GHOST_CELLS) n_loop_faces = nfaces
    
    !=== INTERNAL FACES: ===
    do face_idx = 1, n_loop_faces
		left_cell = face_left_cell(face_idx)
		right_cell = face_right_cell(face_idx)
		
		n = face_normal(:, face_idx)
        ksi = cell_center(:, right_cell) - cell_center(:, left_cell) 
		area = face_area(face_idx)
		
		!=== LEFT/RIGHT STATE: ===
		P_L = Q(1, left_cell); P_R = Q(1, right_cell)
		V2_L = 0.d0; V2_R = 0.d0
		do d = 1, dim
			V_L(d) = Q(d+1, left_cell); V_R(d) = Q(d+1, right_cell)
			V2_L = V2_L + V_L(d)**2; V2_R = V2_R + V_R(d)**2
		end do
		T_L = Q(dim+2, left_cell); T_R = Q(dim+2, right_cell)
		
		Ro_L = P_L/(R_gas*T_L); Ro_R = P_R/(R_gas*T_R)
		H_L = cp*T_L + 0.5d0*V2_L; H_R = cp*T_R + 0.5d0*V2_R
				
		!=== ROE AVERAGING: ===
		Ro_tilda = dsqrt(Ro_L*Ro_R)
		sqrtRoL = dsqrt(Ro_L); sqrtRoR = dsqrt(Ro_R)
		invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)
		
		V2_tilda = 0.d0
		Vn_tilda = 0.d0
		do d = 1, dim
			V_tilda(d) = (sqrtRoL*V_L(d) + sqrtRoR*V_R(d))*invSumSqrt
			V2_tilda = V2_tilda + V_tilda(d)**2
			Vn_tilda = Vn_tilda + V_tilda(d)*n(d)
		end do
		
		H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
		A_tilda = dsqrt((k-1.d0)*(H_tilda - 0.5d0*V2_tilda))
		T_tilda = A_tilda**2/(k*R_gas)
		E_tilda = cv*T_tilda + 0.5d0*V2_tilda
		P_tilda = (k - 1.d0)/k*(H_tilda - 0.5d0*V2_tilda)*Ro_tilda
		
		!=== PHYSICAL PROPERTIES: ===
		mu = MU_air(T_tilda)
		lambda = lambda_air(T_tilda, Pr, cp)
		
		!=== SPECTRAL RADIUS: ===
		lambda_face(face_idx) = dabs(Vn_tilda) + A_tilda
		if (MODEL == 2) then
			lambda_face(face_idx) = lambda_face(face_idx) + 2.d0*mu/(Ro_tilda*dabs(dot_product(n, ksi)))
        end if
        lambda_face(face_idx) = lambda_face(face_idx)*area
        
        !=== JACOBIAN MATRIX: ===
        if (USE_CONSERVATIVE_VARS) then
			!=== CONSERVATIVE VARIABLES: ===
			call get_Au_matrix(dim, k, V_tilda, n, Ro_tilda, E_tilda, P_tilda, Jac)
			Jacobian_face(:, :, face_idx) = Jac
			
			if (MODEL == 2) then
				RL = norm2(ksi)
				call get_Au_visc_matrix(dim, cp, mu, lambda, RL, V_tilda, Ro_tilda, T_tilda, A_tilda, H_tilda, Jac)
				Jacobian_face(:, :, face_idx) = Jacobian_face(:, :, face_idx) + Jac
			end if
			
			Jacobian_face(:, :, face_idx) = Jacobian_face(:, :, face_idx)*area
        else
			!=== PRIMITIVE VARIABLES: ===
			call get_Aq_matrix(dim, R_gas, k, V_tilda, n, Ro_tilda, T_tilda, Jac)
			Jacobian_face(:, :, face_idx) = Jac
			
			if (MODEL == 2) then
				RL = norm2(ksi)
				call get_Aq_visc_matrix(dim, cp, mu, lambda, RL, V_tilda, Ro_tilda, T_tilda, A_tilda, Jac)
				Jacobian_face(:, :, face_idx) = Jacobian_face(:, :, face_idx) + Jac
			end if
			
			Jacobian_face(:, :, face_idx) = Jacobian_face(:, :, face_idx)*area
        end if
    end do
    
        
    !=== EXTERNAL FACES: ===
    if (.not. USE_GHOST_CELLS) then
		do face_idx = nfaces-nbfaces+1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			n = face_normal(:, face_idx)
			ksi = face_center(:, face_idx) - cell_center(:, left_cell) 
			area = face_area(face_idx)
		
			!=== BOUNDARY FACE STATE: ===
			P_R = Q(1, right_cell)
			V2_R = 0.d0
			Vn_R = 0.d0
			do d = 1, dim
				V_R(d) = Q(d+1, right_cell)
				V2_R = V2_R + V_R(d)**2
				Vn_R = Vn_R + V_R(d)*n(d)
			end do
			T_R = Q(dim+2, right_cell)
			A_R = dsqrt(k*R_gas*T_R)
			Ro_R = P_R/(R_gas*T_R)
			H_R = cp*T_R + 0.5d0*V2_R
			E_R = cv*T_R + 0.5d0*V2_R
			
			!=== PHYSICAL PROPERTIES: ===
			mu = MU_air(T_R)
			lambda = lambda_air(T_R, Pr, cp)
			
			!=== SPECTRAL RADIUS: ===
			lambda_face(face_idx) = dabs(Vn_R) + A_R
			if (MODEL == 2) then
				lambda_face(face_idx) = lambda_face(face_idx) + 2.d0*mu/(Ro_R*dabs(dot_product(n, ksi)))
			end if
			lambda_face(face_idx) = lambda_face(face_idx)*area
			
			!=== JACOBIAN MATRIX: ===
			if (USE_CONSERVATIVE_VARS) then
				!=== CONSERVATIVE VARIABLES: ===
				call get_Au_matrix(dim, k, V_R, n, Ro_R, E_R, P_R, Jac)
				Jacobian_face(:, :, face_idx) = Jac
				
				if (MODEL == 2) then
					RL = norm2(ksi)
					call get_Au_visc_matrix(dim, cp, mu, lambda, RL, V_R, Ro_R, T_R, A_R, H_R, Jac)
					Jacobian_face(:, :, face_idx) = Jacobian_face(:, :, face_idx) + Jac
				end if
				
				Jacobian_face(:, :, face_idx) = Jacobian_face(:, :, face_idx)*area
			else
				!=== PRIMITIVE VARIABLES: ===
				call get_Aq_matrix(dim, R_gas, k, V_R, n, Ro_R, T_R, Jac)
				Jacobian_face(:, :, face_idx) = Jac
				
				if (MODEL == 2) then
					RL = norm2(ksi)
					call get_Aq_visc_matrix(dim, cp, mu, lambda, RL, V_R, Ro_R, T_R, A_R, Jac)
					Jacobian_face(:, :, face_idx) = Jacobian_face(:, :, face_idx) + Jac
				end if
				
				Jacobian_face(:, :, face_idx) = Jacobian_face(:, :, face_idx)*area
			end if
		end do
    end if
end subroutine

!=======================================================================
!========= SUBROUTINE FOR COMPUTING DIAGONAL BLOCKS OF MATRIX ==========
!=======================================================================
pure subroutine blusgs_compute_diagonal_blocks(dim, ncells, weight, cell_faces_ptr, cell_faces,&
											   face_left_cell, face_right_cell,&
											   cell_volume, pseudo_time,&
											   lambda_face, Jacobian_face, Jacobian_diag)
	implicit none
	real(8), intent(in) :: weight
	integer, intent(in) :: dim, ncells
	integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
									   face_left_cell(:), face_right_cell(:)
	real(8), intent(in), contiguous :: cell_volume(:), pseudo_time(:) 
	real(8), intent(in), contiguous :: lambda_face(:), Jacobian_face(:, :, :)
	real(8), intent(inout), contiguous :: Jacobian_diag(:, :, :)
	
	integer :: cell_idx, face_idx
	integer :: pos1, pos2, f_ptr, j
	integer :: left_cell, right_cell
	
    do cell_idx = 1, ncells
		Jacobian_diag(:, :, cell_idx) = 0.d0
		
		pos1 = cell_faces_ptr(cell_idx)
		pos2 = cell_faces_ptr(cell_idx+1)-1
        
        do f_ptr = pos1, pos2
            face_idx = cell_faces(f_ptr)
            left_cell = face_left_cell(face_idx)
            right_cell = face_right_cell(face_idx)
                        
            if (left_cell == cell_idx) then
				Jacobian_diag(:, :, cell_idx) = Jacobian_diag(:, :, cell_idx) + weight*Jacobian_face(:, :, face_idx)
            else
				Jacobian_diag(:, :, cell_idx) = Jacobian_diag(:, :, cell_idx) - weight*Jacobian_face(:, :, face_idx)
            end if
         
            
            do j = 1, dim+2
				Jacobian_diag(j, j, cell_idx) = Jacobian_diag(j, j, cell_idx) + weight*lambda_face(face_idx)
            end do              
        end do
        
        do j = 1, dim+2
			Jacobian_diag(j, j, cell_idx) = Jacobian_diag(j, j, cell_idx) + cell_volume(cell_idx)/pseudo_time(cell_idx)
        end do
        
    end do
    
end subroutine

!=======================================================================
!=================== BLUSGS FORWARD SWEEP SUBROUTINE ===================
!=======================================================================
pure subroutine blusgs_forward_sweep(dim, ncells, weight, cell_faces_ptr, cell_faces,&
								face_left_cell, face_right_cell,&
								lambda_face, Jacobian_face, Jacobian_diag,&
								Fluxes, fluxes_sign, dQ)
	implicit none
	real(8), intent(in) :: weight, fluxes_sign
	integer, intent(in) :: dim, ncells
	integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
									   face_left_cell(:), face_right_cell(:)
	real(8), intent(in), contiguous :: lambda_face(:), Jacobian_face(:, :, :), Jacobian_diag(:, :, :)
	real(8), intent(in), contiguous :: Fluxes(:, :)
	real(8), intent(inout), contiguous :: dQ(:)
	
	integer :: cell_idx, neighbor_idx, face_idx
    real(8) :: sum_L(dim+2), RHS(dim+2), A_ij_minus(dim+2, dim+2), LU(dim+2, dim+2), lambda_ij
    integer :: pos1, pos2, f_ptr, j, r, c, off
    integer :: left_cell, right_cell
    
    do cell_idx = 1, ncells
		sum_L = 0.d0
		pos1 = cell_faces_ptr(cell_idx)
		pos2 = cell_faces_ptr(cell_idx+1)-1
		
		do f_ptr = pos1, pos2
			face_idx = cell_faces(f_ptr)
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if (left_cell == cell_idx) then
                neighbor_idx = right_cell
                A_ij_minus = Jacobian_face(:, :, face_idx)
            else
                neighbor_idx = left_cell
                A_ij_minus = -Jacobian_face(:, :, face_idx)
            end if
            
            !=== LOWER DIAGONAL BLOCKS: ===
            if (neighbor_idx <= ncells .and. neighbor_idx < cell_idx) then
				lambda_ij = lambda_face(face_idx)  
				
				do j = 1, dim+2
					A_ij_minus(j, j) = A_ij_minus(j, j) - lambda_ij
				end do									
				
				off = (dim+2)*(neighbor_idx - 1)
				do c = 1, dim+2
					do r = 1, dim+2
						sum_L(r) = sum_L(r) + A_ij_minus(r, c)*dQ(c + off)
					end do
				end do
										
            end if
			
		end do
		
		off = (dim+2)*(cell_idx - 1)
		sum_L = weight*sum_L
		RHS = -fluxes_sign*Fluxes(:, cell_idx) - sum_L
        call solve_system(dim+2, Jacobian_diag(:, :, cell_idx), RHS, dQ(1+off:dim+2+off), LU)
    end do
end subroutine

!=======================================================================
!=================== BLUSGS BACKWARD SWEEP SUBROUTINE ==================
!=======================================================================
pure subroutine blusgs_backward_sweep(dim, ncells, weight, cell_faces_ptr, cell_faces,&
								 face_left_cell, face_right_cell,&
								 lambda_face, Jacobian_face, Jacobian_diag, dQ)
	implicit none
	real(8), intent(in) :: weight
	integer, intent(in) :: dim, ncells
	integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
									   face_left_cell(:), face_right_cell(:)
	real(8), intent(in), contiguous :: lambda_face(:), Jacobian_face(:, :, :), Jacobian_diag(:, :, :)
	real(8), intent(inout), contiguous :: dQ(:)
	
	integer :: cell_idx, neighbor_idx, face_idx
    real(8) :: sum_U(dim+2), temp(dim+2), A_ij_minus(dim+2, dim+2), LU(dim+2, dim+2), lambda_ij
    integer :: pos1, pos2, f_ptr, j, r, c, off
    integer :: left_cell, right_cell
    
    
    do cell_idx = ncells, 1, -1
		sum_U = 0.d0
		pos1 = cell_faces_ptr(cell_idx)
		pos2 = cell_faces_ptr(cell_idx+1)-1
		
		do f_ptr = pos1, pos2
			face_idx = cell_faces(f_ptr)
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if (left_cell == cell_idx) then
                neighbor_idx = right_cell
                A_ij_minus = Jacobian_face(:, :, face_idx)
            else
                neighbor_idx = left_cell
                A_ij_minus = -Jacobian_face(:, :, face_idx)
            end if
            
            !=== UPPER DIAGONAL BLOCKS: ===
            if (neighbor_idx <= ncells .and. neighbor_idx > cell_idx) then
				lambda_ij = lambda_face(face_idx)  
				
				do j = 1, dim+2
					A_ij_minus(j, j) = A_ij_minus(j, j) - lambda_ij
				end do									
				
				off = (dim+2)*(neighbor_idx-1)
				do c = 1, dim+2
					do r = 1, dim+2
						sum_U(r) = sum_U(r) + A_ij_minus(r, c)*dQ(c + off)
					end do
				end do
										
            end if
			
		end do
		
		off = (dim+2)*(cell_idx-1)
		sum_U = weight*sum_U
        call solve_system(dim+2, Jacobian_diag(:, :, cell_idx), sum_U, temp, LU)
        dQ(1+off:dim+2+off) = dQ(1+off:dim+2+off) - temp
    end do
end subroutine
		
pure subroutine solve_system(n, A, b, x, LU)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n, n)
    real(8), intent(in) :: b(n)
    real(8), intent(out) :: x(n)
    real(8), intent(inout) :: LU(n, n)
    
    integer :: i, j, k
    
    LU = A
    
    do i = 1, n
        do j = i, n
            do k = 1, i-1
                LU(i,j) = LU(i,j) - LU(i,k)*LU(k,j)
            end do
        end do
        do j = i+1, n
            do k = 1, i-1
                LU(j,i) = LU(j,i) - LU(j,k)*LU(k,i)
            end do
            LU(j,i) = LU(j,i)/LU(i,i)
        end do
    end do
    
    x = b
    do i = 1, n
        do j = 1, i-1
            x(i) = x(i) - LU(i,j)*x(j)
        end do
    end do
    
    do i = n, 1, -1
        do j = i+1, n
            x(i) = x(i) - LU(i,j)*x(j)
        end do
        x(i) = x(i)/LU(i,i)
    end do
end subroutine

!=======================================================================
!================= BLUSGS MAIN SETUP&APPLY SUBROUTINE ==================
!=======================================================================
pure subroutine blusgs_setup(dim, ncells, nfaces, nbfaces,&
							 face_left_cell, face_right_cell,&
							 cell_faces, cell_faces_ptr,&
							 face_area, face_normal, face_center,&
							 cell_center, cell_volume, pseudo_time,&
							 k, R_gas, cv, cp, Pr, weight,&
							 Q, lambda_face, Jacobian_face, Jacobian_diag,&
							 MODEL, USE_GHOST_CELLS, USE_CONSERVATIVE_VARS)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces, nbfaces
	integer, intent(in) :: MODEL
	logical, intent(in) :: USE_GHOST_CELLS, USE_CONSERVATIVE_VARS
	real(8), intent(in) :: k, R_gas, cp, cv, Pr, weight
	integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
									   cell_faces(:), cell_faces_ptr(:)
	real(8), intent(in), contiguous :: face_area(:), face_normal(:, :),&
									   face_center(:, :), cell_center(:, :),&
									   cell_volume(:), pseudo_time(:), Q(:, :)
	real(8), intent(inout), contiguous :: lambda_face(:), Jacobian_face(:, :, :), Jacobian_diag(:, :, :)
	
	call blusgs_compute_face_values(dim, nfaces, nbfaces,&
									face_left_cell, face_right_cell,&
									face_area, face_normal, face_center,&
									cell_center, k, R_gas, cv, cp, Pr,&
									Q, lambda_face, Jacobian_face,&
									MODEL, USE_GHOST_CELLS, USE_CONSERVATIVE_VARS)
	call blusgs_compute_diagonal_blocks(dim, ncells, weight, cell_faces_ptr, cell_faces,&
										face_left_cell, face_right_cell,&
										cell_volume, pseudo_time,&
										lambda_face, Jacobian_face, Jacobian_diag)									  						 
end subroutine
							  	
pure subroutine blusgs_apply(dim, ncells, weight, fluxes_sign,&
							 face_left_cell, face_right_cell,&
							 cell_faces, cell_faces_ptr,&
							 dQ, Fluxes, lambda_face, Jacobian_face, Jacobian_diag)
	implicit none
	integer, intent(in) :: dim, ncells
	real(8), intent(in) :: weight, fluxes_sign
	integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
									   cell_faces(:), cell_faces_ptr(:)
	real(8), intent(in), contiguous :: Fluxes(:, :), lambda_face(:),&
									   Jacobian_face(:, :, :), Jacobian_diag(:, :, :)
	real(8), intent(inout), contiguous :: dQ(:)
								 
	call blusgs_forward_sweep(dim, ncells, weight, cell_faces_ptr, cell_faces,&
							  face_left_cell, face_right_cell,&
							  lambda_face, Jacobian_face, Jacobian_diag,&
							  Fluxes, fluxes_sign, dQ)
	
	call blusgs_backward_sweep(dim, ncells, weight, cell_faces_ptr, cell_faces,&
							   face_left_cell, face_right_cell,&
							   lambda_face, Jacobian_face, Jacobian_diag, dQ)
end subroutine 



end module
