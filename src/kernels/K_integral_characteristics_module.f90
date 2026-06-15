module integral_characteristics_module
use physical_properties_module
implicit none

contains
!=======================================================================
!========= SUBOUTINE FOR COMPUTING TOTAL FORCE VECTOR ON PATCH =========
!=======================================================================
	subroutine compute_force_patch(face_indices, cell_indices, b_indices,&
								   dim, ncells, face_area, face_normal,&
								   Q, gradQ, force, USE_GHOST_CELLS)
		implicit none
		integer, intent(in) :: dim, ncells
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in), contiguous :: face_indices(:), cell_indices(:), b_indices(:)
		real(8), intent(in), contiguous :: face_area(:), face_normal(:, :)
		real(8), intent(in), contiguous :: Q(:, :), gradQ(:, :)
		real(8), intent(inout) :: force(:)
		
		integer :: f, d, v, off, face_idx, cell_idx, b_idx
		real(8) :: df(dim), tmp(dim)
		real(8) :: P_b, T_b, grad_V_b(dim*dim), mu
		real(8) :: grad_V(dim, dim), V_grad(dim, dim), divV
		
		force = 0.d0
		if (.not. USE_GHOST_CELLS) then
			do f = 1, size(face_indices)
				face_idx = face_indices(f)
				b_idx = b_indices(f)
				
				!=== BOUNDARY VALUES: ===
				P_b = Q(1, ncells+b_idx)
				off = dim
				grad_V_b = gradQ(off+1:off+dim*dim, ncells+b_idx)
				T_b = Q(dim+2, ncells+b_idx)
				mu = MU_air(T_b)
				
				!=== VELOCITY GRADIENT TENSOR: ===
				divV = 0.d0
				do d = 1, dim
					do v = 1, dim
						grad_V(v, d) = grad_V_b(d + dim*(v - 1))
						V_grad(v, d) = grad_V_b(v + dim*(d - 1))
					end do
					divV = divV + grad_V(d, d)
				end do
				
				!=== FORCE CALCULATION: ===				
				tmp = 0.d0
				do v = 1, dim 
					do d = 1, dim
						tmp(d) = tmp(d) + grad_V(d, v)*face_normal(v, face_idx) + V_grad(d, v)*face_normal(v, face_idx)
					end do
				end do
				
				df = (-P_b*face_normal(:, face_idx) + mu*(tmp - (2.d0/3.d0)*divV*face_normal(:, face_idx)))*face_area(face_idx)
				force = force + df
			end do
			
		else
			do f = 1, size(face_indices)
				face_idx = face_indices(f)
				cell_idx = cell_indices(f)
				b_idx = b_indices(f)
				
				!=== BOUNDARY VALUES: ===
				P_b = 0.5d0*(Q(1, cell_idx) + Q(1, ncells+b_idx))
				off = dim
				grad_V_b = 0.5d0*(gradQ(off+1:off+dim*dim, cell_idx) +&
								  gradQ(off+1:off+dim*dim, ncells+b_idx))
				T_b = 0.5d0*(Q(dim+2, cell_idx) + Q(dim+2, ncells+b_idx))
				mu = MU_air(T_b)
				
				!=== VELOCITY GRADIENT TENSOR: ===
				divV = 0.d0
				do d = 1, dim
					do v = 1, dim
						grad_V(v, d) = grad_V_b(d + dim*(v - 1))
						V_grad(v, d) = grad_V_b(v + dim*(d - 1))
					end do
					divV = divV + grad_V(d, d)
				end do
				
				!=== FORCE CALCULATION: ===				
				tmp = 0.d0
				do v = 1, dim 
					do d = 1, dim
						tmp(d) = tmp(d) + grad_V(d, v)*face_normal(v, face_idx) + V_grad(d, v)*face_normal(v, face_idx)
					end do
				end do
				
				df = (-P_b*face_normal(:, face_idx) + mu*(tmp - (2.d0/3.d0)*divV*face_normal(:, face_idx)))*face_area(face_idx)
				force = force + df
			end do
			
			
		end if

	end subroutine

!=======================================================================
!======= SUBROUTINE FOR COMPUTING WALL SHEAR STERSS & HEAT FLUX ========
!=======================================================================
	subroutine compute_wfluxes_patch(filename, dim, ncells, Pr, cp,&
									 face_indices, cell_indices, b_indices,&
									 face_center, face_normal,&
									 Q, gradQ, USE_GHOST_CELLS)
		implicit none
		character(len=*), intent(in) :: filename
		integer, intent(in) :: dim, ncells
		real(8), intent(in) :: Pr, cp
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in), contiguous :: face_indices(:), cell_indices(:), b_indices(:)
		real(8), intent(in), contiguous :: face_center(:, :), face_normal(:, :)
		real(8), intent(in), contiguous :: Q(:, :), gradQ(:, :)
		
		integer :: iunit
		integer :: f, d, v, face_idx, b_idx, cell_idx
		real(8) :: mu, lambda
		real(8) :: P_b, T_b, grad_V_b(dim*dim), grad_T_b(dim)
		real(8) :: grad_V(dim, dim), V_grad(dim, dim), divV
		real(8) :: total_stress(dim), shear_stress(dim), heat_flux(dim)
		
		character(len=16) :: dim_str
		character(len=256) :: fmt_data, fmt_header
		character(len=256) :: data, header
		
		open(newunit=iunit, file=filename, status='replace', action='write')
		!=== FILE HEADER: ===
		write(dim_str, '(I0)') dim 
		fmt_header = '(' // trim(dim_str) // 'A18, ' // trim(dim_str) // 'A18, ' // trim(dim_str) // 'A18, A18' // ')'
		fmt_data = '(' // trim(dim_str) // 'ES18.9, ' // trim(dim_str) // 'ES18.9, ' // trim(dim_str) // 'ES18.9, ES18.9' // ')'
		if (dim == 3) then
			write(header, fmt_header) 'X', 'Y', 'Z', 'Tau_X', 'Tau_Y', 'Tau_Z', 'q_X', 'q_Y', 'q_Z', 'q_n'
		else
			write(header, fmt_header) 'X', 'Y', 'Tau_X', 'Tau_Y', 'q_X', 'q_Y', 'q_n'
		end if
		
		write(iunit, '(a)') trim(header)
		
		!=== FILE DATA: ==
		if (USE_GHOST_CELLS) then
			do f = 1, size(face_indices)
				face_idx = face_indices(f)
				cell_idx = cell_indices(f)
				b_idx = b_indices(f)
				
				!=== BOUNDARY VALUES: ===
				P_b = 0.5d0*(Q(1, cell_idx) + Q(1, ncells+b_idx))
				T_b = 0.5d0*(Q(dim+2, cell_idx) + Q(dim+2, ncells+b_idx))
				grad_V_b = 0.5d0*(gradQ(1+dim:dim+dim*dim, cell_idx) +&
								  gradQ(1+dim:dim+dim*dim, ncells+b_idx))
				grad_T_b = 0.5d0*(gradQ(1+dim*(dim+1):dim+dim*(dim+1), cell_idx) +&
								  gradQ(1+dim*(dim+1):dim+dim*(dim+1), ncells+b_idx))
				mu = MU_air(T_b)
				lambda = lambda_air(T_b, Pr, cp)
				
				!=== VELOCITY GRADIENT TENSOR: ===
				divV = 0.d0
				do d = 1, dim
					do v = 1, dim
						grad_V(v, d) = grad_V_b(d + dim*(v - 1))
						V_grad(v, d) = grad_V_b(v + dim*(d - 1))
					end do
					divV = divV + grad_V(d, d)
				end do
				
				!=== TOTAL STRESS CALCULATION: ===				
				total_stress = 0.d0
				do v = 1, dim 
					do d = 1, dim
						total_stress(d) = total_stress(d) + grad_V(d, v)*face_normal(v, face_idx)&
														  + V_grad(d, v)*face_normal(v, face_idx)
					end do
				end do
				
				total_stress = -P_b*face_normal(:, face_idx) + mu*(total_stress - (2.d0/3.d0)*divV*face_normal(:, face_idx))
				shear_stress = total_stress - dot_product(total_stress, face_normal(:, face_idx))*face_normal(:, face_idx)
				heat_flux = -lambda*grad_T_b
				
				write(data, fmt_data) face_center(:, face_idx), shear_stress(:),&
				                      heat_flux(:), dot_product(heat_flux, face_normal(:, face_idx))
				write(iunit, '(a)') trim(data)
			end do
		else
			do f = 1, size(face_indices)
				face_idx = face_indices(f)
				b_idx = b_indices(f)
				
				!=== BOUNDARY VALUES: ===
				P_b = Q(1, ncells+b_idx)
				T_b = Q(dim+2, ncells+b_idx)
				grad_V_b = gradQ(dim+1:dim+dim*dim, ncells+b_idx)
				grad_T_b = gradQ(1+dim*(dim+1):dim+dim*(dim+1), ncells+b_idx)
				mu = MU_air(T_b)
				lambda = lambda_air(T_b, Pr, cp)
				
				!=== VELOCITY GRADIENT TENSOR: ===
				divV = 0.d0
				do d = 1, dim
					do v = 1, dim
						grad_V(v, d) = grad_V_b(d + dim*(v - 1))
						V_grad(v, d) = grad_V_b(v + dim*(d - 1))
					end do
					divV = divV + grad_V(d, d)
				end do
				
				!=== TOTAL STRESS CALCULATION: ===				
				total_stress = 0.d0
				do v = 1, dim 
					do d = 1, dim
						total_stress(d) = total_stress(d) + grad_V(d, v)*face_normal(v, face_idx)&
														  + V_grad(d, v)*face_normal(v, face_idx)
					end do
				end do
				
				total_stress = -P_b*face_normal(:, face_idx) + mu*(total_stress - (2.d0/3.d0)*divV*face_normal(:, face_idx))
				shear_stress = total_stress - dot_product(total_stress, face_normal(:, face_idx))*face_normal(:, face_idx)
				heat_flux = -lambda*grad_T_b
				
				
				write(data, fmt_data) face_center(:, face_idx), shear_stress(:),&
								      heat_flux(:), dot_product(heat_flux, face_normal(:, face_idx))
				write(iunit, '(a)') trim(data)
			end do	
		end if
			
		close(iunit)
	end subroutine

end module
