module bc_apply_ghost_cells_module
use bc_ghost_kernels_module
implicit none
contains
!=======================================================================
!====================== FIXED VALUE LOOP SUBROUTINE ====================
!=======================================================================
pure subroutine apply_fixed_value_ghost_loop(dim, ncells, nfaces, start_vidx, end_vidx, update_grad,&
										     fixed_values, values_ptr, values_grad_ptr,&
										     b_indices, cell_indices, face_indices,&
										     cell_center, face_center, face_normal)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces, start_vidx, end_vidx
	logical, intent(in) :: update_grad
	integer, intent(in), contiguous :: b_indices(:), cell_indices(:), face_indices(:)
	real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :), face_normal(:, :)
	real(8), intent(in) :: fixed_values(:)
	real(8), intent(inout), contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		
	integer :: i, d, off
	integer :: face_idx, cell_idx, b_idx
	real(8) :: n(dim), ksi(dim), dist

	
	do i = 1, nfaces					
		b_idx = b_indices(i)
		cell_idx = cell_indices(i)
		face_idx = face_indices(i)
		
		n = face_normal(:, face_idx)
		ksi = face_center(:, face_idx) - cell_center(:, cell_idx)
		dist = dot_product(ksi, n)
				
		do d = start_vidx, end_vidx
			off = dim*(d-1)
			call apply_fixed_value_ghost_bc(dim, n,&
											values_ptr(d, cell_idx), values_grad_ptr(1+off:dim+off, cell_idx),&
											fixed_values(d),&
											values_ptr(d, ncells+b_idx), values_grad_ptr(1+off:dim+off, ncells+b_idx),&
											dist, update_grad)
		end do		
	end do
	
end subroutine

!=======================================================================
!====================== EXTRAPOLATION LOOP SUBROUTINE ==================
!=======================================================================
pure subroutine apply_extrapolation0_ghost_loop(dim, ncells, nfaces, start_vidx, end_vidx, update_grad,&
											    values_ptr, values_grad_ptr,&
											    b_indices, cell_indices, face_indices)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces, start_vidx, end_vidx
	logical, intent(in) :: update_grad
	integer, intent(in), contiguous :: b_indices(:), cell_indices(:), face_indices(:)
	real(8), intent(inout), contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		
	integer :: i, d, off
	integer :: cell_idx, b_idx

	do i = 1, nfaces	
		cell_idx = cell_indices(i)
		b_idx = b_indices(i)
		
		do d = start_vidx, end_vidx
			off = dim*(d-1)
			call apply_extrapolation0_ghost_bc(dim,&
								 values_ptr(d, cell_idx), values_grad_ptr(1+off:dim+off, cell_idx),&
								 values_ptr(d, ncells+b_idx), values_grad_ptr(1+off:dim+off, ncells+b_idx), update_grad)
		end do
	end do
	
end subroutine

!=======================================================================
!====================== FIXED GRADIENT LOOP SUBROUTINE =================
!=======================================================================
pure subroutine apply_fixed_gradient_ghost_loop(dim, ncells, nfaces, start_vidx, end_vidx, update_grad,&
										  fixed_ngrads, values_ptr, values_grad_ptr,&
									      b_indices, cell_indices, face_indices,&
									      cell_center, face_center, face_normal)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces, start_vidx, end_vidx
	logical, intent(in) :: update_grad
	integer, intent(in), contiguous :: b_indices(:), cell_indices(:), face_indices(:)
	real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :), face_normal(:, :)
	real(8), intent(in) :: fixed_ngrads(:)
	real(8), intent(inout), contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		
	integer :: i, d, off
	integer :: face_idx, cell_idx, b_idx
	real(8) :: n(dim), ksi(dim), dist

	
	do i = 1, nfaces					
		b_idx = b_indices(i)
		cell_idx = cell_indices(i)
		face_idx = face_indices(i)
		
		n = face_normal(:, face_idx)
		ksi = face_center(:, face_idx) - cell_center(:, cell_idx)
		dist = dot_product(ksi, n)
				
		do d = start_vidx, end_vidx
			off = dim*(d-1)
			call apply_fixed_gradient_ghost_bc(dim, n,&
										 values_ptr(d, cell_idx), values_grad_ptr(1+off:dim+off, cell_idx),&
										 fixed_ngrads(d),&
										 values_ptr(d, ncells+b_idx), values_grad_ptr(1+off:dim+off, ncells+b_idx),&
										 dist, update_grad)
		end do		
	end do
	
end subroutine

!=======================================================================
!====================== SPLIP VECTOR LOOP SUBROUTINE ===================
!=======================================================================
pure subroutine apply_slip_vector_ghost_loop(dim, ncells, nfaces, start_vidx, end_vidx, update_grad,&
									   values_ptr, values_grad_ptr,&
									   b_indices, cell_indices, face_indices,&
									   cell_center, face_center, face_normal)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces, start_vidx, end_vidx
	logical, intent(in) :: update_grad
	integer, intent(in), contiguous :: b_indices(:), cell_indices(:), face_indices(:)
	real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :), face_normal(:, :)
	real(8), intent(inout), contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
	
	integer :: i
	integer :: face_idx, cell_idx, b_idx
	integer :: grad_start_vidx, grad_end_vidx
	real(8) :: n(dim), ksi(dim), dist
	
	grad_start_vidx = 1 + dim*(start_vidx-1)
	grad_end_vidx = dim + dim*(end_vidx-1)
	
	do i = 1, nfaces
		b_idx = b_indices(i)
		cell_idx = cell_indices(i)
		face_idx = face_indices(i)
		
		n = face_normal(:, face_idx)
		ksi = face_center(:, face_idx) - cell_center(:, cell_idx)
		dist = dot_product(ksi, n)
		
		call apply_slip_vector_ghost_bc(dim, n,&
								  values_ptr(start_vidx:end_vidx, cell_idx), values_grad_ptr(grad_start_vidx:grad_end_vidx, cell_idx),&
								  values_ptr(start_vidx:end_vidx, ncells+b_idx), values_grad_ptr(grad_start_vidx:grad_end_vidx, ncells+b_idx),&
								  dist, update_grad)
	end do

end subroutine





!=======================================================================
!=========== PRIMITIVE VARIABLES SPECIFIC BC LOOP SUBROUTINE ===========
!=======================================================================
pure subroutine apply_subsonic_inlet_Q_ghost_loop(dim, ncells, nfaces, update_grad,&
												 fixed_values, values_ptr, values_grad_ptr,&
												 b_indices, cell_indices, face_indices,&
												 cell_center, face_center, face_normal)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces
	logical, intent(in) :: update_grad
	integer, intent(in), contiguous :: b_indices(:), cell_indices(:), face_indices(:)
	real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :), face_normal(:, :)
	real(8), intent(in) :: fixed_values(:)								![velocity_vector, temperature]
	real(8), intent(inout), contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		
	integer :: i, d, off
	integer :: face_idx, cell_idx, b_idx
	real(8) :: n(dim), ksi(dim), dist

	do i = 1, nfaces					
		b_idx = b_indices(i)
		cell_idx = cell_indices(i)
		face_idx = face_indices(i)
		
		n = face_normal(:, face_idx)
		ksi = face_center(:, face_idx) - cell_center(:, cell_idx)
		dist = dot_product(ksi, n)
		
		!pressure:
		call apply_extrapolation0_ghost_bc(dim,&
									 values_ptr(1, cell_idx), values_grad_ptr(1:dim, cell_idx),&
									 values_ptr(1, ncells+b_idx), values_grad_ptr(1:dim, ncells+b_idx), update_grad)
										 
		!velocity&temperature:						  	
		do d = 2, dim+2
			off = dim*(d-1)
			call apply_fixed_value_ghost_bc(dim, n,&
											values_ptr(d, cell_idx), values_grad_ptr(1+off:dim+off, cell_idx),&
											fixed_values(d-1),&
											values_ptr(d, ncells+b_idx), values_grad_ptr(1+off:dim+off, ncells+b_idx),&
											dist, update_grad)
		end do
	end do
				
				
end subroutine

pure subroutine apply_subsonic_outlet_Q_ghost_loop(dim, ncells, nfaces, update_grad,&
											 fixed_values, values_ptr, values_grad_ptr,&
											 b_indices, cell_indices, face_indices,&
											 cell_center, face_center, face_normal)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces
	logical, intent(in) :: update_grad
	integer, intent(in), contiguous :: b_indices(:), cell_indices(:), face_indices(:)
	real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :), face_normal(:, :)
	real(8), intent(in) :: fixed_values(:)								![pressure]
	real(8), intent(inout), contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		
	integer :: i, d, off
	integer :: face_idx, cell_idx, b_idx
	real(8) :: n(dim), ksi(dim), dist

	do i = 1, nfaces					
		b_idx = b_indices(i)
		cell_idx = cell_indices(i)
		face_idx = face_indices(i)
		
		n = face_normal(:, face_idx)
		ksi = face_center(:, face_idx) - cell_center(:, cell_idx)
		dist = dot_product(ksi, n)
		
		!pressure:
		call apply_fixed_value_ghost_bc(dim, n,&
								  values_ptr(1, cell_idx), values_grad_ptr(1:dim, cell_idx),&
								  fixed_values(1),&
								  values_ptr(1, ncells+b_idx), values_grad_ptr(1:dim, ncells+b_idx),&
								  dist, update_grad)
									  
										 
		!velocity&temperature:						  	
		do d = 2, dim+2
			off = dim*(d-1)
			call apply_extrapolation0_ghost_bc(dim,&
										 values_ptr(d, cell_idx), values_grad_ptr(1+off:dim+off, cell_idx),&
										 values_ptr(d, ncells+b_idx), values_grad_ptr(1+off:dim+off, ncells+b_idx), update_grad)
		end do
	end do
				
				
end subroutine

pure subroutine apply_symm_Q_ghost_loop(dim, ncells, nfaces, update_grad,&
								  values_ptr, values_grad_ptr,&
								  b_indices, cell_indices, face_indices,&
						   		  cell_center, face_center, face_normal)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces
	logical, intent(in) :: update_grad
	integer, intent(in), contiguous :: b_indices(:), cell_indices(:), face_indices(:)
	real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :), face_normal(:, :)
	real(8), intent(inout), contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		
	integer :: i, off
	integer :: face_idx, cell_idx, b_idx
	real(8) :: n(dim), ksi(dim), dist

	do i = 1, nfaces					
		b_idx = b_indices(i)
		cell_idx = cell_indices(i)
		face_idx = face_indices(i)
		
		n = face_normal(:, face_idx)
		ksi = face_center(:, face_idx) - cell_center(:, cell_idx)
		dist = dot_product(ksi, n)
		
		!pressure:
		call apply_extrapolation0_ghost_bc(dim,&
									 values_ptr(1, cell_idx), values_grad_ptr(1:dim, cell_idx),&
									 values_ptr(1, ncells+b_idx), values_grad_ptr(1:dim, ncells+b_idx), update_grad)	
									  					 
		!velocity:
		off = dim*dim
		call apply_slip_vector_ghost_bc(dim, n,&
								  values_ptr(2:dim+1, cell_idx), values_grad_ptr(1+dim:dim+off, cell_idx),&
								  values_ptr(2:dim+1, ncells+b_idx), values_grad_ptr(1+dim:dim+off, ncells+b_idx),&
								  dist, update_grad)
							  	
		!temperature:
		off = dim*(dim + 1)
		call apply_fixed_gradient_ghost_bc(dim, n,&
									 values_ptr(dim+2, cell_idx), values_grad_ptr(1+off:dim+off, cell_idx), 0.d0,&
									 values_ptr(dim+2, ncells+b_idx), values_grad_ptr(1+off:dim+off, ncells+b_idx),&
									 dist, update_grad)
										 
	end do
				
				
end subroutine

pure subroutine apply_isothermal_noslip_wall_Q_ghost_loop(dim, ncells, nfaces, update_grad,&
													fixed_values, values_ptr, values_grad_ptr,&
													b_indices, cell_indices, face_indices,&
													cell_center, face_center, face_normal)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces
	logical, intent(in) :: update_grad
	integer, intent(in), contiguous :: b_indices(:), cell_indices(:), face_indices(:)
	real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :), face_normal(:, :)
	real(8), intent(in) :: fixed_values(:)								![velocity_vector, temperature]
	real(8), intent(inout), contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		
	integer :: i, d, off
	integer :: face_idx, cell_idx, b_idx
	real(8) :: n(dim), ksi(dim), dist

	do i = 1, nfaces					
		b_idx = b_indices(i)
		cell_idx = cell_indices(i)
		face_idx = face_indices(i)
		
		n = face_normal(:, face_idx)
		ksi = face_center(:, face_idx) - cell_center(:, cell_idx)
		dist = dot_product(ksi, n)
		
		!pressure:
		call apply_extrapolation0_ghost_bc(dim,&
									 values_ptr(1, cell_idx), values_grad_ptr(1:dim, cell_idx),&
									 values_ptr(1, ncells+b_idx), values_grad_ptr(1:dim, ncells+b_idx), update_grad)
										 
		!velocity&temperature:						  	
		do d = 2, dim+2
			off = dim*(d-1)
			call apply_fixed_value_ghost_bc(dim, n,&
									  values_ptr(d, cell_idx), values_grad_ptr(1+off:dim+off, cell_idx),&
									  fixed_values(d-1),&
									  values_ptr(d, ncells+b_idx), values_grad_ptr(1+off:dim+off, ncells+b_idx),&
									  dist, update_grad)
		end do
		
		!temperature correction:
		if (values_ptr(dim+2, ncells+b_idx) < 1e-10) values_ptr(dim+2, ncells+b_idx) = 0.1d0*fixed_values(dim+1)
	end do
				
				
end subroutine

pure subroutine apply_fixed_flux_noslip_wall_Q_ghost_loop(dim, ncells, nfaces, update_grad,&
													fixed_values, values_ptr, values_grad_ptr,&
													b_indices, cell_indices, face_indices,&
													cell_center, face_center, face_normal)
	implicit none
	integer, intent(in) :: dim, ncells, nfaces
	logical, intent(in) :: update_grad
	integer, intent(in), contiguous :: b_indices(:), cell_indices(:), face_indices(:)
	real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :), face_normal(:, :)
	real(8), intent(in) :: fixed_values(:)								![velocity_vector, heat_flux]
	real(8), intent(inout), contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		
	integer :: i, d, off
	integer :: face_idx, cell_idx, b_idx
	real(8) :: n(dim), ksi(dim), dist

	do i = 1, nfaces					
		b_idx = b_indices(i)
		cell_idx = cell_indices(i)
		face_idx = face_indices(i)
		
		n = face_normal(:, face_idx)
		ksi = face_center(:, face_idx) - cell_center(:, cell_idx)
		dist = dot_product(ksi, n)
		
		!pressure:
		call apply_extrapolation0_ghost_bc(dim,&
									 values_ptr(1, cell_idx), values_grad_ptr(1:dim, cell_idx),&
									 values_ptr(1, ncells+b_idx), values_grad_ptr(1:dim, ncells+b_idx), update_grad)
										 
		!velocity:						  	
		do d = 2, dim+1
			off = dim*(d-1)
			call apply_fixed_value_ghost_bc(dim, n,&
									  values_ptr(d, cell_idx), values_grad_ptr(1+off:dim+off, cell_idx),&
									  fixed_values(d-1),&
									  values_ptr(d, ncells+b_idx), values_grad_ptr(1+off:dim+off, ncells+b_idx),&
									  dist, update_grad)
		end do
		
		!temperature:
		off = dim*(dim+1)
		call apply_fixed_gradient_ghost_bc(dim, n,&
									 values_ptr(dim+2, cell_idx), values_grad_ptr(1+off:dim+off, cell_idx),&
									 fixed_values(dim+1),&
									 values_ptr(dim+2, ncells+b_idx), values_grad_ptr(1+off:dim+off, ncells+b_idx),&
									 dist, update_grad)
	end do
				
				
end subroutine

end module
