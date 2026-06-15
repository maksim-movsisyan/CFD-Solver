module bc_ghost_kernels_module
implicit none
contains
!=======================================================================
!============= FIXED VALUE BOUNDARY CONDITION SUBROUTINES ==============
!=======================================================================
	pure subroutine apply_fixed_value_ghost_bc(dim, n, val_c, grad_c, val_b, val_ghost, grad_ghost, dist, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim), grad_c(dim)
		real(8), intent(in) :: val_c, val_b, dist
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: val_ghost, grad_ghost(dim)
		
		integer :: d
		real(8) :: b_grad(dim), tmp
		
		!=== FI_ghost = 2*FI_fixed - FI_cell: ===
		val_ghost = 2.d0*val_b - val_c
		
		!=== Grad(FI)_boundary = t*Grad(FI)_cell_t + n*(FI_fixed - FI_cell)/dist ===
		!=== Grad(FI)_ghost = 2*Grad(FI)_boundary - Grad(FI)_cell ===
		if (update_grad) then
			tmp = 0.d0
			do d = 1, dim
				tmp = tmp + grad_c(d)*n(d)
			end do
			b_grad = (grad_c - n*tmp) + n*(val_b - val_c)/dist
			
			grad_ghost = 2.d0*b_grad - grad_c
		end if
	end subroutine
	
!=======================================================================
!============ FIXED GRADIENT BOUNDARY CONDITION SUBROUTINES ============
!=======================================================================
	pure subroutine apply_fixed_gradient_ghost_bc(dim, n, val_c, grad_c, grad_bn, val_ghost, grad_ghost, dist, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim), grad_c(dim)
		real(8), intent(in) :: val_c, grad_bn, dist
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: val_ghost, grad_ghost(dim)
		
		integer :: d
		real(8) :: tmp, b_val, b_grad(dim)
		
		!=== FI_boundary = FI_cell + GRAD(FI)_b_n*(r_cb, n): ===
		!=== FI_ghost = 2*FI_boundary - FI_cell: ===
		b_val = val_c + dist*grad_bn
		val_ghost = 2.d0*b_val - val_c
		
		!=== Grad(FI)_boundary = t*Grad(FI)_cell_t + n*GRAD(FI)_n_fixed ===
		!=== Grad(FI)_ghost = 2*Grad(FI)_boundary - Grad(FI)_cell ===
		if (update_grad) then
			tmp = 0.d0
			do d = 1, dim
				tmp = tmp + grad_c(d)*n(d)
			end do
			b_grad = (grad_c - n*tmp) + n*grad_bn
			
			grad_ghost = 2.d0*b_grad - grad_c
		end if
	end subroutine

!=======================================================================
!============ EXTRAPOLATION BOUNDARY CONDITION SUBROUTINES =============
!=======================================================================	
	pure subroutine apply_extrapolation0_ghost_bc(dim, val_c, grad_c, val_ghost, grad_ghost, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: grad_c(dim)
		real(8), intent(in) :: val_c
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: val_ghost, grad_ghost(dim)
				
		
		!=== FI_ghost = FI_cell: ===
		val_ghost = val_c
		
		!=== GRAD(FI)_ghost = GRAD(FI)_cell: ===
		if (update_grad) then
			grad_ghost = grad_c
		end if
	end subroutine
	
	pure subroutine apply_extrapolation1_ghost_bc(dim, n, r_cf, val_c, grad_c, val_ghost, grad_ghost, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim), r_cf(dim), grad_c(dim)
		real(8), intent(in) :: val_c
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: val_ghost, grad_ghost(dim)
		
		!=== FI_boundary = FI_cell + (GRAD(FI)_cell, r_cf): ===
		val_ghost = val_c + 2.d0*dot_product(grad_c, r_cf)
		
		!=== GRAD(FI)_boundary = GRAD(FI)_cell: ===
		if (update_grad) then
			grad_ghost = grad_c
		end if
	end subroutine

!=======================================================================
!============= SLIP VECTOR BOUNDARY CONDITION SUBROUTINES ==============
!=======================================================================	
	pure subroutine apply_slip_vector_ghost_bc(dim, n, vec_c, grad_c, vec_ghost, grad_ghost, dist, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim), vec_c(dim), grad_c(dim*dim)
		real(8), intent(in) :: dist
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: vec_ghost(dim), grad_ghost(dim*dim)
		
		real(8) :: vec_n, g_n, b_grad(dim*dim)
		integer :: d, off
		
		!=== VEC_n_boundary = 0, VEC_t_boundary = VEC_t_cell: ===
		vec_n = 0.d0
		do d = 1, dim
			vec_n = vec_n + vec_c(d)*n(d)
		end do
		vec_ghost = vec_c - 2.d0*vec_n*n
		
		!=== GRAD(VEC)_boundary = t*GRAD(VEC)_t_cell + n*0.0 ===
		
		if (update_grad) then
			do d = 1, dim
				off = (d-1)*dim
				g_n = dot_product(grad_c(off+1:off+dim), n)
				b_grad(off+1: off+dim) = (grad_c(off+1:off+dim) - n*g_n) !+ n*(b_vec(d) - vec_c(d))/dist
			end do
			
			grad_ghost = 2.d0*b_grad - grad_c
		end if
	end subroutine
	
end module
