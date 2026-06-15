module bc_kernels_module
implicit none
contains
!=======================================================================
!============= FIXED VALUE BOUNDARY CONDITION SUBROUTINES ==============
!=======================================================================
	pure subroutine apply_fixed_value_bc(dim, n, val_c, grad_c, val_b, b_val, b_grad, dist, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim), grad_c(dim)
		real(8), intent(in) :: val_c, val_b, dist
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: b_val, b_grad(dim)
		
		integer :: d
		real(8) :: tmp
		
		!=== FI_boundary = FI_fixed: ===
		b_val = val_b
		
		!=== GRAD(FI)_t_boundary = GRAD(FI)_t_cell ===
		!=== GRAD(FI)_n_boundary = n*(FI_fixed - FI_cell)/dist: ===
		if (update_grad) then
			tmp = 0.d0
			do d = 1, dim
				tmp = tmp + grad_c(d)*n(d)
			end do
			b_grad = (grad_c - n*tmp) + n*(val_b - val_c)/dist
		end if
	end subroutine
	
!=======================================================================
!============ FIXED GRADIENT BOUNDARY CONDITION SUBROUTINES ============
!=======================================================================
	pure subroutine apply_fixed_gradient_bc(dim, n, val_c, grad_c, grad_bn, b_val, b_grad, dist, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim), grad_c(dim)
		real(8), intent(in) :: val_c, grad_bn, dist
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: b_val, b_grad(dim)
		
		integer :: d
		real(8) :: tmp
		
		!=== FI_boundary = FI_cell + GRAD(FI)_b_n*(r_cb, n): ===
		b_val = val_c + dist*grad_bn
		
		!=== GRAD(FI)_t_boundary = GRAD(FI)_t_cell ===
		!=== GRAD(FI)_n_boundary = n*GRAD(FI)_n_fixed: ===
		if (update_grad) then
			tmp = 0.d0
			do d = 1, dim
				tmp = tmp + grad_c(d)*n(d)
			end do
			b_grad = (grad_c - n*tmp) + n*grad_bn
		end if
	end subroutine

!=======================================================================
!============ EXTRAPOLATION BOUNDARY CONDITION SUBROUTINES =============
!=======================================================================	
	pure subroutine apply_extrapolation0_bc(dim, val_c, grad_c, b_val, b_grad, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: grad_c(dim)
		real(8), intent(in) :: val_c
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: b_val, b_grad(dim)
		
		
		!=== FI_boundary = FI_cell: ===
		b_val = val_c
		
		!=== GRAD(FI)_boundary = GRAD(FI)_cell: ===
		if (update_grad) then
			b_grad = grad_c
		end if
	end subroutine
	
	pure subroutine apply_extrapolation1_bc(dim, n, r_cf, val_c, grad_c, b_val, b_grad, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim), r_cf(dim), grad_c(dim)
		real(8), intent(in) :: val_c
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: b_val, b_grad(dim)
		
		!=== FI_boundary = FI_cell + (GRAD(FI)_cell, r_cf): ===
		b_val = val_c + dot_product(grad_c, r_cf)
		
		!=== GRAD(FI)_boundary = GRAD(FI)_cell: ===
		if (update_grad) then
			b_grad = grad_c
		end if
	end subroutine

!=======================================================================
!============= SLIP VECTOR BOUNDARY CONDITION SUBROUTINES ==============
!=======================================================================	
	pure subroutine apply_slip_vector_bc(dim, n, vec_c, grad_c, b_vec, b_grad, dist, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim), vec_c(dim), grad_c(dim*dim)
		real(8), intent(in) :: dist
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: b_vec(dim), b_grad(dim*dim)
		
		real(8) :: vec_n, g_n
		integer :: d, off
		
		!=== VEC_n_boundary = 0, VEC_t_boundary = VEC_t_cell: ===
		vec_n = 0.d0
		do d = 1, dim
			vec_n = vec_n + vec_c(d)*n(d)
		end do
		b_vec = vec_c - vec_n*n
		
		!=== GRAD(VEC)_t_boundary = GRAD(VEC)_t_cell ===
		!=== GRAD(VEC)_n_boundary = 0.0*n*(VEC_boundary - VEC_cell)/dist: ===
		if (update_grad) then
			do d = 1, dim
				off = (d-1)*dim
				g_n = dot_product(grad_c(off+1:off+dim), n)
				b_grad(off+1: off+dim) = (grad_c(off+1:off+dim) - n*g_n) !+ n*(b_vec(d) - vec_c(d))/dist
			end do
		end if
	end subroutine
	
!=======================================================================
!================ AUX GRADIENTS MANIPULATION SUBROUTINES ===============
!=======================================================================
	pure subroutine copy_grad_from_cell(dim, grad_c, b_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: grad_c(dim)
		real(8), intent(inout) :: b_grad(dim)
		
		b_grad = grad_c
	end subroutine
        
	pure subroutine copy_grad_t_from_cell(dim, n, grad_c, b_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim), grad_c(dim)
		real(8), intent(inout) :: b_grad(dim)
	
		integer :: d
		real(8) :: tmp
		
		!== tangential component: ==
		tmp = 0.d0
		do d = 1, dim
			tmp = tmp + grad_c(d)*n(d)
		end do
		
		do d = 1, dim
			b_grad(d) = grad_c(d) - tmp*n(d)
		end do 
	end subroutine
        
	pure subroutine add_grad_n(dim, n, dist, val_c, val_b, b_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: n(dim)
		real(8), intent(in) :: dist, val_c, val_b
		real(8), intent(inout) :: b_grad(dim)
	
		integer :: d
		!== normal component: ==
		do d = 1, dim
			b_grad(d) = b_grad(d) + n(d)*(val_b - val_c)/dist
		end do
	end subroutine




!=======================================================================
!============= ADDITIONAL BOUNDARY CONDITION SUBROUTINES ===============
!=======================================================================
	pure subroutine apply_reiman_pressure_bc(dim, prs_c, grad_prs_c, vel_c, tmp_c, n, k, R_gas, b_prs, b_prs_grad, update_grad)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: prs_c, grad_prs_c(dim), vel_c(dim), tmp_c, n(dim)
		real(8), intent(in) :: k, R_gas
		logical, intent(in) :: update_grad
		real(8), intent(inout) :: b_prs, b_prs_grad(dim)
		
		real(8) :: v_n
		real(8) :: a_c, ro_c
		
		v_n = dot_product(-n, vel_c)
					
		if (V_n > 0.d0) then
			a_c = dsqrt(k*R_gas*tmp_c)
			b_prs = prs_c*(1.d0 - (k - 1.d0)/2.d0 * v_n/a_c)**(2.d0*k/(k - 1.d0))
		else
			ro_c = prs_c/(R_gas*tmp_c)
			b_prs = prs_c + (k + 1.d0)/4.d0*ro_c*v_n**2 + dsqrt((1.d0/16.d0)*(k + 1.d0)**2*ro_c**2*v_n**4 + prs_c*ro_c*v_n**2*k)
		end if
		
		!=== GRAD(FI)_boundary = GRAD(FI)_cell: ===
		if (update_grad) then
			b_prs_grad = grad_prs_c
		end if
	
	end subroutine
	
end module
