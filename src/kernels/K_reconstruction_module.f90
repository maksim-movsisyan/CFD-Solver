module reconstruction_module
implicit none

contains
!=======================================================================
!================= COUPLED VARIABLES RECONSTRUCTION ====================
!=======================================================================
    pure subroutine reconstruct_iface_states_coupled(dim, n_comp, face_idx, ORDER, LIMITER,&
													left_cell, right_cell,&
													face_center, cell_center,&
													Q, gradQ, W_L, W_R)
        integer, intent(in) :: dim, n_comp, face_idx, LIMITER, ORDER
        integer, intent(in) :: left_cell, right_cell
        real(8), intent(in), contiguous :: face_center(:, :), cell_center(:, :)
        real(8), intent(in), contiguous :: Q(:,:)
        real(8), intent(in), contiguous :: gradQ(:,:)
        real(8), intent(out) :: W_L(n_comp), W_R(n_comp)

        integer :: c, d
        real(8) :: r_cf_L(3), r_cf_R(3), ksi(3)
        real(8) :: val_L, val_R, grad_L(3), grad_R(3)
        real(8) :: LL, RR

        r_cf_L(1:dim) = face_center(1:dim, face_idx) - cell_center(1:dim, left_cell)
        r_cf_R(1:dim) = face_center(1:dim, face_idx) - cell_center(1:dim, right_cell)
        
        ksi(1:dim) = cell_center(1:dim, right_cell) - cell_center(1:dim, left_cell)

        do c = 1, n_comp
            val_L = Q(c, left_cell)
            val_R = Q(c, right_cell)

			do d = 1, dim
				grad_L(d) = gradQ(d + (c-1)*dim, left_cell)
				grad_R(d) = gradQ(d + (c-1)*dim, right_cell)
			end do
				
            select case(ORDER)
            case(1)
                W_L(c) = val_L
                W_R(c) = val_R
                
            case(2)             
                LL = val_R - 2.0d0*dot_product(grad_L(1:dim), ksi(1:dim))
                RR = val_L + 2.0d0*dot_product(grad_R(1:dim), ksi(1:dim))

                call muscl2(LL, val_L, val_R, W_L(c), LIMITER)
                call muscl2(RR, val_R, val_L, W_R(c), LIMITER)
            
            case(3)
				LL = val_R - 2.0d0*dot_product(grad_L(1:dim), ksi(1:dim))
                RR = val_L + 2.0d0*dot_product(grad_R(1:dim), ksi(1:dim))

                call weno3(LL, val_L, val_R, W_L(c))
                call weno3(RR, val_R, val_L, W_R(c))
                
            end select
        end do
    end subroutine



!=======================================================================
!============================== TVD LIMMITERS ==========================
!=======================================================================
	pure real(8) function psi_VA(r, LIMITER)
		implicit none
		integer, intent(in) :: LIMITER
		real(8), intent(in) :: r
		real(8) :: temp
		if (r < 0.d0) then
			temp = 0.d0
		else if (r >= 0.d0) then
			if (LIMITER == 1) then
				temp = min(r, 1.d0) !minmod
			else if (LIMITER == 2) then
				temp = (r**2 + r)/(r**2 + 1.d0) !van Albada
			else if (LIMITER == 3) then
				temp = 2.d0*r/(r + 1.d0) !van Leer
			else if (LIMITER == 4) then
				temp = max(min(2.d0*r, 1.d0), min(r, 2.d0)) !Superbee
			else if (LIMITER == 5) then
				temp = 0.5d0*(r + 1.d0)*min(min(1.d0, 4.d0*r/(r + 1.d0)), min(1.d0, 4.d0/(r + 1.d0)))
			end if
		end if
		psi_VA = temp
	end function

!=======================================================================
!=================== MUSCL 2ND ORDER RECONSTRUCTION ====================
!=======================================================================
	pure subroutine muscl2(a_mm, a_m, a, reconstructed, LIMITER)
		implicit none
		real(8), parameter :: eps = 1e-10
		integer, intent(in) :: LIMITER
		real(8), intent(in) :: a_mm, a_m, a
		real(8), intent(out) :: reconstructed
		real(8) :: r_l, delta_minus, delta_plus
		
		delta_minus = a_m - a_mm
		delta_plus  = a - a_m
		
		if (dabs(delta_minus) > eps .AND. dabs(delta_plus) > eps) then
			r_l = delta_plus/delta_minus  
			reconstructed = a_m + 0.5d0*psi_VA(r_l, LIMITER)*delta_minus
			return
		else
			reconstructed = a_m
			return
		end if
	end subroutine


!=======================================================================
!=================== WENO 3RD ORDER RECONSTRUCTION =====================
!=======================================================================
	pure subroutine weno3(a_m, a, a_p, reconstructed)
		implicit none
		real(8), intent(in) :: a_m, a, a_p
		real(8), intent(out) :: reconstructed
		real(8), parameter :: d0 = 3.d0/4.d0, d1 = 1.d0/4.d0, eps = 1e-6
		real(8) :: beta0, beta1, alfa0, alfa1, alfa_sum, omega0, omega1
		real(8) :: a0, a1
		
		a0 = 0.5d0*(a + a_p)
		a1 = -0.5d0*a_m + 1.5d0*a
		
		beta0 = (a - a_p)**2
		beta1 = (a_m - a)**2
		
		alfa0 = d0/(eps + beta0)**2
		alfa1 = d1/(eps + beta1)**2
		
		alfa_sum = alfa0 + alfa1
		
		omega0 = alfa0/alfa_sum
		omega1 = alfa1/alfa_sum
		
		reconstructed = omega0*a0 + omega1*a1
	end subroutine


!=======================================================================
!=================== WENO 5TH ORDER RECONSTRUCTION =====================
!=======================================================================
	pure subroutine weno5(a_mm, a_m, a, a_p, a_pp, reconstructed)
		implicit none
		real(8), intent(in) :: a_mm, a_m, a, a_p, a_pp
		real(8), intent(out) :: reconstructed
		real(8), parameter :: d0 = 5.d0/16.d0, d1 = 10.d0/16.d0, d2 = 1.d0/16.d0, eps = 1e-6
		real(8) :: beta0, beta1, beta2, alfa0, alfa1, alfa2, alfa_sum, omega0, omega1, omega2
		real(8) :: a0, a1, a2
		
		a0 = (3.d0/8.d0)*a + (3.d0/4.d0)*a_p - (1.d0/8.d0)*a_pp
		a1 = -(1.d0/8.d0)*a_m + (3.d0/4.d0)*a + (3.d0/8.d0)*a_p
		a2 = (3.d0/8.d0)*a_mm - (5.d0/4.d0)*a_m + (15.d0/8.d0)*a
		
		beta0 = (13.d0/12.d0)*(a - 2.d0*a_p + a_pp)**2 + (1.d0/4.d0)*(3.d0*a - 4.d0*a_p + a_pp)**2
		beta1 = (13.d0/12.d0)*(a_m - 2.d0*a + a_p)**2 + (1.d0/4.d0)*(a_m - a_p)**2
		beta2 = (13.d0/12.d0)*(a_mm - 2.d0*a_m + a)**2 + (1.d0/4.d0)*(a_mm - 4.d0*a_m + 3.d0*a)**2
		
		alfa0 = d0/(eps + beta0)**2
		alfa1 = d1/(eps + beta1)**2
		alfa2 = d2/(eps + beta2)**2
		
		alfa_sum = alfa0 + alfa1 + alfa2
		
		omega0 = alfa0/alfa_sum
		omega1 = alfa1/alfa_sum
		omega2 = alfa2/alfa_sum
		
		reconstructed = omega0*a0 + omega1*a1 + omega2*a2

	end subroutine


!=======================================================================
!=================== WENO 7TH ORDER RECONSTRUCTION =====================
!=======================================================================
	pure subroutine weno7(a_mmm, a_mm, a_m, a, a_p, a_pp, a_ppp, reconstructed)
		implicit none
		real(8), intent(in) :: a_mmm, a_mm, a_m, a, a_p, a_pp, a_ppp
		real(8), intent(out) :: reconstructed
		real(8), parameter :: d0 = 7.d0/64.d0, d1 = 35.d0/64.d0, d2 = 21.d0/64.d0, d3 = 1.d0/64.d0, eps = 1e-6
		real(8) :: beta0, beta1, beta2, beta3, alfa0, alfa1, alfa2, alfa3, alfa_sum, omega0, omega1, omega2, omega3
		real(8) :: a0, a1, a2, a3
		
		a0 = (5.d0/16.d0)*a + (15.d0/16.d0)*a_p - (5.d0/16.d0)*a_pp + (1.d0/16.d0)*a_ppp
		a1 = -(1.d0/16.d0)*a_m + (9.d0/16.d0)*a + (9.d0/16.d0)*a_p - (1.d0/16.d0)*a_pp
		a2 = (1.d0/16.d0)*a_mm - (5.d0/16.d0)*a_m + (15.d0/16.d0)*a + (5.d0/16.d0)*a_p
		a3 = -(5.d0/16.d0)*a_mmm + (21.d0/16.d0)*a_mm - (35.d0/16.d0)*a_m + (35.d0/16.d0)*a
		
		beta0 = (1.d0/64.d0)*(-15.d0*a + 25.d0*a_p - 13.d0*a_pp + 3.d0*a_ppp)**2 +&
				(13.d0/12.d0)*(2.d0*a - 5.d0*a_p + 4.d0*a_pp - a_ppp)**2 +&
				(781.d0/720.d0)*(-a + 3.d0*a_p - 3.d0*a_pp + a_ppp)**2
				
		beta1 = (1.d0/64.d0)*(-3.d0*a_m - 3.d0*a + 7.d0*a_p - a_pp)**2 +&
				(13.d0/12.d0)*(a_m - 2.d0*a + a_p)**2. +&
				(781.d0/720.d0)*(-a_m + 3.d0*a - 3.d0*a_p + a_pp)**2
				
		beta2 = (1.d0/64.d0)*(a_mm - 7.d0*a_m + 3.d0*a + 3.d0*a_p)**2 +&
				(13.d0/12.d0)*(a_m - 2.d0*a + a_p)**2 +&
				(781.d0/720.d0)*(-a_mm + 3.d0*a_m - 3.d0*a + a_p)**2
				
		beta3 = (1.d0/64.d0)*(-3.d0*a_mmm + 13.d0*a_mm - 25.d0*a_m + 15.d0*a)**2 +&
				(13.d0/12.d0)*(-a_mmm + 4.d0*a_mm - 5.d0*a_m + 2.d0*a)**2 +&
				(781.d0/720.d0)*(-a_mmm + 3.d0*a_mm - 3.d0*a_m + a)**2
		
		
		
		alfa0 = d0/(eps + beta0)**2
		alfa1 = d1/(eps + beta1)**2
		alfa2 = d2/(eps + beta2)**2
		alfa3 = d3/(eps + beta3)**2
		
		alfa_sum = alfa0 + alfa1 + alfa2 + alfa3
		
		omega0 = alfa0/alfa_sum
		omega1 = alfa1/alfa_sum
		omega2 = alfa2/alfa_sum
		omega3 = alfa3/alfa_sum
		
		reconstructed = omega0*a0 + omega1*a1 + omega2*a2 + omega3*a3

	end subroutine


!=======================================================================
!=================== WENO 9TH ORDER RECONSTRUCTION =====================
!=======================================================================
	pure subroutine weno9(a_mmmm, a_mmm, a_mm, a_m, a, a_p, a_pp, a_ppp, a_pppp, reconstructed)
		implicit none
		real(8), intent(in) :: a_mmmm, a_mmm, a_mm, a_m, a, a_p, a_pp, a_ppp, a_pppp
		real(8), intent(out) :: reconstructed
		real(8), parameter :: d0 = 9.d0/256.d0, d1 = 21.d0/64.d0, d2 = 63.d0/128.d0, &
							  d3 = 9.d0/64.d0, d4 = 1.d0/256.d0, eps = 1e-6
		real(8) :: beta0, beta1, beta2, beta3, beta4
		real(8) :: alfa0, alfa1, alfa2, alfa3, alfa4, alfa_sum
		real(8) :: omega0, omega1, omega2, omega3, omega4
		real(8) :: a0, a1, a2, a3, a4

		a0 = (35.d0/128.d0)*a + (35.d0/32.d0)*a_p - (35.d0/64.d0)*a_pp + &
		(7.d0/32.d0)*a_ppp - (5.d0/128.d0)*a_pppp

		a1 = -(5.d0/128.d0)*a_m + (15.d0/32.d0)*a + (45.d0/64.d0)*a_p - &
		(5.d0/32.d0)*a_pp + (3.d0/128.d0)*a_ppp

		a2 = (3.d0/128.d0)*a_mm - (5.d0/32.d0)*a_m + (45.d0/64.d0)*a + &
		(15.d0/32.d0)*a_p - (5.d0/128.d0)*a_pp

		a3 = -(5.d0/128.d0)*a_mmm + (7.d0/32.d0)*a_mm - (35.d0/64.d0)*a_m + &
		(35.d0/32.d0)*a + (35.d0/128.d0)*a_p

		a4 = (35.d0/128.d0)*a_mmmm - (45.d0/32.d0)*a_mmm + (189.d0/64.d0)*a_mm - &
		(105.d0/32.d0)*a_m + (315.d0/128.d0)*a

		beta0 = ((-35.d0*a + 70.d0*a_p - 56.d0*a_pp + 26.d0*a_ppp - 5.d0*a_pppp)**2)/256.d0 + &
				((4613.d0*a - 13772.d0*a_p + 15198.d0*a_pp - 7532.d0*a_ppp + 1493.d0*a_pppp)**2)/2246400.d0 + &
				(781.d0/2880.d0)*(-5.d0*a + 18.d0*a_p - 24.d0*a_pp + 14.d0*a_ppp - 3.d0*a_pppp)**2 + &
				(1421461.d0/1310400.d0)*(a - 4.d0*a_p + 6.d0*a_pp - 4.d0*a_ppp + a_pppp)**2

		beta1 = ((-5.d0*a_m - 10.d0*a + 20.d0*a_p - 6.d0*a_pp + a_ppp)**2)/256.d0 + &
				((1493.d0*a_m - 2852.d0*a + 1158.d0*a_p + 268.d0*a_pp - 67.d0*a_ppp)**2)/2246400.d0 + &
				(781.d0/2880.d0)*(-3.d0*a_m + 10.d0*a - 12.d0*a_p + 6.d0*a_pp - a_ppp)**2 + &
				(1421461.d0/1310400.d0)*(a_m - 4.d0*a + 6.d0*a_p - 4.d0*a_pp + a_ppp)**2

		beta2 = ((a_mm - 10.d0*a_m + 10.d0*a_p - a_pp)**2)/256.d0 + &
				((-67.d0*a_mm + 1828.d0*a_m - 3522.d0*a + 1828.d0*a_p - 67.d0*a_pp)**2)/2246400.d0 + &
				(781.d0/2880.d0)*(-a_mm + 2.d0*a_m - 2.d0*a_p + a_pp)**2 + &
				(1421461.d0/1310400.d0)*(a_mm - 4.d0*a_m + 6.d0*a - 4.d0*a_p + a_pp)**2

		beta3 = ((-a_mmm + 6.d0*a_mm - 20.d0*a_m + 10.d0*a + 5.d0*a_p)**2)/256.d0 + &
				((-67.d0*a_mmm + 268.d0*a_mm + 1158.d0*a_m - 2852.d0*a + 1493.d0*a_p)**2)/2246400.d0 + &
				(781.d0/2880.d0)*(a_mmm - 6.d0*a_mm + 12.d0*a_m - 10.d0*a + 3.d0*a_p)**2 + &
				(1421461.d0/1310400.d0)*(a_mmm - 4.d0*a_mm + 6.d0*a_m - 4.d0*a + a_p)**2

		beta4 = ((5.d0*a_mmmm - 26.d0*a_mmm + 56.d0*a_mm - 70.d0*a_m + 35.d0*a)**2)/256.d0 + &
				((1493.d0*a_mmmm - 7532.d0*a_mmm + 15198.d0*a_mm - 13772.d0*a_m + 4613.d0*a)**2)/2246400.d0 + &
				(781.d0/2880.d0)*(3.d0*a_mmmm - 14.d0*a_mmm + 24.d0*a_mm - 18.d0*a_m + 5.d0*a)**2 + &
				(1421461.d0/1310400.d0)*(a_mmmm - 4.d0*a_mmm + 6.d0*a_mm - 4.d0*a_m + a)**2

		alfa0 = d0/(eps + beta0)**2
		alfa1 = d1/(eps + beta1)**2
		alfa2 = d2/(eps + beta2)**2
		alfa3 = d3/(eps + beta3)**2
		alfa4 = d4/(eps + beta4)**2

		alfa_sum = alfa0 + alfa1 + alfa2 + alfa3 + alfa4

		omega0 = alfa0/alfa_sum
		omega1 = alfa1/alfa_sum
		omega2 = alfa2/alfa_sum
		omega3 = alfa3/alfa_sum
		omega4 = alfa4/alfa_sum

		reconstructed = omega0*a0 + omega1*a1 + omega2*a2 + omega3*a3 + omega4*a4
	end subroutine

end module
