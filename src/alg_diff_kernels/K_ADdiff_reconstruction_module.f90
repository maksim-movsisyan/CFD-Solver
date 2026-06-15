module ADdiff_reconstruction_module
implicit none

contains

!  Differentiation of reconstruct_iface_states_coupled in forward (tangent) mode:
!   variations   of useful results: w_l w_r
!   with respect to varying inputs: w_l q w_r gradq
!=======================================================================
!================= COUPLED VARIABLES RECONSTRUCTION ====================
!=======================================================================
  SUBROUTINE RECONSTRUCT_IFACE_STATES_COUPLED_D(dim, n_comp, face_idx, &
&   order, limiter, left_cell, right_cell, face_center, cell_center, q, &
&   qd, gradq, gradqd, w_l, w_ld, w_r, w_rd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim, n_comp, face_idx, limiter, order
    INTEGER, INTENT(IN) :: left_cell, right_cell
    REAL*8, INTENT(IN), CONTIGUOUS :: face_center(:, :), cell_center(:, &
&   :)
    REAL*8, INTENT(IN), CONTIGUOUS :: q(:, :)
    REAL*8, INTENT(IN), CONTIGUOUS :: qd(:, :)
    REAL*8, INTENT(IN), CONTIGUOUS :: gradq(:, :)
    REAL*8, INTENT(IN), CONTIGUOUS :: gradqd(:, :)
    REAL*8, INTENT(OUT) :: w_l(n_comp), w_r(n_comp)
    REAL*8, INTENT(OUT) :: w_ld(n_comp), w_rd(n_comp)
    INTEGER :: c, d
    REAL*8 :: r_cf_l(3), r_cf_r(3), ksi(3)
    REAL*8 :: val_l, val_r, grad_l(3), grad_r(3)
    REAL*8 :: val_ld, val_rd, grad_ld(3), grad_rd(3)
    REAL*8 :: ll, rr, dot_prct
    REAL*8 :: lld, rrd, dot_prctd
    r_cf_l(1:dim) = face_center(1:dim, face_idx) - cell_center(1:dim, &
&     left_cell)
    r_cf_r(1:dim) = face_center(1:dim, face_idx) - cell_center(1:dim, &
&     right_cell)
    ksi(1:dim) = cell_center(1:dim, right_cell) - cell_center(1:dim, &
&     left_cell)
    grad_ld = 0.0_8
    grad_rd = 0.0_8
    DO c=1,n_comp
      val_ld = qd(c, left_cell)
      val_l = q(c, left_cell)
      val_rd = qd(c, right_cell)
      val_r = q(c, right_cell)
      DO d=1,dim
        grad_ld(d) = gradqd(d+(c-1)*dim, left_cell)
        grad_l(d) = gradq(d+(c-1)*dim, left_cell)
        grad_rd(d) = gradqd(d+(c-1)*dim, right_cell)
        grad_r(d) = gradq(d+(c-1)*dim, right_cell)
      END DO
      SELECT CASE (order)
      CASE (1)
        w_ld(c) = val_ld
        w_l(c) = val_l
        w_rd(c) = val_rd
        w_r(c) = val_r
      CASE (2)
        dot_prct = 0.d0
        dot_prctd = 0.0_8
        DO d=1,dim
          dot_prctd = dot_prctd + ksi(d)*grad_ld(d)
          dot_prct = dot_prct + grad_l(d)*ksi(d)
        END DO
        lld = val_rd - 2.0d0*dot_prctd
        ll = val_r - 2.0d0*dot_prct
        dot_prct = 0.d0
        dot_prctd = 0.0_8
        DO d=1,dim
          dot_prctd = dot_prctd + ksi(d)*grad_rd(d)
          dot_prct = dot_prct + grad_r(d)*ksi(d)
        END DO
        rrd = val_ld + 2.0d0*dot_prctd
        rr = val_l + 2.0d0*dot_prct
        CALL MUSCL2_D(ll, lld, val_l, val_ld, val_r, val_rd, w_l(c), &
&               w_ld(c), limiter)
        CALL MUSCL2_D(rr, rrd, val_r, val_rd, val_l, val_ld, w_r(c), &
&               w_rd(c), limiter)
      CASE (3)
        dot_prct = 0.d0
        dot_prctd = 0.0_8
        DO d=1,dim
          dot_prctd = dot_prctd + ksi(d)*grad_ld(d)
          dot_prct = dot_prct + grad_l(d)*ksi(d)
        END DO
        lld = val_rd - 2.0d0*dot_prctd
        ll = val_r - 2.0d0*dot_prct
        dot_prct = 0.d0
        dot_prctd = 0.0_8
        DO d=1,dim
          dot_prctd = dot_prctd + ksi(d)*grad_rd(d)
          dot_prct = dot_prct + grad_r(d)*ksi(d)
        END DO
        rrd = val_ld + 2.0d0*dot_prctd
        rr = val_l + 2.0d0*dot_prct
        CALL WENO3_D(ll, lld, val_l, val_ld, val_r, val_rd, w_l(c), w_ld&
&              (c))
        CALL WENO3_D(rr, rrd, val_r, val_rd, val_l, val_ld, w_r(c), w_rd&
&              (c))
      END SELECT
    END DO
  END SUBROUTINE RECONSTRUCT_IFACE_STATES_COUPLED_D

!  Differentiation of psi_va in forward (tangent) mode:
!   variations   of useful results: psi_va
!   with respect to varying inputs: r
!=======================================================================
!============================== TVD LIMMITERS ==========================
!=======================================================================
  REAL*8 FUNCTION PSI_VA_D(r, rd, limiter, psi_va)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: limiter
    REAL*8, INTENT(IN) :: r
    REAL*8, INTENT(IN) :: rd
    REAL*8 :: temp
    REAL*8 :: tempd
    INTRINSIC MIN
    INTRINSIC MAX
    REAL*8 :: x1
    REAL*8 :: x1d
    REAL*8 :: y1
    REAL*8 :: y1d
    DOUBLE PRECISION :: x2
    DOUBLE PRECISION :: x2d
    DOUBLE PRECISION :: y2
    DOUBLE PRECISION :: y2d
    DOUBLE PRECISION :: min1
    DOUBLE PRECISION :: min1d
    REAL*8 :: temp0
    REAL*8 :: psi_va
    IF (r .LT. 0.d0) THEN
      temp = 0.d0
      tempd = 0.0_8
    ELSE IF (r .GE. 0.d0) THEN
      IF (limiter .EQ. 1) THEN
        IF (r .GT. 1.d0) THEN
          temp = 1.d0
          tempd = 0.0_8
        ELSE
          tempd = rd
          temp = r
        END IF
      ELSE IF (limiter .EQ. 2) THEN
!van Albada
        temp0 = (r*r+r)/(r*r+1.d0)
        tempd = (2*r-temp0*2*r+1.0)*rd/(r**2+1.d0)
        temp = temp0
      ELSE IF (limiter .EQ. 3) THEN
!van Leer
        temp0 = r/(r+1.d0)
        tempd = 2.d0*(1.0-temp0)*rd/(r+1.d0)
        temp = 2.d0*temp0
      ELSE IF (limiter .EQ. 4) THEN
        IF (2.d0*r .GT. 1.d0) THEN
          x1 = 1.d0
          x1d = 0.0_8
        ELSE
          x1d = 2.d0*rd
          x1 = 2.d0*r
        END IF
        IF (r .GT. 2.d0) THEN
          y1 = 2.d0
          y1d = 0.0_8
        ELSE
          y1d = rd
          y1 = r
        END IF
        IF (x1 .LT. y1) THEN
          tempd = y1d
          temp = y1
        ELSE
          tempd = x1d
          temp = x1
        END IF
      ELSE IF (limiter .EQ. 5) THEN
        IF (1.d0 .GT. 4.d0*r/(r+1.d0)) THEN
          temp0 = r/(r+1.d0)
          x2d = 4.d0*(1.0-temp0)*rd/(r+1.d0)
          x2 = 4.d0*temp0
        ELSE
          x2 = 1.d0
          x2d = 0.D0
        END IF
        IF (1.d0 .GT. 4.d0/(r+1.d0)) THEN
          y2d = -(4.d0*rd/(r+1.d0)**2)
          y2 = 4.d0/(r+1.d0)
        ELSE
          y2 = 1.d0
          y2d = 0.D0
        END IF
        IF (x2 .GT. y2) THEN
          min1d = y2d
          min1 = y2
        ELSE
          min1d = x2d
          min1 = x2
        END IF
        tempd = 0.5d0*(min1*rd+(r+1.d0)*min1d)
        temp = 0.5d0*(r+1.d0)*min1
      ELSE
        tempd = 0.0_8
      END IF
    ELSE
      tempd = 0.0_8
    END IF
    psi_va_d = tempd
    psi_va = temp
  END FUNCTION PSI_VA_D

!  Differentiation of muscl2 in forward (tangent) mode:
!   variations   of useful results: reconstructed
!   with respect to varying inputs: a_mm a_m a
!=======================================================================
!=================== MUSCL 2ND ORDER RECONSTRUCTION ====================
!=======================================================================
  SUBROUTINE MUSCL2_D(a_mm, a_mmd, a_m, a_md, a, ad, reconstructed, &
&   reconstructedd, limiter)
    IMPLICIT NONE
    REAL*8, PARAMETER :: eps=1e-10
    INTEGER, INTENT(IN) :: limiter
    REAL*8, INTENT(IN) :: a_mm, a_m, a
    REAL*8, INTENT(IN) :: a_mmd, a_md, ad
    REAL*8, INTENT(OUT) :: reconstructed
    REAL*8, INTENT(OUT) :: reconstructedd
    REAL*8 :: r_l, delta_minus, delta_plus
    REAL*8 :: r_ld, delta_minusd, delta_plusd
    INTRINSIC DABS
    DOUBLE PRECISION :: dabs0
    DOUBLE PRECISION :: dabs1
    REAL*8 :: result1
    REAL*8 :: result1d
    delta_minusd = a_md - a_mmd
    delta_minus = a_m - a_mm
    delta_plusd = ad - a_md
    delta_plus = a - a_m
    IF (delta_minus .GE. 0.) THEN
      dabs0 = delta_minus
    ELSE
      dabs0 = -delta_minus
    END IF
    IF (delta_plus .GE. 0.) THEN
      dabs1 = delta_plus
    ELSE
      dabs1 = -delta_plus
    END IF
    IF (dabs0 .GT. eps .AND. dabs1 .GT. eps) THEN
      r_ld = (delta_plusd-delta_plus*delta_minusd/delta_minus)/&
&       delta_minus
      r_l = delta_plus/delta_minus
      result1d = PSI_VA_D(r_l, r_ld, limiter, result1)
      reconstructedd = a_md + 0.5d0*(delta_minus*result1d+result1*&
&       delta_minusd)
      reconstructed = a_m + 0.5d0*result1*delta_minus
      RETURN
    ELSE
      reconstructedd = a_md
      reconstructed = a_m
      RETURN
    END IF
  END SUBROUTINE MUSCL2_D

!  Differentiation of weno3 in forward (tangent) mode:
!   variations   of useful results: reconstructed
!   with respect to varying inputs: a_m a_p a
!=======================================================================
!=================== WENO 3RD ORDER RECONSTRUCTION =====================
!=======================================================================
  SUBROUTINE WENO3_D(a_m, a_md, a, ad, a_p, a_pd, reconstructed, &
&   reconstructedd)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: a_m, a, a_p
    REAL*8, INTENT(IN) :: a_md, ad, a_pd
    REAL*8, INTENT(OUT) :: reconstructed
    REAL*8, INTENT(OUT) :: reconstructedd
    REAL*8, PARAMETER :: d0=3.d0/4.d0, d1=1.d0/4.d0, eps=1e-6
    REAL*8 :: beta0, beta1, alfa0, alfa1, alfa_sum, omega0, omega1
    REAL*8 :: beta0d, beta1d, alfa0d, alfa1d, alfa_sumd, omega0d, &
&   omega1d
    REAL*8 :: a0, a1
    REAL*8 :: a0d, a1d
    REAL*8 :: temp
    a0d = 0.5d0*(ad+a_pd)
    a0 = 0.5d0*(a+a_p)
    a1d = 1.5d0*ad - 0.5d0*a_md
    a1 = -(0.5d0*a_m) + 1.5d0*a
    beta0d = 2*(a-a_p)*(ad-a_pd)
    beta0 = (a-a_p)**2
    beta1d = 2*(a_m-a)*(a_md-ad)
    beta1 = (a_m-a)**2
    temp = d0/((eps+beta0)*(eps+beta0))
    alfa0d = -(temp*2*beta0d/(eps+beta0))
    alfa0 = temp
    temp = d1/((eps+beta1)*(eps+beta1))
    alfa1d = -(temp*2*beta1d/(eps+beta1))
    alfa1 = temp
    alfa_sumd = alfa0d + alfa1d
    alfa_sum = alfa0 + alfa1
    omega0d = (alfa0d-alfa0*alfa_sumd/alfa_sum)/alfa_sum
    omega0 = alfa0/alfa_sum
    omega1d = (alfa1d-alfa1*alfa_sumd/alfa_sum)/alfa_sum
    omega1 = alfa1/alfa_sum
    reconstructedd = a0*omega0d + omega0*a0d + a1*omega1d + omega1*a1d
    reconstructed = omega0*a0 + omega1*a1
  END SUBROUTINE WENO3_D

end module
