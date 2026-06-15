module ADdiff_bc_kernels_module
implicit none
contains

!  Differentiation of apply_fixed_value_bc in forward (tangent) mode:
!   variations   of useful results: b_grad
!   with respect to varying inputs: grad_c b_grad val_c
!=======================================================================
!============= FIXED VALUE BOUNDARY CONDITION SUBROUTINES ==============
!=======================================================================
  PURE SUBROUTINE APPLY_FIXED_VALUE_BC_D(dim, n, val_cd, grad_cd, &
&   b_gradd, dist, update_grad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: n(dim), grad_cd(dim)
    REAL*8, INTENT(IN) :: dist
    REAL*8, INTENT(IN) :: val_cd
    LOGICAL, INTENT(IN) :: update_grad
    REAL*8, INTENT(INOUT) :: b_gradd(dim)
    INTEGER :: d
    REAL*8 :: tmpd

    IF (update_grad) THEN
      tmpd = 0.0_8
      DO d=1,dim
        tmpd = tmpd + n(d)*grad_cd(d)
      END DO
      b_gradd = grad_cd - n*tmpd - n*val_cd/dist
    END IF
  END SUBROUTINE

!  Differentiation of apply_extrapolation0_bc in forward (tangent) mode:
!   variations   of useful results: b_val b_grad
!   with respect to varying inputs: grad_c b_grad val_c
!=======================================================================
!============ EXTRAPOLATION BOUNDARY CONDITION SUBROUTINES =============
!=======================================================================	
  PURE SUBROUTINE APPLY_EXTRAPOLATION0_BC_D(dim, val_cd, grad_cd, &
&   b_vald, b_gradd, update_grad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: grad_cd(dim)
    REAL*8, INTENT(IN) :: val_cd
    LOGICAL, INTENT(IN) :: update_grad
    REAL*8, INTENT(INOUT) :: b_vald, b_gradd(dim)
    
    b_vald = val_cd

    IF (update_grad) THEN
      b_gradd = grad_cd
    END IF
  END SUBROUTINE

!  Differentiation of apply_slip_vector_bc in forward (tangent) mode:
!   variations   of useful results: b_grad b_vec
!   with respect to varying inputs: grad_c b_grad vec_c
!=======================================================================
!============= SLIP VECTOR BOUNDARY CONDITION SUBROUTINES ==============
!=======================================================================	
  PURE SUBROUTINE APPLY_SLIP_VECTOR_BC_D(dim, n, vec_cd, &
&   grad_cd, b_vecd, b_gradd, dist, update_grad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: n(dim)
    REAL*8, INTENT(IN) :: vec_cd(dim), grad_cd(dim*dim)
    REAL*8, INTENT(IN) :: dist
    LOGICAL, INTENT(IN) :: update_grad
    REAL*8, INTENT(INOUT) :: b_vecd(dim), b_gradd(dim*dim)
    REAL*8 :: vec_nd, g_nd
    INTEGER :: d, v, off

    vec_nd = 0.0_8
    DO d=1,dim
      vec_nd = vec_nd + n(d)*vec_cd(d)
    END DO
    b_vecd = vec_cd - n*vec_nd

    IF (update_grad) THEN
      DO d=1,dim
        off = (d-1)*dim
        g_nd = 0.0_8
        DO v=1,dim
          g_nd = g_nd + n(v)*grad_cd(off+v)
        END DO
        b_gradd(off+1:off+dim) = grad_cd(off+1:off+dim) - n*g_nd
      END DO
    END IF
  END SUBROUTINE

!  Differentiation of copy_grad_t_from_cell in forward (tangent) mode:
!   variations   of useful results: b_grad
!   with respect to varying inputs: grad_c b_grad
!=======================================================================
!================ AUX GRADIENTS MANIPULATION SUBROUTINES ===============
!=======================================================================
  PURE SUBROUTINE COPY_GRAD_T_FROM_CELL_D(dim, n, grad_cd, b_gradd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: n(dim)
    REAL*8, INTENT(IN) :: grad_cd(dim)
    REAL*8, INTENT(INOUT) :: b_gradd(dim)
    INTEGER :: d
    REAL*8 :: tmpd

    tmpd = 0.0_8
    DO d=1,dim
      tmpd = tmpd + n(d)*grad_cd(d)
    END DO
    DO d=1,dim
      b_gradd(d) = grad_cd(d) - n(d)*tmpd
    END DO
  END SUBROUTINE

!  Differentiation of apply_fixed_gradient_bc in forward (tangent) mode:
!   variations   of useful results: b_val b_grad
!   with respect to varying inputs: grad_c b_grad val_c
!=======================================================================
!============ FIXED GRADIENT BOUNDARY CONDITION SUBROUTINES ============
!=======================================================================
  PURE SUBROUTINE APPLY_FIXED_GRADIENT_BC_D(dim, n, val_cd, grad_cd, b_vald, b_gradd, update_grad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: n(dim)
    REAL*8, INTENT(IN) :: grad_cd(dim)
    REAL*8, INTENT(IN) :: val_cd
    LOGICAL, INTENT(IN) :: update_grad
    REAL*8, INTENT(INOUT) :: b_vald, b_gradd(dim)
    INTEGER :: d
    REAL*8 :: tmpd

    b_vald = val_cd

    IF (update_grad) THEN
      tmpd = 0.0_8
      DO d=1,dim
        tmpd = tmpd + n(d)*grad_cd(d)
      END DO
      b_gradd = grad_cd - n*tmpd
    END IF
  END SUBROUTINE

!  Differentiation of apply_reiman_pressure_bc in forward (tangent) mode:
!   variations   of useful results: b_prs_grad b_prs
!   with respect to varying inputs: vel_c tmp_c grad_prs_c b_prs_grad
!                prs_c
!=======================================================================
!============= ADDITIONAL BOUNDARY CONDITION SUBROUTINES ===============
!=======================================================================
  PURE SUBROUTINE APPLY_REIMAN_PRESSURE_BC_D(dim, prs_c, prs_cd, grad_prs_c, &
&   grad_prs_cd, vel_c, vel_cd, tmp_c, tmp_cd, n, k, r_gas, &
&   b_prsd, b_prs_gradd, update_grad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: prs_c, grad_prs_c(dim), vel_c(dim), tmp_c, n(&
&   dim)
    REAL*8, INTENT(IN) :: prs_cd, grad_prs_cd(dim), vel_cd(dim), tmp_cd
    REAL*8, INTENT(IN) :: k, r_gas
    LOGICAL, INTENT(IN) :: update_grad
    REAL*8, INTENT(INOUT) :: b_prsd, b_prs_gradd(dim)
    REAL*8 :: v_n
    REAL*8 :: v_nd
    REAL*8 :: a_c, ro_c
    REAL*8 :: a_cd, ro_cd
    INTEGER :: d
    INTRINSIC DSQRT
    REAL*8 :: arg1
    REAL*8 :: arg1d
    DOUBLE PRECISION :: arg10
    DOUBLE PRECISION :: arg10d
    DOUBLE PRECISION :: result1
    DOUBLE PRECISION :: result1d
    DOUBLE PRECISION :: temp
    REAL*8 :: temp0
    REAL*8 :: temp1
    DOUBLE PRECISION :: temp2
    DOUBLE PRECISION :: tempd
    v_n = 0.d0
    v_nd = 0.0_8
    DO d=1,dim
      v_nd = v_nd - n(d)*vel_cd(d)
      v_n = v_n - n(d)*vel_c(d)
    END DO
    IF (v_n .GT. 0.d0) THEN
      arg1d = k*r_gas*tmp_cd
      arg1 = k*r_gas*tmp_c
      temp = DSQRT(arg1)
      IF (arg1 .EQ. 0.0) THEN
        a_cd = 0.0_8
      ELSE
        a_cd = arg1d/(2.D0*DSQRT(arg1))
      END IF
      a_c = temp
      temp0 = 2.d0*k/(k-1.d0)
      temp1 = v_n/(2.d0*a_c)
      temp = -((k-1.d0)*temp1) + 1.d0
      temp2 = temp**temp0
      IF (temp .LE. 0.0 .AND. (temp0 .EQ. 0.0 .OR. temp0 .NE. INT(temp0)&
&         )) THEN
        tempd = 0.D0
      ELSE
        tempd = -(temp0*temp**(temp0-1)*(k-1.d0)*(v_nd-temp1*2.d0*a_cd)/&
&         (2.d0*a_c))
      END IF
      b_prsd = temp2*prs_cd + prs_c*tempd
    ELSE
      temp1 = prs_c/(r_gas*tmp_c)
      ro_cd = (prs_cd-temp1*r_gas*tmp_cd)/(r_gas*tmp_c)
      ro_c = temp1
      temp1 = v_n**4
      temp0 = ro_c*ro_c/16.d0
      arg10d = (k+1.d0)**2*(temp1*2*ro_c*ro_cd/16.d0+temp0*4*v_n**3*v_nd&
&       ) + k*(v_n**2*(ro_c*prs_cd+prs_c*ro_cd)+prs_c*ro_c*2*v_n*v_nd)
      arg10 = (k+1.d0)*(k+1.d0)*(temp0*temp1) + k*(prs_c*ro_c*(v_n*v_n))
      temp2 = DSQRT(arg10)
      IF (arg10 .EQ. 0.0) THEN
        result1d = 0.D0
      ELSE
        result1d = arg10d/(2.D0*DSQRT(arg10))
      END IF
      result1 = temp2
      b_prsd = prs_cd + (k+1.d0)*(v_n**2*ro_cd/4.d0+ro_c*2*v_n*v_nd/4.d0&
&       ) + result1d
    END IF
!=== GRAD(FI)_boundary = GRAD(FI)_cell: ===
    IF (update_grad) THEN
      b_prs_gradd = grad_prs_cd
    END IF
  END SUBROUTINE

end module
