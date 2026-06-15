module ADdiff_bc_ghost_kernels_module
implicit none
contains
!  Differentiation of apply_fixed_value_ghost_bc in forward (tangent) mode:
!   variations   of useful results: grad_ghost val_ghost
!   with respect to varying inputs: grad_c val_c grad_ghost
!   RW status of diff variables: grad_c:in val_c:in grad_ghost:in-out
!                val_ghost:out
!=======================================================================
!============= FIXED VALUE BOUNDARY CONDITION SUBROUTINES ==============
!=======================================================================
  PURE SUBROUTINE APPLY_FIXED_VALUE_GHOST_BC_D(dim, n, val_cd, grad_cd, val_ghostd, grad_ghostd, dist, update_grad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: n(dim), grad_cd(dim)
    REAL*8, INTENT(IN) :: val_cd, dist
    LOGICAL, INTENT(IN) :: update_grad
    REAL*8, INTENT(INOUT) :: val_ghostd, grad_ghostd(dim)
    INTEGER :: d
    REAL*8 :: b_gradd(dim), tmpd
!=== FI_ghost = 2*FI_fixed - FI_cell: ===
    val_ghostd = -val_cd

!=== Grad(FI)_boundary = t*Grad(FI)_cell_t + n*(FI_fixed - FI_cell)/dist ===
!=== Grad(FI)_ghost = 2*Grad(FI)_boundary - Grad(FI)_cell ===
    IF (update_grad) THEN
      tmpd = 0.0_8
      DO d=1,dim
        tmpd = tmpd + n(d)*grad_cd(d)
      END DO
      b_gradd = grad_cd - n*tmpd - n*val_cd/dist
      grad_ghostd = 2.d0*b_gradd - grad_cd

    END IF
  END SUBROUTINE APPLY_FIXED_VALUE_GHOST_BC_D



!  Differentiation of apply_fixed_gradient_ghost_bc in forward (tangent) mode:
!   variations   of useful results: grad_ghost val_ghost
!   with respect to varying inputs: grad_c val_c grad_ghost
!   RW status of diff variables: grad_c:in val_c:in grad_ghost:in-out
!                val_ghost:out
!=======================================================================
!============ FIXED GRADIENT BOUNDARY CONDITION SUBROUTINES ============
!=======================================================================
  PURE SUBROUTINE APPLY_FIXED_GRADIENT_GHOST_BC_D(dim, n, val_cd, grad_cd, val_ghostd, grad_ghostd, update_grad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: n(dim), grad_cd(dim)
    REAL*8, INTENT(IN) :: val_cd
    LOGICAL, INTENT(IN) :: update_grad
    REAL*8, INTENT(INOUT) :: val_ghostd, grad_ghostd(dim)
    INTEGER :: d
    REAL*8 :: tmpd, b_vald, b_gradd(dim)
!=== FI_boundary = FI_cell + GRAD(FI)_b_n*(r_cb, n): ===
!=== FI_ghost = 2*FI_boundary - FI_cell: ===
    b_vald = val_cd
    val_ghostd = 2.d0*b_vald - val_cd

!=== Grad(FI)_boundary = t*Grad(FI)_cell_t + n*GRAD(FI)_n_fixed ===
!=== Grad(FI)_ghost = 2*Grad(FI)_boundary - Grad(FI)_cell ===
    IF (update_grad) THEN
      tmpd = 0.0_8
      DO d=1,dim
        tmpd = tmpd + n(d)*grad_cd(d)
      END DO
      b_gradd = grad_cd - n*tmpd
      grad_ghostd = 2.d0*b_gradd - grad_cd
    END IF
  END SUBROUTINE APPLY_FIXED_GRADIENT_GHOST_BC_D



!  Differentiation of apply_extrapolation0_ghost_bc in forward (tangent) mode:
!   variations   of useful results: grad_ghost val_ghost
!   with respect to varying inputs: grad_c val_c grad_ghost
!   RW status of diff variables: grad_c:in val_c:in grad_ghost:in-out
!                val_ghost:out
!=======================================================================
!============ EXTRAPOLATION BOUNDARY CONDITION SUBROUTINES =============
!=======================================================================	
  PURE SUBROUTINE APPLY_EXTRAPOLATION0_GHOST_BC_D(dim, val_cd, grad_cd, val_ghostd, grad_ghostd, update_grad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: grad_cd(dim)
    REAL*8, INTENT(IN) :: val_cd
    LOGICAL, INTENT(IN) :: update_grad
    REAL*8, INTENT(INOUT) :: val_ghostd, grad_ghostd(dim)
!=== FI_ghost = FI_cell: ===
    val_ghostd = val_cd
!=== GRAD(FI)_ghost = GRAD(FI)_cell: ===
    IF (update_grad) THEN
      grad_ghostd = grad_cd
    END IF
  END SUBROUTINE APPLY_EXTRAPOLATION0_GHOST_BC_D




!  Differentiation of apply_slip_vector_ghost_bc in forward (tangent) mode:
!   variations   of useful results: grad_ghost vec_ghost
!   with respect to varying inputs: grad_c vec_c grad_ghost
!   RW status of diff variables: grad_c:in vec_c:in grad_ghost:in-out
!                vec_ghost:out
!=======================================================================
!============= SLIP VECTOR BOUNDARY CONDITION SUBROUTINES ==============
!=======================================================================	
  PURE SUBROUTINE APPLY_SLIP_VECTOR_GHOST_BC_D(dim, n, vec_cd, grad_cd, vec_ghostd, grad_ghostd, update_grad)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: n(dim), vec_cd(dim), grad_cd(dim*dim)
    LOGICAL, INTENT(IN) :: update_grad
    REAL*8, INTENT(INOUT) :: vec_ghostd(dim), grad_ghostd(dim*dim)
    REAL*8 :: vec_nd, g_nd, b_gradd(dim*dim)
    INTEGER :: d, v, off
!=== VEC_n_boundary = 0, VEC_t_boundary = VEC_t_cell: ===
    vec_nd = 0.0_8
    DO d=1,dim
      vec_nd = vec_nd + n(d)*vec_cd(d)
    END DO
    vec_ghostd = vec_cd - n*2.d0*vec_nd
!=== GRAD(VEC)_boundary = t*GRAD(VEC)_t_cell + n*0.0 ===
    IF (update_grad) THEN
      b_gradd = 0.0_8
      DO d=1,dim
        off = (d-1)*dim
        g_nd = 0.0_8
        DO v=1,dim
          g_nd = g_nd + n(v)*grad_cd(off+v)
        END DO
        b_gradd(off+1:off+dim) = grad_cd(off+1:off+dim) - n*g_nd
      END DO
      grad_ghostd = 2.d0*b_gradd - grad_cd
    END IF
  END SUBROUTINE APPLY_SLIP_VECTOR_GHOST_BC_D

end module
