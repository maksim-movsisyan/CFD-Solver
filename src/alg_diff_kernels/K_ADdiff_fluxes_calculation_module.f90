module ADdiff_fluxes_calculation_module
use ADdiff_reiman_solver_module
use ADdiff_reconstruction_module
use ADdiff_physical_properties_module
implicit none

contains

!  Differentiation of compute_inviscid_fluxes_cmprs in forward (tangent) mode:
!   variations   of useful results: fluxes
!   with respect to varying inputs: q gradq fluxes
!   RW status of diff variables: q:in gradq:in fluxes:in-out
!=======================================================================
!=========== CMPRS INVISCID FLUXES CALCLUCATION SUBROUTINE: ============
!=======================================================================
  SUBROUTINE COMPUTE_INVISCID_FLUXES_CMPRS_D(dim, nfaces, nbfaces, &
&   face_left_cell, face_right_cell, face_normal, face_area, face_center&
&   , cell_center, q, qd, gradq, gradqd, fluxesd, k, r_gas, cp, &
&   cv, order, limiter, reiman_solver_type, use_ghost_cells)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim, nfaces, nbfaces
    INTEGER, INTENT(IN) :: order, limiter, reiman_solver_type
    LOGICAL, INTENT(IN) :: use_ghost_cells
    INTEGER, INTENT(IN), CONTIGUOUS :: face_left_cell(:), &
&   face_right_cell(:)
    REAL*8, INTENT(IN), CONTIGUOUS :: face_normal(:, :), face_area(:), &
&   face_center(:, :), cell_center(:, :)
    REAL*8, INTENT(IN), CONTIGUOUS :: q(:, :), gradq(:, :)
    REAL*8, INTENT(IN), CONTIGUOUS :: qd(:, :), gradqd(:, :)
    REAL*8, INTENT(IN) :: k, r_gas, cp, cv
    REAL*8, INTENT(INOUT), CONTIGUOUS :: fluxesd(:, :)
    INTEGER :: n_loop_faces
    INTEGER :: face_idx
    INTEGER :: left_cell, right_cell
    REAL*8 :: n(dim), area
    REAL*8 :: var_l(dim+2), var_r(dim+2)
    REAL*8 :: var_ld(dim+2), var_rd(dim+2)
    REAL*8 :: face_fluxesd(dim+2)
    REAL*8 :: p, v(dim), t, ro, h, v_n, v2
    REAL*8 :: pd, vd(dim), td, rod, hd, v_nd, v2d
    INTEGER :: d
    REAL*8 :: temp
    n_loop_faces = nfaces - nbfaces
    IF (use_ghost_cells) THEN
      n_loop_faces = nfaces
    END IF
    var_ld = 0.0_8
    var_rd = 0.0_8
      
	!=== INTERNAL FACES: ===
	SELECT CASE(reiman_solver_type)
		CASE(1)
#			define REIMAN_SOLVER_DNAME HLL_REIMAN_SOLVER_D
#			include "ADdiff_inviscid_flux_kernel.inc"
#			undef REIMAN_SOLVER_DNAME
		CASE(2)
#			define REIMAN_SOLVER_DNAME HLLC_REIMAN_SOLVER_D
#			include "ADdiff_inviscid_flux_kernel.inc"
#			undef REIMAN_SOLVER_DNAME
		CASE(3)
#			define REIMAN_SOLVER_DNAME ROE_REIMAN_SOLVER_D
#			include "ADdiff_inviscid_flux_kernel.inc"
#			undef REIMAN_SOLVER_DNAME
		CASE(4)
#			define REIMAN_SOLVER_DNAME AUSM_REIMAN_SOLVER_D
#			include "ADdiff_inviscid_flux_kernel.inc"
#			undef REIMAN_SOLVER_DNAME
		CASE(5)
#			define REIMAN_SOLVER_DNAME AUSMPLUS_REIMAN_SOLVER_D
#			include "ADdiff_inviscid_flux_kernel.inc"
#			undef REIMAN_SOLVER_DNAME
	END SELECT
!=== EXTERNAL FACES: ===
    IF (.NOT.use_ghost_cells) THEN
      vd = 0.0_8
!=== UPDATING FLUXES USING BOUNDARY FACES VALUES: ===
      DO face_idx=nfaces-nbfaces+1,nfaces
        left_cell = face_left_cell(face_idx)
        right_cell = face_right_cell(face_idx)
        n = face_normal(:, face_idx)
        area = face_area(face_idx)
        v2 = 0.d0
        v_n = 0.d0
        pd = qd(1, right_cell)
        p = q(1, right_cell)
        v2d = 0.0_8
        v_nd = 0.0_8
        DO d=1,dim
          vd(d) = qd(1+d, right_cell)
          v(d) = q(1+d, right_cell)
          v2d = v2d + 2*v(d)*vd(d)
          v2 = v2 + v(d)**2
          v_nd = v_nd + n(d)*vd(d)
          v_n = v_n + v(d)*n(d)
        END DO
        td = qd(dim+2, right_cell)
        t = q(dim+2, right_cell)
        temp = p/(r_gas*t)
        rod = (pd-temp*r_gas*td)/(r_gas*t)
        ro = temp
        hd = cp*td + 0.5d0*v2d
        h = cp*t + 0.5d0*v2
        face_fluxesd(1) = v_n*rod + ro*v_nd
        DO d=1,dim
          face_fluxesd(1+d) = ro*v_n*vd(d) + v(d)*(v_n*rod+ro*v_nd) + n(&
&           d)*pd
        END DO
        face_fluxesd(dim+2) = v_n*(h*rod+ro*hd) + ro*h*v_nd
!=== FLUXES UPDATING: ===
        fluxesd(:, left_cell) = fluxesd(:, left_cell) + area*&
&         face_fluxesd
      END DO
    END IF
  END SUBROUTINE

!  Differentiation of compute_viscous_fluxes_cmprs in forward (tangent) mode:
!   variations   of useful results: fluxes
!   with respect to varying inputs: q gradq fluxes
!   RW status of diff variables: q:in gradq:in fluxes:in-out
!=======================================================================
!=========== CMPRS VISCOUS FLUXES CALCLUCATION SUBROUTINE: =============
!=======================================================================
  SUBROUTINE COMPUTE_VISCOUS_FLUXES_CMPRS_D(dim, nfaces, nbfaces, &
&   face_left_cell, face_right_cell, face_normal, face_area, face_weight&
&   , face_center, cell_center, q, qd, gradq, gradqd, fluxesd, k&
&   , r_gas, cp, cv, pr, use_ghost_cells)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim, nfaces, nbfaces
    LOGICAL, INTENT(IN) :: use_ghost_cells
    INTEGER, INTENT(IN), CONTIGUOUS :: face_left_cell(:), &
&   face_right_cell(:)
    REAL*8, INTENT(IN), CONTIGUOUS :: face_normal(:, :), face_area(:), &
&   face_weight(:), face_center(:, :), cell_center(:, :)
    REAL*8, INTENT(IN), CONTIGUOUS :: q(:, :), gradq(:, :)
    REAL*8, INTENT(IN), CONTIGUOUS :: qd(:, :), gradqd(:, :)
    REAL*8, INTENT(IN) :: k, r_gas, cp, cv, pr
    REAL*8, INTENT(INOUT), CONTIGUOUS :: fluxesd(:, :)
    INTEGER :: n_loop_faces
    INTEGER :: face_idx, d, v
    INTEGER :: left_cell, right_cell
    INTEGER :: offset, v_offset
    REAL*8 :: w, n(dim), ksi(dim), area, dist
    REAL*8 :: q_f(dim+2), gradq_f(dim*(dim+2))
    REAL*8 :: q_fd(dim+2), gradq_fd(dim*(dim+2))
    REAL*8 :: divv, lambda, mu
    REAL*8 :: divvd, lambdad, mud
    REAL*8 :: qn_f, taun_f(dim)
    REAL*8 :: qn_fd, taun_fd(dim)
    REAL*8 :: face_fluxesd(dim+2)
    REAL*8 :: dot_prct
    REAL*8 :: dot_prctd
    INTRINSIC DSQRT
    REAL*8 :: temp
    n_loop_faces = nfaces - nbfaces
    IF (use_ghost_cells) THEN
      n_loop_faces = nfaces
    END IF
!=== INTERNAL FACES: ===
    DO face_idx=1,n_loop_faces
      left_cell = face_left_cell(face_idx)
      right_cell = face_right_cell(face_idx)
      w = face_weight(face_idx)
      n = face_normal(:, face_idx)
      area = face_area(face_idx)
      ksi = cell_center(:, right_cell) - cell_center(:, left_cell)
      dist = 0.d0
      DO d=1,dim
        dist = dist + ksi(d)**2
      END DO
      dist = DSQRT(dist)
      ksi = ksi/dist
!=== LINEAR INTERPOLATION ON FACE: ===
      q_fd = w*qd(:, right_cell) + (1.d0-w)*qd(:, left_cell)
      q_f = w*q(:, right_cell) + (1.d0-w)*q(:, left_cell)
      gradq_fd = w*gradqd(:, right_cell) + (1.d0-w)*gradqd(:, left_cell)
      gradq_f = w*gradq(:, right_cell) + (1.d0-w)*gradq(:, left_cell)
!=== PHYSICAL TRANSPORT PROPERTIES: ===
      mud = MU_AIR_D(q_f(dim+2), q_fd(dim+2), mu)
      lambdad = LAMBDA_AIR_D(q_f(dim+2), q_fd(dim+2), pr, cp, lambda)
!=== ENERGY DIFFUSION: ===
      offset = dim*(dim+1)
      dot_prct = 0.d0
      dot_prctd = 0.0_8
      DO d=1,dim
        dot_prctd = dot_prctd + (n(d)-ksi(d))*gradq_fd(offset+d)
        dot_prct = dot_prct + gradq_f(offset+d)*(n(d)-ksi(d))
      END DO
      temp = (q(dim+2, right_cell)-q(dim+2, left_cell))/dist + dot_prct
      qn_fd = -(temp*lambdad+lambda*((qd(dim+2, right_cell)-qd(dim+2, &
&       left_cell))/dist+dot_prctd))
      qn_f = -(lambda*temp)
!=== MOMENTUM DIFFUSION: ===
      divv = 0.d0
      divvd = 0.0_8
      DO d=1,dim
        offset = d*dim
        dot_prct = 0.d0
        dot_prctd = 0.0_8
        DO v=1,dim
          dot_prctd = dot_prctd + (n(v)-ksi(v))*gradq_fd(offset+v)
          dot_prct = dot_prct + gradq_f(offset+v)*(n(v)-ksi(v))
        END DO
        taun_fd(d) = (qd(d+1, right_cell)-qd(d+1, left_cell))/dist + &
&         dot_prctd
        taun_f(d) = (q(d+1, right_cell)-q(d+1, left_cell))/dist + &
&         dot_prct
        DO v=1,dim
          v_offset = v*dim
          taun_fd(d) = taun_fd(d) + n(v)*gradq_fd(d+v_offset)
          taun_f(d) = taun_f(d) + gradq_f(d+v_offset)*n(v)
        END DO
        divvd = divvd + gradq_fd(d+offset)
        divv = divv + gradq_f(d+offset)
      END DO
      DO d=1,dim
        temp = taun_f(d) - 2.0d0*n(d)*divv/3.0d0
        taun_fd(d) = temp*mud + mu*(taun_fd(d)-n(d)*2.0d0*divvd/3.0d0)
        taun_f(d) = mu*temp
      END DO
!=== DIFFUSION FLUXES: ===
      face_fluxesd = 0.0_8
      DO d=1,dim
        face_fluxesd(d+1) = -taun_fd(d)
      END DO
      dot_prctd = 0.0_8
      DO v=1,dim
        dot_prctd = dot_prctd + q_f(v+1)*taun_fd(v) + taun_f(v)*q_fd(v+1&
&         )
      END DO
      face_fluxesd(dim+2) = qn_fd - dot_prctd
!=== FLUXES UPDATING: ===
      fluxesd(:, left_cell) = fluxesd(:, left_cell) + area*face_fluxesd
      fluxesd(:, right_cell) = fluxesd(:, right_cell) - area*&
&       face_fluxesd
    END DO
!=== EXTERNAL FACES: ===
    IF (.NOT.use_ghost_cells) THEN
!=== UPDATING FLUXES USING BOUNDARY FACES VALUES: ===
      DO face_idx=nfaces-nbfaces+1,nfaces
        left_cell = face_left_cell(face_idx)
        right_cell = face_right_cell(face_idx)
        n = face_normal(:, face_idx)
        area = face_area(face_idx)
!=== FACE VALUES: ===
        q_fd = qd(:, right_cell)
        q_f = q(:, right_cell)
        gradq_fd = gradqd(:, right_cell)
        gradq_f = gradq(:, right_cell)
!=== PHYSICAL TRANSPORT PROPERTIES: ===
        mud = MU_AIR_D(q_f(dim+2), q_fd(dim+2), mu)
        lambdad = LAMBDA_AIR_D(q_f(dim+2), q_fd(dim+2), pr, cp, lambda)
!=== ENERGY DIFFUSION: ===
        offset = dim*(dim+1)
        dot_prct = 0.d0
        dot_prctd = 0.0_8
        DO v=1,dim
          dot_prctd = dot_prctd + n(v)*gradq_fd(offset+v)
          dot_prct = dot_prct + gradq_f(offset+v)*n(v)
        END DO
        qn_fd = -(dot_prct*lambdad+lambda*dot_prctd)
        qn_f = -(lambda*dot_prct)
!=== MOMENTUM DIFFUSION: ===
        taun_f = 0.d0
        divv = 0.d0
        taun_fd = 0.0_8
        divvd = 0.0_8
        DO d=1,dim
          offset = d*dim
          DO v=1,dim
            v_offset = v*dim
            taun_fd(d) = taun_fd(d) + n(v)*(gradq_fd(v+offset)+gradq_fd(&
&             d+v_offset))
            taun_f(d) = taun_f(d) + (gradq_f(v+offset)+gradq_f(d+&
&             v_offset))*n(v)
          END DO
          divvd = divvd + gradq_fd(d+offset)
          divv = divv + gradq_f(d+offset)
        END DO
        DO d=1,dim
          temp = taun_f(d) - 2.d0*n(d)*divv/3.d0
          taun_fd(d) = temp*mud + mu*(taun_fd(d)-n(d)*2.d0*divvd/3.d0)
          taun_f(d) = mu*temp
        END DO
!=== DIFFUSION FLUXES: ===
        face_fluxesd = 0.0_8
        DO d=1,dim
          face_fluxesd(d+1) = -taun_fd(d)
        END DO
        dot_prctd = 0.0_8
        DO v=1,dim
          dot_prctd = dot_prctd + q_f(v+1)*taun_fd(v) + taun_f(v)*q_fd(v&
&           +1)
        END DO
        face_fluxesd(dim+2) = qn_fd - dot_prctd
!=== FLUXES UPDATING: ===
        fluxesd(:, left_cell) = fluxesd(:, left_cell) + area*&
&         face_fluxesd
      END DO
    END IF
  END SUBROUTINE COMPUTE_VISCOUS_FLUXES_CMPRS_D

end module
