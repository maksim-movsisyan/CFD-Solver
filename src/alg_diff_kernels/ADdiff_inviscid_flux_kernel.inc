DO face_idx=1,n_loop_faces
  left_cell = face_left_cell(face_idx)
  right_cell = face_right_cell(face_idx)
  n = face_normal(:, face_idx)
  area = face_area(face_idx)
!=== VARIABLES RECONSTRUCTION: ===
  CALL RECONSTRUCT_IFACE_STATES_COUPLED_D(dim, dim + 2, face_idx, &
&                                       order, limiter, left_cell, &
&                                       right_cell, face_center, &
&                                       cell_center, q, qd, gradq, &
&                                       gradqd, var_l, var_ld, var_r, &
&                                       var_rd)
!=== REIMAN SOLVER FLUXES: ===
  face_fluxesd = 0.0_8
  CALL REIMAN_SOLVER_DNAME(dim, k, r_gas, cp, cv, var_l(1), var_ld(1&
&                        ), var_l(2), var_ld(2), var_l(3), var_ld(3), &
&                        var_l(4), var_ld(4), var_l(dim+2), var_ld(dim+2&
&                        ), var_r(1), var_rd(1), var_r(2), var_rd(2), &
&                        var_r(3), var_rd(3), var_r(4), var_rd(4), var_r&
&                        (dim+2), var_rd(dim+2), n, &
&                        face_fluxesd)
!=== FLUXES UPDATING: ===
  fluxesd(:, left_cell) = fluxesd(:, left_cell) + area*face_fluxesd
  fluxesd(:, right_cell) = fluxesd(:, right_cell) - area*face_fluxesd
END DO
