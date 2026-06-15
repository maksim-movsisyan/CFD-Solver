do face_idx = 1, n_loop_faces
    left_cell = face_left_cell(face_idx)
	right_cell = face_right_cell(face_idx)
	
	n = face_normal(:, face_idx)
	area = face_area(face_idx)
    
    !=== VARIABLES RECONSTRUCTION: ===
	call reconstruct_iface_states_coupled(dim, dim+2, face_idx, ORDER, LIMITER,&
										  left_cell, right_cell,&
										  face_center, cell_center,&
										  Q, gradQ, Var_L, Var_R)
    
    !=== REIMAN SOLVER FLUXES: ===
	call REIMAN_SOLVER_NAME(dim, k, R_gas, cp, cv, Var_L(1), Var_L(2), Var_L(3), Var_L(4), Var_L(dim+2),&
							Var_R(1), Var_R(2), Var_R(3), Var_R(4), Var_R(dim+2), n, face_fluxes)

    !=== FLUXES UPDATING: ===
	Fluxes(:, left_cell) = Fluxes(:, left_cell) + face_fluxes*area
	Fluxes(:, right_cell) = Fluxes(:, right_cell) - face_fluxes*area
end do
