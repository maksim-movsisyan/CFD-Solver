do face_idx = 1, nfaces - nbfaces
	L = face_left_cell(face_idx)
	R = face_right_cell(face_idx)
	
	n = face_normal(:, face_idx)
	area = face_area(face_idx)
	g = face_weight(face_idx)
	
	!=== STABILITY OPERATOR: ===
	call STABILITY_OPERATOR_NAME(dim, k, R_gas, cp, cv, Q(1, L), Q(2, L), Q(3, L), Q(4, L), Q(dim+2, L),&
														Q(1, R), Q(2, R), Q(3, R), Q(4, R), Q(dim+2, R), g, n, Mtrx_L, Mtrx_R)
	
	!=== JACOBIAN UPDATING: ===
	values(:, :, map_LL(face_idx)) = values(:, :, map_LL(face_idx)) + Mtrx_L*area
	values(:, :, map_LR(face_idx)) = values(:, :, map_LR(face_idx)) + Mtrx_R*area
	values(:, :, map_RL(face_idx)) = values(:, :, map_RL(face_idx)) - Mtrx_L*area
	values(:, :, map_RR(face_idx)) = values(:, :, map_RR(face_idx)) - Mtrx_R*area
	
end do
