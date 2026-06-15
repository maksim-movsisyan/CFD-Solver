module metric3D_computation_module
implicit none
contains
!=======================================================================
!================ SUBROUTINE FOR COMPUTING 2D METRIC  ================== 
!=======================================================================
	pure subroutine compute_3d_face_centroid_and_area(nfaces, dim, node_coords,&
													  face_nodes_ptr, face_nodes,&
													  face_area, face_center)
		implicit none
		integer, intent(in) :: nfaces, dim
		integer, intent(in), contiguous :: face_nodes_ptr(:), face_nodes(:)
		real(8), intent(in), contiguous :: node_coords(:, :)
		real(8), intent(inout), contiguous :: face_area(:), face_center(:,:)
		
		integer :: face_idx
		integer :: pos1, pos2, nnodes
		real(8) :: r_g(dim), r(dim), r1(dim), r2(dim)
		real(8) :: dr1(dim), dr2(dim), tmp(dim)
		real(8) :: total_area, tri_area
		real(8) :: tri_center(dim), r_center_sum(dim)
		integer :: node_idx, node1_idx, node2_idx
		integer :: i
		
		
		do face_idx = 1, nfaces
			pos1 = face_nodes_ptr(face_idx)
			pos2 = face_nodes_ptr(face_idx + 1) - 1
			nnodes = pos2 - pos1  + 1
			
			!=== GEOMETRIC FACE CENTER ====
			r_g = 0.d0
			do i = pos1, pos2
				node_idx = face_nodes(i)
				r_g = r_g + node_coords(:, node_idx)
			end do 
			r_g = r_g/nnodes
			
			!=== TRIANGLE CENTROIDS AND AREAS ===
			!=== FACE AREA = SUM(TRIANGEL_AREAS) ===
			!=== FACE CENTROID = AREA_AVG(TRIANGEL CENTROIDS) ===
			total_area = 0.0d0
			r_center_sum = 0.0d0
			do i = 1, nnodes
				node1_idx = face_nodes(pos1 + i - 1)
				if (i == nnodes) then
					node2_idx = face_nodes(pos1)
				else
					node2_idx = face_nodes(pos1 + i)
				end if
				
				r1 = node_coords(:, node1_idx)
				r2 = node_coords(:, node2_idx)
				
				dr1 = r2 - r1
				dr2 = r_g - r1
				tmp = cross_product_3d(dr1, dr2)
				
				tri_area = 0.5d0*Norm2(tmp)
				tri_center = (1.0d0/3.0d0)*(r1 + r2 + r_g)
				
				total_area = total_area + tri_area
				r_center_sum = r_center_sum + tri_area*tri_center
			end do
			
			face_area(face_idx) = total_area
			face_center(:, face_idx) = r_center_sum/total_area	
			
		end do
	end subroutine

	pure subroutine compute_3d_face_normal(nfaces, dim, node_coords,&
										   face_nodes_ptr, face_nodes,&
										   face_normal)
		implicit none
		integer, intent(in) :: nfaces, dim
		integer, intent(in), contiguous :: face_nodes_ptr(:), face_nodes(:)
		real(8), intent(in), contiguous :: node_coords(:, :)
		real(8), intent(inout), contiguous :: face_normal(:,:)
		
		integer :: face_idx
		integer :: pos1, pos2, nnodes
		integer :: node1_idx, node2_idx
		real(8) :: r1(3), r2(3), normal(3)
		integer :: i
		
		do face_idx = 1, nfaces
			pos1 = face_nodes_ptr(face_idx)
			pos2 = face_nodes_ptr(face_idx + 1) - 1
			nnodes = pos2 - pos1 + 1
			
			normal = 0.d0
			
			do i = 1, nnodes
				node1_idx = face_nodes(pos1 + i - 1)
				if (i == nnodes) then
					node2_idx = face_nodes(pos1)
				else
					node2_idx = face_nodes(pos1 + i)
				end if
				
				r1 = node_coords(:, node1_idx)
				r2 = node_coords(:, node2_idx)
				
				normal(1) = normal(1) + (r1(2) - r2(2)) * (r1(3) + r2(3))
				normal(2) = normal(2) + (r1(3) - r2(3)) * (r1(1) + r2(1))
				normal(3) = normal(3) + (r1(1) - r2(1)) * (r1(2) + r2(2))
			end do
			face_normal(:, face_idx) = normal/norm2(normal)
		end do
	end subroutine
	  
	pure subroutine compute_3d_order_face_normal_direction(nfaces, nbfaces, dim,&
														   face_left_cell, face_right_cell,&
														   cell_center, face_center,&
												           face_normal)
			implicit none
			integer, intent(in) :: nfaces, nbfaces, dim
			integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:)
			real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :)
			real(8), intent(inout), contiguous :: face_normal(:,:)
			
			integer :: face_idx
			integer :: left_cell, right_cell
			real(8) :: direction_indicator, r_lr(dim)
			
			do face_idx = 1, nfaces
				!=== FACE VECTOR DIRECTION: FROM LEFT CELL TO RIGHT CELL ===
				!=== FOR BOUNDARY FACES FACE VECTOR DIRECTE OUTWARD ===
				left_cell = face_left_cell(face_idx)
				right_cell = face_right_cell(face_idx)
				if (face_idx <= nfaces - nbfaces) then
					r_lr(:) = cell_center(:, right_cell) - cell_center(:, left_cell)
				else 
					r_lr(:) = face_center(:, face_idx) - cell_center(:, left_cell)
				end if
				
				direction_indicator = DOT_PRODUCT(face_normal(:, face_idx), r_lr)
				if (direction_indicator < 0.d0) then
					face_normal(:, face_idx) = -face_normal(:, face_idx)
				end if
			end do
		end subroutine

	pure subroutine compute_3d_face_weight(nfaces, nbfaces, dim,&
										   face_left_cell, face_right_cell,&
										   face_center, cell_center,&
										   face_normal, face_weight)
		implicit none
		integer, intent(in) :: nfaces, nbfaces, dim
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:)
		real(8), intent(in), contiguous :: face_center(:, :), cell_center(:, :),&
										   face_normal(:, :)
		real(8), intent(inout), contiguous :: face_weight(:)
		
		integer :: face_idx
		integer :: left_cell, right_cell
		real(8) :: Lf, Rf
		real(8) :: e_f(dim), d_Cf(dim), d_fF(dim)
				
		face_weight = 0.5d0
		do face_idx = 1, nfaces - nbfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			e_f = face_normal(:, face_idx)
			
			d_Cf = face_center(:, face_idx) - cell_center(:, left_cell)
			d_fF = cell_center(:, right_cell) - face_center(:, face_idx)
			
			Lf = DOT_PRODUCT(e_f, d_Cf)
			Rf = DOT_PRODUCT(e_f, d_fF)
			
			!===FI(FACE_CENTROID) = G*RIGHT_CELL_VAL + (1 - G)*LEFT_CELL_VAL====
			face_weight(face_idx) = Lf/(Lf + Rf)	
		end do
	end subroutine

	pure subroutine compute_3d_cell_volume_and_centroid(ncells, dim,&
													    cell_nodes_ptr, cell_nodes,&
													    cell_faces_ptr, cell_faces,&
													    node_coords, face_center,&
													    face_normal, face_area,&
													    cell_volume, cell_center)
		implicit none
		integer, intent(in) :: ncells, dim	
		integer, intent(in), contiguous :: cell_nodes_ptr(:), cell_nodes(:),&
										   cell_faces_ptr(:), cell_faces(:)	
		real(8), intent(in), contiguous :: node_coords(:, :), face_center(:, :),&
										   face_normal(:, :), face_area(:)
		real(8), intent(inout), contiguous :: cell_volume(:), cell_center(:, :)			    
		
		integer :: cell_idx, face_idx, node_idx													
		integer :: pos1n, pos2n
		integer :: pos1f, pos2f
		real(8) :: r_g(dim), r_gP(dim), S_f(dim), d_Gf(dim)
		real(8) :: total_volume, pyrmd_volume
		real(8) :: pyrmd_center(dim), r_center_sum(dim)
		integer :: i, nnodes, nfaces
		
		
		do cell_idx = 1, ncells
			pos1n = cell_nodes_ptr(cell_idx)
			pos2n = cell_nodes_ptr(cell_idx + 1) - 1
			pos1f = cell_faces_ptr(cell_idx)
			pos2f = cell_faces_ptr(cell_idx + 1) - 1
			
			nnodes = pos2n - pos1n + 1
			nfaces = pos2f - pos1f + 1
			
			!=== GEOMETRIC CELL CENTER ====
			r_g = 0.d0
			do i = pos1n, pos2n
				node_idx = cell_nodes(i)
				r_g = r_g + node_coords(:, node_idx)
			end do 
			r_g = r_g/nnodes
			
			!=== PYRAMID CENTROIDS AND AREAS ===
			!=== CELL VOLUME = SUM(PYRAMID AREAS) ===
			!=== CELL CENTROID = VOLUME_AVG(PYRAMID CENTROIDS) ===
			total_volume = 0.0d0
			r_center_sum = 0.0d0
			do i = 1, nfaces
				face_idx = cell_faces(pos1f + i - 1)
				S_f = face_normal(:, face_idx)*face_area(face_idx)
				d_Gf = face_center(:, face_idx) - r_g
				
				pyrmd_center = 0.75d0*face_center(:, face_idx) + 0.25d0*r_g
				pyrmd_volume = (1.d0/3.d0)*dabs(DOT_PRODUCT(S_f, d_Gf))
				
				total_volume = total_volume + pyrmd_volume
				r_center_sum = r_center_sum + pyrmd_volume*pyrmd_center
				
			end do
			
			cell_volume(cell_idx) = total_volume
			cell_center(:, cell_idx) = r_center_sum/total_volume	
			
		end do		

	end subroutine 
	
	pure subroutine compute_3d_ghost_cell_centroid(dim, nfaces, nbfaces,&
												   face_left_cell, face_right_cell,&
												   face_normal, face_center, cell_center)
		implicit none
		integer, intent(in) :: dim, nfaces, nbfaces
		integer, intent(in) :: face_left_cell(:), face_right_cell(:)
		real(8), intent(in) :: face_normal(:, :), face_center(:, :)
		real(8), intent(inout) :: cell_center(:, :)
		
		real(8) :: n(dim), C_L(dim), C_f(dim), dr(dim), d_normal
		integer :: face_idx
		integer :: left_cell, right_cell
		
		!=== BOUNDARY FACES: ===
		do face_idx = nfaces - nbfaces + 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			n = face_normal(:, face_idx)
			C_L = cell_center(:, left_cell)
			C_f = face_center(:, face_idx)
			
			dr = C_f - C_L
			
			d_normal = dot_product(dr, n)
			
			cell_center(:, right_cell) = C_L + 2.d0*d_normal*n
		end do
	end subroutine
		
		
	pure subroutine compute_metric_3d(ncells, nfaces, nbfaces, dim,&
									  face_left_cell, face_right_cell,&
									  face_nodes_ptr, face_nodes,&
									  cell_nodes_ptr, cell_nodes,&
									  cell_faces_ptr, cell_faces,&
									  node_coords,&
									  face_center, face_normal,&
									  face_area, face_weight,&
									  cell_center, cell_volume)
		implicit none
		integer, intent(in) :: ncells, nfaces, nbfaces, dim
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
										   face_nodes_ptr(:), face_nodes(:),&
										   cell_nodes_ptr(:), cell_nodes(:),&
										   cell_faces_ptr(:), cell_faces(:)
		real(8), intent(in), contiguous :: node_coords(:, :)
		real(8), intent(inout), contiguous :: face_center(:, :), face_normal(:, :),&
											  face_area(:), face_weight(:),&
											  cell_center(:, :), cell_volume(:)
		
		
		call compute_3d_face_centroid_and_area(nfaces, dim, node_coords,&
											   face_nodes_ptr, face_nodes,&
											   face_area, face_center)
		call compute_3d_face_normal(nfaces, dim, node_coords,&
									face_nodes_ptr, face_nodes,&
									face_normal)
		call compute_3d_cell_volume_and_centroid(ncells, dim,&
											     cell_nodes_ptr, cell_nodes,&
												 cell_faces_ptr, cell_faces,&
												 node_coords, face_center,&
												 face_normal, face_area,&
												 cell_volume, cell_center)
		call compute_3d_order_face_normal_direction(nfaces, nbfaces, dim,&
													face_left_cell, face_right_cell,&
													cell_center, face_center,&
												    face_normal)
		call compute_3d_face_weight(nfaces, nbfaces, dim,&
									face_left_cell, face_right_cell,&
									face_center, cell_center,&
									face_normal, face_weight)
		call compute_3d_ghost_cell_centroid(dim, nfaces, nbfaces,&
										    face_left_cell, face_right_cell,&
										    face_normal, face_center, cell_center)
	end subroutine
	
	
!=======================================================================
!======================== AUXILIARY SUBROUTINES ========================
!=======================================================================
	pure function cross_product_3d(a, b) result(cross)
		implicit none
		real(8), intent(in) :: a(3), b(3)
		real(8) :: cross(3)
		
		cross(1) = a(2)*b(3) - a(3)*b(2)
		cross(2) = a(3)*b(1) - a(1)*b(3)
		cross(3) = a(1)*b(2) - a(2)*b(1)
	end function 








end module
