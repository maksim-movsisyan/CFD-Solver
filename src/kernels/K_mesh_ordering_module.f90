module mesh_ordering_module
use mesh_connectivity_module
use array_sorting_module
use metric3D_computation_module, only: cross_product_3d
implicit none
contains
!=======================================================================
!====== SUBROTUINE FOR LEFT/RIGHT CELLS RERODERING (LEFT->REAL) ========
!=======================================================================
	subroutine reorder_left_right_cells(nfaces, ncells,&
									    face_left_cell, face_right_cell,&
									    face_nodes_ptr, face_nodes)
		implicit none
		integer, intent(in) :: nfaces, ncells
		integer, intent(in) :: face_nodes_ptr(:)
		integer, intent(inout) :: face_left_cell(:), face_right_cell(:),&
								  face_nodes(:)
		
		integer face_idx
		integer :: left_cell, right_cell
		integer :: j, num_nodes, pos1, pos2
		
		!=== LEFT CELL ALWAYS REAL, LEFT CELL IDX < RIGHT_CELL IDX ===
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if (left_cell == 0 .or. left_cell > ncells) then
				call swap_int(face_left_cell(face_idx), face_right_cell(face_idx))
				
				pos1 = face_nodes_ptr(face_idx)
				pos2 = face_nodes_ptr(face_idx+1) - 1
				num_nodes = pos2 - pos1 + 1
				do j = 1, num_nodes/2
					call swap_int(face_nodes(pos1 + j - 1), face_nodes(pos2 - j + 1))
				end do
				
			else if ((left_cell > right_cell) .and. (right_cell .ne. 0)) then
				call swap_int(face_left_cell(face_idx), face_right_cell(face_idx))
				
				pos1 = face_nodes_ptr(face_idx)
				pos2 = face_nodes_ptr(face_idx+1) - 1
				num_nodes = pos2 - pos1 + 1
				do j = 1, num_nodes/2
					call swap_int(face_nodes(pos1 + j - 1), face_nodes(pos2 - j + 1))
				end do
			end if
		end do
	end subroutine
		
!=======================================================================
!================== SUBROUTINE FOR RCM CELLS REORDERING ================
!=======================================================================
	subroutine reorder_cells_rcm(ncells, nfaces, face_left_cell, face_right_cell,&
								 face_nodes_ptr, face_nodes, cell_type)
		implicit none
		integer, intent(in) :: ncells, nfaces
		integer, intent(in) :: face_nodes_ptr(:)
		integer, intent(inout) :: face_left_cell(:), face_right_cell(:),&
								  face_nodes(:), cell_type(:)
		
		integer :: cell_idx, face_idx
		integer :: head, tail, cur_idx, j, k, curr, neighbor, n_neighs, tmp_c
		integer :: left_cell, right_cell, start_node, tmp_val
		integer, allocatable :: neighs(:)
		integer, allocatable :: c2c_ptr(:), c2c_neighs(:)
		integer, allocatable :: degree(:), new_to_old(:), old_to_new(:), queue(:)
		logical, allocatable :: visited(:)

		!=== AUX CELL-TO-CELL CONNECTIVITY: ===
		call cell_neighbor_connectivity(ncells, nfaces,&
										face_left_cell, face_right_cell,&
										c2c_ptr, c2c_neighs)			
		allocate(degree(ncells), new_to_old(ncells), old_to_new(ncells),&
				 queue(ncells), visited(ncells))
		
		!=== DEGREE (=NUMBER OF NEIGHBORS): ===
		do cell_idx = 1, ncells
			degree(cell_idx) = c2c_ptr(cell_idx+1) - c2c_ptr(cell_idx)
		end do
		
		!=== RCM ALGHORYTHM: ===
		visited = .false.; cur_idx = 1
		do cell_idx = 1, ncells
			if (visited(cell_idx)) cycle
			
			!=== FINDING CELL WITH MIN(DEGREE): ===
			start_node = cell_idx
			do j = cell_idx + 1, ncells
				if (.not. visited(j)) then
					if (degree(j) < degree(start_node)) start_node = j
				end if
				if (degree(start_node) <= 1) exit 
			end do

			head = 1; tail = 1
			queue(tail) = start_node
			visited(start_node) = .true.
			
			do while (head <= tail)
				curr = queue(head); head = head + 1
				new_to_old(cur_idx) = curr
				cur_idx = cur_idx + 1
				
				!=== FINDING ALL NIGHBORS OF CUR CELL: ===
				n_neighs = c2c_ptr(curr+1) - c2c_ptr(curr)
				if (n_neighs > 0) then
					allocate(neighs(n_neighs))
					neighs = c2c_neighs(c2c_ptr(curr):c2c_ptr(curr+1)-1)
					
					!=== BUBBLE SORT NEIGHBORS BY DEGREE: ===
					do j = 1, n_neighs - 1
						do k = 1, n_neighs - j
							if (degree(neighs(k)) > degree(neighs(k+1))) then
								tmp_val = neighs(k)
								neighs(k) = neighs(k+1)
								neighs(k+1) = tmp_val
							end if
						end do
					end do
					
					do j = 1, n_neighs
						neighbor = neighs(j)
						if (.not. visited(neighbor)) then
							visited(neighbor) = .true.
							tail = tail + 1
							queue(tail) = neighbor
						end if
					end do
					deallocate(neighs)
				end if
			end do
		end do

		!=== REVERSING: ===
		do cell_idx = 1, ncells
			queue(cell_idx) = new_to_old(ncells - cell_idx + 1)
		end do
		new_to_old = queue

		!=== OLD TO NEW MAP: ===
		do cell_idx = 1, ncells
			old_to_new(new_to_old(cell_idx)) = cell_idx
		end do

		!=== UPDATING DATA: ===
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if (left_cell > 0 .and. left_cell <= ncells)  face_left_cell(face_idx) = old_to_new(left_cell)
			if (right_cell > 0 .and. right_cell <= ncells) face_right_cell(face_idx) = old_to_new(right_cell)
		end do
		
		!=== MAKING SURE THAT LEFT_CELLS .NE. 0: ===
		call reorder_left_right_cells(nfaces, ncells,&
									  face_left_cell, face_right_cell,&
									  face_nodes_ptr, face_nodes)
			
		do cell_idx = 1, ncells
			queue(cell_idx) = cell_type(new_to_old(cell_idx))
		end do
		cell_type = queue

		deallocate(queue, visited, degree, new_to_old, old_to_new, c2c_ptr, c2c_neighs)
	end subroutine

!=======================================================================
!==================== SUBROUTINE FOR FACES REORDERING ==================
!=======================================================================
	subroutine reorder_faces(ncells, nfaces, nbfaces, face_left_cell, face_right_cell,&
							 face_zone, face_type, face_nodes_ptr, face_nodes)
		implicit none
		integer, intent(in) :: ncells, nfaces
		integer, intent(inout) :: nbfaces
		integer, allocatable, intent(inout) :: face_left_cell(:), face_right_cell(:),&
											   face_zone(:), face_type(:),&
											   face_nodes_ptr(:), face_nodes(:)
		
		integer :: face_idx
		integer :: left_cell, right_cell
		integer :: n_int, n_ext, idx, pos1, pos2, pos1_old, pos2_old
		integer, allocatable :: old_to_new(:), new_to_old(:),& 
								new_left(:), new_right(:),&
								new_zone(:), new_type(:),&
								new_face_nodes(:), new_face_nodes_ptr(:)

		allocate(old_to_new(nfaces), new_to_old(nfaces))
		n_int = 0; n_ext = 0
		
		!=== FINDING NUMBER OF INTERNAL & EXTERNAL FACES: ===
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if ((left_cell > 0 .and. left_cell <= ncells) .and.&
				(right_cell > 0 .and. right_cell <= ncells)) then
				n_int = n_int + 1
			else
				n_ext = n_ext + 1
			end if
		end do

		nbfaces = n_ext
		
		!=== REORDERING FACES: [INTERNAL BLOCK, EXTERNAL BLOCK] ===
		n_int = 0
		n_ext = nfaces - nbfaces
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if ((left_cell > 0 .and. left_cell <= ncells) .and.&
				(right_cell > 0 .and. right_cell <= ncells)) then
				n_int = n_int + 1
				old_to_new(face_idx) = n_int
				new_to_old(n_int) = face_idx
			else
				n_ext = n_ext + 1
				old_to_new(face_idx) = n_ext
				new_to_old(n_ext) = face_idx
			end if
		end do
		
		!=== REORDERING FACE DATA: ===
		allocate(new_left(nfaces), new_right(nfaces),&
				 new_zone(nfaces), new_type(nfaces))
		do face_idx = 1, nfaces
			new_left(old_to_new(face_idx)) = face_left_cell(face_idx)
			new_right(old_to_new(face_idx)) = face_right_cell(face_idx)
			new_zone(old_to_new(face_idx)) = face_zone(face_idx)
			new_type(old_to_new(face_idx)) = face_type(face_idx)
		end do
		call move_alloc(new_left, face_left_cell)
		call move_alloc(new_right, face_right_cell)
		call move_alloc(new_zone, face_zone)
		call move_alloc(new_type, face_type)
		
		allocate(new_face_nodes(size(face_nodes)),&
				 new_face_nodes_ptr(nfaces+1))
		
		new_face_nodes_ptr(1) = 1
		do face_idx = 1, nfaces
			new_face_nodes_ptr(face_idx+1) = new_face_nodes_ptr(face_idx) + &
										face_nodes_ptr(new_to_old(face_idx) + 1) - face_nodes_ptr(new_to_old(face_idx))				
		end do
		
		do face_idx = 1, nfaces
			pos1 = new_face_nodes_ptr(face_idx)
			pos2 = new_face_nodes_ptr(face_idx+1) - 1
			
			pos1_old = face_nodes_ptr(new_to_old(face_idx))
			pos2_old = face_nodes_ptr(new_to_old(face_idx) + 1) - 1
			
			new_face_nodes(pos1:pos2) = face_nodes(pos1_old:pos2_old)
		end do
		call move_alloc(new_face_nodes_ptr, face_nodes_ptr)
		call move_alloc(new_face_nodes, face_nodes)
		
		
		deallocate(old_to_new, new_to_old)
	end subroutine
		
!=======================================================================
!=============== SUBROUTINE FOR ADVANCED FACES REORDERING ==============
!=======================================================================
	subroutine reorder_faces_advanced(ncells, nfaces, nbfaces,&
									  face_left_cell, face_right_cell,&
									  face_zone, face_type,&
									  face_nodes_ptr, face_nodes)
		implicit none
		integer, intent(in) :: ncells, nfaces
		integer, intent(inout) :: nbfaces
		integer, allocatable, intent(inout) :: face_left_cell(:), face_right_cell(:),&
											   face_zone(:), face_type(:),&
											   face_nodes_ptr(:), face_nodes(:)
		
		integer :: face_idx, cell_idx, num_nodes
		integer :: left_cell, right_cell, target_cell
		integer :: j, f, n_int, n_ext
		integer :: idx, pos, pos1, pos2, pos1_old, pos2_old
		integer, allocatable :: face_counts(:), offsets(:), old_to_new(:), new_to_old(:)
		integer, allocatable :: new_left(:), new_right(:),&
								new_zone(:), new_type(:),&
								new_face_nodes(:), new_face_nodes_ptr(:)


		allocate(old_to_new(nfaces), new_to_old(nfaces))
		allocate(face_counts(ncells), source=0)

		!=== FINDING ALL INTERNAL & EXTERNAL FACES: ===
		n_int = 0
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			if ((left_cell > 0 .and. left_cell <= ncells) .and.&
				(right_cell > 0 .and. right_cell <= ncells)) then
				n_int = n_int + 1
				face_counts(min(left_cell, right_cell)) = face_counts(min(left_cell, right_cell)) + 1
			end if
		end do
		nbfaces = nfaces - n_int

		!=== CREATING OFFSETS IN INTERNAL FACES BLOCK: ===
		allocate(offsets(ncells + 1))
		offsets(1) = 1
		do cell_idx = 1, ncells
			offsets(cell_idx+1) = offsets(cell_idx) + face_counts(cell_idx)
		end do

		!=== FILLIN OLD<->TO<->NEW MAP: ===
		n_ext = n_int + 1
		face_counts = 0
		
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if ((left_cell > 0 .and. left_cell <= ncells) .and.&
				(right_cell > 0 .and. right_cell <= ncells)) then
				
																		!=== WARNING!!!! ===
				!if (left_cell > right_cell) then
				!	call swap_int(face_left_cell(face_idx), face_right_cell(face_idx))
				!	pos1 = face_nodes_ptr(face_idx)
				!	pos2 = face_nodes_ptr(face_idx+1) - 1
				!	num_nodes = pos2 - pos1 + 1
				!	do j = 1, num_nodes/2
				!		call swap_int(face_nodes(pos1 + j - 1), face_nodes(pos2 - j + 1))
				!	end do
				!end if
			
				target_cell = min(left_cell, right_cell)
				pos = offsets(target_cell) + face_counts(target_cell)
				old_to_new(face_idx) = pos
				new_to_old(pos) = face_idx
				face_counts(target_cell) = face_counts(target_cell) + 1
			else
				old_to_new(face_idx) = n_ext
				new_to_old(n_ext) = face_idx
				n_ext = n_ext + 1
			end if
		end do
		
		
		!=== REORDERING DATA: ===
		allocate(new_left(nfaces), new_right(nfaces),&
				 new_zone(nfaces), new_type(nfaces))
		do face_idx = 1, nfaces
			new_left(old_to_new(face_idx)) = face_left_cell(face_idx)
			new_right(old_to_new(face_idx)) = face_right_cell(face_idx)
			new_zone(old_to_new(face_idx)) = face_zone(face_idx)
			new_type(old_to_new(face_idx)) = face_type(face_idx)
		end do
		call move_alloc(new_left, face_left_cell)
		call move_alloc(new_right, face_right_cell)
		call move_alloc(new_zone, face_zone)
		call move_alloc(new_type, face_type)
		
		allocate(new_face_nodes(size(face_nodes)),&
				 new_face_nodes_ptr(nfaces+1))
		
		new_face_nodes_ptr(1) = 1
		do face_idx = 1, nfaces
			new_face_nodes_ptr(face_idx+1) = new_face_nodes_ptr(face_idx) + &
										face_nodes_ptr(new_to_old(face_idx) + 1) - face_nodes_ptr(new_to_old(face_idx))				
		end do
		
		do face_idx = 1, nfaces
			pos1 = new_face_nodes_ptr(face_idx)
			pos2 = new_face_nodes_ptr(face_idx+1) - 1
			
			pos1_old = face_nodes_ptr(new_to_old(face_idx))
			pos2_old = face_nodes_ptr(new_to_old(face_idx) + 1) - 1
			
			new_face_nodes(pos1:pos2) = face_nodes(pos1_old:pos2_old)
		end do
		call move_alloc(new_face_nodes_ptr, face_nodes_ptr)
		call move_alloc(new_face_nodes, face_nodes)
		
		
		deallocate(old_to_new, new_to_old, face_counts, offsets)
	
	end subroutine
		
!=======================================================================
!===================== NODES REORDERING SUBROUTINES ====================
!=======================================================================
	subroutine order_cell_nodes_2d(ncells, node_coords, cell_nodes_ptr, cell_nodes)
		implicit none
		real(8), parameter :: eps = 1.0d-12
		integer, intent(in) :: ncells
		real(8), intent(in) :: node_coords(:, :)
		integer, intent(in) :: cell_nodes_ptr(:)
		integer, intent(inout) :: cell_nodes(:)
		
		integer :: cell_idx
		integer :: pos1, pos2
		integer :: ordered_nodes(40)
		integer :: i, j, num_nodes, node_idx
		
		real(8) :: angles(40)
		real(8) :: center(2)
		real(8) :: dx, dy 
		
			
		do cell_idx = 1, ncells
			pos1 = cell_nodes_ptr(cell_idx)
			pos2 = cell_nodes_ptr(cell_idx+1) - 1
			
			num_nodes = pos2 - pos1 + 1
			
			!== GEOMETRIC CENTER: ==
			center = 0.0d0
			do i = 1, num_nodes
				node_idx = cell_nodes(pos1 + i - 1) 
				center = center + node_coords(:, node_idx)
				ordered_nodes(i) = cell_nodes(pos1 + i - 1)
			end do
			center = center/num_nodes
			
			do i = 1, num_nodes
				node_idx = cell_nodes(pos1 + i - 1) 
				dx = node_coords(1, node_idx) - center(1)
				dy = node_coords(2, node_idx) - center(2)
				angles(i) = atan2(dy, dx)
				if (angles(i) < 0.0d0) angles(i) = angles(i) + 2.0d0*3.141592653589793d0
			end do
				
			do i = 1, num_nodes - 1
				do j = i + 1, num_nodes
					if (angles(i) > angles(j)) then
						call swap_real(angles(i), angles(j))
						call swap_int(ordered_nodes(i), ordered_nodes(j))
					end if
				end do
			end do
			
			cell_nodes(pos1:pos2) = ordered_nodes(1:num_nodes)
		end do
	end subroutine
		
	subroutine order_face_nodes_3d(nfaces, node_coords, face_nodes_ptr, face_nodes)
		implicit none
		real(8), parameter :: eps = 1.0d-12
		integer, intent(in) :: nfaces
		real(8), intent(in) :: node_coords(:, :)
		integer, intent(in) :: face_nodes_ptr(:)
		integer, intent(inout) :: face_nodes(:)
		
		integer :: face_idx
		integer :: pos1, pos2
		integer :: ordered_nodes(40)
		integer :: i, j, num_nodes, node_idx, next_idx
		
		real(8) :: angles(40)
		real(8) :: r_i(3), r_j(3), r_vec(3), center(3), normal(3), basis_x(3), basis_y(3)
		real(8) :: x, y 
		
		
		do face_idx = 1, nfaces
			pos1 = face_nodes_ptr(face_idx)
			pos2 = face_nodes_ptr(face_idx+1) - 1
			
			num_nodes = pos2 - pos1 + 1
			
			!== GEOMETRIC CENTER: ==
			center = 0.0d0
			do i = 1, num_nodes
				node_idx = face_nodes(pos1 + i - 1) 
				center = center + node_coords(:, node_idx)
			end do
			center = center/num_nodes
			
			!== UNIT NORMAL VECTOR: ==
			normal = 0.0d0
			do i = 1, num_nodes
				node_idx = face_nodes(pos1 + i - 1) 
				if (i == num_nodes) then
					next_idx = face_nodes(pos1)
				else
					next_idx = face_nodes(pos1 + i)
				end if
				
				r_i = node_coords(:, node_idx)
				r_j = node_coords(:, next_idx)
				
				normal(1) = normal(1) + (r_i(2) - r_j(2))*(r_i(3) + r_j(3))
				normal(2) = normal(2) + (r_i(3) - r_j(3))*(r_i(1) + r_j(1))
				normal(3) = normal(3) + (r_i(1) - r_j(1))*(r_i(2) + r_j(2))
			end do
			normal = normal/norm2(normal)
			
			!== LOCAL 2D BASIS: ==
			basis_x = node_coords(:, face_nodes(pos1)) - center
			basis_x = basis_x - dot_product(basis_x, normal)*normal
			basis_x = basis_x/norm2(basis_x)
			basis_y = cross_product_3d(basis_x, normal)
			basis_y = basis_y/sqrt(sum(basis_y**2))
			
			do i = 1, num_nodes
				node_idx = face_nodes(pos1 + i - 1)
				
				r_vec = node_coords(:, node_idx) - center
				
				x = dot_product(r_vec, basis_x)
				y = dot_product(r_vec, basis_y)
				
				angles(i) = atan2(y, x)
				if (angles(i) < 0.0d0) angles(i) = angles(i) + 2.0d0 * 3.141592653589793d0
				
				ordered_nodes(i) = face_nodes(pos1 + i - 1)
			end do
				
			do i = 1, num_nodes - 1
				do j = i + 1, num_nodes
					if (angles(i) > angles(j)) then
						call swap_real(angles(i), angles(j))
						call swap_int(ordered_nodes(i), ordered_nodes(j))
					end if
				end do
			end do
			
			
			face_nodes(pos1:pos2) = ordered_nodes(1:num_nodes)
		end do
	end subroutine
		
!=======================================================================
!================ BOUNDARY FACES ORDERING SUBROUTINE ===================
!=======================================================================
	subroutine order_boundary_faces_by_zone(dim, ncells, nfaces, nbfaces,&
                                            face_left_cell, face_right_cell,&
                                            face_zone, face_type, face_bidx,&
                                            face_nodes_ptr, face_nodes,&
                                            face_normal, face_area, face_center)
        implicit none
        integer, intent(in) :: dim, ncells, nfaces, nbfaces
        integer, intent(inout) :: face_left_cell(:), face_right_cell(:)
        integer, intent(inout) :: face_zone(:), face_type(:)
        integer, intent(inout), optional :: face_bidx(:)
        integer, intent(inout) :: face_nodes_ptr(:), face_nodes(:)
        real(8), optional, intent(inout) :: face_normal(:,:), face_area(:), face_center(:,:)

        integer :: i, j, b_start, n_nodes, pos_old, pos_new
        integer, allocatable :: perm(:), temp_int(:), new_nodes(:)
        real(8), allocatable :: temp_real_2d(:,:), temp_real_1d(:)
        integer :: zone_i, zone_j, tmp
        
        b_start = nfaces - nbfaces + 1
        allocate(perm(nbfaces))
        do i = 1, nbfaces
            perm(i) = b_start + i - 1
        end do

        !=== BUBLE SORT: ===
        do i = 1, nbfaces - 1
            do j = i + 1, nbfaces
                if (face_zone(perm(j)) < face_zone(perm(i))) then
                    tmp = perm(i)
                    perm(i) = perm(j)
                    perm(j) = tmp
                end if
            end do
        end do

        !=== REORDER DATA: ===
        allocate(temp_int(nbfaces))
        allocate(temp_real_1d(nbfaces))
        allocate(temp_real_2d(dim, nbfaces))

        !=== REODER LEFT/RIGHT CELLS: ===
        do i = 1, nbfaces
            temp_int(i) = face_left_cell(perm(i))
        end do
        face_left_cell(b_start:nfaces) = temp_int
        
        do i = 1, nbfaces
            temp_int(i) = face_right_cell(perm(i))
        end do
        face_right_cell(b_start:nfaces) = temp_int

        !=== REODER FACE ZONES/TYPES: ===
        do i = 1, nbfaces
            temp_int(i) = face_zone(perm(i))
        end do
        face_zone(b_start:nfaces) = temp_int

        do i = 1, nbfaces
            temp_int(i) = face_type(perm(i))
        end do
        face_type(b_start:nfaces) = temp_int

        !=== REORDER METRIC: ===
        if (present(face_normal)) then
			do i = 1, nbfaces
				temp_real_2d(:, i) = face_normal(:, perm(i))
			end do
			face_normal(:, b_start:nfaces) = temp_real_2d
		end if
		
		if (present(face_area)) then
			do i = 1, nbfaces
				temp_real_1d(i) = face_area(perm(i))
			end do
			face_area(b_start:nfaces) = temp_real_1d
		end if

		if (present(face_center)) then
			do i = 1, nbfaces
				temp_real_2d(:, i) = face_center(:, perm(i))
			end do
			face_center(:, b_start:nfaces) = temp_real_2d
		end if

		!=== REORDER NODES: ===
        n_nodes = size(face_nodes)
        allocate(new_nodes(n_nodes))
        new_nodes(1:face_nodes_ptr(b_start)-1) = face_nodes(1:face_nodes_ptr(b_start)-1)
        
        !=== BOUNDARY FACES NODES: ===
        pos_new = face_nodes_ptr(b_start)
        do i = 1, nbfaces
            pos_old = face_nodes_ptr(perm(i))
            n_nodes = face_nodes_ptr(perm(i)+1) - pos_old
            new_nodes(pos_new:pos_new + n_nodes - 1) = face_nodes(pos_old:pos_old + n_nodes - 1)
            
            temp_int(i) = n_nodes
            pos_new = pos_new + n_nodes
        end do
        
        do i = 1, nbfaces
            face_nodes_ptr(b_start + i) = face_nodes_ptr(b_start + i - 1) + temp_int(i)
        end do
        face_nodes = new_nodes

        !=== UPDATE FACE_BIDX/GHOST CELL IDX: ===
        if (present(face_bidx)) then
			do i = 1, nbfaces
				face_bidx(b_start + i - 1) = i
				face_right_cell(b_start + i - 1) = ncells + i
			end do
		end if

        deallocate(perm, temp_int, temp_real_1d, temp_real_2d, new_nodes)
    end subroutine


end module
