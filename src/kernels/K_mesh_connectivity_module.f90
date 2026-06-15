module mesh_connectivity_module
implicit none
contains
!=======================================================================
!================ FACE-NODE CONNECTIVITY SUBROUTINE ====================
!=======================================================================
	subroutine face_node_connectivity(nfaces, nnodes_per_face, face_nodes_ptr, face_nodes)
		implicit none
		integer, intent(in) :: nfaces
		integer, intent(in) :: nnodes_per_face(:)
		integer, allocatable, intent(inout) :: face_nodes_ptr(:), face_nodes(:)
		
		integer :: face_idx
		integer :: total_nodes
		
		total_nodes = sum(nnodes_per_face)
		if (allocated(face_nodes_ptr)) deallocate(face_nodes_ptr)
		if (allocated(face_nodes)) deallocate(face_nodes)
		
		allocate(face_nodes_ptr(nfaces+1))
		allocate(face_nodes(total_nodes))
				
		face_nodes_ptr(1) = 1
		do face_idx = 1, nfaces
			face_nodes_ptr(face_idx+1) = face_nodes_ptr(face_idx) + nnodes_per_face(face_idx)
		end do
	end subroutine
		
!=======================================================================
!============ CELL-CELL(NEIGHBORS) CONNECTIVITY SUBROUTINE =============
!=======================================================================	
	subroutine cell_neighbor_connectivity(ncells, nfaces,&
										  face_left_cell, face_right_cell,&
										  cell_neighbors_ptr, cell_neighbors)
		implicit none
		integer, intent(in) :: ncells, nfaces
		integer, intent(in) :: face_left_cell(:), face_right_cell(:)   
		integer, allocatable, intent(out) :: cell_neighbors_ptr(:), cell_neighbors(:)
		
		integer :: face_idx, cell_idx
		integer :: left_cell, right_cell
		integer :: curr_pos
		integer, allocatable :: count(:)

		
		if (allocated(cell_neighbors_ptr)) deallocate(cell_neighbors_ptr)
		allocate(cell_neighbors_ptr(ncells+1))
		
		allocate(count(ncells))
		
		cell_neighbors_ptr = 0
		count = 0

		!=== FINDING NUMBER OF NEIGHBORS FOR EACH CELL: ===
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			if ((left_cell > 0 .and. left_cell <= ncells) .and.&
				(right_cell > 0 .and. right_cell <= ncells)) then
				
				cell_neighbors_ptr(left_cell + 1) = cell_neighbors_ptr(left_cell + 1) + 1
				cell_neighbors_ptr(right_cell + 1) = cell_neighbors_ptr(right_cell + 1) + 1
			end if
		end do

		!=== CEATING CELL-NEIGHBORS PTR: ===
		cell_neighbors_ptr(1) = 1
		do cell_idx = 1, ncells
			cell_neighbors_ptr(cell_idx+1) = cell_neighbors_ptr(cell_idx) + cell_neighbors_ptr(cell_idx+1)
		end do
		
		
		!=== CREATING CELL-NEIGHBORS CONNECTIVITY: ===
		if (allocated(cell_neighbors)) deallocate(cell_neighbors)
		allocate(cell_neighbors(cell_neighbors_ptr(ncells+1) - 1))
		
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if ((left_cell > 0 .and. left_cell <= ncells) .and.&
				(right_cell > 0 .and. right_cell <= ncells)) then
				!=== RIGHT AS NEIGHBOR OF LEFT: ===
				curr_pos = cell_neighbors_ptr(left_cell) + count(left_cell)
				cell_neighbors(curr_pos) = right_cell
				count(left_cell) = count(left_cell) + 1
				
				!=== LEFT AS NEIGHBOR OF RIGHT: ===
				curr_pos = cell_neighbors_ptr(right_cell) + count(right_cell)
				cell_neighbors(curr_pos) = left_cell
				count(right_cell) = count(right_cell) + 1
			end if
		end do

		deallocate(count)
	end subroutine
	
	subroutine cell_neighbor_connectivity_advanced(ncells, nfaces,&
												   face_left_cell, face_right_cell,&
												   cell_faces_ptr, cell_faces,&
												   cell_neighbors_ptr, cell_neighbors)
		implicit none
		integer, intent(in) :: ncells, nfaces
		integer, intent(in) :: face_left_cell(:), face_right_cell(:)   
		integer, intent(in) :: cell_faces_ptr(:), cell_faces(:)
		integer, allocatable, intent(out) :: cell_neighbors_ptr(:), cell_neighbors(:)
		
		integer :: face_idx, cell_idx
		integer :: left_cell, right_cell, neighbor_cell
		integer :: curr_pos, nneights
		integer :: pos1, pos2, f_ptr, i
		integer :: cell_neighbor_idx(200)
		logical, allocatable :: visited(:)

		
		if (allocated(cell_neighbors_ptr)) deallocate(cell_neighbors_ptr)
		allocate(cell_neighbors_ptr(ncells+1))
		
		allocate(visited(ncells))
		
		cell_neighbors_ptr = 0
		visited = .false.
		cell_neighbor_idx = -1
		
		!=== ADVANCED FINDING NUMBER OF NEIGHBORS FOR EACH CELL: ===
		do cell_idx = 1, ncells
			nneights = 0
			
			pos1 = cell_faces_ptr(cell_idx)
			pos2 = cell_faces_ptr(cell_idx+1)-1
			
			do f_ptr = pos1, pos2
				face_idx = cell_faces(f_ptr)
				left_cell = face_left_cell(face_idx)
				right_cell = face_right_cell(face_idx) 
				
				neighbor_cell = left_cell
				if (left_cell == cell_idx) neighbor_cell = right_cell
				
				if (neighbor_cell <= 0 .or. neighbor_cell > ncells) cycle
				
				if (.not. visited(neighbor_cell)) then
					nneights = nneights + 1
					visited(neighbor_cell) = .true.
					cell_neighbor_idx(nneights) = neighbor_cell
				end if
			end do
			
			cell_neighbors_ptr(cell_idx + 1) = nneights
			
			do i = 1, nneights
				visited(cell_neighbor_idx(i)) = .false.
			end do
		end do
	

		!=== CEATING CELL-NEIGHBORS PTR: ===
		cell_neighbors_ptr(1) = 1
		do cell_idx = 1, ncells
			cell_neighbors_ptr(cell_idx+1) = cell_neighbors_ptr(cell_idx) + cell_neighbors_ptr(cell_idx+1)
		end do
		
		
		!=== CREATING CELL-NEIGHBORS CONNECTIVITY: ===
		if (allocated(cell_neighbors)) deallocate(cell_neighbors)
		allocate(cell_neighbors(cell_neighbors_ptr(ncells+1) - 1))
		
		do cell_idx = 1, ncells
			nneights = 0
			
			pos1 = cell_faces_ptr(cell_idx)
			pos2 = cell_faces_ptr(cell_idx+1)-1
			
			do f_ptr = pos1, pos2
				face_idx = cell_faces(f_ptr)
				left_cell = face_left_cell(face_idx)
				right_cell = face_right_cell(face_idx) 
				
				neighbor_cell = left_cell
				if (left_cell == cell_idx) neighbor_cell = right_cell
				
				if (neighbor_cell <= 0 .or. neighbor_cell > ncells) cycle
				
				if (.not. visited(neighbor_cell)) then
					nneights = nneights + 1
					visited(neighbor_cell) = .true.
					cell_neighbor_idx(nneights) = neighbor_cell
				end if
			end do
			
			pos1 = cell_neighbors_ptr(cell_idx) - 1
			
			do i = 1, nneights
				cell_neighbors(pos1+i) = cell_neighbor_idx(i)
				visited(cell_neighbor_idx(i)) = .false.
			end do
		end do
		
		deallocate(visited)
	end subroutine
			
!=======================================================================
!================= CELL-FACE CONNECTIVITY SUBROUTINE ===================
!=======================================================================
	subroutine cell_face_connectivity(ncells, nfaces,&
									  face_left_cell, face_right_cell,&
									  cell_faces_ptr, cell_faces)
		implicit none
		integer, intent(in) :: ncells, nfaces
		integer, intent(in) :: face_left_cell(:), face_right_cell(:)
		integer, allocatable, intent(inout) :: cell_faces_ptr(:), cell_faces(:)
		
		integer, allocatable :: nfaces_per_cell(:)
		integer :: face_idx, cell_idx, idx
		integer :: left_cell, right_cell
		integer :: cell_faces_num
		
		!=== FINDING NUMBER OF FACES PER EACH CELL: ===
		allocate(nfaces_per_cell(ncells))
		nfaces_per_cell = 0
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if (left_cell > 0 .and. left_cell <= ncells) then
				nfaces_per_cell(left_cell) = nfaces_per_cell(left_cell) + 1
			end if
			
			if (right_cell > 0 .and. right_cell <= ncells) then
				nfaces_per_cell(right_cell) = nfaces_per_cell(right_cell) + 1
			end if
		end do
		
		cell_faces_num = sum(nfaces_per_cell)
		if (allocated(cell_faces_ptr)) deallocate(cell_faces_ptr)
		if (allocated(cell_faces)) deallocate(cell_faces)
		allocate(cell_faces_ptr(ncells+1))
		allocate(cell_faces(cell_faces_num))
		
		!=== CEATING CELL-FACES PTR: ===
		cell_faces_ptr(1) = 1
		do cell_idx = 1, ncells
			cell_faces_ptr(cell_idx+1) = cell_faces_ptr(cell_idx) + nfaces_per_cell(cell_idx)
		end do
		
		!=== CREATING CELL-FACES CONNECTIVITY: ===
		nfaces_per_cell = 0
		do face_idx = 1, nfaces
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if (left_cell > 0 .and. left_cell <= ncells) then
				idx = cell_faces_ptr(left_cell) + nfaces_per_cell(left_cell)
				cell_faces(idx) = face_idx
				nfaces_per_cell(left_cell) = nfaces_per_cell(left_cell) + 1
			end if
			
			if (right_cell > 0 .and. right_cell <= ncells) then
				idx = cell_faces_ptr(right_cell) + nfaces_per_cell(right_cell)
				cell_faces(idx) = face_idx
				nfaces_per_cell(right_cell) = nfaces_per_cell(right_cell) + 1
			end if
		end do		
		deallocate(nfaces_per_cell)
	end subroutine
		
!=======================================================================
!================= CELL-NODE CONNECTIVITY SUBROUTINE ===================
!=======================================================================
	subroutine cell_node_connectivity(nnodes, ncells,&
									  face_nodes_ptr, face_nodes,&
									  cell_faces_ptr, cell_faces,&
									  cell_nodes_ptr, cell_nodes)
		implicit none
		integer, intent(in) :: nnodes, ncells
		integer, intent(in) :: face_nodes_ptr(:), face_nodes(:),&
							   cell_faces_ptr(:), cell_faces(:)
		integer, allocatable, intent(inout) :: cell_nodes_ptr(:), cell_nodes(:)
		
		integer :: cell_idx, face_idx, node_idx
		integer :: f, n, count
		integer, allocatable :: marker(:)
		
		if (allocated(cell_nodes_ptr)) deallocate(cell_nodes_ptr)
	    allocate(cell_nodes_ptr(ncells+1))
		
		!=== CEATING CELL-NODES PTR: ===
		allocate(marker(nnodes)); marker = 0
		
		cell_nodes_ptr(1) = 1
		do cell_idx = 1, ncells
			count = 0
			do f = cell_faces_ptr(cell_idx), cell_faces_ptr(cell_idx+1) - 1
				face_idx = cell_faces(f)
				do n = face_nodes_ptr(face_idx), face_nodes_ptr(face_idx+1) - 1
					node_idx = face_nodes(n)
					if (marker(node_idx) /= cell_idx) then
						count = count + 1
						marker(node_idx) = cell_idx
					end if
				end do
			end do
			cell_nodes_ptr(cell_idx+1) = cell_nodes_ptr(cell_idx) + count
		end do
		
		!=== CREATING CELL-NODES CONNECTIVITY: ===
		if (allocated(cell_nodes)) deallocate(cell_nodes)
		allocate(cell_nodes(cell_nodes_ptr(ncells + 1) - 1))
		marker = 0
		
		do cell_idx = 1, ncells
			count = cell_nodes_ptr(cell_idx)
			do f = cell_faces_ptr(cell_idx), cell_faces_ptr(cell_idx+1) - 1
				face_idx = cell_faces(f)
				do n = face_nodes_ptr(face_idx), face_nodes_ptr(face_idx+1) - 1
					node_idx = face_nodes(n)
					if (marker(node_idx) /= cell_idx) then
						cell_nodes(count) = node_idx
						count = count + 1
						marker(node_idx) = cell_idx
					end if
				end do
			end do
		end do
		
		deallocate(marker)
	end subroutine

!=======================================================================
!================= NODE-CELL CONNECTIVITY SUBROUTINE ===================
!=======================================================================
	subroutine node_cell_connectivity(ncells, nnodes,&
									  cell_nodes_ptr, cell_nodes,&
									  node_cells_ptr, node_cells)
		implicit none
		integer, intent(in) :: ncells, nnodes
		integer, intent(in) :: cell_nodes_ptr(:), cell_nodes(:)
		integer, allocatable, intent(inout) :: node_cells_ptr(:), node_cells(:)
		
		integer :: cell_idx, node_idx
		integer :: total_cells
		integer :: i, pos, pos1, pos2
		integer, allocatable :: counter(:)
		
		!=== CEATING NODE-CELLS PTR: ===
		if (allocated(node_cells_ptr)) deallocate(node_cells_ptr)
		allocate(node_cells_ptr(nnodes+1))
		node_cells_ptr = 0
		
		do cell_idx = 1, ncells
			pos1 = cell_nodes_ptr(cell_idx)
			pos2 = cell_nodes_ptr(cell_idx+1)-1
			
			do i = pos1, pos2
				node_idx = cell_nodes(i)
				node_cells_ptr(node_idx + 1) = node_cells_ptr(node_idx + 1) + 1
			end do
		end do
		
		node_cells_ptr(1) = 1
		do node_idx = 1, nnodes
			node_cells_ptr(node_idx + 1) = node_cells_ptr(node_idx) + node_cells_ptr(node_idx + 1)
		end do
		total_cells = node_cells_ptr(nnodes+1)-1
		
		!=== CREATING NODE-CELLS CONNECTIVITY: ===
		if (allocated(node_cells)) deallocate(node_cells)
		allocate(node_cells(total_cells))
		allocate(counter(nnodes), source=0)
		do cell_idx = 1, ncells
			pos1 = cell_nodes_ptr(cell_idx)
			pos2 = cell_nodes_ptr(cell_idx+1)-1
			
			do i = pos1, pos2
				node_idx = cell_nodes(i)
				pos = node_cells_ptr(node_idx) + counter(node_idx)
				node_cells(pos) = cell_idx
				counter(node_idx) = counter(node_idx) + 1
			end do
		end do

		deallocate(counter)
	end subroutine
	
!=======================================================================
!============== BOUNDARY FACES CONNECTIVITY SUBROUTINE =================
!=======================================================================
	subroutine boundary_faces_connectivity(nfaces, nbfaces, face_bidx)
		implicit none
		integer, intent(in) :: nfaces, nbfaces
		integer, allocatable, intent(inout) :: face_bidx(:)
		
		integer :: face_idx, counter
		
		if (allocated(face_bidx)) deallocate(face_bidx)
		allocate(face_bidx(nfaces))
		face_bidx = -1
		
		counter = 0
		do face_idx = nfaces - nbfaces + 1, nfaces
			counter = counter + 1
			face_bidx(face_idx) = counter
		end do
		
	end subroutine
		
		


end module
