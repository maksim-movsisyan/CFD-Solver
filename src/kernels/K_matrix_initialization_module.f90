module matrix_initialization_module
use array_sorting_module
implicit none
contains
!=======================================================================
!================ SPARCE 1ST ORDER MATRIX INITIALIZATION ===============
!=======================================================================
	pure subroutine spr_mtrx_initialize1(ncells, nfaces, nbfaces,&
										 cell_faces_ptr, cell_faces,&
										 face_left_cell, face_right_cell,&
										 n, nnz, col_indices, row_ptr, diag_indices)
		implicit none
		integer, intent(in) :: ncells, nfaces, nbfaces
		integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
										   face_left_cell(:), face_right_cell(:)
		integer, intent(inout) :: n, nnz
		integer, allocatable, intent(inout) :: col_indices(:), row_ptr(:), diag_indices(:)
		
		integer :: cell_id, neighbor_id, face_id
		integer :: left_cell, right_cell
		integer :: pos1, pos2, current_pos
		integer :: i, j, k, row, col
		integer, allocatable :: neighbor_count(:), nnz_per_row(:)
		logical, allocatable :: visited(:)
		integer :: visited_idx(100)
		integer :: counter
		integer, allocatable :: temp_col_indices(:)
		integer :: row_size, max_row_size
		
		allocate(neighbor_count(ncells), visited(ncells))
		
		!=== COUNING NUMBER OF NEIGHBORS FOR EACH CELL: ===
		visited = .false.
		do cell_id = 1, ncells
			visited(cell_id) = .true.
			neighbor_count(cell_id) = 0
			
			counter = 1
			visited_idx(counter) = cell_id
			
			pos1 = cell_faces_ptr(cell_id)
			pos2 = cell_faces_ptr(cell_id+1) - 1
			
			do i = pos1, pos2
				face_id = cell_faces(i)
				left_cell = face_left_cell(face_id)
				right_cell = face_right_cell(face_id)
				
				if (left_cell == cell_id) then
					neighbor_id = right_cell
				else
					neighbor_id = left_cell
				end if
				
				if (neighbor_id > ncells) cycle
				
				if (.not. visited(neighbor_id)) then
					visited(neighbor_id) = .true.
					neighbor_count(cell_id) = neighbor_count(cell_id) + 1
					
					counter = counter + 1
					visited_idx(counter) = neighbor_id
				end if
			end do
			
			
			do i = 1, counter
				visited(visited_idx(i)) = .false.
			end do
				
		end do
		
		!=== CSR MATRIX INITIALIZATION: ===
		n = ncells
		allocate(nnz_per_row(n))
		nnz_per_row(:) = 0
		
		do cell_id = 1, ncells
			nnz_per_row(cell_id) = 1 + neighbor_count(cell_id)
		end do
		
		nnz = sum(nnz_per_row)
		allocate(col_indices(nnz),&
				 diag_indices(n), row_ptr(n + 1))
				 
		col_indices = 0
		diag_indices = 0
		row_ptr = 0
		
		row_ptr(1) = 1
		do i = 2, n + 1
			row_ptr(i) = row_ptr(i - 1) + nnz_per_row(i - 1)
		end do
		
		!=== INDICES ARRAYS FILLING: ===
		visited = .false.
		do cell_id = 1, ncells  
			row = cell_id
			col = cell_id
			current_pos = row_ptr(row)
			
			col_indices(current_pos) = col
			current_pos = current_pos + 1
			
			visited(cell_id) = .true.
			
			counter = 1
			visited_idx(counter) = cell_id
			
			pos1 = cell_faces_ptr(cell_id)
			pos2 = cell_faces_ptr(cell_id+1) - 1
				
			do i = pos1, pos2
				face_id = cell_faces(i)
				left_cell = face_left_cell(face_id)
				right_cell = face_right_cell(face_id)
				
				if (left_cell == cell_id) then
					neighbor_id = right_cell
				else
					neighbor_id = left_cell
				end if
				
				if (neighbor_id > ncells) cycle
				
				if (.not. visited(neighbor_id)) then
					visited(neighbor_id) = .true.
					col = neighbor_id
					col_indices(current_pos) = col
					current_pos = current_pos + 1
					
					counter = counter + 1
					visited_idx(counter) = neighbor_id
				end if
			end do
			
			do i = 1, counter
				visited(visited_idx(i)) = .false.
			end do
			
		end do
		deallocate(neighbor_count, nnz_per_row, visited)
		
		!=== SORTING COL ARRAYS: ===
		max_row_size = 0
		do i = 1, n
			row_size = row_ptr(i + 1) - row_ptr(i)
			max_row_size = max(max_row_size, row_size)
		end do
		
		allocate(temp_col_indices(max_row_size))
		
		do i = 1, n
			pos1 = row_ptr(i)
			pos2 = row_ptr(i + 1) - 1
			
			row_size = pos2 - pos1 + 1
			
			temp_col_indices(1:row_size) = col_indices(pos1:pos2)
			call hybrid_sort(temp_col_indices(1:row_size), row_size)
			col_indices(pos1:pos2) = temp_col_indices(1:row_size)
		end do
		deallocate(temp_col_indices)
		
		!=== FINDING DIAG INDICES: ===
		do i = 1, n
			pos1 = row_ptr(i)
			pos2 = row_ptr(i + 1) - 1
			
			do j = pos1, pos2
				if (col_indices(j) == i) then
					diag_indices(i) = j
					exit
				end if
			end do
		end do
	
	end subroutine

!=======================================================================
!================== CSR MATRIX INITIALIZATION ==========================
!=======================================================================
	pure subroutine csr_mtrx_initialize1(ncells, nfaces, nbfaces,&
										 cell_faces_ptr, cell_faces,&
										 face_left_cell, face_right_cell,&
										 n, nnz, values, col_indices, row_ptr, diag_indices)
		implicit none
		integer, intent(in) :: ncells, nfaces, nbfaces
		integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
										   face_left_cell(:), face_right_cell(:)
		integer, intent(inout) :: n, nnz
		real(8), allocatable, intent(inout) :: values(:)
		integer, allocatable, intent(inout) :: col_indices(:), row_ptr(:), diag_indices(:)
		
		call spr_mtrx_initialize1(ncells, nfaces, nbfaces,&
								  cell_faces_ptr, cell_faces,&
								  face_left_cell, face_right_cell,&
								  n, nnz, col_indices, row_ptr, diag_indices)
		allocate(values(nnz))
		values = 0.d0
	end subroutine

	pure subroutine bcsr_mtrx_initialize1(ncells, nfaces, nbfaces, bs,&
										  cell_faces_ptr, cell_faces,&
										  face_left_cell, face_right_cell,&
										  n, nnz, values, col_indices, row_ptr, diag_indices)
		implicit none
		integer, intent(in) :: ncells, nfaces, nbfaces, bs
		integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
										   face_left_cell(:), face_right_cell(:)
		integer, intent(inout) :: n, nnz
		real(8), allocatable, intent(inout) :: values(:, :, :)
		integer, allocatable, intent(inout) :: col_indices(:), row_ptr(:), diag_indices(:)
		
		call spr_mtrx_initialize1(ncells, nfaces, nbfaces,&
								  cell_faces_ptr, cell_faces,&
								  face_left_cell, face_right_cell,&
								  n, nnz, col_indices, row_ptr, diag_indices)
		allocate(values(bs, bs, nnz))
		values = 0.d0
	end subroutine
end module
