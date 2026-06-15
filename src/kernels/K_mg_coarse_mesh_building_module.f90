module mg_coarse_mesh_building_module
use mesh_connectivity_module
use mesh_ordering_module
implicit none
!=======================================================================
!===================== AUXILIARY DATA STRUCTURES =======================
!=======================================================================
type face_candidate
	integer :: left_coarse, right_coarse
	integer :: face_zone, face_type
	real(8) :: center(3) = 0.d0 
	real(8) :: face_area  
	real(8) :: normal(3) = 0.d0
	real(8) :: line_constant
	integer, allocatable :: original_faces(:)              
	integer :: n_original_faces   
	integer, allocatable :: fine_nodes(:)
	integer :: n_fine_nodes
end type


contains

!=======================================================================
!===================== COARSE MESH FACES BUILDING ======================
!=======================================================================
	!face candidates allocation subroutines:
	subroutine collect_face_candidates(dim, fncells, fnfaces, fnbfaces,&
									   fface_left_cell, fface_right_cell,&
									   fface_zone, fface_type,&
								  	   fface_nodes_ptr, fface_nodes,&
									   fface_center, fface_normal, fface_area,&
									   fine_to_coarse,&
									   candidates, num_candidates)
		implicit none
		integer, intent(in) :: dim, fncells, fnfaces, fnbfaces
		integer, intent(in) :: fface_left_cell(:), fface_right_cell(:)
		integer, intent(in) :: fface_zone(:), fface_type(:)
		integer, intent(in) :: fface_nodes_ptr(:), fface_nodes(:)
		integer, intent(in) :: fine_to_coarse(:)
		real(8), intent(in) :: fface_center(:, :), fface_normal(:, :), fface_area(:)
		
		type(face_candidate), allocatable, intent(inout) :: candidates(:)
		integer, intent(inout) :: num_candidates
		
		integer :: fface_idx, i
		integer :: fleft_cell, fright_cell
		integer :: cleft_cell, cright_cell, temp_ccell
		real(8) :: normal(dim)

		!=== FINE MESH FACES LOOP: ===
		num_candidates = 0
		do fface_idx = 1, fnfaces
			fleft_cell = fface_left_cell(fface_idx)
			fright_cell = fface_right_cell(fface_idx)
			
			cleft_cell = 0; cright_cell = 0
			
			if (fleft_cell > 0 .and. fleft_cell <= fncells) then
				cleft_cell = fine_to_coarse(fleft_cell)
			end if
			
			if (fright_cell > 0 .and. fright_cell <= fncells) then
				cright_cell = fine_to_coarse(fright_cell)
			end if			
			
			!=== INTERNAL AGGLOMERATE FACE => MUST BE DELETED ===
			if (cleft_cell == cright_cell) cycle
			
			num_candidates = num_candidates + 1
		end do
		
		if (allocated(candidates)) deallocate(candidates)
		allocate(candidates(num_candidates))
		num_candidates = 0
		
		!=== FINE MESH FACES LOOP: ===
		do fface_idx = 1, fnfaces
			fleft_cell = fface_left_cell(fface_idx)
			fright_cell = fface_right_cell(fface_idx)
			
			cleft_cell = 0; cright_cell = 0
			
			if (fleft_cell > 0 .and. fleft_cell <= fncells) then
				cleft_cell = fine_to_coarse(fleft_cell)
			end if
			
			if (fright_cell > 0 .and. fright_cell <= fncells) then
				cright_cell = fine_to_coarse(fright_cell)
			end if			
			
			!=== INTERNAL AGGLOMERATE FACE => MUST BE DELETED ===
			if (cleft_cell == cright_cell) cycle
			
			normal = fface_normal(:, fface_idx)
			
			if ((cleft_cell > cright_cell .and. cright_cell .ne. 0) .or. cleft_cell == 0) then
				temp_ccell = cleft_cell
				cleft_cell = cright_cell
				cright_cell = temp_ccell
				normal = -normal
			end if
		
			num_candidates = num_candidates + 1
			candidates(num_candidates)%left_coarse = cleft_cell
			candidates(num_candidates)%right_coarse = cright_cell
			candidates(num_candidates)%face_zone = fface_zone(fface_idx)
			candidates(num_candidates)%face_type = fface_type(fface_idx)
			candidates(num_candidates)%center(1:dim) = fface_center(:, fface_idx)
			candidates(num_candidates)%face_area = fface_area(fface_idx)
			candidates(num_candidates)%normal(1:dim) = normal
			
			candidates(num_candidates)%line_constant = dot_product(normal, candidates(num_candidates)%center(1:dim))
			
			
			if (.not. allocated(candidates(num_candidates)%original_faces)) allocate(candidates(num_candidates)%original_faces(1))
			candidates(num_candidates)%original_faces(1) = fface_idx
			candidates(num_candidates)%n_original_faces = 1
			
			candidates(num_candidates)%n_fine_nodes = fface_nodes_ptr(fface_idx+1) - fface_nodes_ptr(fface_idx)
			if (.not. allocated(candidates(num_candidates)%fine_nodes)) &
				allocate(candidates(num_candidates)%fine_nodes(candidates(num_candidates)%n_fine_nodes))
			do i = fface_nodes_ptr(fface_idx), fface_nodes_ptr(fface_idx+1)-1
				candidates(num_candidates)%fine_nodes(i - fface_nodes_ptr(fface_idx) + 1) = fface_nodes(i)
			end do
		end do			
	end subroutine
	
	!aux merging subroutines:
	logical function can_merge_faces(face1, face2, tolerance, dim)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: tolerance
		type(face_candidate), intent(in) :: face1, face2
		
		real(8) :: dot_prct, line_constant1, line_constant2
		
		can_merge_faces = .false.
		
		if (face1%left_coarse /= face2%left_coarse) return
		if (face1%right_coarse /= face2%right_coarse) return
		
		if (face1%face_zone /= face2%face_zone) return
		if (face1%face_type /= face2%face_type) return
		
		dot_prct = dot_product(face1%normal, face2%normal)
		if (dabs(1.d0 - dot_prct) > tolerance) return
		
		line_constant1 = dot_product(face1%normal, face2%center)
		if (dabs(line_constant1 - face1%line_constant) > tolerance) return
		
		line_constant2 =  dot_product(face2%normal, face1%center)
		if (dabs(line_constant2 - face2%line_constant) > tolerance) return
		
		if (dim == 2) then
			if (.not. have_same_node(face1, face2)) return
		else
			if (.not. have_common_edge(face1, face2)) return
		end if
		
		can_merge_faces = .true.
	end function
				
	logical function have_same_node(face1, face2)
		implicit none
		type(face_candidate), intent(in) :: face1, face2
		
		integer :: i
		
		have_same_node = .false.
		
		do i = 1, face1%n_fine_nodes
			if (any(face2%fine_nodes == face1%fine_nodes(i))) then
				have_same_node = .true.
				return
			end if
		end do

	end function	
	
	logical function have_common_edge(face1, face2)
		implicit none
		type(face_candidate), intent(in) :: face1, face2
		
		integer :: i, j, ncommon
		
		ncommon = 0
		do i = 1, face1%n_fine_nodes
			do j = 1, face2%n_fine_nodes
				if (face1%fine_nodes(i) == face2%fine_nodes(j)) then
					ncommon = ncommon + 1
					exit
				end if
			end do
		end do
		
		have_common_edge = (ncommon >= 2)
	end function
	
	
	subroutine merge_candidates(dim, face1, face2)
		implicit none
		integer, intent(in) :: dim
		type(face_candidate), intent(inout) :: face1
		type(face_candidate), intent(in) :: face2
			   
		face1%center = face1%center*face1%face_area + face2%center*face2%face_area
		face1%face_area = face1%face_area + face2%face_area	
		face1%center = face1%center/face1%face_area
					
		call merge_arrays(face1%original_faces, face1%n_original_faces,&
						  face2%original_faces, face2%n_original_faces)
		
		!call merge_nodes(face1, face2)		
		call merge_formative_nodes(dim, face1, face2)			  
	end subroutine
	
	subroutine merge_arrays(array1, n1, array2, n2)
		implicit none
		integer, allocatable, intent(inout) :: array1(:)
		integer, intent(inout) :: n1
		integer, intent(in) :: array2(:)
		integer, intent(in) :: n2
		
		integer, allocatable :: temp_array(:)
		integer :: new_size
		
		new_size = n1 + n2
		allocate(temp_array(new_size))
		
		temp_array(1:n1) = array1(1:n1)
		temp_array(n1+1:new_size) = array2(1:n2)
		
		call move_alloc(temp_array, array1)
		n1 = new_size
	end subroutine 
	
	subroutine merge_nodes(face1, face2)
		implicit none
		type(face_candidate), intent(inout) :: face1
		type(face_candidate), intent(in) :: face2
		
		integer, allocatable :: temp_nodes(:)
		integer :: n_total, n_unique
		logical :: dup
		integer :: i, j
		
		n_total = face1%n_fine_nodes + face2%n_fine_nodes
		allocate(temp_nodes(n_total))
		
		temp_nodes(1:face1%n_fine_nodes) = face1%fine_nodes(1:face1%n_fine_nodes)
		temp_nodes(face1%n_fine_nodes+1:n_total) = face2%fine_nodes(1:face2%n_fine_nodes)
		
		n_unique = 0
		do i = 1, n_total
			dup = .false.
			do j = 1, n_unique
				if (temp_nodes(i) == temp_nodes(j)) then
					dup = .true.
					exit
				end if
			end do
			if (.not. dup) then
				n_unique = n_unique + 1
				temp_nodes(n_unique) = temp_nodes(i)
			end if
		end do
				
		if (allocated(face1%fine_nodes)) deallocate(face1%fine_nodes)
		allocate(face1%fine_nodes(n_unique))
		face1%fine_nodes(1:n_unique) = temp_nodes(1:n_unique)
		face1%n_fine_nodes = n_unique
		
		deallocate(temp_nodes)
	end subroutine
	
	subroutine merge_formative_nodes(dim, face1, face2)
		implicit none
		integer, intent(in) :: dim
		type(face_candidate), intent(inout) :: face1
		type(face_candidate), intent(in) :: face2
		
		integer, allocatable :: temp_nodes(:)
		integer :: n_total, n_unique
		logical :: dup
		integer :: i, j, k, common_node, node1, node2
		logical :: found
		
		n_total = face1%n_fine_nodes + face2%n_fine_nodes
		allocate(temp_nodes(n_total))
		
		temp_nodes(1:face1%n_fine_nodes) = face1%fine_nodes(1:face1%n_fine_nodes)
		temp_nodes(face1%n_fine_nodes+1:n_total) = face2%fine_nodes(1:face2%n_fine_nodes)
		
		n_unique = 0
		do i = 1, n_total
			dup = .false.
			do j = 1, n_unique
				if (temp_nodes(i) == temp_nodes(j)) then
					dup = .true.
					exit
				end if
			end do
			if (.not. dup) then
				n_unique = n_unique + 1
				temp_nodes(n_unique) = temp_nodes(i)
			end if
		end do
		
		if (dim == 2 .and. n_unique > 2) then
			common_node = 0
			do i = 1, face1%n_fine_nodes
				do j = 1, face2%n_fine_nodes
					if (face1%fine_nodes(i) == face2%fine_nodes(j)) then
						common_node = face1%fine_nodes(i)
						exit
					end if
				end do
				if (common_node /= 0) exit
			end do
			
			node1 = 0; node2 = 0
			do i = 1, n_unique
				if (temp_nodes(i) /= common_node) then
					if (node1 == 0) then
						node1 = temp_nodes(i)
					else
						node2 = temp_nodes(i)
					end if
				end if
			end do
			temp_nodes(1) = node1
			temp_nodes(2) = node2
			n_unique = 2
		end if
		
		if (allocated(face1%fine_nodes)) deallocate(face1%fine_nodes)
		allocate(face1%fine_nodes(n_unique))
		face1%fine_nodes(1:n_unique) = temp_nodes(1:n_unique)
		face1%n_fine_nodes = n_unique
		
		deallocate(temp_nodes)
	end subroutine

	
	!aux sorting subrouitnes:
	recursive subroutine qsort_perm(perm, candidates, low, high)
		implicit none
        integer, intent(inout) :: perm(:)
        type(face_candidate), intent(in) :: candidates(:)
        integer, intent(in) :: low, high
        integer :: i, j, pivot_idx, temp
        type(face_candidate) :: pivot_val

        if (low < high) then
            pivot_idx = perm((low + high)/2)
            pivot_val = candidates(pivot_idx)
            i = low
            j = high
            do
                do while (is_less(candidates(perm(i)), pivot_val))
                    i = i + 1
                end do
                do while (is_less(pivot_val, candidates(perm(j))))
                    j = j - 1
                end do
                if (i >= j) exit
                temp = perm(i)
                perm(i) = perm(j)
                perm(j) = temp
                i = i + 1
                j = j - 1
            end do
            call qsort_perm(perm, candidates, low, j)
            call qsort_perm(perm, candidates, j + 1, high)
        end if
    end subroutine
	
	logical function is_less(face1, face2)
		implicit none
		type(face_candidate), intent(in) :: face1, face2
		if (face1%left_coarse < face2%left_coarse) then
			is_less = .true.
		else if (face1%left_coarse == face2%left_coarse) then
			is_less = (face1%right_coarse < face2%right_coarse)
		else
			is_less = .false.
		end if
	end function
	
	!face candidates processing subroutines: 
	subroutine process_face_candidates(dim, tolerance, candidates, num_candidates, cnfaces)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: tolerance
		integer, intent(in) :: num_candidates
		type(face_candidate), allocatable, intent(inout) :: candidates(:)
		integer, intent(inout) :: cnfaces

		integer, allocatable :: perm(:)
		integer :: i, g_start, g_end, n_group
		integer :: left_key, right_key
		type(face_candidate), allocatable :: coarse_faces(:)

		!=== QUICK SORT OF CANDIDATES BY LEFT/RIGHT CELLS: ===
		allocate(perm(num_candidates))
		do i = 1, num_candidates
			perm(i) = i
		end do
		call qsort_perm(perm, candidates, 1, num_candidates)

		!=== FACE GROUPS LOOP: ===
		allocate(coarse_faces(num_candidates))
		cnfaces = 0
		
		i = 1
		do while (i <= num_candidates)
			g_start = i
			left_key = candidates(perm(i))%left_coarse
			right_key = candidates(perm(i))%right_coarse
			i = i + 1
			do while (i <= num_candidates)
				if (candidates(perm(i))%left_coarse /= left_key .or. &
					candidates(perm(i))%right_coarse /= right_key) exit
				i = i + 1
			end do
			g_end = i - 1
			n_group = g_end - g_start + 1

			!=== PROCESSIGN FACE GROUPD [LEFT/RIGHT]: ===
			call process_group(dim, tolerance, &
							   candidates, perm(g_start:g_end), n_group, &
							   cnfaces, coarse_faces)
		end do
		
		if (allocated(candidates)) deallocate(candidates)
		allocate(candidates(cnfaces))
		candidates(1:cnfaces) = coarse_faces(1:cnfaces)

		deallocate(perm, coarse_faces)
	end subroutine
	
	subroutine process_group(dim, tolerance, all_candidates,&
							 group_idx, n_group,&
							 cnfaces, out_faces)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: tolerance
		type(face_candidate), intent(in) :: all_candidates(:)
		integer, intent(in) :: group_idx(:)
		integer, intent(in) :: n_group
		integer, intent(inout) :: cnfaces
		type(face_candidate), intent(inout) :: out_faces(:)

		integer :: i, j
		logical :: merged
		type(face_candidate), allocatable :: temp_faces(:)
		integer :: n_temp

		allocate(temp_faces(n_group))
		n_temp = 0

		do i = 1, n_group
			merged = .false.
			
			do j = 1, n_temp
				if (can_merge_faces(temp_faces(j), all_candidates(group_idx(i)), tolerance, dim)) then
					call merge_candidates(dim, temp_faces(j), all_candidates(group_idx(i)))
					merged = .true.
					exit
				end if
			end do
			if (.not. merged) then
				n_temp = n_temp + 1
				temp_faces(n_temp) = all_candidates(group_idx(i))
			end if
		end do
		
		do i = 1, n_temp
			cnfaces = cnfaces + 1
			out_faces(cnfaces) = temp_faces(i)
		end do

		deallocate(temp_faces)
	end subroutine
	

	!faces reordering subroutine: 
	subroutine reorder_faces_n_data(dim, ncells, nfaces, nbfaces, face_left_cell, face_right_cell,&
								    face_zone, face_type, face_nodes_ptr, face_nodes,&
								    face_normal, face_area, face_center)
		implicit none
		integer, intent(in) :: dim, ncells, nfaces
		integer, intent(inout) :: nbfaces
		integer, allocatable, intent(inout) :: face_left_cell(:), face_right_cell(:),&
											   face_zone(:), face_type(:),&
											   face_nodes_ptr(:), face_nodes(:)
		real(8), allocatable, intent(inout) :: face_normal(:, :), face_area(:), face_center(:, :)
		
		integer :: face_idx
		integer :: left_cell, right_cell
		integer :: n_int, n_ext, idx, pos1, pos2, pos1_old, pos2_old
		integer, allocatable :: old_to_new(:), new_to_old(:),& 
								new_left(:), new_right(:),&
								new_zone(:), new_type(:),&
								new_face_nodes(:), new_face_nodes_ptr(:)
		real(8), allocatable :: new_face_normal(:, :), new_face_area(:), new_face_center(:, :)

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
				 new_zone(nfaces), new_type(nfaces),&
				 new_face_normal(dim, nfaces), new_face_area(nfaces),&
				 new_face_center(dim, nfaces))
		do face_idx = 1, nfaces
			new_left(old_to_new(face_idx)) = face_left_cell(face_idx)
			new_right(old_to_new(face_idx)) = face_right_cell(face_idx)
			new_zone(old_to_new(face_idx)) = face_zone(face_idx)
			new_type(old_to_new(face_idx)) = face_type(face_idx)
			new_face_normal(:, old_to_new(face_idx)) = face_normal(:, face_idx)
			new_face_center(:, old_to_new(face_idx)) = face_center(:, face_idx)
			new_face_area(old_to_new(face_idx)) = face_area(face_idx)
		end do
		call move_alloc(new_left, face_left_cell)
		call move_alloc(new_right, face_right_cell)
		call move_alloc(new_zone, face_zone)
		call move_alloc(new_type, face_type)
		call move_alloc(new_face_normal, face_normal)
		call move_alloc(new_face_center, face_center)
		call move_alloc(new_face_area, face_area)
		
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
		
	
	!coarse faces buiding: 
	subroutine mg_build_coarse_faces(dim, fncells, fnfaces, fnbfaces,&
									 fface_left_cell, fface_right_cell,&
									 fface_zone, fface_type,&
									 fface_nodes_ptr, fface_nodes,&
									 fface_center, fface_normal, fface_area,&
									 fine_to_coarse,&
									 cncells, cnfaces, cnbfaces,&
									 cface_center, cface_normal, cface_area,&
									 cface_zone, cface_type, cface_bidx,&
									 cface_left_cell, cface_right_cell,&
									 cface_nodes_ptr, cface_nodes)
									 
		implicit none
		integer, intent(in) :: dim, fncells, fnfaces, fnbfaces
		integer, intent(in) :: fface_left_cell(:), fface_right_cell(:)
		integer, intent(in) :: fface_zone(:), fface_type(:)
		integer, intent(in) :: fface_nodes_ptr(:), fface_nodes(:)
		integer, intent(in) :: fine_to_coarse(:)
		real(8), intent(in) :: fface_center(:, :), fface_normal(:, :), fface_area(:)
		
		integer, intent(in) :: cncells
		integer, intent(inout) :: cnfaces, cnbfaces
		integer, allocatable, intent(inout) :: cface_zone(:), cface_type(:), cface_bidx(:)
		integer, allocatable, intent(inout) :: cface_left_cell(:), cface_right_cell(:)
		integer, allocatable, intent(inout) :: cface_nodes_ptr(:), cface_nodes(:)
		real(8), allocatable, intent(inout) :: cface_center(:, :), cface_normal(:, :), cface_area(:)
		
		
		type(face_candidate), allocatable :: candidates(:)
		integer :: num_candidates
		integer :: i
		
		
		!=== FINDING FACE CANDIDATES: ===
		call collect_face_candidates(dim, fncells, fnfaces, fnbfaces,&
									 fface_left_cell, fface_right_cell,&
									 fface_zone, fface_type,&
									 fface_nodes_ptr, fface_nodes,&
									 fface_center, fface_normal, fface_area,&
									 fine_to_coarse,&
									 candidates, num_candidates)
		
		!=== MERGING FACE CANDIDATES: ===
		call process_face_candidates(dim, 0.001d0, candidates, num_candidates, cnfaces)
		
		!=== BUILDING COARSE FACES: ===
		allocate(cface_left_cell(cnfaces), cface_right_cell(cnfaces))
		allocate(cface_zone(cnfaces), cface_type(cnfaces), cface_bidx(cnfaces))
		allocate(cface_area(cnfaces))
		allocate(cface_center(dim, cnfaces), cface_normal(dim, cnfaces))
		allocate(cface_nodes_ptr(cnfaces + 1))
		cface_nodes_ptr(1) = 1
		
		do i = 1, cnfaces
			cface_left_cell(i) = candidates(i)%left_coarse
			cface_right_cell(i) = candidates(i)%right_coarse
			cface_zone(i) = candidates(i)%face_zone
			cface_type(i) = candidates(i)%face_type
			
			cface_area(i) = candidates(i)%face_area
			cface_center(:, i) = candidates(i)%center(1:dim)
			cface_normal(:, i) = candidates(i)%normal(1:dim)
			
			cface_nodes_ptr(i+1) = cface_nodes_ptr(i) + candidates(i)%n_fine_nodes
		end do
	
		allocate(cface_nodes(cface_nodes_ptr(cnfaces+1) - 1))
		
		do i = 1, cnfaces
			cface_nodes(cface_nodes_ptr(i):cface_nodes_ptr(i+1)-1) =&
				candidates(i)%fine_nodes(1:candidates(i)%n_fine_nodes)
		end do
		
		call reorder_faces_n_data(dim, cncells, cnfaces, cnbfaces, cface_left_cell, cface_right_cell,&
								  cface_zone, cface_type, cface_nodes_ptr, cface_nodes,&
								  cface_normal, cface_area, cface_center)
		
		call order_boundary_faces_by_zone(dim, cncells, cnfaces, cnbfaces,&
                                          cface_left_cell, cface_right_cell,&
                                          cface_zone, cface_type, cface_bidx,&
                                          cface_nodes_ptr, cface_nodes,&
                                          cface_normal, cface_area, cface_center)
                                         
		cface_bidx = -1
		do i = cnfaces - cnbfaces + 1, cnfaces
			cface_bidx(i) = i - (cnfaces - cnbfaces)
			cface_right_cell(i) = cncells+cface_bidx(i)
		end do
		
		deallocate(candidates)
	end subroutine
	

!=======================================================================
!=============== COARSE MESH CELLS GEOMETRY BUILDING ===================
!=======================================================================
	subroutine mg_build_coarse_cells_geometry(dim, fcell_volume, fcell_center,&
											  fine_cells_ptr, fine_cells,&
											  coarse_ncells, ccell_volume, ccell_center)
		implicit none
		integer, intent(in) :: dim
		integer, intent(in) :: coarse_ncells
		real(8), intent(in) :: fcell_volume(:), fcell_center(:, :)
		integer, intent(in) :: fine_cells_ptr(:), fine_cells(:)
		real(8), allocatable, intent(inout) :: ccell_volume(:), ccell_center(:, :)
		
		integer :: i
		integer :: ccell_idx, fcell_idx
		real(8) :: centroid(dim), total_volume
		
		
		if (allocated(ccell_volume)) deallocate(ccell_volume)
		if (allocated(ccell_center)) deallocate(ccell_center)
		
		allocate(ccell_volume(coarse_ncells), ccell_center(dim, coarse_ncells))
		ccell_volume = 0.d0; ccell_center = 0.d0
		
		
		!=== COARSE CELLS LOOP: ===
		do ccell_idx = 1, coarse_ncells
			total_volume = 0.d0
			centroid = 0.d0
			
			do i = fine_cells_ptr(ccell_idx), fine_cells_ptr(ccell_idx+1)-1
				fcell_idx = fine_cells(i)
				
				total_volume = total_volume + fcell_volume(fcell_idx)
				centroid = centroid + fcell_volume(fcell_idx)*fcell_center(:, fcell_idx)
			end do
			
			ccell_volume(ccell_idx) =  total_volume
			ccell_center(:, ccell_idx) = centroid/total_volume			
		end do	
				
	end subroutine
	

!=======================================================================
!===================== COARSE MESH NODES BUILDING ======================
!=======================================================================	
	subroutine mg_build_coarse_nodes(fnodes, fnode_coords,&
									 cnfaces, cface_nodes_ptr, cface_nodes,&
									 cnnodes, cnode_coords)
		implicit none
		integer, intent(in) :: fnodes
		real(8), intent(in) :: fnode_coords(:,:)
		integer, intent(in) :: cnfaces
		integer, intent(in) :: cface_nodes_ptr(:)
		integer, intent(inout) :: cface_nodes(:)
		integer, intent(inout) :: cnnodes
		real(8), allocatable, intent(inout) :: cnode_coords(:,:)
		
		integer, allocatable :: fnode_to_cnode(:)
		integer :: i, fnode_idx, node_count
		
		!=== FINE TO COARSE MAP: ===
		allocate(fnode_to_cnode(fnodes))
		fnode_to_cnode = 0
		node_count = 0
		
		!=== ALL NODES LOOP: ===
		do i = 1, size(cface_nodes)
			fnode_idx = cface_nodes(i)
			
			if (fnode_to_cnode(fnode_idx) == 0) then
				node_count = node_count + 1
				fnode_to_cnode(fnode_idx) = node_count
			end if
			
			cface_nodes(i) = fnode_to_cnode(fnode_idx)
		end do
		
		cnnodes = node_count
		
		!=== COARSE MESH NODE COORDINATES: ===
		allocate(cnode_coords(size(fnode_coords, 1), cnnodes))
		
		do fnode_idx = 1, fnodes
			if (fnode_to_cnode(fnode_idx) > 0) then
				i = fnode_to_cnode(fnode_idx)
				cnode_coords(:, i) = fnode_coords(:, fnode_idx)
			end if
		end do
		
		deallocate(fnode_to_cnode)
		
	end subroutine
	

!=======================================================================
!===================== COARSE MESH FACES BUILDING ======================
!=======================================================================	
	subroutine mg_build_coarse_mesh(dim, fncells, fnfaces, fnbfaces, fnnodes,&
									fface_left_cell, fface_right_cell,&
									fface_zone, fface_type,&
									fface_nodes_ptr, fface_nodes,&
									fface_center, fface_normal, fface_area,&
									fnode_coords, fcell_volume, fcell_center,&
									fine_to_coarse, fine_cells_ptr, fine_cells,&
									cncells, cnfaces, cnbfaces, cnnodes,&
									cface_center, cface_normal, cface_area, cface_weight,&
									ccell_volume, ccell_center, cnode_coords,&
									cface_zone, cface_type, cface_bidx,&
									cface_left_cell, cface_right_cell,&
									cface_nodes_ptr, cface_nodes,&
									ccell_nodes_ptr, ccell_nodes,&
									ccell_faces_ptr, ccell_faces, ccell_type)
									 
		implicit none
		integer, intent(in) :: dim, fncells, fnfaces, fnbfaces, fnnodes
		integer, intent(in) :: fface_left_cell(:), fface_right_cell(:)
		integer, intent(in) :: fface_zone(:), fface_type(:)
		integer, intent(in) :: fface_nodes_ptr(:), fface_nodes(:)
		integer, intent(in) :: fine_to_coarse(:), fine_cells_ptr(:), fine_cells(:)
		real(8), intent(in) :: fface_center(:, :), fface_normal(:, :), fface_area(:),&
							   fnode_coords(:, :), fcell_volume(:), fcell_center(:, :)
		
		integer, intent(in) :: cncells
		integer, intent(inout) :: cnfaces, cnbfaces, cnnodes
		integer, allocatable, intent(inout) :: cface_zone(:), cface_type(:), cface_bidx(:)
		integer, allocatable, intent(inout) :: cface_left_cell(:), cface_right_cell(:)
		integer, allocatable, intent(inout) :: cface_nodes_ptr(:), cface_nodes(:),&
											   ccell_nodes_ptr(:), ccell_nodes(:),&
											   ccell_faces_ptr(:), ccell_faces(:), ccell_type(:)
		real(8), allocatable, intent(inout) :: cface_center(:, :), cface_normal(:, :),&
											   cface_area(:), cnode_coords(:, :),&
											   ccell_volume(:), ccell_center(:, :),&
											   cface_weight(:)
		
		!=== COARSE FACES: ===
		call mg_build_coarse_faces(dim, fncells, fnfaces, fnbfaces,&
								   fface_left_cell, fface_right_cell,&
								   fface_zone, fface_type,&
								   fface_nodes_ptr, fface_nodes,&
								   fface_center, fface_normal, fface_area,&
								   fine_to_coarse,&
								   cncells, cnfaces, cnbfaces,&
								   cface_center, cface_normal, cface_area,&
								   cface_zone, cface_type, cface_bidx,&
								   cface_left_cell, cface_right_cell,&
								   cface_nodes_ptr, cface_nodes)
						   
		!=== COARSE NODES: ===
		call mg_build_coarse_nodes(fnnodes, fnode_coords,&
								   cnfaces, cface_nodes_ptr, cface_nodes,&
								   cnnodes, cnode_coords)
		
		!=== COARSE CELLS: ===
		call mg_build_coarse_cells_geometry(dim, fcell_volume, fcell_center,&
										    fine_cells_ptr, fine_cells,&
										    cncells, ccell_volume, ccell_center)
										    								    
		call cell_face_connectivity(cncells, cnfaces,&
									cface_left_cell, cface_right_cell,&
									ccell_faces_ptr, ccell_faces)	
									
		call cell_node_connectivity(cnnodes, cncells,&
									cface_nodes_ptr, cface_nodes,&
									ccell_faces_ptr, ccell_faces,&
									ccell_nodes_ptr, ccell_nodes)
									
		if (dim == 2) then
			call order_cell_nodes_2d(cncells, cnode_coords,&
									 ccell_nodes_ptr, ccell_nodes)	
		end if						
		allocate(ccell_type(cncells))
		ccell_type = 3
		
		
		allocate(cface_weight(cnfaces))
		cface_weight = 0.5d0
		
	end subroutine












end module
