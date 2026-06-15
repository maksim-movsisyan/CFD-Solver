module mg_aglomeration_module
implicit none
contains
!=======================================================================
!============== AGGLOMERATION INFO GENERATION SUBROUTINE ===============
!=======================================================================
	subroutine mg_fill_agglomeation_info(fine_ncells, coarse_ncells, fine_to_coarse,&
										 fine_cells_ptr, fine_cells)
		implicit none
		integer, intent(in) :: fine_ncells, coarse_ncells
		integer, intent(in) :: fine_to_coarse(:)
		integer, allocatable, intent(inout) :: fine_cells_ptr(:), fine_cells(:)
		
		integer, allocatable :: nfine_per_coarse(:)
		integer :: cell_idx, coarse_idx
		integer :: i, pos
		
		allocate(nfine_per_coarse(coarse_ncells))
		nfine_per_coarse = 0
		
		!=== FINDING NUMBER OF FINE CELLS PER EACH COARSE CELL: ===
		do cell_idx = 1, fine_ncells
			coarse_idx = fine_to_coarse(cell_idx)
			nfine_per_coarse(coarse_idx) = nfine_per_coarse(coarse_idx) + 1
		end do 
		
		!=== FILLING FINE_CELLS_PTR: ===
		if (allocated(fine_cells_ptr)) deallocate(fine_cells_ptr)
		allocate(fine_cells_ptr(coarse_ncells + 1))
		fine_cells_ptr(1) = 1
		do i = 1, coarse_ncells
			fine_cells_ptr(i+1) = fine_cells_ptr(i) + nfine_per_coarse(i)
		end do
		
		!=== FILLING FINE_CELLS: ===
		if (allocated(fine_cells)) deallocate(fine_cells)
		allocate(fine_cells(fine_cells_ptr(coarse_ncells + 1) - 1))
		nfine_per_coarse = 0
		do cell_idx = 1, fine_ncells
			coarse_idx = fine_to_coarse(cell_idx)
			pos = fine_cells_ptr(coarse_idx) + nfine_per_coarse(coarse_idx)
			fine_cells(pos) = cell_idx
			nfine_per_coarse(coarse_idx) = nfine_per_coarse(coarse_idx) + 1
		end do
		
		deallocate(nfine_per_coarse)
	end subroutine


!=======================================================================
!================== FRONTAL AGGLOMETARION SUBROUTINE ===================
!=======================================================================
	subroutine mg_perform_frontal_agglomeration(dim, ncells, nfaces, nbfaces,&
												face_left_cell, face_right_cell,&
												cell_faces_ptr, cell_faces,&
												face_area, cell_center,&
												target_fine_per_coarse,&
												coarse_ncells, fine_to_coarse)
		implicit none
		integer, intent(in) :: target_fine_per_coarse
		integer, intent(in) :: dim, ncells, nfaces, nbfaces
		integer, intent(in) :: face_left_cell(:), face_right_cell(:)
		integer, intent(in) :: cell_faces_ptr(:), cell_faces(:)
		real(8), intent(in) :: face_area(:), cell_center(:, :)
		integer, intent(inout) :: coarse_ncells
		integer, allocatable, intent(inout) :: fine_to_coarse(:)
		
		integer :: face_idx, f_idx, cell_idx, c_idx, j, k
		integer :: i, neighbor, nb_of_nb, best_neighbor
		integer :: current_seed, cells_in_agg
		real(8) :: max_quality, current_quality, total_interface_area 
		integer, allocatable :: cell_status(:) 							!=== 0: FREE, -1: IN QUEUE, 1: AGGREGATED ===
		integer, allocatable :: queue(:), agg_cells(:)
		integer :: head, tail, next_search_idx
		
		!optional params:
		!1) direction of agglomeration
		logical, parameter :: use_direction = .false.
		real(8) :: pref_dir(3) = [1.0d0, 0.0d0, 0.0d0]
		real(8) :: beta = 1.25d0
		real(8) :: dx(dim), proj
			
		!=== INITIALIZATION ===
		if (allocated(fine_to_coarse)) deallocate(fine_to_coarse)
		allocate(cell_status(ncells), fine_to_coarse(ncells))
		allocate(queue(ncells), agg_cells(target_fine_per_coarse))
		
		cell_status = 0
		fine_to_coarse = 0
		coarse_ncells = 0
		head = 1; tail = 0
		next_search_idx = 1
				
		!=== WAVE FRONT INITIALIZATION: ===
		block
			integer, allocatable :: internal_neigh_count(:), head_buckets(:), next_cell(:)
			logical, allocatable :: visited(:)
			integer :: cell_neighbor_idx(200)
			integer :: max_neigh, n_neigh
			
			allocate(internal_neigh_count(ncells))
			allocate(visited(ncells))
			internal_neigh_count = 0
			visited = .false.
			max_neigh = 0
            
			!=== FINDING ALL NEIGHBORS FOR EACH BOUND CELL: ====
			do f_idx = nfaces - nbfaces + 1, nfaces
				c_idx = face_left_cell(f_idx)
				if (cell_status(c_idx) == 0) then
                    n_neigh = 0
                  
                    do k = cell_faces_ptr(c_idx), cell_faces_ptr(c_idx+1)-1
                        face_idx = cell_faces(k)
                        neighbor = face_left_cell(face_idx)
                        if (neighbor == c_idx) neighbor = face_right_cell(face_idx)
                        
                        if (neighbor <= 0 .or. neighbor > ncells) cycle
                        
                        if (.not. visited(neighbor)) then
							n_neigh = n_neigh + 1
							visited(neighbor) = .true.
							cell_neighbor_idx(n_neigh) = neighbor
						end if
                    end do
                    
					
                    internal_neigh_count(c_idx) = n_neigh
                    max_neigh = max(max_neigh, n_neigh)
                    cell_status(c_idx) = -2
                    
                     do i = 1, n_neigh
						visited(cell_neighbor_idx(i)) = .false.
					end do
				end if
			end do
			
			!=== COUNTING SORT OF BOUNDARY CELLS: ===
			allocate(head_buckets(0:max_neigh))
            allocate(next_cell(ncells))
            head_buckets = 0; next_cell = 0
            
            do c_idx = 1, ncells
                if (cell_status(c_idx) == -2) then
                    n_neigh = internal_neigh_count(c_idx)
                    next_cell(c_idx) = head_buckets(n_neigh)
                    head_buckets(n_neigh) = c_idx
                end if
            end do
            
			!=== QUEUE FILLING: ===
			do i = 0, max_neigh
                c_idx = head_buckets(i)
                do while (c_idx > 0)
                    tail = tail + 1
                    queue(tail) = c_idx
                    cell_status(c_idx) = -1
                    c_idx = next_cell(c_idx)
                end do
			end do
			
			deallocate(internal_neigh_count, head_buckets, visited, next_cell)
		end block
		
		!=== MAIN LOOP: ===
		do while (next_search_idx <= ncells .or. head <= tail)
			if (head <= tail) then
				current_seed = queue(head)
				head = head + 1
			else
				current_seed = -1
				do i = next_search_idx, ncells
					if (cell_status(i) <= 0) then
						current_seed = i
						next_search_idx = i + 1
						exit
					end if
				end do
				if (current_seed == -1) exit
			end if
			
			
			if (cell_status(current_seed) == 1) cycle
			
			!=== START NEW AGGLOMERATE: ===
			coarse_ncells = coarse_ncells + 1
			fine_to_coarse(current_seed) = coarse_ncells
			cell_status(current_seed) = 1
			cells_in_agg = 1
            agg_cells(1) = current_seed
			
			!=== FINDING BEST NEIGBORS: ===
			do while (cells_in_agg < target_fine_per_coarse)
				best_neighbor = -1
				max_quality = -1.0e30 
				
				!=== ALL AGLOMERATE NEIGHBORS LOOP: ===
                do j = 1, cells_in_agg
                    cell_idx = agg_cells(j)
                    
					do k = cell_faces_ptr(cell_idx), cell_faces_ptr(cell_idx+1)-1
						face_idx = cell_faces(k)
						neighbor = face_left_cell(face_idx)
						if (neighbor == cell_idx) neighbor = face_right_cell(face_idx)
						
						if (neighbor > 0 .and. neighbor <= ncells) then
                            if (cell_status(neighbor) /= 1) then
                            
								!=== FINDING NEIGHBOR QUALITY: ===
								total_interface_area = 0.0d0
								do i = cell_faces_ptr(neighbor), cell_faces_ptr(neighbor+1)-1
									f_idx = cell_faces(i)
									nb_of_nb = face_left_cell(f_idx)
									if (nb_of_nb == neighbor) nb_of_nb = face_right_cell(f_idx)
									
									if (nb_of_nb > 0 .and. nb_of_nb <= ncells) then
										if (fine_to_coarse(nb_of_nb) == coarse_ncells) then
											total_interface_area = total_interface_area + face_area(f_idx)
										end if
									end if
								end do
                            
								current_quality = compute_quality(dim, agg_cells, cells_in_agg, neighbor,&
                                                                  cell_center, total_interface_area)
								
								!optional feature:
								!1) direction of agglomeration
								if (use_direction .and. cells_in_agg == 1) then
									dx = cell_center(:,neighbor) - cell_center(:,current_seed)
									dx = dx/norm2(dx)
									proj = dabs(dot_product(dx, pref_dir(1:dim)))
									current_quality = current_quality*(1.0d0 + beta*proj)
								end if

                                if (current_quality > max_quality) then
									max_quality = current_quality
									best_neighbor = neighbor
								end if
                            end if
                        end if
					end do
				end do
				
				if (best_neighbor /= -1) then
					cells_in_agg = cells_in_agg + 1
                    agg_cells(cells_in_agg) = best_neighbor
					fine_to_coarse(best_neighbor) = coarse_ncells
					cell_status(best_neighbor) = 1
					
					!=== ADDING NEIGHBORS OF SEED NEIGBOR TO QUEUE: ===
					do k = cell_faces_ptr(best_neighbor), cell_faces_ptr(best_neighbor+1)-1
						face_idx = cell_faces(k)
						neighbor = face_left_cell(face_idx)
						if (neighbor == best_neighbor) neighbor = face_right_cell(face_idx)
						if (neighbor > 0 .and. neighbor <= ncells) then
							if (cell_status(neighbor) == 0) then
								tail = tail + 1
								queue(tail) = neighbor
								cell_status(neighbor) = -1
							end if
						end if
					end do
				else
					exit
				end if
			end do
			
		end do
		
		
		!=== POST-PROCESSING (REMOVE SINGLETONS): ===
		block
			integer, allocatable :: cells_per_agg(:)
			integer :: c_idx, f_idx, best_agg, neighbor_agg
			real(8) :: max_q, q
			
			allocate(cells_per_agg(coarse_ncells))
			cells_per_agg = 0
			do i = 1, ncells
				cells_per_agg(fine_to_coarse(i)) = cells_per_agg(fine_to_coarse(i)) + 1
			end do
			
			do c_idx = 1, ncells
				if (cells_per_agg(fine_to_coarse(c_idx)) == 1) then
					best_agg = -1
					max_q = -1.0e30
					
					!=== FINDING BEST AGGLOMERATE TO ADD: ===
					do k = cell_faces_ptr(c_idx), cell_faces_ptr(c_idx+1)-1
						face_idx = cell_faces(k)
						neighbor = face_left_cell(face_idx)
						if (neighbor == c_idx) neighbor = face_right_cell(face_idx)
						
						if (neighbor > 0 .and. neighbor <= ncells) then
							neighbor_agg = fine_to_coarse(neighbor)
							
							if (neighbor_agg /= fine_to_coarse(c_idx) .and. &
                                cells_per_agg(neighbor_agg) > 1) then
                                
								total_interface_area = 0.0d0
								do i = cell_faces_ptr(c_idx), cell_faces_ptr(c_idx+1)-1
									f_idx = cell_faces(i)
									nb_of_nb = face_left_cell(f_idx)
									if (nb_of_nb == c_idx) nb_of_nb = face_right_cell(f_idx)
									
									if (nb_of_nb > 0 .and. nb_of_nb <= ncells) then
										if (fine_to_coarse(nb_of_nb) == neighbor_agg) then
											total_interface_area = total_interface_area + face_area(f_idx)
										end if
									end if
								end do
								
								q = total_interface_area/(norm2(cell_center(:,c_idx) - cell_center(:,neighbor)))
								
								if (q > max_q) then
									max_q = q
									best_agg = neighbor_agg
								end if
							end if
						end if
					end do
					
					if (best_agg /= -1) then
						fine_to_coarse(c_idx) = best_agg
					end if
				end if
			end do
			deallocate(cells_per_agg)
		end block
		
		!=== RENUMBERING: REMOVE EMPTY COARSE CELLS ===
        block
            integer, allocatable :: old_to_new(:)
            integer :: n_active, old_id
            
            allocate(old_to_new(coarse_ncells))
            old_to_new = 0
            n_active = 0
            
            do i = 1, ncells
                old_id = fine_to_coarse(i)
                if (old_to_new(old_id) == 0) then
                    n_active = n_active + 1
                    old_to_new(old_id) = n_active
                end if
            end do
            
            do i = 1, ncells
                fine_to_coarse(i) = old_to_new(fine_to_coarse(i))
            end do
            
            coarse_ncells = n_active
            
            deallocate(old_to_new)
            print*, 'Agglomeration finished. Coarse cells:', coarse_ncells
        end block


		deallocate(cell_status, queue, agg_cells)
		
	end subroutine


!=======================================================================
!==================== COMPUTE AGGLOMERATION QUALITY ====================
!=======================================================================
	function compute_quality(dim, agg_cells, cells_in_agg, potential_neighbor,&
						     cell_center, face_area) result(q)
		implicit none
		integer, intent(in) :: dim
		integer, intent(in) :: agg_cells(:), cells_in_agg, potential_neighbor
		real(8), intent(in) :: cell_center(:,:)
		real(8), intent(in) :: face_area       
		real(8) :: q
		
		real(8) :: agg_centroid(dim), dist_sq
		integer :: i, d
		
		agg_centroid = 0.0d0
		do i = 1, cells_in_agg
			agg_centroid(:) = agg_centroid(:) + cell_center(:, agg_cells(i))
		end do
		agg_centroid = agg_centroid/real(cells_in_agg, 8)
		
		dist_sq = 0.0d0
		do d = 1, dim
			dist_sq = dist_sq + (cell_center(d, potential_neighbor) - agg_centroid(d))**2
		end do
		
		q = face_area/(dist_sq)
		
	end function


!=======================================================================
!==================== MAIN AGGLOMETARION SUBROUTINE ====================
!=======================================================================
	subroutine mg_perform_agglomeration(dim, ncells, nfaces, nbfaces,&
										face_left_cell, face_right_cell,&
										cell_faces_ptr, cell_faces,&
										face_area, cell_center,&
										target_fine_per_coarse,&
										coarse_ncells, fine_to_coarse,&
										fine_cells_ptr, fine_cells) 
		implicit none
		integer, intent(in) :: target_fine_per_coarse
		integer, intent(in) :: dim, ncells, nfaces, nbfaces
		integer, intent(in) :: face_left_cell(:), face_right_cell(:)
		integer, intent(in) :: celL_faces_ptr(:), cell_faces(:)
		real(8), intent(in) :: face_area(:), cell_center(:, :)
		integer, intent(inout) :: coarse_ncells
		integer, allocatable, intent(inout) :: fine_to_coarse(:), fine_cells_ptr(:), fine_cells(:)
		
		
		call mg_perform_frontal_agglomeration(dim, ncells, nfaces, nbfaces,&
											  face_left_cell, face_right_cell,&
											  cell_faces_ptr, cell_faces,&
											  face_area, cell_center,&
											  target_fine_per_coarse,&
											  coarse_ncells, fine_to_coarse)
		call mg_fill_agglomeation_info(ncells, coarse_ncells, fine_to_coarse,&
									   fine_cells_ptr, fine_cells)
		
	end subroutine






end module
