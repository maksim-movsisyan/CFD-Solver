module gradient_computation_module
implicit none
contains

!=======================================================================
!============ GREEN GAUSS GRADIENT COMPUTATION SUBROUTINE ==============
!=======================================================================
	!=== FIELD SUBS: ===
	pure subroutine compute_grad_scalarfield_gg(Field, Field_grad,&
										        face_left_cell, face_right_cell,&
										        face_normal, face_area, cell_volume, face_weight,&
										        ncells, nfaces, nbfaces, dim, USE_GHOST_CELLS)
		implicit none
		integer, intent(in) :: ncells, nfaces, nbfaces, dim
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:),&
										   cell_volume(:), face_weight(:)
		real(8), intent(in), contiguous :: Field(:)
		real(8), intent(inout), contiguous :: Field_grad(:, :)
		
		integer :: n_loop_faces
		integer :: i, left_cell, right_cell
		real(8) :: w, avg_face_value, val(dim), inv_vol
		
		Field_grad = 0.d0
		
		n_loop_faces = nfaces - nbfaces
		if (USE_GHOST_CELLS) n_loop_faces = nfaces
		
		!=== INTERNAL FACES: ===
		do i = 1, n_loop_faces
			left_cell = face_left_cell(i)
			right_cell = face_right_cell(i)
			w = face_weight(i)
			
			avg_face_value = w*Field(right_cell) + (1.d0 - w)*Field(left_cell)
			val = avg_face_value*face_normal(:, i)*face_area(i)
			
			Field_grad(:, left_cell) = Field_grad(:, left_cell) + val
			Field_grad(:, right_cell) = Field_grad(:, right_cell) - val
		end do
		
		if (.not. USE_GHOST_CELLS) then
			!=== UPDATING ACCORDING BOUNDARY FACES VALUES: ===
			do i = nfaces - nbfaces + 1, nfaces
				left_cell = face_left_cell(i)
				right_cell = face_right_cell(i)
							
				avg_face_value = Field(right_cell)
				val = avg_face_value*face_normal(:, i)*face_area(i)
				
				Field_grad(:, left_cell) = Field_grad(:, left_cell) + val
			end do
		end if
		
		!=== DEVIDING BY CELL VOLUME: ===
		do i = 1, ncells
			inv_vol = 1.d0/cell_volume(i)
			Field_grad(:, i) = Field_grad(:, i)*inv_vol
		end do
	end subroutine

	pure subroutine compute_grad_vectorfield_gg(Field, Field_grad,&
										        face_left_cell, face_right_cell,&
										        face_normal, face_area, cell_volume, face_weight,&
										        ncells, nfaces, nbfaces, dim, nvars, USE_GHOST_CELLS)
		implicit none
		integer, intent(in) :: ncells, nfaces, nbfaces, dim, nvars
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:),&
										   cell_volume(:), face_weight(:)
		real(8), intent(in), contiguous :: Field(:, :)
		real(8), intent(inout), contiguous :: Field_grad(:, :)
		
		integer :: n_loop_faces
		integer :: i, j, d, offset
		integer :: left_cell, right_cell
		real(8) :: w, avg_face_value(nvars), val, inv_vol
		
		Field_grad = 0.d0
		
		n_loop_faces = nfaces - nbfaces
		if (USE_GHOST_CELLS) n_loop_faces = nfaces
		
		!=== INTERNAL FACES: ===
		do i = 1, n_loop_faces
			left_cell = face_left_cell(i)
			right_cell = face_right_cell(i)
			w = face_weight(i)
			
			avg_face_value = w*Field(:, right_cell) + (1.d0 - w)*Field(:, left_cell)
			
			do j = 1, nvars
				offset = (j-1)*dim
				do d = 1, dim
					val = avg_face_value(j)*face_normal(d, i)*face_area(i)
					Field_grad(d + offset, left_cell) = Field_grad(d + offset, left_cell) + val
					Field_grad(d + offset, right_cell) = Field_grad(d + offset, right_cell) - val
				end do
			end do
		end do
		
		if (.not. USE_GHOST_CELLS) then
			!=== UPDATING ACCORDING BOUNDARY FACES VALUES: ===
			do i = nfaces - nbfaces + 1, nfaces
				left_cell = face_left_cell(i)
				right_cell = face_right_cell(i)	
				
				avg_face_value = Field(:, right_cell)
				
				do j = 1, nvars
					offset = (j-1)*dim
					do d = 1, dim
						val = avg_face_value(j)*face_normal(d, i)*face_area(i)
						Field_grad(d + offset, left_cell) = Field_grad(d + offset, left_cell) + val
					end do
				end do
			end do
		end if
		
		!=== DEVIDING BY CELL VOLUME: ===
		do i = 1, ncells
			inv_vol = 1.d0/cell_volume(i)
			Field_grad(:, i) = Field_grad(:, i)*inv_vol
		end do
	end subroutine

	!=== CELL SUBS: ===
	pure subroutine compute_grad_scalarcell_gg(Field, Cell_grad, cell_idx, ncells, dim,&
										       face_left_cell, face_right_cell,&
										       cell_faces, cell_faces_ptr,&
										       face_normal, face_area,&
										       cell_volume, face_weight, USE_GHOST_CELLS)
		implicit none
		integer, intent(in) :: cell_idx, ncells, dim
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
										   cell_faces(:), cell_faces_ptr(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:),&
										   cell_volume(:), face_weight(:)
		real(8), intent(in), contiguous :: Field(:)
		real(8), intent(inout), contiguous :: Cell_grad(:)
		
		integer :: i, left_cell, right_cell, face_idx
		integer :: pos1, pos2
		real(8) :: w, avg_face_value, val(dim), inv_vol
		
		Cell_grad = 0.d0
		
		pos1 = cell_faces_ptr(cell_idx)
		pos2 = cell_faces_ptr(cell_idx + 1)
		
		do i = pos1, pos2 - 1
			face_idx = cell_faces(i)
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if (right_cell <= ncells .OR. USE_GHOST_CELLS) then
				w = face_weight(face_idx)
				
				avg_face_value = w*Field(right_cell) + (1.d0 - w)*Field(left_cell)
				val = avg_face_value*face_normal(:, face_idx)*face_area(face_idx)
				
				if (left_cell == cell_idx) then
					Cell_grad(:) = Cell_grad(:) + val
				else	
					Cell_grad(:) = Cell_grad(:) - val
				end if
							
			else				
				avg_face_value = Field(right_cell)
				val = avg_face_value*face_normal(:, face_idx)*face_area(face_idx)
				Cell_grad(:) = Cell_grad(:) + val
			end if
		end do
		
		inv_vol = 1.d0/cell_volume(cell_idx)
		Cell_grad(:) = Cell_grad(:)*inv_vol
	end subroutine

	pure subroutine compute_grad_vectorcell_gg(Field, Cell_grad, cell_idx,&
											   ncells, dim, nvars,&
										       face_left_cell, face_right_cell,&
										       cell_faces, cell_faces_ptr,&
										       face_normal, face_area,&
										       cell_volume, face_weight, USE_GHOST_CELLS)
		implicit none
		integer, intent(in) :: cell_idx, ncells, dim, nvars
		logical, intent(in) :: USE_GHOST_CELLS
		integer, intent(in), contiguous :: face_left_cell(:), face_right_cell(:),&
										   cell_faces(:), cell_faces_ptr(:)
		real(8), intent(in), contiguous :: face_normal(:, :), face_area(:),&
										   cell_volume(:), face_weight(:)
		real(8), intent(in), contiguous :: Field(:, :)
		real(8), intent(inout), contiguous :: Cell_grad(:)
		
		integer :: i, d, j, left_cell, right_cell, b_idx, face_idx
		integer :: pos1, pos2, offset
		real(8) :: w, avg_face_value(nvars), val, inv_vol
		
		Cell_grad = 0.d0
		
		pos1 = cell_faces_ptr(cell_idx)
		pos2 = cell_faces_ptr(cell_idx + 1)
		
		do i = pos1, pos2 - 1
			face_idx = cell_faces(i)
			left_cell = face_left_cell(face_idx)
			right_cell = face_right_cell(face_idx)
			
			if (right_cell <= ncells .or. USE_GHOST_CELLS) then
				w = face_weight(face_idx)
				
				avg_face_value = w*Field(:, right_cell) + (1.d0 - w)*Field(:, left_cell)
				
				if (left_cell == cell_idx) then
					do j = 1, nvars
						offset = (j-1)*dim
						do d = 1, dim
							val = avg_face_value(j)*face_normal(d, face_idx)*face_area(face_idx)
							Cell_grad(d + offset) = Cell_grad(d + offset) + val
						end do
					end do
				else	
					do j = 1, nvars
						offset = (j-1)*dim
						do d = 1, dim
							val = avg_face_value(j)*face_normal(d, face_idx)*face_area(face_idx)
							Cell_grad(d + offset) = Cell_grad(d + offset) - val
						end do
					end do
				end if
							
			else				
				avg_face_value = Field(:, right_cell)

				do j = 1, nvars
					offset = (j-1)*dim
					do d = 1, dim
						val = avg_face_value(j)*face_normal(d, face_idx)*face_area(face_idx)
						Cell_grad(d + offset) = Cell_grad(d + offset) + val
					end do
				end do
			end if
		end do
		
		inv_vol = 1.d0/cell_volume(cell_idx)
		Cell_grad(:) = Cell_grad(:)*inv_vol
	end subroutine

!=======================================================================
!============ LEAST SQUARES GRADIENT COMPUTATION SUBROUTINE ============
!=======================================================================
	!=== FIELD SUBS: ===
	pure subroutine precompute_lsq_weights(dim, ncells, cell_faces_ptr, cell_faces,&
                                           face_left_cell, face_right_cell,&
                                           cell_center, face_center, lsq_w)
        implicit none
        integer, intent(in) :: dim, ncells
        integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
                                           face_left_cell(:), face_right_cell(:)
        real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :)
        real(8), intent(inout), contiguous :: lsq_w(:, :) 				! (dim, size(cell_faces))
        
        real(8) :: A(dim, dim), invA(dim, dim), delta_r(dim, 50)
        real(8) :: center_j(3), dr(3), weights(50)
        real(8) :: r2, w, det, invDet, eps, val
        integer :: i, m, n, cell_idx, face_idx, neighbor, pos1, pos2, nfaces
        
        eps = 1.0d-16

        do cell_idx = 1, ncells
            A = 0.0d0
            pos1 = cell_faces_ptr(cell_idx)
            pos2 = cell_faces_ptr(cell_idx+1)-1
            nfaces = pos2 - pos1 + 1
            
            do i = 1, nfaces
                face_idx = cell_faces(pos1 + i - 1)
                
                if (face_left_cell(face_idx) == cell_idx) then
                    neighbor = face_right_cell(face_idx)
                else
                    neighbor = face_left_cell(face_idx)
                end if
                
                if (neighbor > ncells) then
                    center_j(1:dim) = face_center(1:dim, face_idx)
                else
                    center_j(1:dim) = cell_center(1:dim, neighbor)
                end if
            
                dr(1:dim) = center_j(1:dim) - cell_center(1:dim, cell_idx)
                r2 = dot_product(dr(1:dim), dr(1:dim)) + eps
                w = 1.0d0/r2
                
                delta_r(1:dim, i) = dr(1:dim)
                weights(i) = w
                
                do n = 1, dim
                    do m = 1, dim
                        A(m, n) = A(m, n) + w*dr(m)*dr(n)
                    end do
                end do
            end do
            
            invA = 0.0d0
            if (dim == 2) then
                det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
                if (abs(det) > eps) then
                    invDet = 1.0d0/det
                    invA(1,1) =  A(2,2)*invDet
                    invA(1,2) = -A(1,2)*invDet
                    invA(2,1) = -A(2,1)*invDet
                    invA(2,2) =  A(1,1)*invDet
                end if
            else if (dim == 3) then
                det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - &
                      A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + &
                      A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
                if (abs(det) > eps) then
                    invDet = 1.0d0/det
                    invA(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2))*invDet
                    invA(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3))*invDet
                    invA(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2))*invDet
                    invA(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3))*invDet
                    invA(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1))*invDet
                    invA(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3))*invDet
                    invA(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1))*invDet
                    invA(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2))*invDet
                    invA(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))*invDet
                end if
            end if
            
            do i = 1, nfaces
                w = weights(i)
                do m = 1, dim
                    val = 0.0d0
                    do n = 1, dim
                        val = val + invA(m, n)*(w*delta_r(n, i))
                    end do
                    lsq_w(m, pos1 + i - 1) = val
                end do
            end do
        end do
    end subroutine
	
	pure subroutine compute_grad_scalarfield_lsq(Field, Field_grad,&
											     ncells, dim, cell_faces_ptr, cell_faces,&
                                                 face_left_cell, face_right_cell, lsq_w)
        implicit none
        integer, intent(in) :: ncells, dim
        integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
										   face_left_cell(:), face_right_cell(:)
        real(8), intent(in), contiguous :: Field(:)
        real(8), intent(in), contiguous :: lsq_w(:, :) 					! (dim, size(cell_faces))
        real(8), intent(inout), contiguous :: Field_grad(:, :)

        integer :: cell_idx, f, face_idx, neighbor, pos1, pos2
        real(8) :: dphi, weight(3)

        Field_grad = 0.0d0

        do cell_idx = 1, ncells
            pos1 = cell_faces_ptr(cell_idx)
            pos2 = cell_faces_ptr(cell_idx+1) - 1

            do f = pos1, pos2
                face_idx = cell_faces(f)
                
                if (face_left_cell(face_idx) == cell_idx) then
                    neighbor = face_right_cell(face_idx)
                else
                    neighbor = face_left_cell(face_idx)
                end if

                weight(1:dim) = lsq_w(1:dim, f)

				dphi = Field(neighbor) - Field(cell_idx)
				Field_grad(:, cell_idx) = Field_grad(:, cell_idx) + weight(1:dim)*dphi
				
            end do
        end do
    end subroutine
    
	pure subroutine compute_grad_vectorfield_lsq(Field, Field_grad,&
											     ncells, nvars, dim,&
											     cell_faces_ptr, cell_faces,&
                                                 face_left_cell, face_right_cell, lsq_w)
        implicit none
        integer, intent(in)  :: ncells, nvars, dim
        integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
										   face_left_cell(:), face_right_cell(:)
        real(8), intent(in), contiguous :: Field(:, :)
        real(8), intent(in), contiguous :: lsq_w(:, :) 					! (dim, size(cell_faces))
        real(8), intent(inout), contiguous :: Field_grad(:, :)

        integer :: cell_idx, f, face_idx, neighbor, v, d, pos1, pos2, offset
        real(8) :: dphi, weight(3)

        Field_grad = 0.0d0

        do cell_idx = 1, ncells
            pos1 = cell_faces_ptr(cell_idx)
            pos2 = cell_faces_ptr(cell_idx+1) - 1

            do f = pos1, pos2
                face_idx = cell_faces(f)
                
                if (face_left_cell(face_idx) == cell_idx) then
                    neighbor = face_right_cell(face_idx)
                else
                    neighbor = face_left_cell(face_idx)
                end if

                weight(1:dim) = lsq_w(1:dim, f)

                do v = 1, nvars
                    dphi = Field(v, neighbor) - Field(v, cell_idx)

                    offset = (v-1)*dim
                    do d = 1, dim
                        Field_grad(d + offset, cell_idx) = Field_grad(d + offset, cell_idx) + weight(d)*dphi
                    end do
                end do
            end do
        end do
    end subroutine
    
    !=== CELL SUBS: ===
	pure subroutine compute_grad_scalarcell_lsq(Field, Cell_grad, cell_idx, ncells, dim,&
												cell_faces_ptr, cell_faces,& 
												face_left_cell, face_right_cell,&
												cell_center, face_center)
		implicit none
		integer, intent(in) :: cell_idx, ncells, dim
		integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
										   face_left_cell(:), face_right_cell(:)
		real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :)
		real(8), intent(in), contiguous :: Field(:)
		real(8), intent(inout), contiguous :: Cell_grad(:)
		
		real(8), parameter :: epsilon = 1.0d-12
		integer :: pos1, pos2
		integer :: i, m, n, n_points
		integer :: nfaces, face_idx, neighbor
		
		real(8) :: center_k(dim), center_j(dim), delta(dim), value_diff
		real(8) :: A(dim, dim), b(dim), det, weight, r
		logical :: is_boundary
				
		pos1 = cell_faces_ptr(cell_idx)
		pos2 = cell_faces_ptr(cell_idx + 1) - 1
		nfaces = pos2 - pos1 + 1
		center_k = cell_center(:, cell_idx)
		
		A = 0.0d0
		b = 0.0d0
		n_points = 0
		Cell_grad = 0.d0
		do i = 1, nfaces
			face_idx = cell_faces(pos1 + i - 1)
			
			if (face_left_cell(face_idx) == cell_idx) then
				neighbor = face_right_cell(face_idx)
			else
				neighbor = face_left_cell(face_idx)
			end if
			
			is_boundary = (neighbor > ncells)
			
			if (is_boundary) then
				center_j = face_center(:, face_idx)
			else
				center_j = cell_center(:, neighbor)
			end if
			
			value_diff = Field(neighbor) - Field(cell_idx)
			delta = center_j - center_k
			r = norm2(delta)
			
			if (r > epsilon) then
				weight = 1.0d0/(r*r)
				
				do n = 1, dim
					b(n) = b(n) + weight*delta(n)*value_diff
					do m = 1, n
						A(n, m) = A(n, m) + weight*delta(n)*delta(m)
					end do
				end do
				
				n_points = n_points + 1
			end if
		end do
		
		do n = 1, dim
			do m = n+1, dim
				A(n, m) = A(m, n)
			end do
		end do
		
		select case (dim)
			case (2)
				det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
				if (abs(det) > epsilon .and. n_points >= 2) then
					Cell_grad(1) = (A(2,2)*b(1) - A(1,2)*b(2))/det
					Cell_grad(2) = (-A(2,1)*b(1) + A(1,1)*b(2))/det
				end if
				
			case (3)
				det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
					- A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
					+ A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
				
				if (abs(det) > epsilon .and. n_points >= 3) then
					Cell_grad(1) = (b(1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) +&
								  A(1,2)*(A(2,3)*b(3) - b(2)*A(3,3)) +&
								  A(1,3)*(b(2)*A(3,2) - A(2,2)*b(3)))/det
					
					Cell_grad(2) = (A(1,1)*(b(2)*A(3,3) - A(2,3)*b(3)) +&
								  b(1)*(A(2,3)*A(3,1) - A(2,1)*A(3,3)) +&
								  A(1,3)*(A(2,1)*b(3) - b(2)*A(3,1)))/det
					
					Cell_grad(3) = (A(1,1)*(A(2,2)*b(3) - b(2)*A(3,2)) +&
								  A(1,2)*(b(2)*A(3,1) - A(2,1)*b(3)) +&
								  b(1)*(A(2,1)*A(3,2) - A(2,2)*A(3,1)))/det
				end if
		end select
	end subroutine

	pure subroutine compute_grad_vectorcell_lsq(Field, Cell_grad, cell_idx, ncells, dim, nvars,&
                                                cell_faces_ptr, cell_faces,& 
                                                face_left_cell, face_right_cell,&
                                                cell_center, face_center)
        implicit none
        integer, intent(in) :: cell_idx, ncells, dim, nvars
        integer, intent(in), contiguous :: cell_faces_ptr(:), cell_faces(:),&
                                           face_left_cell(:), face_right_cell(:)
        real(8), intent(in), contiguous :: cell_center(:, :), face_center(:, :)
        real(8), intent(in), contiguous :: Field(:, :)
        real(8), intent(inout), contiguous :: Cell_grad(:)
        
        real(8), parameter :: eps = 1.0d-14
        integer :: pos1, pos2, i, m, n, face_idx, neighbor, v, offset
        real(8) :: center_k(3), center_j(3), delta(3), r2, weight, det, invDet, val
        real(8) :: A(3, 3), invA(3, 3), b(3, nvars), dphi(nvars)

        Cell_grad = 0.0d0
        A = 0.0d0
        b = 0.0d0
        center_k(1:dim) = cell_center(1:dim, cell_idx)
        
        pos1 = cell_faces_ptr(cell_idx)
        pos2 = cell_faces_ptr(cell_idx + 1) - 1
        
        do i = pos1, pos2
            face_idx = cell_faces(i)
            
            if (face_left_cell(face_idx) == cell_idx) then
                neighbor = face_right_cell(face_idx)
            else
                neighbor = face_left_cell(face_idx)
            end if
            
            if (neighbor > ncells) then
                center_j(1:dim) = face_center(1:dim, face_idx)
            else
                center_j(1:dim) = cell_center(1:dim, neighbor)
            end if
            
            
            do v = 1, nvars
				dphi(v) = Field(v, neighbor) - Field(v, cell_idx)
			end do
            delta(1:dim) = center_j(1:dim) - center_k(1:dim)
            r2 = dot_product(delta(1:dim), delta(1:dim)) + eps
            weight = 1.0d0/r2
            
            do n = 1, dim
                do m = 1, dim
                    A(m, n) = A(m, n) + weight*delta(m)*delta(n)
                end do
                do v = 1, nvars
                    b(n, v) = b(n, v) + weight*delta(n)*dphi(v)
                end do
            end do
        end do
        
        invA = 0.0d0
        det = 0.0d0
        if (dim == 2) then
            det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
            if (abs(det) > eps) then
                invDet = 1.0d0/det
                invA(1,1) =  A(2,2)*invDet
                invA(1,2) = -A(1,2)*invDet
                invA(2,1) = -A(2,1)*invDet
                invA(2,2) =  A(1,1)*invDet
            end if
        else if (dim == 3) then
            det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - &
                  A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + &
                  A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
            if (abs(det) > eps) then
                invDet = 1.0d0/det
                invA(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2))*invDet
                invA(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3))*invDet
                invA(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2))*invDet
                invA(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3))*invDet
                invA(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1))*invDet
                invA(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3))*invDet
                invA(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1))*invDet
                invA(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2))*invDet
                invA(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1))*invDet
            end if
        end if
        
        if (abs(det) > eps) then
            do v = 1, nvars
                offset = (v-1)*dim
                do m = 1, dim
                    val = 0.0d0
                    do n = 1, dim
                        val = val + invA(m, n)*b(n, v)
                    end do
                    Cell_grad(m + offset) = val
                end do
            end do
        end if
    end subroutine









end module
