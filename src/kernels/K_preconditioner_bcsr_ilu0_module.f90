module preconditioner_bcsr_ilu0_module
use linear_operator_dense_module
implicit none
contains
!=======================================================================
!================ BCSR ILU(0) PRECONDITIONER SUBROUTINES ===============
!=======================================================================
	pure subroutine apply_bcsr_ilu0_preconditioner(n, bs, values, col_indices,&
											  diag_indices, row_ptr, r, z)
		implicit none
		integer, intent(in) :: n, bs
		integer, intent(in), contiguous :: col_indices(:), diag_indices(:), row_ptr(:)
        real(8), intent(in), contiguous :: values(:, :, :), r(:)
        real(8), intent(inout), contiguous :: z(:)
       
        integer :: i, j, r_idx, c_idx
        integer :: j_start, j_diag, j_end
        integer :: x_ptr, y_ptr
        real(8) :: temp_vec(bs), val
		
		do i = 1, n
			j_start = row_ptr(i)
			j_diag = diag_indices(i)
			y_ptr = (i-1)*bs
			
			z(1+y_ptr:bs+y_ptr) = r(1+y_ptr:bs+y_ptr)
			
			do j = j_start, j_diag - 1
				x_ptr = (col_indices(j)-1)*bs
				
				temp_vec = 0.d0
				do c_idx = 1, bs
					val = z(x_ptr + c_idx)
					do r_idx = 1, bs
						temp_vec(r_idx) = temp_vec(r_idx) + values(r_idx, c_idx, j)*val
					end do
				end do
				
				z(1+y_ptr:bs+y_ptr) = z(1+y_ptr:bs+y_ptr) - temp_vec
			end do
		end do
		
		do i = n, 1, -1 
			j_diag = diag_indices(i)
			j_end = row_ptr(i + 1) - 1
			y_ptr = (i-1)*bs
			    
			do j = j_diag + 1, j_end
				x_ptr = (col_indices(j)-1)*bs
				
				temp_vec = 0.d0
				do c_idx = 1, bs
					val = z(x_ptr + c_idx)
					do r_idx = 1, bs
						temp_vec(r_idx) = temp_vec(r_idx) + values(r_idx, c_idx, j)*val
					end do
				end do	
				
				z(1+y_ptr:bs+y_ptr) = z(1+y_ptr:bs+y_ptr) - temp_vec
			end do
			
			do r_idx = 1, bs
				temp_vec(r_idx) = z(y_ptr + r_idx)
				z(y_ptr + r_idx) = 0.d0
			end do
			
			do c_idx = 1, bs
				val = temp_vec(c_idx)
				do r_idx = 1, bs
					z(y_ptr + r_idx) = z(y_ptr + r_idx) + values(r_idx, c_idx, j_diag)*val
				end do
			end do				
		end do
    end subroutine
    
	pure subroutine update_bcsr_ilu0_preconditioner(n, bs, values, col_indices,&
												    diag_indices, row_ptr, iw)
		!=== values INITIALLY CONTAINS MATRIX VALUES 				 ===
		!=== NO OTHER INPUTS NEEDED									 ===
		implicit none
		real(8), parameter :: threshold = 1e-12
		integer, intent(in) :: n, bs
		integer, intent(in), contiguous :: col_indices(:), diag_indices(:), row_ptr(:)
		integer, intent(inout), contiguous :: iw(:)
        real(8), intent(inout), contiguous :: values(:, :, :)
				
		integer :: i, j, k, j1, j2, jrow, jj, jw
		integer :: c_idx, r_idx, l_idx
		real(8) :: tl, temp_block(bs, bs), work_block(bs, bs)
		    
        iw = 0 
		do k = 1, n
			j1 = row_ptr(k)
			j2 = row_ptr(k + 1) - 1
			
			do j = j1, j2
				iw(col_indices(j)) = j
			end do
			
			j = j1
			do while (j <= j2)
				jrow = col_indices(j)
				if (jrow >= k) exit  
				                
                temp_block = 0.d0
                do l_idx = 1, bs  
					do c_idx = 1, bs      
						do r_idx = 1, bs  
							temp_block(r_idx, l_idx) = temp_block(r_idx, l_idx) + &
								values(r_idx, c_idx, j)*values(c_idx, l_idx, diag_indices(jrow))
						end do
					end do
				end do
                
                values(:,:,j) =  temp_block(:, :)
                               
				
				do jj = diag_indices(jrow) + 1, row_ptr(jrow+1) - 1
					jw = iw(col_indices(jj))
					if (jw /= 0) then
					
						work_block = 0.d0
						do l_idx = 1, bs  
							do c_idx = 1, bs      
								do r_idx = 1, bs  
									work_block(r_idx, l_idx) = work_block(r_idx, l_idx) + &
										temp_block(r_idx, c_idx)*values(c_idx, l_idx, jj)
								end do
							end do
						end do
				
                        values(:,:,jw) = values(:,:,jw) - work_block(:, :)
                                       
					end if
				end do
				j = j + 1
			end do
			
			
			if (bs == 3) then
				call invert_3x3(values(:,:,j), temp_block, threshold)
			else if (bs == 4) then
				call invert_4x4(values(:,:,j), temp_block, threshold)
			else
				call invert_general(bs, values(:,:,j), temp_block, threshold)
			end if
			
			values(:,:,j) = temp_block
			
			do i = j1, j2
				iw(col_indices(i)) = 0
			end do
		end do
		
			
    end subroutine
    
end module
