module preconditioner_csr_ilu0_module
implicit none
contains
!=======================================================================
!=============== CSR ILU(0)	 PRECONDITIONER SUBROUTINES ================
!=======================================================================
	pure subroutine apply_csr_ilu0_preconditioner(n, values, col_indices,&
												  diag_indices, row_ptr, r, z)
		implicit none
		integer, intent(in) :: n
		integer, intent(in), contiguous :: col_indices(:), diag_indices(:), row_ptr(:)
        real(8), intent(in), contiguous :: values(:), r(:)
        real(8), intent(inout), contiguous :: z(:)
        
        integer :: i, j, j1, j2
        real(8) :: res
		
		!=== FORWARD SWEEP: ===
		do i = 1, n
			j1 = row_ptr(i)          
			j2 = diag_indices(i) - 1
			res = r(i)              
			do j = j1, j2
				res = res - values(j)*z(col_indices(j))
			end do
			z(i) = res 
		end do
		
		!=== BACKWARD SWEEP: ===
		do i = n, 1, -1           
			j1 = diag_indices(i) + 1      
			j2 = row_ptr(i + 1) - 1  
			res = z(i)  
			do j = j1, j2
				res = res - values(j)*z(col_indices(j))
			end do
			z(i) = res*values(diag_indices(i))
		end do
    end subroutine
	
	pure subroutine update_csr_ilu0_preconditioner(n, values, col_indices,&
												   diag_indices, row_ptr, iw)
		!=== values INITIALLY CONTAINS MATRIX VALUES 				 ===
		!=== NO OTHER INPUTS NEEDED									 ===
		implicit none
		integer, intent(in) :: n
		integer, intent(in), contiguous :: col_indices(:), diag_indices(:), row_ptr(:)
		integer, intent(inout), contiguous :: iw(:)
        real(8), intent(inout), contiguous :: values(:)
				
		integer :: i, j, k, j1, j2, jrow, jj, jw
		real(8) :: tl
		     
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
				
				tl = values(j)*values(diag_indices(jrow))
				values(j) = tl
				
				do jj = diag_indices(jrow) + 1, row_ptr(jrow+1) - 1
					jw = iw(col_indices(jj))
					if (jw /= 0) values(jw) = values(jw) - tl*values(jj)
				end do
				j = j + 1
			end do
			
			values(j) = 1.0d0/values(j)
			
			do i = j1, j2
				iw(col_indices(i)) = 0
			end do
		end do
			
    end subroutine
    
end module
