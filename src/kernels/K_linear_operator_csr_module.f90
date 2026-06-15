module linear_operator_csr_module
implicit none
contains
!=======================================================================
!===================== CSR MATRIX LINEAR OPERATIONS ====================
!=======================================================================
	pure subroutine csr_mtrx_matvec(n, values, col_indices, row_ptr, x, y)
		implicit none
		integer, intent(in) :: n
		integer, intent(in), contiguous :: col_indices(:), row_ptr(:)
		real(8), intent(in), contiguous :: values(:), x(:)
		real(8), intent(inout), contiguous :: y(:)
		
		integer :: i, k, k1, k2
				
		do i = 1, n
			k1 = row_ptr(i)
			k2 = row_ptr(i + 1) - 1
			
			y(i) = 0.d0
			do k = k1, k2
				y(i) = y(i) + values(k)*x(col_indices(k))
			end do
		end do
	end subroutine
	
	pure subroutine csr_mtrx_get_diagonal(n, values, diag_indices, diag)
		implicit none
		integer, intent(in) :: n
		integer, intent(in), contiguous :: diag_indices(:)
		real(8), intent(in), contiguous :: values(:)
		real(8), intent(inout), contiguous :: diag(:)
		
		integer :: i
		
		do i = 1, n
			diag(i) = values(diag_indices(i))
		end do
	end subroutine
	
end module
