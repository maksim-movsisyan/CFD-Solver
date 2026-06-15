module linear_operator_bcsr_module
implicit none
contains
!=======================================================================
!===================== BCSR MATRIX LINEAR OPERATIONS ===================
!=======================================================================	
	pure subroutine bcsr_mtrx_matvec(n, bs, values, col_indices, row_ptr, x, y)
		implicit none
		integer, intent(in) :: n, bs
		integer, intent(in), contiguous :: col_indices(:), row_ptr(:)
		real(8), intent(in), contiguous :: values(:, :, :), x(:)
		real(8), intent(inout), contiguous :: y(:)
		
		integer :: i, k, row, col, r, c
		integer :: b_start, b_end, x_ptr, y_ptr
		real(8) :: val
				
		do i = 1, n
			b_start = row_ptr(i)
			b_end = row_ptr(i+1) - 1
			y_ptr = (i-1)*bs
			
			y(y_ptr + 1:y_ptr + bs) = 0.0d0
			do k = b_start, b_end
				x_ptr = (col_indices(k)-1)*bs
				
				do c = 1, bs
					val = x(x_ptr + c)
					do r = 1, bs
						y(y_ptr + r) = y(y_ptr + r) + values(r, c, k)*val
					end do
				end do
			end do
		end do
	end subroutine
	
	pure subroutine bcsr_mtrx_get_diagonal(n, values, diag_indices, diag)
		implicit none
		integer, intent(in) :: n
		integer, intent(in), contiguous :: diag_indices(:)
		real(8), intent(in), contiguous :: values(:, :, :)
		real(8), intent(inout), contiguous :: diag(:, :, :)
		
		integer :: i
		
		do i = 1, n
			diag(:, :, i) = values(:, :, diag_indices(i))
		end do
	end subroutine
	
	pure subroutine bcsr_mtrx_get_main_diagonal(n, bs, values, diag_indices, diag)
		implicit none
		integer, intent(in) :: n, bs
		integer, intent(in), contiguous :: diag_indices(:)
		real(8), intent(in), contiguous :: values(:, :, :)
		real(8), intent(inout), contiguous :: diag(:)
		
		integer :: i, j, offset
		
		do i = 1, n
			offset = bs*(i - 1)
			
			do j = 1, bs
				diag(j + offset) = values(j, j, diag_indices(i))
			end do			
		end do
		
	end subroutine
	
end module
