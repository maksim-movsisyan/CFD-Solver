module preconditioner_block_jacobi_module
use linear_operator_dense_module
implicit none
contains
!=======================================================================
!============= BLOCK JACOBI PRECONDITIONER SUBROUTINES =================
!=======================================================================
	pure subroutine apply_block_jacobi_preconditioner(n, bs, invDIAG, r, z)
		implicit none
		integer, intent(in) :: n, bs
        real(8), intent(in), contiguous :: invDIAG(:, :, :), r(:)
        real(8), intent(inout), contiguous :: z(:)
        
        integer :: i, y_ptr, c_idx, r_idx
        real(8) :: val
              
		do i = 1, n
			y_ptr = (i-1)*bs
			z(1 + y_ptr:bs + y_ptr) = 0.d0
			do c_idx = 1, bs
				val = r(y_ptr + c_idx)
				do r_idx = 1, bs
					z(y_ptr + r_idx) = z(y_ptr + r_idx) + invDIAG(r_idx, c_idx, i)*val
				end do
			end do
		end do
    end subroutine
	
	pure subroutine update_block_jacobi_preconditioner1(n, bs, invDIAG)
		!=== invDIAG INITIALLY CONTAINS NON INVERSED DIAGONAL BLOCKS ===
		!=== NO OTHER INPUTS NEEDED									 ===
		implicit none
		real(8), parameter :: threshold = 1e-12
        integer, intent(in) :: n, bs
        real(8), intent(inout), contiguous :: invDIAG(:, :, :)
		
		real(8) :: tmp_block(bs, bs)
		integer :: i
        
        if (bs == 3) then
			do i = 1, n
				call invert_3x3(invDIAG(:, :, i), tmp_block, threshold)
				invDIAG(:, :, i) = tmp_block
			end do
			
        else if (bs == 4) then
			do i = 1, n
				call invert_4x4(invDIAG(:, :, i), tmp_block, threshold)
				invDIAG(:, :, i) = tmp_block
			end do
			
        else
			do i = 1, n
				call invert_general(bs, invDIAG(:, :, i), tmp_block, threshold)
				invDIAG(:, :, i) = tmp_block
			end do
			
        end if
    end subroutine
    
    pure subroutine update_block_jacobi_preconditioner2(n, bs, values, diag_indices, invDIAG)
		!=== invDIAG INITIALLY NON INITIALIZER						 ===
		!=== NEED BCSR MATRIX FORMAT DATA FOR INITIALIZATION		 ===
		implicit none
		real(8), parameter :: threshold = 1e-12
        integer, intent(in) :: n, bs
        integer, intent(in), contiguous :: diag_indices(:)							
        real(8), intent(in), contiguous :: values(:, :, :)
        real(8), intent(inout), contiguous :: invDIAG(:, :, :)
		
		integer :: i, diag_pos
                
        if (bs == 3) then
			do i = 1, n
				diag_pos = diag_indices(i)
				call invert_3x3(values(:, :, diag_pos), invDIAG(:, :, i), threshold)
			end do
        else if (bs == 4) then
			do i = 1, n
				diag_pos = diag_indices(i)
				call invert_4x4(values(:, :, diag_pos), invDIAG(:, :, i), threshold)
			end do
        else
			do i = 1, n
				diag_pos = diag_indices(i)
				call invert_general(bs, values(:, :, diag_pos), invDIAG(:, :, i), threshold)
			end do
        end if
    end subroutine
    
end module
