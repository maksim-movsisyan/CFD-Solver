module mg_operators_module
implicit none
contains
!=======================================================================
!============= INTERPOLATION: FINE -> COARSE SUBROUTINE ================
!=======================================================================
	!scalar field:
	subroutine mg_restriction_scalar_cell(cncells, ccell_volume, fcell_volume,&
										  fine_cells_ptr, fine_cells,&
										  fine_field, coarse_field)
		implicit none
		integer, intent(in) :: cncells
		integer, intent(in), contiguous :: fine_cells_ptr(:), fine_cells(:)
		real(8), intent(in), contiguous :: ccell_volume(:), fcell_volume(:)
		real(8), intent(in), contiguous :: fine_field(:)
		real(8), intent(inout), contiguous :: coarse_field(:)
		
		integer :: ccell_idx, fcell_idx
		integer :: pos1, pos2, pos
		real(8) :: weight, weighted_sum
		
		do ccell_idx = 1, cncells			
			pos1 = fine_cells_ptr(ccell_idx)
			pos2 = fine_cells_ptr(ccell_idx + 1) - 1
			
			weighted_sum = 0.d0
			do pos = pos1, pos2
				fcell_idx = fine_cells(pos)
				weight = fcell_volume(fcell_idx)
				weighted_sum = weighted_sum + weight*fine_field(fcell_idx)
			end do
			
			coarse_field(ccell_idx) = weighted_sum/ccell_volume(ccell_idx)
			
		end do

	end subroutine

	!vector field:
	subroutine mg_restriction_vector_cell(cncells, nvars,&
										  ccell_volume, fcell_volume,&
										  fine_cells_ptr, fine_cells,&
										  fine_field, coarse_field)
		implicit none
		integer, intent(in) :: cncells, nvars
		integer, intent(in), contiguous :: fine_cells_ptr(:), fine_cells(:)
		real(8), intent(in), contiguous :: ccell_volume(:), fcell_volume(:)
		real(8), intent(in), contiguous :: fine_field(:, :)
		real(8), intent(inout), contiguous :: coarse_field(:, :)
		
		integer :: v
		integer :: ccell_idx, fcell_idx
		integer :: pos1, pos2, pos
		real(8) :: weight, weighted_sum(nvars), inv_vol
		
		do ccell_idx = 1, cncells			
			pos1 = fine_cells_ptr(ccell_idx)
			pos2 = fine_cells_ptr(ccell_idx + 1) - 1
			
			weighted_sum = 0.d0
			do pos = pos1, pos2
				fcell_idx = fine_cells(pos)
				weight = fcell_volume(fcell_idx)
				do v = 1, nvars
					weighted_sum(v) = weighted_sum(v) + weight*fine_field(v, fcell_idx)
				end do
			end do
			
			inv_vol = 1.d0/ccell_volume(ccell_idx)
			do v = 1, nvars
				coarse_field(v, ccell_idx) = weighted_sum(v)*inv_vol
			end do
		end do

	end subroutine

!=======================================================================
!======== FIRST ORDER EXTRAPOLATION: COARSE -> FINE SUBROUTINE =========
!=======================================================================
	!scalar field:
	subroutine mg_prolongation1_scalar_cell(fncells, fine_to_coarse,&
											coarse_field, fine_field)
		implicit none
		integer, intent(in) :: fncells
		integer, intent(in), contiguous :: fine_to_coarse(:)
		real(8), intent(in), contiguous :: coarse_field(:)
		real(8), intent(inout), contiguous :: fine_field(:)
		
		integer :: fcell_idx, ccell_idx
		
		do fcell_idx = 1, fncells
			ccell_idx = fine_to_coarse(fcell_idx)
			fine_field(fcell_idx) = coarse_field(ccell_idx)
		end do
	end subroutine
	
	!vectors field:
	subroutine mg_prolongation1_vector_cell(fncells, nvars, fine_to_coarse,&
											coarse_field, fine_field)
		implicit none
		integer, intent(in) :: fncells, nvars
		integer, intent(in), contiguous :: fine_to_coarse(:)
		real(8), intent(in), contiguous :: coarse_field(:, :)
		real(8), intent(inout), contiguous :: fine_field(:, :)
		
		integer :: v
		integer :: fcell_idx, ccell_idx
		
		do fcell_idx = 1, fncells
			ccell_idx = fine_to_coarse(fcell_idx)
			
			do v = 1, nvars
				fine_field(v, fcell_idx) = coarse_field(v, ccell_idx)
			end do
			
		end do
	end subroutine
	
!=======================================================================
!======== SECOND ORDER EXTRAPOLATION: COARSE -> FINE SUBROUTINE ========
!=======================================================================
	!scalar field:
	subroutine mg_prolongation2_scalar_cell(dim, fncells, fine_to_coarse,&
											fcell_center, ccell_center,&
											coarse_field, coarse_gradient,&
											fine_field)
		implicit none
		integer, intent(in) :: dim, fncells
		integer, intent(in), contiguous :: fine_to_coarse(:)
		real(8), intent(in), contiguous :: fcell_center(:, :), ccell_center(:, :)
		real(8), intent(in), contiguous :: coarse_field(:), coarse_gradient(:, :)
		real(8), intent(inout), contiguous :: fine_field(:)
		
		integer :: fcell_idx, ccell_idx
		real(8) :: dist(dim), grad(dim)
		
		do fcell_idx = 1, fncells
			ccell_idx = fine_to_coarse(fcell_idx)
			dist = ccell_center(:, ccell_idx) - fcell_center(:, fcell_idx)
			grad = coarse_gradient(:, ccell_idx)
			
			
			fine_field(fcell_idx) = coarse_field(ccell_idx) + dot_product(grad, dist)
		end do
	end subroutine
	
	!vector field:
	subroutine mg_prolongation2_vector_cell(dim, fncells, nvars, fine_to_coarse,&
											fcell_center, ccell_center,&
											coarse_field, coarse_gradient,&
											fine_field)
		implicit none
		integer, intent(in) :: dim, fncells, nvars
		integer, intent(in), contiguous :: fine_to_coarse(:)
		real(8), intent(in), contiguous :: fcell_center(:, :), ccell_center(:, :)
		real(8), intent(in), contiguous :: coarse_field(:, :), coarse_gradient(:, :)
		real(8), intent(inout), contiguous :: fine_field(:, :)
		
		integer :: v, off
		integer :: fcell_idx, ccell_idx
		real(8) :: dist(dim), grad(dim)
		
		do fcell_idx = 1, fncells
			ccell_idx = fine_to_coarse(fcell_idx)
			dist = ccell_center(:, ccell_idx) - fcell_center(:, fcell_idx)
			do v = 1, nvars
				off = (v - 1)*dim
				grad = coarse_gradient(1+off:dim+off, ccell_idx)
				fine_field(v, fcell_idx) = coarse_field(v, ccell_idx) + dot_product(grad, dist)
			end do
		end do
	end subroutine
		
end module
