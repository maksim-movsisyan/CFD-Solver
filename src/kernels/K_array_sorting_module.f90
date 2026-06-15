module array_sorting_module
implicit none
contains
!=======================================================================
!============ HYBRID SORT OF INTEGER ARRAY SUBROUITNE ==================
!=======================================================================
	pure subroutine hybrid_sort(arr, n)
		implicit none
		integer, intent(inout) :: arr(n)
		integer, intent(in) :: n
		
		integer, parameter :: THRESHOLD = 16
		
		if (n <= THRESHOLD) then
			call insertion_sort(arr, n)
		else
			call quick_sort(arr, 1, n)
		end if
	end subroutine

	pure subroutine insertion_sort(arr, n)
		implicit none
		integer, intent(inout) :: arr(n)
		integer, intent(in) :: n
		
		integer :: i, j, temp
		
		do i = 2, n
			temp = arr(i)
			j = i - 1
			do while (j >= 1)
				if (arr(j) <= temp) exit
				arr(j + 1) = arr(j)
				j = j - 1
			end do
			arr(j + 1) = temp
		end do
	end subroutine

	pure recursive subroutine quick_sort(arr, left, right)
		implicit none
		integer, intent(inout) :: arr(*)
		integer, intent(in) :: left, right
		integer :: i, j, pivot, temp
		
		if (left < right) then
			pivot = arr((left + right)/2)
			i = left - 1
			j = right + 1
			
			do
				do
					i = i + 1
					if (arr(i) >= pivot) exit
				end do
				do
					j = j - 1
					if (arr(j) <= pivot) exit
				end do
				if (i >= j) exit
				temp = arr(i)
				arr(i) = arr(j)
				arr(j) = temp
			end do
			
			call quick_sort(arr, left, j)
			call quick_sort(arr, j + 1, right)
		end if
	end subroutine
	
	pure recursive subroutine quick_sort_real(angles, indices, n)
		implicit none
		real(8), intent(inout) :: angles(:)
		integer, intent(inout) :: indices(:)
		integer, intent(in) :: n
		real(8) :: pivot
		integer :: i, j, temp_idx
		real(8) :: temp_angle
		
		if (n <= 1) return
		
		pivot = angles(n/2)
		i = 1
		j = n
		
		do
			do while (angles(i) < pivot)
				i = i + 1
			end do
			do while (angles(j) > pivot)
				j = j - 1
			end do
			if (i >= j) exit
			
			temp_angle = angles(i)
			angles(i) = angles(j)
			angles(j) = temp_angle
			
			temp_idx = indices(i)
			indices(i) = indices(j)
			indices(j) = temp_idx
			
			i = i + 1
			j = j - 1
		end do
	
		call quick_sort_real(angles(1:i-1), indices(1:i-1), i-1)
		call quick_sort_real(angles(j+1:n), indices(j+1:n), n-j)
	end subroutine
	
!=======================================================================
!============ BINARTY SEARCH IN INTEGER ARRAY SUBROUITNE ===============
!=======================================================================
	pure function binary_search(array, low_bound, high_bound, target_val) result(val)
		implicit none
		integer, intent(in) :: array(:)
		integer, intent(in) :: low_bound, high_bound, target_val
		integer :: val
		
		integer :: low, mid, high
		
		val = -1
		low = low_bound
		high = high_bound
		
		do while (low <= high)
			mid = low + (high - low)/2
			if (array(mid) == target_val) then
				val = mid
				return
			else if (array(mid) < target_val) then
				low = mid + 1
			else
				high = mid - 1
			end if
		end do
	end function

!=======================================================================
!========================= SWAP SUBROUTINES ============================
!=======================================================================
	pure subroutine swap_real(a, b)
		real(8), intent(inout) :: a, b
		real(8) :: temp
		temp = a
		a = b
		b = temp
	end subroutine

	pure subroutine swap_int(a, b)
		integer, intent(inout) :: a, b
		integer :: temp
		temp = a
		a = b
		b = temp
	end subroutine
	

end module
