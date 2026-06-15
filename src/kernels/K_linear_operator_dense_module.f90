module linear_operator_dense_module
implicit none
contains
!=======================================================================
!=================== DENSE MATRIX LINEAR OPERATIONS ====================
!=======================================================================
	pure subroutine dense_mtrx_matvec(n, values, x, y)
		implicit none
		integer, intent(in) :: n
		real(8), intent(in), contiguous :: values(:, :), x(:)
		real(8), intent(inout), contiguous :: y(:)
		
		integer :: r, c
				
		y = 0.d0
		do c = 1, n
			do r = 1, n
				y(r) = y(r) + values(r, c)*x(c)
			end do
		end do
		
	end subroutine
	
	pure subroutine dense_mtrx_matmat(n, values1, values2, values3)
		implicit none
		integer, intent(in) :: n
		real(8), intent(in), contiguous :: values1(:, :), values2(:, :)
		real(8), intent(inout), contiguous :: values3(:, :)
		
		real(8) :: tmp
		integer :: r, c, k
		
		values3 = 0.d0
		do c = 1, n
			do k = 1, n
				tmp = values2(k, c)
				do r = 1, n
					values3(r, c) = values3(r, c) + values1(r, k)*tmp
				end do
			end do
		end do
	end subroutine
	
	pure subroutine dense_mtrx_get_diagonal(n, values, diag)
		implicit none
		integer, intent(in) :: n
		real(8), intent(in), contiguous :: values(:, :)
		real(8), intent(inout), contiguous :: diag(:)
		
		integer :: i
		
		do i = 1, n
			diag(i) = values(i, i)
		end do
	end subroutine
	
	pure subroutine invert_3x3(A, A_inv, threshold)
		implicit none
		real(8), intent(in) :: threshold
		real(8), intent(in) :: A(3,3)
		real(8), intent(inout) :: A_inv(3,3)
		
		real(8) :: det, inv_det
		real(8) :: m11, m12, m13, m21, m22, m23, m31, m32, m33
		
		!=== MATRIX ELEMENTS ===
		m11 = A(1,1); m12 = A(1,2); m13 = A(1,3)
		m21 = A(2,1); m22 = A(2,2); m23 = A(2,3)
		m31 = A(3,1); m32 = A(3,2); m33 = A(3,3)
		
		!=== MATRIX DETERMINANT ===
		det = m11*(m22*m33 - m23*m32) -&
			  m12*(m21*m33 - m23*m31) +&
			  m13*(m21*m32 - m22*m31)
		
		if (abs(det) < threshold) then
			A_inv = 0.0d0
			if (abs(m11) > threshold) A_inv(1,1) = 1.0d0/m11
			if (abs(m22) > threshold) A_inv(2,2) = 1.0d0/m22
			if (abs(m33) > threshold) A_inv(3,3) = 1.0d0/m33
			return
		end if
		
		inv_det = 1.0d0/det
		
		A_inv(1,1) =  (m22*m33 - m23*m32)*inv_det
		A_inv(1,2) = -(m12*m33 - m13*m32)*inv_det
		A_inv(1,3) =  (m12*m23 - m13*m22)*inv_det
		
		A_inv(2,1) = -(m21*m33 - m23*m31)*inv_det
		A_inv(2,2) =  (m11*m33 - m13*m31)*inv_det
		A_inv(2,3) = -(m11*m23 - m13*m21)*inv_det
		
		A_inv(3,1) =  (m21*m32 - m22*m31)*inv_det
		A_inv(3,2) = -(m11*m32 - m12*m31)*inv_det
		A_inv(3,3) =  (m11*m22 - m12*m21)*inv_det
	end subroutine
	
	pure subroutine invert_4x4(A, A_inv, threshold)
		implicit none
		real(8), intent(in) :: threshold
		real(8), intent(in) :: A(4,4)
		real(8), intent(inout) :: A_inv(4,4)
		
		real(8) :: det, inv_det
		real(8) :: m11, m12, m13, m14, m21, m22, m23, m24
		real(8) :: m31, m32, m33, m34, m41, m42, m43, m44
		real(8) :: t11, t12, t13, t14, t21, t22, t23, t24
		real(8) :: t31, t32, t33, t34, t41, t42, t43, t44

		m11 = A(1,1); m12 = A(1,2); m13 = A(1,3); m14 = A(1,4)
		m21 = A(2,1); m22 = A(2,2); m23 = A(2,3); m24 = A(2,4)
		m31 = A(3,1); m32 = A(3,2); m33 = A(3,3); m34 = A(3,4)
		m41 = A(4,1); m42 = A(4,2); m43 = A(4,3); m44 = A(4,4)

		t11 = m22*(m33*m44 - m34*m43) - m23*(m32*m44 - m34*m42) + m24*(m32*m43 - m33*m42)
		t12 = m21*(m33*m44 - m34*m43) - m23*(m31*m44 - m34*m41) + m24*(m31*m43 - m33*m41)
		t13 = m21*(m32*m44 - m34*m42) - m22*(m31*m44 - m34*m41) + m24*(m31*m42 - m32*m41)
		t14 = m21*(m32*m43 - m33*m42) - m22*(m31*m43 - m33*m41) + m23*(m31*m42 - m32*m41)

		det = m11*t11 - m12*t12 + m13*t13 - m14*t14

		if (abs(det) < threshold) then
			A_inv = 0.0d0
			if (abs(m11) > threshold) A_inv(1,1) = 1.0d0/m11
			if (abs(m22) > threshold) A_inv(2,2) = 1.0d0/m22
			if (abs(m33) > threshold) A_inv(3,3) = 1.0d0/m33
			if (abs(m44) > threshold) A_inv(4,4) = 1.0d0/m44
			return
		end if

		inv_det = 1.0d0/det

		t11 = m22*(m33*m44 - m34*m43) - m23*(m32*m44 - m34*m42) + m24*(m32*m43 - m33*m42)
		t12 = -(m21*(m33*m44 - m34*m43) - m23*(m31*m44 - m34*m41) + m24*(m31*m43 - m33*m41))
		t13 = m21*(m32*m44 - m34*m42) - m22*(m31*m44 - m34*m41) + m24*(m31*m42 - m32*m41)
		t14 = -(m21*(m32*m43 - m33*m42) - m22*(m31*m43 - m33*m41) + m23*(m31*m42 - m32*m41))

		t21 = -(m12*(m33*m44 - m34*m43) - m13*(m32*m44 - m34*m42) + m14*(m32*m43 - m33*m42))
		t22 = m11*(m33*m44 - m34*m43) - m13*(m31*m44 - m34*m41) + m14*(m31*m43 - m33*m41)
		t23 = -(m11*(m32*m44 - m34*m42) - m12*(m31*m44 - m34*m41) + m14*(m31*m42 - m32*m41))
		t24 = m11*(m32*m43 - m33*m42) - m12*(m31*m43 - m33*m41) + m13*(m31*m42 - m32*m41)

		t31 = m12*(m23*m44 - m24*m43) - m13*(m22*m44 - m24*m42) + m14*(m22*m43 - m23*m42)
		t32 = -(m11*(m23*m44 - m24*m43) - m13*(m21*m44 - m24*m41) + m14*(m21*m43 - m23*m41))
		t33 = m11*(m22*m44 - m24*m42) - m12*(m21*m44 - m24*m41) + m14*(m21*m42 - m22*m41)
		t34 = -(m11*(m22*m43 - m23*m42) - m12*(m21*m43 - m23*m41) + m13*(m21*m42 - m22*m41))

		t41 = -(m12*(m23*m34 - m24*m33) - m13*(m22*m34 - m24*m32) + m14*(m22*m33 - m23*m32))
		t42 = m11*(m23*m34 - m24*m33) - m13*(m21*m34 - m24*m31) + m14*(m21*m33 - m23*m31)
		t43 = -(m11*(m22*m34 - m24*m32) - m12*(m21*m34 - m24*m31) + m14*(m21*m32 - m22*m31))
		t44 = m11*(m22*m33 - m23*m32) - m12*(m21*m33 - m23*m31) + m13*(m21*m32 - m22*m31)

		A_inv(1,1) =  t11*inv_det
		A_inv(1,2) =  t21*inv_det
		A_inv(1,3) =  t31*inv_det
		A_inv(1,4) =  t41*inv_det

		A_inv(2,1) =  t12*inv_det
		A_inv(2,2) =  t22*inv_det
		A_inv(2,3) =  t32*inv_det
		A_inv(2,4) =  t42*inv_det

		A_inv(3,1) =  t13*inv_det
		A_inv(3,2) =  t23*inv_det
		A_inv(3,3) =  t33*inv_det
		A_inv(3,4) =  t43*inv_det

		A_inv(4,1) =  t14*inv_det
		A_inv(4,2) =  t24*inv_det
		A_inv(4,3) =  t34*inv_det
		A_inv(4,4) =  t44*inv_det

	end subroutine
	
	pure subroutine invert_general(n, A, A_inv, threshold)
		implicit none
		integer, intent(in) :: n
		real(8), intent(in) :: threshold
		real(8), intent(in) :: A(n,n)
		real(8), intent(inout) :: A_inv(n,n)
		
		real(8) :: Aug(n, 2*n), temp
		integer :: i, j, k, pivot
		
		Aug(:, 1:n) = A
		Aug(:, n+1:2*n) = 0.0d0
		do i = 1, n
			Aug(i, n+i) = 1.0d0
		end do
		
		do k = 1, n
			pivot = k
			do i = k+1, n
				if (abs(Aug(i,k)) > abs(Aug(pivot,k))) then
					pivot = i
				end if
			end do
			
			if (abs(Aug(pivot,k)) < threshold) then
				A_inv = 0.0d0
				do i = 1, n
					if (abs(A(i,i)) > threshold) then
						A_inv(i,i) = 1.0d0/A(i,i)
					end if
				end do
				return
			end if
			
			if (pivot /= k) then
				do j = k, 2*n
					temp = Aug(k,j)
					Aug(k,j) = Aug(pivot,j)
					Aug(pivot,j) = temp
				end do
			end if
			
			temp = Aug(k,k)
			do j = k, 2*n
				Aug(k,j) = Aug(k,j)/temp
			end do
			
			do i = 1, n
				if (i /= k) then
					temp = Aug(i,k)
					do j = k, 2*n
						Aug(i,j) = Aug(i,j) - temp*Aug(k,j)
					end do
				end if
			end do
		end do
		
		A_inv = Aug(:, n+1:2*n)
	end subroutine

end module
