module string_utils_module
implicit none
contains
!=======================================================================
!================== PRECOND PARAMS PARSING SUBROUTINE ==================
!=======================================================================
subroutine get_key_value_real(params, key, value, found)
    character(len=*), intent(in) :: params  							! "omega=0.8; threshold=1e-10"
    character(len=*), intent(in) :: key     							! "omega"
    real(8), intent(out) :: value
    logical, intent(out) :: found
    
    integer :: p1, p2, p_eq
    character(len=len_trim(params)) :: temp_params
    character(len=len_trim(key)) :: temp_key
    
    found = .false.
    if (len_trim(params) == 0) return
    

    temp_params = adjustl(params)
    temp_key = trim(adjustl(key))
    
    p1 = index(temp_params, temp_key // "=")
    if (p1 == 0) return
    
    p_eq = p1 + len(temp_key)
    
    p2 = index(temp_params(p_eq+1:), ";")
    if (p2 == 0) then
        p2 = len_trim(temp_params)
    else
        p2 = p_eq + p2
    end if
    
    read(temp_params(p_eq+1:p2), *, iostat=p1) value
    if (p1 == 0) found = .true.
end subroutine

!=======================================================================
!=================== ADVANCED STRIP/TRIM FUNCTION ======================
!=======================================================================
function strip(str) result(res)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: res
    integer :: i, j, n

    n = len_trim(str)
    i = 1
    do while (i <= n .and. is_whitespace(str(i:i)))
        i = i + 1
    end do
    j = n
    do while (j >= 1 .and. is_whitespace(str(j:j)))
        j = j - 1
    end do
    if (i > j) then
        res = ''
    else
        res = str(i:j)
    end if
end function

function is_whitespace(ch) result(w)
    character(len=1), intent(in) :: ch
    logical :: w
    integer :: code
    code = iachar(ch)
    w = (code == 32) .or. (code == 9) .or. (code == 10) .or. (code == 13) .or. (code == 160)
end function

!=======================================================================
!=================== MESH PARSING AUX SUBROUTINES ======================
!=======================================================================
	integer function hex_to_int(hex_str) result(int_val)
		implicit none
		character(len=*), intent(in) :: hex_str
		character(len=len(hex_str)) :: temp_str
		integer :: i, digit, base, strlen
		
		int_val = 0
		temp_str = trim(adjustl(hex_str))
		strlen = len_trim(temp_str)
		
		do i = 1, strlen
			select case (temp_str(i:i))
			case ('0':'9')
				digit = ichar(temp_str(i:i)) - ichar('0')
			case ('a':'f')
				digit = ichar(temp_str(i:i)) - ichar('a') + 10
			case ('A':'F')  
				digit = ichar(temp_str(i:i)) - ichar('A') + 10
			case default
				digit = 0
			end select
			int_val = int_val*16 + digit
		end do
	end function

	function int_to_hex(n) result(hex_str)
		integer, intent(in) :: n
		character(len=:), allocatable :: hex_str
		character(len=100) :: temp

		write(temp, '(Z0)') n
		hex_str = trim(temp)
	end function
 
	integer function get_num_tokens(line) result(num_tokens)
		implicit none
		character(len=*), intent(in) :: line
		integer :: i, n
		logical :: in_token
		
		num_tokens = 0
		in_token = .false.
		n = len_trim(line)
		
		if (n == 0) return
		
		do i = 1, n
			if (line(i:i) /= ' ') then
				if (.not. in_token) then
					num_tokens = num_tokens + 1
					in_token = .true.
				end if
			else
				in_token = .false.
			end if
		end do
	end function

	subroutine get_token(line, idx, token)
		implicit none
		character(len=*), intent(in) :: line
		integer, intent(in) :: idx
		character(len=:), allocatable, intent(out) :: token
		
		integer :: i, n, token_count, start_pos, end_pos
		logical :: in_token, found
		
		token = ''
		n = len_trim(line)
		
		if (idx <= 0 .or. n == 0) return
		
		in_token = .false.
		token_count = 0
		start_pos = 0
		end_pos = 0
		found = .false.
		
		do i = 1, n
			if (line(i:i) /= ' ') then
				if (.not. in_token) then
					token_count = token_count + 1
					start_pos = i
					in_token = .true.
					
					if (token_count == idx) then
						found = .true.
					end if
				end if
			else
				if (in_token) then
					end_pos = i - 1
					in_token = .false.
					
					if (found) then
						exit
					end if
				end if
			end if
		end do
		
		if (in_token .and. token_count == idx) then
			end_pos = n
			found = .true.
		end if
		
		if (found .and. start_pos > 0 .and. end_pos >= start_pos) then
			token = line(start_pos:end_pos)
		end if
	end subroutine
		
end module
