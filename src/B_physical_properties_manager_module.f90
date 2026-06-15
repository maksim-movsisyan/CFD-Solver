module physical_properties_manager_module
use string_utils_module
implicit none
!=======================================================================
!============== PHYSICAL PROPERIES MANAGER DATA STRUCTURE ==============
!=======================================================================
type, public :: phys_prop_manager_t 
    real(8) :: k = 1.4d0
    real(8) :: R_gas = 287.0d0
    real(8) :: Pr = 0.71d0
    real(8) :: cp, cv
contains
    procedure :: initialize => ppm_read_input_file
end type

contains
!=======================================================================
!================== PHYSICAL PROPERIES MANAGER METHODS =================
!=======================================================================
	subroutine ppm_read_input_file(this, filename)
		implicit none
		class(phys_prop_manager_t), intent(inout) :: this
		character(len=*), intent(in) :: filename
		
		integer :: iunit, ios, pos
        character(len=128) :: data_line, line, key, value
		
		
		open(newunit=iunit, file=filename, status='old', action='read')
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            
            line = adjustl(line)
            if (line(1:1) == '!' .or. line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            data_line = line
            pos = index(line, '!')
            if (pos > 0) data_line = line(:pos-1)
            
            pos = index(line, '#')
            if (pos > 0) data_line = line(:pos-1)
            !read(data_line, *) key, value
            
            pos = index(data_line, '=')
            if (pos <= 0) cycle
            
            key = adjustl(data_line(:pos-1))
            key = strip(key)
            
            value = adjustl(data_line(pos+1:))
            value = strip(value)
            
            
            select case (trim(adjustl(key)))
            case ('PR', 'Pr', 'Prandtl', 'PRANDTL')
                read(value, *) this%Pr  
			case ('k', 'K', 'gamma', 'GAMMA')
				read(value, *) this%k
			case ('R_gas', 'R_GAS', 'R0')
				read(value, *) this%R_gas
            end select
        end do
        close(iunit)
        
        
        this%cp = this%k*this%R_gas/(this%k - 1.d0)
        this%cv = this%cp/this%k
        
	end subroutine

end module
