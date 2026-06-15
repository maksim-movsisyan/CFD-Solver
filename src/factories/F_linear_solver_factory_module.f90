module linear_solver_factory_module
use linear_solver_module
use preconditioner_module
use string_utils_module
implicit none

!=======================================================================
!============== LINEAR SOLVER CFG FILE DATA STRUCTURE ==================
!=======================================================================
type :: linear_solver_config_t
	character(len=32) :: solver_type = 'GMRES'
	
	!=== GENERAL PARAMETERS ===
	integer :: max_iter = 1000
	real(8) :: rel_tol = 1.0d-3
	real(8) :: abs_tol = 1.0d-6
	integer :: verbosity = 0
	character(len=128) :: preconditioner_type = 'NONE'
	
	!=== GMRES PARAMETERS ===
	integer :: restart = 30
	logical :: use_no_restart = .false.
	
	!=== JACOBI PARAMETERS ===
	real(8) :: omega = 1.0d0
	
	contains
	procedure :: read_cfg => linear_sover_read_cfg
end type


!=======================================================================
!=============== LINEAR SOLVER FACTORY DATA STRUCTURE ==================
!=======================================================================
type :: linear_solver_factory_t
contains
	procedure :: create => factory_create_linear_solver
end type

type(linear_solver_factory_t) :: linear_solver_factory
contains
!=======================================================================
!==================== LINEAR SOLVER CFG FILE METHODS ===================
!=======================================================================
	subroutine linear_sover_read_cfg(this, filename)
		implicit none
		class(linear_solver_config_t), intent(inout) :: this
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
            case ('SOLVER_TYPE', 'solver_type')
                this%solver_type = trim(adjustl(value))
            case ('MAX_ITER', 'max_iter')
                read(value, *) this%max_iter
            case ('REL_TOL', 'rel_tol')
                read(value, *) this%rel_tol
            case ('ABS_TOL', 'abs_tol')
                read(value, *) this%abs_tol
            case ('VERBOSITY', 'verbosity')
                read(value, *) this%verbosity
            case ('RESTART', 'restart')
                read(value, *) this%restart
            case ('OMEGA', 'omega')
                read(value, *) this%omega
            case ('PRECONDITIONER_TYPE', 'preconditioner_type',&
				  'PRECONDITIONER', 'PRECOND', 'preconditioner', 'precond')
                this%preconditioner_type = trim(adjustl(value))
            case ('USE_NO_RETSTART', 'use_no_restart')
                read(value, *) this%use_no_restart
            end select
        end do
        
        close(iunit)
    end subroutine


!=======================================================================
!================== LINEAR SOLVER FACTORY METHODSa =====================
!=======================================================================
	function create_linear_solver(config) result(solver_ptr)
		implicit none
        type(linear_solver_config_t), intent(in) :: config
        class(linear_solver_t), pointer :: solver_ptr
        
        nullify(solver_ptr)
        
        select case (trim(adjustl(config%solver_type)))
        case ('GMRES', 'gmres')
            allocate(gmres_solver_t::solver_ptr)
			select type (solver_ptr)
			type is (gmres_solver_t)
				solver_ptr%restart = config%restart
				solver_ptr%rel_tol = config%rel_tol
				solver_ptr%abs_tol = config%abs_tol
				solver_ptr%max_iter = config%max_iter
				solver_ptr%verbosity = config%verbosity
			end select
            
        case ('BiCGSTAB', 'BICGSTAB', 'bicgstab')
            allocate(bicgstab_solver_t::solver_ptr)
            select type (solver_ptr)
            type is (bicgstab_solver_t)
                solver_ptr%rel_tol = config%rel_tol
                solver_ptr%abs_tol = config%abs_tol
                solver_ptr%max_iter = config%max_iter
                solver_ptr%verbosity = config%verbosity
            end select
            
        case ('JACOBI', 'JACOBI_SOLVER', 'Jacobi', 'jacobi')
            allocate(jacobi_solver_t::solver_ptr)
            select type (solver_ptr)
            type is (jacobi_solver_t)
                solver_ptr%omega = config%omega
                solver_ptr%rel_tol = config%rel_tol
                solver_ptr%abs_tol = config%abs_tol
                solver_ptr%max_iter = config%max_iter
                solver_ptr%verbosity = config%verbosity
            end select
            
        case default
            write(*,*) 'ERROR: Unknown linear solver type: ', trim(config%solver_type)
            error stop 'Linear solver factory error'
        end select
        
        
        if (trim(config%preconditioner_type) /= 'NONE') then
			call attach_preconditioner(solver_ptr, config%preconditioner_type)
		end if        
    end function
	
	subroutine attach_preconditioner(solver, prec_type)
		implicit none
		class(linear_solver_t), intent(inout) :: solver
		character(len=*), intent(in) :: prec_type
		
		class(preconditioner_t), pointer :: prec_ptr
		
		nullify(prec_ptr)
		
		select case (trim(prec_type))
		case ('CSR_ILU0', 'csr_ilu0', 'CSR-ILU0', 'csr-ilu0',&
			  'CSR_ILU(0)', 'csr_ilu(0)', 'CSR-ILU(0)', 'csr-ilu(0)')
			allocate(csr_ilu0_preconditioner_t :: prec_ptr)
			
		case ('BCSR_ILU0', 'bcsr_ilu0', 'BCSR-ILU0', 'bcsr-ilu0',&
			  'BCSR_ILU(0)', 'bcsr_ilu(0)', 'BCSR-ILU(0)', 'bcsr-ilu(0)')
			allocate(bcsr_ilu0_preconditioner_t :: prec_ptr)
				
		case('BLOCK-JACOBI', 'JACOBI', 'Block-Jacobi', 'block-jacobi',&
			 'BLOCK_JACOBI', 'jacobi', 'Block_Jacobi', 'block_jacobi')
			allocate(block_jacobi_preconditioner_t :: prec_ptr)
		
		case('BLUSGS', 'blusgs', 'BlockLUGSG', 'blockLUSGS',&
			 'BLOCKLUSGS')
			allocate(blusgs_preconditioner_t :: prec_ptr)
		
		case('GMG', 'gmg', 'GeometricMultiGrid', 'MultiGrid',&
			 'MULTIGRID')
			allocate(gmg_preconditioner_t :: prec_ptr)
					 
		case default
			write(*,*) 'Warning: Unknown preconditioner type: ', trim(prec_type)
			return
		end select
		
		call solver%set_preconditioner(prec_ptr)
	end subroutine
	
	function factory_create_linear_solver(this, lsover_cfg_filename) result(solver_ptr)
		implicit none
        class(linear_solver_factory_t), intent(in) :: this
        character(len=*), optional, intent(in) :: lsover_cfg_filename
        class(linear_solver_t), pointer :: solver_ptr
        
        character(len=128) :: filename
        type(linear_solver_config_t) :: config
        
        filename = 'data\input\linear_solver_settings.txt'
        if (present(lsover_cfg_filename)) filename = lsover_cfg_filename
        
        call config%read_cfg(filename)
		
        solver_ptr => create_linear_solver(config)
    end function
    
    
end module
