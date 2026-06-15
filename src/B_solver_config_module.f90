module solver_config_module
use string_utils_module
implicit none

!=======================================================================
!===================== SOLVER CFG FILE DATA STRUCTURE ==================
!=======================================================================
type, abstract :: solver_cfg_t

	contains
	procedure(solver_cfg_read_interface), deferred :: read_cfg
end type

abstract interface
	subroutine solver_cfg_read_interface(this, filename)
		import solver_cfg_t
		class(solver_cfg_t), intent(inout) :: this
		character(len=*), intent(in) :: filename
	end subroutine
end interface


!=======================================================================
!================= DENSITY-BASED CFG FILE DATA STRUCTURE ===============
!=======================================================================
type, extends(solver_cfg_t) :: db_solver_config_t	
	character(len=128) :: GRADIENT_TYPE = 'GG'							!=== GRADIENTS TYPE ===
	real(8) :: CFL = 1.d0												!=== COURANT NUMBER ===
	
	integer :: NITER = 1000												!=== NUMBER OF NONLINEAR ITERATIONS ===
	real(8) :: TOLERANCE = 1e-8											!=== NON LINEAR ITERATIONS TOLERANCE ===
	
	integer :: MODEL = 1												!=== MODEL: EULER/NAVIER-STOKES ===
	integer :: SCHEME = 1												!=== SCHEME ===
	integer :: ORDER = 1												!=== ORDER ===
	integer :: LIMITER = 1												!=== TVD LIMITER ===
	
	integer :: REIMAN_SOLVER_TYPE = 1									!=== REIMAN SOLVER TYPE ===
	integer :: STABILITY_OPERATOR_TYPE = 1								!=== STABILITY OPERATOR TYPE ===
	
	integer :: FOS_NITER = 0											!=== FIRST ORDER SMOOTHING NUMBER OF ITERATION ===
	real(8) :: FOS_CFL = 0.001d0										!=== FIRST ORDER SMOOTHING COURANT NUMBER ===
	
	integer :: INITIALIZATION_TYPE = 1									!=== INITIALIZATION TYPE ===
	
	integer :: SAVE_NITER = 1000										!=== SOLUTION SAVE FREQUENCY ===
	integer :: MAS_NITER = 1											!=== MATRIX ASSEMBLING STEP ===
	integer :: PAS_NITER = 1											!=== PRECONDITIONER ASSEMBLING STEP ===
	
	logical :: IS_PRINT = .true.										!=== PRING LOG FLAG ===			
	integer :: PRINT_NITER = 1											!=== PRING LOG FREQUENCY ===
	logical :: IS_LOG = .true.											!=== WRITE LOG FLAG ===	
	integer :: LOG_NITER = 1											!=== WRITE LOG FREQUENCY ===
	
	logical :: USE_GHOST_CELLS = .false.								!=== GHOST CELLS FLAG ===
	logical :: USE_CONSERVATIVE_VARS = .false.							!=== CONSERVATIVER VARIABLES FLAG ===
	
	contains
	procedure :: read_cfg => db_read_solver_cfg
end type

!=======================================================================
!=================== MULTIGRID CFG FILE DATA STRUCTURE =================
!=======================================================================
type, extends(solver_cfg_t) :: multigrid_config_t	
	logical :: output_mesh = .false.									!=== OUTPUT MESH FLAG ===
	integer :: target_size = 2											!=== TARGET SIZE OF AGGLOMERATORS ===
	integer :: nlevels = 1												!=== NUMBER OF MG LEVELS ===
	integer :: nu1 = 1, nu2 = 1											!=== NUMBER OF PRE/POST SMOOTHING ITERATIONS ===
	real(8) :: omega = 1.0d0											!=== RELAXATION FACTOR FOR SMOOTHER ===
	integer :: max_iter = 10											!=== MAX ITERATIONS FOR THE COARSET LEVEL SLAE SOLUTIUON ===
	integer :: p_order = 1												!=== ORDER OF PROLONGATION OPERATOR ===
	
	
	contains
	procedure :: read_cfg => mg_read_solver_cfg
end type

contains
!=======================================================================
!===================== DENSITY-BASED CFG FILE METHODS ==================
!=======================================================================
	subroutine db_read_solver_cfg(this, filename)
		implicit none
		class(db_solver_config_t), intent(inout) :: this
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
            case ('GRADIENT_TYPE', 'gradient_type')
                this%GRADIENT_TYPE = trim(adjustl(value))
            case ('CFL', 'cfl')
				read(value, *) this%CFL                    
            case ('NITER', 'niter')
                read(value, *) this%NITER 
            case ('TOLERANCE', 'tolerance')
                read(value, *) this%TOLERANCE 
			case ('SCHEME', 'scheme')
                read(value, *) this%SCHEME 
            case ('MODEL', 'model')
                read(value, *) this%MODEL 
            case ('ORDER', 'order')
                read(value, *) this%ORDER 
            case ('LIMITER', 'limiter')
                read(value, *) this%LIMITER     
            case ('REIMAN_SOLVER_TYPE', 'reiman_solver_type')
                read(value, *) this%REIMAN_SOLVER_TYPE 
            case ('STABILITY_OPERATOR_TYPE', 'stability_operator_type')
                read(value, *) this%STABILITY_OPERATOR_TYPE  
            case ('FOS_NITER', 'fos_niter')
                read(value, *) this%FOS_NITER
            case ('FOS_CFL', 'fos_cfl')
                read(value, *) this%FOS_CFL
            case ('IS_PRINT', 'is_print')
                read(value, *) this%IS_PRINT
            case ('PRINT_NITER', 'print_niter')
                read(value, *) this%PRINT_NITER
            case ('IS_LOG', 'is_log')
                read(value, *) this%IS_LOG
            case ('LOG_NITER', 'log_niter')
                read(value, *) this%LOG_NITER
            case ('INITIALIZATION_TYPE', 'initialization_type')
                read(value, *) this%INITIALIZATION_TYPE
            case ('SAVE_NITER', 'save_niter')
                read(value, *) this%SAVE_NITER
            case ('MAS_NITER', 'mas_niter')
                read(value, *) this%MAS_NITER
            case ('PAS_NITER', 'pas_niter')
                read(value, *) this%PAS_NITER
            case ('USE_GHOST_CELLS', 'use_ghost_cells')
				read(value, *) this%USE_GHOST_CELLS
			case ('USE_CONSERVATIVE_VARS', 'use_conservative_vars')
				read(value, *) this%USE_CONSERVATIVE_VARS
            end select
        end do
        
        close(iunit)
    end subroutine


!=======================================================================
!===================== MULTIGRID CFG FILE METHODS ======================
!=======================================================================
	subroutine mg_read_solver_cfg(this, filename)
		implicit none
		class(multigrid_config_t), intent(inout) :: this
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
            case ('TARGET_SIZE', 'target_size')
                read(value, *) this%target_size
			case ('nu1', 'NU1', 'Nu1')
                read(value, *) this%nu1
			case ('nu2', 'NU2', 'Nu2')
                read(value, *) this%nu2
			case ('MAX_ITER', 'max_iter', 'Max_Iter')
                read(value, *) this%max_iter
			case ('NLEVELS', 'nlevels', 'n_levels', 'N_LEVELS')
                read(value, *) this%nlevels
			case ('OMEGA', 'omega', 'Omega', 'OmEgA')
                read(value, *) this%omega
            case ('OUTPUT_MESH', 'output_mesh')
				read(value, *) this%output_mesh
			case ('P_ORDER', 'p_order', 'PROLONGATION_ORDER', 'prolongation_order')
				read(value, *) this%p_order
            end select
        end do
        
        close(iunit)
    end subroutine












end module
