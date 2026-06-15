module bc_manager_module
use bc_base_module
use bc_wall_module
use bc_symm_module
use bc_inlet_module
use bc_outlet_module
use field_manager_module, only: field_manager_t
use physical_properties_manager_module, only: phys_prop_manager_t
use mesh_module, only: mesh_t
implicit none

!=======================================================================
!=============== BOUNDARY CONDITION MANAGER DATA STRUCTURE =============
!=======================================================================
type :: bc_container_t
	class(bc_base_t), allocatable :: bc
end type
    
    
type :: bc_manager_t
	integer :: num_bc = 0
	type(bc_container_t), allocatable :: bc_list(:)
	
	type(field_manager_t), pointer :: fldsm => null()
	type(phys_prop_manager_t), pointer :: ppm => null()
	
	logical :: initialized = .false.
	integer :: verbosity = -1
	
	integer, allocatable :: zone_to_index(:)							!=== mesh%face_zone(i) => j: local bc index ===
	integer :: max_zone_id = 0
contains
	procedure :: initialize => bc_manager_initialize
	procedure :: read_config_file => bc_manager_read_config
	procedure :: get_ref_values => bc_manager_get_ref_values
	procedure :: update_boundary_values => bc_manager_update_boundary_values
	procedure :: ADupdate_boundary_values => bc_manager_ADupdate_boundary_values
	procedure :: add_boundary_jacobians => bc_manager_add_boundary_jacobians
end type



!=== AUX PARSING DATA STRUCTRURE ===
type :: config_line_t
	character(len=256) :: raw_line
	integer :: zone_id
	character(len=64) :: bc_type_name
	integer :: bc_type_id
	character(len=64) :: bc_subtype_name
	integer :: bc_subtype_id
	character(len=128) :: args_line
end type
    
contains

!=======================================================================
!================= BOUNDARY CONDITION MANAGER METHODS ==================
!=======================================================================
	subroutine bc_manager_initialize(this, mesh, fldsm, ppm, config_file)
		implicit none
        class(bc_manager_t), intent(inout) :: this
        type(mesh_t), target, intent(in) :: mesh
        type(field_manager_t), target, intent(in) :: fldsm
        type(phys_prop_manager_t), target, intent(in) :: ppm
        character(len=*), intent(in) :: config_file
        
        integer :: i, max_zone
        
        if (this%verbosity > 0) then
			print*, ' '
			print*, ' Initializing boundary conditions manager:'
        end if
		
		this%fldsm => fldsm
		this%ppm => ppm

        !=== READING INPUT FILE ===
        call this%read_config_file(config_file, mesh)
        
        max_zone = 0
        do i = 1, this%num_bc
            max_zone = max(max_zone, this%bc_list(i)%bc%zone_id)
        end do
        
        if (max_zone > 0) then
            allocate(this%zone_to_index(max_zone))
            this%zone_to_index = 0
            
            do i = 1, this%num_bc
                this%zone_to_index(this%bc_list(i)%bc%zone_id) = i
            end do
            this%max_zone_id = max_zone
        end if
        
        !=== BC INITIALIZATION ===
        if (this%verbosity > 0) then
			print*, ' BC initialization: '
		end if
		
        do i = 1, this%num_bc
            call this%bc_list(i)%bc%initialize(mesh)
        end do
        
        this%initialized = .true.
        if (this%verbosity > 0) then
			print*, ' Boundary conditions manager initialized with ', this%num_bc, ' boundary zones'
		end if
    end subroutine
    
    subroutine bc_manager_read_config(this, filename, mesh)
		implicit none
        class(bc_manager_t), intent(inout) :: this
        character(len=*), intent(in) :: filename
        type(mesh_t), target, intent(in) :: mesh
		
		integer :: i, iunit, ios, line_count, pos1, pos2
		character(len=256) :: line, trim_line
		type(config_line_t), allocatable :: lines(:)
		real(8), allocatable :: args(:)
		
		
		if (this%verbosity > 0) then
			print*, ' Reading inputfile: ',  filename
        end if
		
		!=== NUMBER OF BOUNDARY CONDITIONS ===
		open(newunit=iunit, file=trim(filename), status='old', action='read')
		line_count = 0
		do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line_count = line_count + 1
        end do
        rewind(iunit)
		
		if (this%verbosity > 0) then
			print*, ' Inputfile contains ',  line_count, ' lines.'
			print*, ' BC info:'
        end if
		
		!=== LINES PARSING ===
		this%num_bc = line_count
		allocate(lines(line_count))
		
		do i = 1, line_count
            read(iunit, '(A)') lines(i)%raw_line
            
            read(lines(i)%raw_line, *) lines(i)%zone_id
									   
			pos1=index(lines(i)%raw_line, '(')
			pos2=index(lines(i)%raw_line, ')')
			lines(i)%bc_type_name = lines(i)%raw_line(pos1+1:pos2-1)
			lines(i)%bc_type_id = bc_type_from_string(lines(i)%bc_type_name)
			
			
			pos1=index(lines(i)%raw_line, '{')
			pos2=index(lines(i)%raw_line, '}')
			if (pos1 > 0 .and. pos2 > pos1) then
				lines(i)%bc_subtype_name = lines(i)%raw_line(pos1+1:pos2-1)
				lines(i)%bc_subtype_id = bc_subtype_from_string(lines(i)%bc_type_id, lines(i)%bc_subtype_name)
			else
				lines(i)%bc_subtype_name='None'
				lines(i)%bc_subtype_id = 0
			end if
			
			pos1=index(lines(i)%raw_line, '[')
			pos2=index(lines(i)%raw_line, ']')
			if (pos1 > 0 .and. pos2 > pos1) then
				lines(i)%args_line = lines(i)%raw_line(pos1+1:pos2-1)
			else
				lines(i)%args_line=''
			end if
			
			
			
			if (this%verbosity > 0) then
				print*, i, ')', ' ZONE ID: ', lines(i)%zone_id,&
						' BC TYPE NAME: ', lines(i)%bc_type_name,&
						' BC TYPE ID: ', lines(i)%bc_type_id,&
						' BC SUBTYPE NAME: ', lines(i)%bc_subtype_name,&
						' BC SUBTYPE ID: ', lines(i)%bc_subtype_id,&
						' ARGS LINE : ', lines(i)%args_line
			end if
			
        end do
		
		!=== BC ALLOCATION ===
		allocate(this%bc_list(line_count))
		
		do i = 1, line_count
			if (len_trim(lines(i)%args_line) > 0) then
				call parse_arguments(lines(i)%args_line, args)
			else
				allocate(args(0))
			end if
			
			select case(lines(i)%bc_type_id)
			case(BC_TYPE_WALL)
				allocate(bc_wall_t :: this%bc_list(i)%bc)
                select type(bc => this%bc_list(i)%bc)
                type is (bc_wall_t)
                    call setup_wall_bc(bc, lines(i), mesh, args, this%verbosity)
                end select
			
			case(BC_TYPE_SYMMETRY)
				allocate(bc_symm_t :: this%bc_list(i)%bc)
                select type(bc => this%bc_list(i)%bc)
                type is (bc_symm_t)
                    call setup_symm_bc(bc, lines(i), mesh, args, this%verbosity)
                end select
                
			case(BC_TYPE_INLET)
				allocate(bc_inlet_t :: this%bc_list(i)%bc)
                select type(bc => this%bc_list(i)%bc)
                type is (bc_inlet_t)
                    call setup_inlet_bc(bc, lines(i), mesh, args, this%verbosity)
                end select
				
			case(BC_TYPE_OUTLET)
				allocate(bc_outlet_t :: this%bc_list(i)%bc)
                select type(bc => this%bc_list(i)%bc)
                type is (bc_outlet_t)
                    call setup_outlet_bc(bc, lines(i), mesh, args, this%verbosity)
                end select
				
			case default
				error stop ' Unknown boundary condition type: '
				
			end select
			
			if (allocated(args)) deallocate(args)
        end do
		
		deallocate(lines)
    end subroutine
    
    subroutine parse_arguments(line, args)
		implicit none
        character(len=*), intent(in) :: line
        real(8), allocatable, intent(out) :: args(:)
        
        character(len=256) :: temp_line, token
        integer :: arg_count, i, pos
        
        temp_line = trim(adjustl(line))
        if (len_trim(temp_line) == 0) then
            allocate(args(0))
            return
        end if
        
        
        arg_count = 1
        do i = 1, len_trim(temp_line)
            if (temp_line(i:i) == ' ') arg_count = arg_count + 1
        end do
        
        allocate(args(arg_count))
        
		read(temp_line, *) args
    end subroutine
    
    !================= CONVERTER STRING TO BC IDNEX ====================
    function bc_type_from_string(type_str) result(bc_type)
		implicit none
        character(len=*), intent(in) :: type_str
        integer :: bc_type
        
        select case(trim(adjustl(type_str)))
        case('wall', 'WALL', 'Wall')
            bc_type = BC_TYPE_WALL
        case('inlet', 'INLET', 'Inlet')
            bc_type = BC_TYPE_INLET
        case('outlet', 'OUTLET', 'Outlet')
            bc_type = BC_TYPE_OUTLET
        case('symmetry', 'SYMMETRY', 'Symmetry')
            bc_type = BC_TYPE_SYMMETRY
        case('farfield', 'FARFIELD', 'Farfield')
            bc_type = BC_TYPE_FARFIELD
        case('pressure_inlet', 'PRESSURE_INLET')
            bc_type = BC_TYPE_PRESSURE_INLET
        case('pressure_outlet', 'PRESSURE_OUTLET')
            bc_type = BC_TYPE_PRESSURE_OUTLET
        case('mass_flow_inlet', 'MASS_FLOW_INLET')
            bc_type = BC_TYPE_MASS_FLOW_INLET
        case('periodic', 'PERIODIC', 'Periodic')
            bc_type = BC_TYPE_PERIODIC
        case default
            error stop ' Unknown boundary condition type string: ' // trim(type_str)
        end select
    end function
    
    function bc_subtype_from_string(bc_type, subtype_str) result(bc_subtype)
		implicit none
        integer, intent(in) :: bc_type
        character(len=*), intent(in) :: subtype_str
        integer :: bc_subtype
        
        select case(bc_type)
        case(BC_TYPE_WALL)
            select case(trim(adjustl(subtype_str)))
            case('no_slip', 'NO_SLIP', 'no-slip', 'NOSLIP', 'noslip')
                bc_subtype = BC_WALL_NO_SLIP
            case('slip', 'SLIP', 'Slip')
                bc_subtype = BC_WALL_SLIP
            case('slip-reiman', 'SLIP_REIMAN', 'Slip-Reiman', 'SLIP-REIMAN', 'Slip_Reiman', 'slip_reiman')
                bc_subtype = BC_WALL_SLIP_REIMAN	
            case('adiabatic', 'ADIABATIC', 'Adiabatic')
                bc_subtype = BC_WALL_ADIABATIC
			case('heat_flux', 'HEAT_FLUX', 'Heat_Flux', 'HEATFLUX', 'HEAT-FLUX', 'heatflux', 'heat-flux')
                bc_subtype = BC_WALL_HEAT_FLUX
            case('isothermal', 'ISOTHERMAL', 'Isothermal')
                bc_subtype = BC_WALL_ISOTHERMAL
            case default
                bc_subtype = BC_WALL_NO_SLIP
            end select
            
        case(BC_TYPE_INLET)
            select case(trim(adjustl(subtype_str)))
            case('subsonic', 'SUBSONIC', 'Subsonic')
                bc_subtype = BC_INLET_SUBSONIC
            case('supersonic', 'SUPERSONIC', 'Supersonic')
                bc_subtype = BC_INLET_SUPERSONIC
            case('total', 'TOTAL', 'Total')
                bc_subtype = BC_INLET_TOTAL
            case default
                bc_subtype = BC_INLET_SUBSONIC
            end select
            
        case(BC_TYPE_OUTLET)
            select case(trim(adjustl(subtype_str)))
            case('subsonic', 'SUBSONIC', 'Subsonic')
                bc_subtype = BC_OUTLET_SUBSONIC
            case('supersonic', 'SUPERSONIC', 'Supersonic')
                bc_subtype = BC_OUTLET_SUPERSONIC
            case default
                bc_subtype = BC_OUTLET_SUBSONIC
            end select
            
        case default
            bc_subtype = 0
        end select
    end function
    
    !======================= BC SETUP SUBROUTINES ======================
    subroutine setup_wall_bc(bc, line, mesh, args, verbosity)
		implicit none
        type(bc_wall_t), intent(inout) :: bc
        type(config_line_t), intent(in) :: line
        type(mesh_t), target, intent(in) :: mesh
        real(8), intent(in) :: args(:)
        integer, intent(in) :: verbosity
        
		bc%verbosity = verbosity
        bc%zone_id = line%zone_id
        bc%bc_name = line%bc_type_name
        bc%bc_type = line%bc_type_id
        bc%bc_subtype = line%bc_subtype_id
        
        select case(line%bc_subtype_id)
        case(BC_WALL_NO_SLIP, BC_WALL_ADIABATIC)
            bc%wall_velocity(1:mesh%dim) = args(1:mesh%dim)
            
            if (verbosity > 0) then
				print*, '		ZONE ', bc%zone_id, ' WALL VELOCITY: ', bc%wall_velocity(1:mesh%dim)
            end if
            
        case(BC_WALL_ISOTHERMAL)
            bc%wall_velocity(1:mesh%dim) = args(1:mesh%dim)
            bc%wall_temperature = args(mesh%dim + 1)
            
            if (verbosity > 0) then
				print*, '		ZONE ', bc%zone_id, ' WALL VELOCITY: ', bc%wall_velocity(1:mesh%dim),&
				                                    ' WALL TEMPERATURE: ', bc%wall_temperature
            end if
            
        case(BC_WALL_HEAT_FLUX)
            bc%wall_velocity(1:mesh%dim) = args(1:mesh%dim)
            bc%wall_heat_flux = args(mesh%dim + 1)
            
            if (verbosity > 0) then
				print*, '		ZONE ', bc%zone_id, ' WALL VELOCITY: ', bc%wall_velocity(1:mesh%dim),&
				                                    ' WALL HEAT FLUX: ', bc%wall_heat_flux
            end if
            
        end select 
    end subroutine
    
    subroutine setup_symm_bc(bc, line, mesh, args, verbosity)
		implicit none
        type(bc_symm_t), intent(inout) :: bc
        type(config_line_t), intent(in) :: line
        type(mesh_t), target, intent(in) :: mesh
        real(8), intent(in) :: args(:)
        integer, intent(in) :: verbosity
        
		bc%verbosity = verbosity
        bc%zone_id = line%zone_id
        bc%bc_name = line%bc_type_name
        bc%bc_type = line%bc_type_id
        bc%bc_subtype = line%bc_subtype_id

    end subroutine
    
    subroutine setup_inlet_bc(bc, line, mesh, args, verbosity)
		implicit none
        type(bc_inlet_t), intent(inout) :: bc
        type(config_line_t), intent(in) :: line
        type(mesh_t), target, intent(in) :: mesh
        real(8), intent(in) :: args(:)
        integer, intent(in) :: verbosity
        
        bc%verbosity = verbosity
        bc%zone_id = line%zone_id
        bc%bc_name = line%bc_type_name
        bc%bc_type = line%bc_type_id
        bc%bc_subtype = line%bc_subtype_id
        
        select case(line%bc_subtype_id)
        case(BC_INLET_SUBSONIC)
            bc%inlet_velocity(1:mesh%dim) = args(1:mesh%dim)
            
            if (size(args) > mesh%dim) then
				bc%inlet_temperature = args(mesh%dim + 1)
            end if
            
            
            if (verbosity > 0 .and. size(args) > mesh%dim) then
				print*, '		ZONE ', bc%zone_id, ' INLET VELOCITY: ', bc%inlet_velocity(1:mesh%dim),&
											        ' INLET TEMPERATURE: ', bc%inlet_temperature
			else if (verbosity > 0) then
				print*, '		ZONE ', bc%zone_id, ' INLET VELOCITY: ', bc%inlet_velocity(1:mesh%dim)
            end if
            
        case(BC_INLET_SUPERSONIC)
			bc%inlet_pressure = args(1)
			bc%inlet_velocity(1:mesh%dim) = args(2:mesh%dim + 1)
			bc%inlet_temperature = args(mesh%dim + 2)
			
			if (verbosity > 0) then
				print*, '		ZONE ', bc%zone_id, ' INLET PRESSURE: ', bc%inlet_pressure,&
												    ' INLET VELOCITY: ', bc%inlet_velocity(1:mesh%dim),&
													' INLET TEMPERATURE: ', bc%inlet_temperature
            end if
			
        end select
    end subroutine

	subroutine setup_outlet_bc(bc, line, mesh, args, verbosity)
		implicit none
        type(bc_outlet_t), intent(inout) :: bc
        type(config_line_t), intent(in) :: line
        type(mesh_t), target, intent(in) :: mesh
        real(8), intent(in) :: args(:)
        integer, intent(in) :: verbosity 
        
        bc%verbosity = verbosity
        bc%zone_id = line%zone_id
        bc%bc_name = line%bc_type_name
        bc%bc_type = line%bc_type_id
        bc%bc_subtype = line%bc_subtype_id
        
        select case(line%bc_subtype_id)
        case(BC_OUTLET_SUBSONIC)
            bc%outlet_pressure = args(1)
            
            if (verbosity > 0) then
				print*, '		ZONE ', bc%zone_id, ' OUTLET PRESSURE: ', bc%outlet_pressure
            end if
            
        end select
                
    end subroutine
    
    
!=======================================================================
!================= BOUNDARY CONDITION MANAGER METHODS 2 ================
!=======================================================================
	subroutine bc_manager_update_boundary_values(this, name, update_grad, USE_GHOST_CELLS, fidx)
		implicit none
		class(bc_manager_t), intent(inout) :: this
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
		logical, intent(in) :: USE_GHOST_CELLS
		integer, optional, intent(in) :: fidx
		
		integer :: i, field_idx
		real(8), pointer, contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		
		if (present(fidx)) then
			field_idx = fidx
		else
			field_idx = this%fldsm%get_idx(name)
		end if
		
        values_ptr => this%fldsm%registry(field_idx)%values
        values_grad_ptr => this%fldsm%registry(field_idx)%grad
		
		if (USE_GHOST_CELLS) then
			do i = 1, this%num_bc
				call this%bc_list(i)%bc%apply_ghosts_patch(values_ptr, values_grad_ptr, name, update_grad)
			end do
		else
			do i = 1, this%num_bc
				call this%bc_list(i)%bc%apply_patch(this%ppm%R_gas, this%ppm%k, values_ptr, values_grad_ptr, name, update_grad)
			end do
		end if
		
	end subroutine
    
    subroutine bc_manager_ADupdate_boundary_values(this, name, update_grad, USE_GHOST_CELLS, fidx, dfidx)
		implicit none
		class(bc_manager_t), intent(inout) :: this
		character(len=*), intent(in) :: name
		logical, intent(in) :: update_grad
		logical, intent(in) :: USE_GHOST_CELLS
		integer, optional, intent(in) :: fidx, dfidx
		
		integer :: i, field_idx, dfield_idx
		real(8), pointer, contiguous :: values_ptr(:, :), values_grad_ptr(:, :)
		real(8), pointer, contiguous :: values_ptrd(:, :), values_grad_ptrd(:, :)
		
		if (present(fidx)) then
			field_idx = fidx
		else
			field_idx = this%fldsm%get_idx(name)
		end if
		
		if (present(dfidx)) then
			dfield_idx = dfidx
		else
			dfield_idx = this%fldsm%get_idx(trim(name//'_eps'))
		end if
		
        values_ptr => this%fldsm%registry(field_idx)%values
        values_grad_ptr => this%fldsm%registry(field_idx)%grad
        
        values_ptrd => this%fldsm%registry(dfield_idx)%values
        values_grad_ptrd => this%fldsm%registry(dfield_idx)%grad
		
		if (USE_GHOST_CELLS) then
			do i = 1, this%num_bc
				call this%bc_list(i)%bc%ADdiff_apply_ghosts_patch(values_ptr, values_ptrd,&
														   values_grad_ptr, values_grad_ptrd, name, update_grad)
			end do
		else
			do i = 1, this%num_bc
				call this%bc_list(i)%bc%ADdiff_apply_patch(this%ppm%R_gas, this%ppm%k, values_ptr, values_ptrd,&
														   values_grad_ptr, values_grad_ptrd, name, update_grad)
			end do
		end if
		
	end subroutine
    
    subroutine bc_manager_add_boundary_jacobians(this, values, map_LB, name, MODEL,&
												 USE_GHOST_CELLS, USE_CONSERVATIVE_VARS, fidx)
		implicit none
		class(bc_manager_t), intent(inout) :: this
		integer, intent(in), contiguous :: map_LB(:)
		integer, intent(in) :: MODEL
		character(len=*), intent(in) :: name
		logical, intent(in) :: USE_GHOST_CELLS, USE_CONSERVATIVE_VARS
		integer, optional, intent(in) :: fidx
		real(8), intent(inout), contiguous :: values(:, :, :)
		
		integer :: i, field_idx
		real(8), pointer, contiguous :: values_ptr(:, :)
		
		if (USE_GHOST_CELLS) return
		
		if (USE_CONSERVATIVE_VARS) then
			print*, 'WARNING IN <bc_manager_add_boundary_jacobians> USE_CONSERVATIVE_VARS FLAG TBD '
			return
		end if
		
		if (present(fidx)) then
			field_idx = fidx
		else
			field_idx = this%fldsm%get_idx(name)
		end if
		
        values_ptr => this%fldsm%registry(field_idx)%values
				
		do i = 1, this%num_bc
			call this%bc_list(i)%bc%apply_jmatrix_inv_patch(this%ppm%R_gas, this%ppm%k, values_ptr, values, map_LB, name)
			if (MODEL == 2) call this%bc_list(i)%bc%apply_jmatrix_visc_patch(this%ppm%R_gas, this%ppm%k, this%ppm%cp, this%ppm%Pr, values_ptr, values, map_LB, name)
		end do
		
	end subroutine
    
    subroutine bc_manager_get_ref_values(this, V_ref, P_ref, T_ref)
		implicit none
		class(bc_manager_t), intent(in) :: this
		real(8), intent(inout) :: V_ref(3), P_ref, T_ref
		
		integer :: i
		
		V_ref = 0.d0
		P_ref = 100000.d0
		T_ref = 300.d0
		
		do i = 1, this%num_bc
			select case(this%bc_list(i)%bc%bc_type)
			case(BC_TYPE_INLET)
				select case(this%bc_list(i)%bc%bc_subtype)
				case(BC_INLET_SUBSONIC)
					select type(bc => this%bc_list(i)%bc)
					type is (bc_inlet_t)
						V_ref = bc%inlet_velocity
						T_ref = bc%inlet_temperature
					end select
				case(BC_INLET_SUPERSONIC)
					select type(bc => this%bc_list(i)%bc)
					type is (bc_inlet_t)
						V_ref = bc%inlet_velocity
						T_ref = bc%inlet_temperature
						P_ref = bc%inlet_pressure
					end select
				end select
			
			case(BC_TYPE_OUTLET)
				select case(this%bc_list(i)%bc%bc_subtype)
				case(BC_OUTLET_SUBSONIC)
					select type(bc => this%bc_list(i)%bc)
					type is (bc_outlet_t)
						P_ref = bc%outlet_pressure
					end select
				end select
			end select
		end do
    end subroutine
    
end module
