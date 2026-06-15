module output_fields_module
use string_utils_module
implicit none
!=======================================================================
!================= FLUENT .DAT VARIABLES INDEXATION ====================
!=======================================================================
integer, parameter :: VAR_PRESSURE = 1
integer, parameter :: VAR_TEMPERATURE = 3
integer, parameter :: VAR_DENSITY = 101
integer, parameter :: VAR_VELOCITY_X = 111
integer, parameter :: VAR_VELOCITY_Y = 112
integer, parameter :: VAR_VELOCITY_Z = 113
integer, parameter :: VAR_MACH = 7 


contains
!=======================================================================
!================= PARAVIEW .VTK OUTPUT SUBROUTINE =====================
!=======================================================================
	subroutine output_fields_vtk(fields, field_names, filename, title,&
								 dim, ncells, nfaces, nnodes,&
								 cell_nodes, cell_nodes_ptr, node_coords)
		implicit none
        real(8), intent(in) :: fields(:, :)
        character(len=32), intent(in) :: field_names(:)
        character(len=*), intent(in) :: filename, title
        integer, intent(in) :: dim, ncells, nfaces, nnodes
        integer, intent(in), contiguous :: cell_nodes(:), cell_nodes_ptr(:)
        real(8), intent(in), contiguous :: node_coords(:, :)
				
		integer :: cell_id, face_id, node_id
		integer :: j, d, total_connectivity_size
		integer :: pos1, pos2
		integer :: iunit, ios
				
        total_connectivity_size = 0
        do cell_id = 1, ncells
			pos1 = cell_nodes_ptr(cell_id)
			pos2 = cell_nodes_ptr(cell_id+1)
            total_connectivity_size = total_connectivity_size + (pos2 - pos1) + 1
        end do
       
        open(newunit=iunit, file=filename, status='replace', action='write') 
        write(iunit, '(a)') '# vtk DataFile Version 3.0'
        write(iunit, '(a)') trim(title)
        write(iunit, '(a)') 'ASCII'
        write(iunit, '(a)') 'DATASET UNSTRUCTURED_GRID'
        
        !=== MESH NODES: ===
        write(iunit, '(a,i0,a)') 'POINTS ', nnodes, ' double'
        if (dim == 2) then
			do node_id = 1, nnodes
				write(iunit, '(3e20.10)') node_coords(:, node_id), 0.0d0
			end do
		else
			do node_id = 1, nnodes
				write(iunit, '(3e20.10)') node_coords(:, node_id)
			end do
		end if
        write(iunit, '()')
        
        !=== MESH CELLS: ===
        write(iunit, '(a,i0,1x,i0)') 'CELLS ', ncells, total_connectivity_size
        do cell_id = 1, ncells
			pos1 = cell_nodes_ptr(cell_id)
			pos2 = cell_nodes_ptr(cell_id+1)
            write(iunit, '(i0)', advance='no') pos2 - pos1
            do j = pos1, pos2-1
                write(iunit, '(1x,i0)', advance='no') cell_nodes(j) - 1
            end do
            write(iunit, '()')
        end do
        write(iunit, '()')
        
        !=== MESH CELL TYPES: ===
        write(iunit, '(a,i0)') 'CELL_TYPES ', ncells
        do cell_id = 1, ncells
            write(iunit, '(i0)') 7
        end do
        write(iunit, '()')
        
        !=== CELL DATA: ===
        write(iunit, '(a,i0)') 'CELL_DATA ', ncells
        write(iunit, '(a)') 'SCALARS cell_id int 1'
        write(iunit, '(a)') 'LOOKUP_TABLE default'
        do cell_id = 1, ncells
            write(iunit, '(i0)') cell_id
        end do
        write(iunit, '()')
        
        !=== OUTPUT FIELDS: ===
        do j = 1, size(field_names)
			write(iunit, '(a)') 'SCALARS '//trim(field_names(j))//' double 1'
			write(iunit, '(a)') 'LOOKUP_TABLE default'
			do cell_id = 1, ncells
				write(iunit, '(f22.15)') fields(j, cell_id)
			end do
			write(iunit, '()')
        end do
        		
        close(iunit)
        
       
        print *, ' '
        print *, 'VTK file: ', filename
        print *, 'nodes: ', nnodes
        print *, 'cells: ', ncells
        print *, 'total connectivity size: ', total_connectivity_size
        print *, ' '
    end subroutine

!=======================================================================
!================== PARAVIEW .VTK INPUT SUBROUTINE =====================
!=======================================================================
	subroutine input_fields_vtk(fields, field_names, filename, dim, ncells)
		implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: dim, ncells
		real(8), intent(inout) :: fields(:, :)
        character(len=32), allocatable, intent(inout) :: field_names(:)
        		
		integer :: cell_id, j
		integer :: iunit, ios
		integer :: nvars
		character(len=128) :: line, aux_line, name
				       
        open(newunit=iunit, file=filename, status='old', action='read') 
        nvars = 0
        !=== READING FILE: ===
        do
			read(iunit, '(a)', iostat=ios) line
			if (ios /= 0) then
				exit
			end if
			
			if (index(line, 'SCALARS') > 0) then
				read(line, *) aux_line, name
				read(iunit, '(a)') aux_line 
				
				if (trim(name) == 'cell_id') cycle
				
				nvars = nvars + 1
				
				if (.not. allocated(field_names)) then
					field_names = [character(len=32) :: name]
				else
					field_names = [character(len=32) :: field_names, name]
				end if
					    
				do cell_id = 1, ncells
					read(iunit, *) fields(nvars, cell_id)
				end do
			end if			
		end do
	
        close(iunit)
        
    end subroutine




!=======================================================================
!=================== FLUENT .DAT OUTPUT SUBROUTINE =====================
!=======================================================================
	subroutine output_fields_dat(fields, field_names, filename, title,&
								 ncells, nfaces, nbfaces,&
								 face_zone, face_type,&
								 face_left_cell, face_right_cell,&
								 USE_GHOST_CELLS)			 
		implicit none
        real(8), intent(in) :: fields(:, :)
        character(len=32), intent(in) :: field_names(:)
        character(len=*), intent(in) :: filename, title
        integer, intent(in) :: ncells, nfaces, nbfaces
        integer, intent(in) :: face_zone(:), face_type(:),&
							   face_left_cell(:), face_right_cell(:)
		logical, intent(in) :: USE_GHOST_CELLS
		
		integer :: iunit, ios
		integer :: cell_idx, face_idx, zone_idx, field_idx, var_idx
		integer :: nzones, counter, current_zone
		integer :: zones_id(1:32), zones_type(1:32), zones_indexes(1:2, 1:32)
		
		!=== BOUNDARY ZONES INFORMATION: ===
		zones_id(:) = 0
		nzones = 0
		do face_idx = nfaces-nbfaces+1, nfaces
			current_zone = face_zone(face_idx)
			if (any(zones_id == current_zone)) cycle
			nzones = nzones + 1
			zones_id(nzones) = current_zone
			zones_type(nzones) = face_type(face_idx)
		end do
		
		zones_indexes(1, :) = huge(1); zones_indexes(2, :) = -huge(1)	
		do zone_idx = 1, nzones
			do face_idx = nfaces-nbfaces+1, nfaces
				current_zone = face_zone(face_idx)
				if (current_zone .ne. zones_id(zone_idx)) cycle
				
				zones_indexes(1, zone_idx) = min(zones_indexes(1, zone_idx), face_idx)
				zones_indexes(2, zone_idx) = max(zones_indexes(2, zone_idx), face_idx)
			end do
		end do
		
		!=== WRITING DATA: ===
		open(newunit=iunit, file=filename, status='replace', action='write')
        write(iunit, '(A)') '(0 "Fluent Data File")'
        write(iunit, '(A)') '(0 "' // title // '")'
        
        do field_idx = 1, size(field_names)
			var_idx = name_to_var(field_names(field_idx))
			
			!=== BOUNDARY FACES DATA: ===
			do zone_idx = 1, nzones
				call write_face_data(field_idx, zones_id(zone_idx), var_idx, zones_indexes(1, zone_idx), zones_indexes(2, zone_idx))
			end do
			
			!=== CELL-CENTERED DATA: ===
			call write_cell_data(field_idx, 2, var_idx)
			
        end do
        
        print *, ' '
        print *, 'DAT file: ', filename
        print *, ' '
        contains
			!=== AUXILIATY SUBROUTINES: ===
			subroutine write_face_data(field_id, zone_id, var_id, start_idx, end_idx)
				integer, intent(in) :: field_id, zone_id, var_id, start_idx, end_idx
								
				!=== OUTPUT FACE DATA: ===
				write(iunit, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
					'(', 300, ' (', var_id, ' ', zone_id, ' ', 1, ' ', 0, ' ', 1, ' ', start_idx, ' ', end_idx, ')'
				
				write(iunit, '(A)') ' ('
				
				do face_idx = start_idx, end_idx
					if (USE_GHOST_CELLS) then
						write(iunit, '(E20.12E3)') 0.5d0*(fields(field_id, face_right_cell(face_idx)) + fields(field_id, face_left_cell(face_idx)))
					else
						write(iunit, '(E20.12E3)') fields(field_id, face_right_cell(face_idx))
					end if
				end do
				
				write(iunit, '(A)') ')'
				write(iunit, '(A)') ')'
							
			end subroutine
			
			subroutine write_cell_data(field_id, zone_id, var_id)
				integer, intent(in) :: field_id, zone_id, var_id
				
				!=== OUTPUT CELL DATA: ===
				write(iunit, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
					'(', 300, ' (', var_id, ' ', zone_id, ' ', 1, ' ', 0, ' ', 1, ' ', 1, ' ', ncells, ')'
				
				write(iunit, '(A)') ' ('
				
				do cell_idx = 1, ncells
					write(iunit, '(E20.12E3)') fields(field_id, cell_idx)
				end do
				
				write(iunit, '(A)') ')'
				write(iunit, '(A)') ')'
							
			end subroutine
			
			function name_to_var(field_name) result(var_idx)
				implicit none
				character(len=*), intent(in) :: field_name
				integer :: var_idx
				
				select case(trim(adjustl(field_name)))
				case('Velocity_x', 'Velocity_X', 'VELOCITY_X',&
					 'velocity_x', 'vel_x', 'Vel_x', 'Vel_X', 'VEL_X')
					var_idx = VAR_VELOCITY_X
				case('Velocity_y', 'Velocity_Y', 'VELOCITY_Y',&
					 'velocity_y', 'vel_y', 'Vel_y', 'Vel_Y', 'VEL_Y')
					var_idx = VAR_VELOCITY_Y
				case('Velocity_z', 'Velocity_Z', 'VELOCITY_Z',&
					 'velocity_z', 'vel_z', 'Vel_z', 'Vel_Z', 'VEL_Z')
					var_idx = VAR_VELOCITY_Z
				case('Pressure', 'pressure', 'PRESSURE', 'PRS', 'prs')
					var_idx = VAR_PRESSURE	
				case('Temperature', 'temperature', 'TEMPERATURE', 'TMP', 'tmp')
					var_idx = VAR_TEMPERATURE	
				case default
					error stop ' Unknown field name: ' // trim(field_name)
				end select
			end function
			
	end subroutine
	
!=======================================================================
!==================== FLUENT .DAT INPUT SUBROUTINE =====================
!=======================================================================
	subroutine input_fields_dat(fields, field_names, filename, ncells)
		implicit none
        character(len=*), intent(in) :: filename
        integer, intent(in) :: ncells
		real(8), intent(inout) :: fields(:, :)
        character(len=32), allocatable, intent(inout) :: field_names(:)
        		
		integer :: cell_id, j
		integer :: iunit, ios, pos1, pos2
		integer :: nvars, zone_idx, var_idx
		character(len=128) :: line, aux_line, info_line
		character(len=:), allocatable :: token
		character(len=32) :: name
				       
        open(newunit=iunit, file=filename, status='old', action='read') 
        nvars = 0
        !=== READING FILE: ===
        do
			read(iunit, '(a)', iostat=ios) line
			if (ios /= 0) then
				exit
			end if
			
			if (index(line, '(300') > 0) then
				pos1 = index(trim(line(2:)), '(') + 1
				pos2 = index(trim(line), ')')
				
				info_line = line(pos1 + 1:pos2 - 1)					
				call get_token(info_line, 2, token)
							
				read(token, *) zone_idx
				if (zone_idx .ne. 2) cycle
					
				nvars = nvars + 1
				call get_token(info_line, 1, token)
				read(token, *) var_idx
				name = var_to_name(var_idx)
				
				if (.not. allocated(field_names)) then
					field_names = [character(len=32) :: name]
				else
					field_names = [character(len=32) :: field_names, name]
				end if
				
				read(iunit, '(a)') aux_line
					    
				do cell_id = 1, ncells
					read(iunit, *) fields(nvars, cell_id)
				end do
			end if			
		end do
	
        close(iunit)
        
        contains
			!=== AUX SUBROUTINES: ===
			function var_to_name(var_idx) result(field_name)
				implicit none
				integer, intent(in) :: var_idx
				character(len=32) :: field_name
				
				select case(var_idx)
				case(VAR_VELOCITY_X)
					field_name = 'Velocity_x' 
				case(VAR_VELOCITY_Y)
					field_name = 'Velocity_y' 
				case(VAR_VELOCITY_Z)
					field_name = 'Velocity_z' 
				case(VAR_PRESSURE)
					field_name = 'Pressure' 
				case(VAR_TEMPERATURE)
					field_name = 'Temperature' 
				case default
					error stop ' Unknown field index: '
				end select
			end function
			
    end subroutine




























!=======================================================================
!===================== AUX MESH OUTPUT SUBROUTINES =====================
!=======================================================================
	subroutine write_mesh_faces_vtk(filename, dim, nnodes, node_coords,&
									nfaces, face_nodes_ptr, face_nodes,&
									face_zone)
		implicit none
		character(len=*), intent(in) :: filename
		integer, intent(in) :: dim, nnodes, nfaces
		real(8), intent(in) :: node_coords(:, :)
		integer, intent(in) :: face_nodes_ptr(:), face_nodes(:)
		integer, optional, intent(in) :: face_zone(:)

		integer :: i, j, iunit, ntotal

		open(newunit=iunit, file=filename, status='replace', action='write')
	
		write(iunit, '(a)') '# vtk DataFile Version 3.0'
		write(iunit, '(a)') 'Coarse mesh output (polygons)'
		write(iunit, '(a)') 'ASCII'
		write(iunit, '(a)') 'DATASET UNSTRUCTURED_GRID'

		!=== POINTS OUTPUT: ===
		write(iunit, '(a, i0, a)') 'POINTS ', nnodes, ' double'
		do i = 1, nnodes
			if (dim == 2) then
				write(iunit, '(3(1x,e23.15))') node_coords(1,i), node_coords(2,i), 0.0d0
			else
				write(iunit, '(3(1x,e23.15))') node_coords(1,i), node_coords(2,i), node_coords(3,i)
			end if
		end do

		!=== TOTAL NODE-FACES CONNECTIVITIS: ===
		ntotal = 0
		do i = 1, nfaces
			ntotal = ntotal + (face_nodes_ptr(i+1) - face_nodes_ptr(i))
		end do

		!=== OUTPUT FACES: ===
		write(iunit, '(a, i0, a, i0)') 'CELLS ', nfaces, ' ', ntotal + nfaces
		do i = 1, nfaces
			write(iunit, '(i0)', advance='no') (face_nodes_ptr(i+1) - face_nodes_ptr(i))
			do j = face_nodes_ptr(i), face_nodes_ptr(i+1)-1
				write(iunit, '(1x,i0)', advance='no') face_nodes(j) - 1
			end do
			write(iunit, *)
		end do

		!=== CELL TYPES: ===
		write(iunit, '(a, i0)') 'CELL_TYPES ', nfaces
		do i = 1, nfaces
			write(iunit, '(i0)') 7
		end do

		if (present(face_zone)) then
			write(iunit, '(a, i0)') 'CELL_DATA ', nfaces
			write(iunit, '(a)') 'SCALARS face_zone int 1'
			write(iunit, '(a)') 'LOOKUP_TABLE default'
			do i = 1, nfaces
				write(iunit, '(i0)') face_zone(i)
			end do
		end if

		close(iunit)
	end subroutine
	
	subroutine write_mesh_cell2d_vtk(filename, nnodes, node_coords,&
									 ncells, cell_nodes_ptr, cell_nodes,&
									 cell_type)
		implicit none
		character(len=*), intent(in) :: filename
		integer, intent(in) :: nnodes, ncells
		real(8), intent(in) :: node_coords(:,:)
		integer, intent(in) :: cell_nodes_ptr(:), cell_nodes(:)
		integer, optional, intent(in) :: cell_type(:)

		integer :: i, j, iunit, ntotal

		open(newunit=iunit, file=filename, status='replace', action='write')
		write(iunit, '(a)') '# vtk DataFile Version 3.0'
		write(iunit, '(a)') '2D coarse mesh'
		write(iunit, '(a)') 'ASCII'
		write(iunit, '(a)') 'DATASET UNSTRUCTURED_GRID'

		!=== OUTPUT POINTS: ===
		write(iunit, '(a, i0, a)') 'POINTS ', nnodes, ' double'
		do i = 1, nnodes
			write(iunit, '(3(1x,e23.15))') node_coords(1,i), node_coords(2,i), 0.0d0
		end do

		!=== TOTAL CELL CONNECTIVITIES: ===
		ntotal = 0
		do i = 1, ncells
			ntotal = ntotal + (cell_nodes_ptr(i+1) - cell_nodes_ptr(i))
		end do
		
		!=== OUTPUT CELLS: ===
		write(iunit, '(a, i0, a, i0)') 'CELLS ', ncells, ' ', ntotal + ncells
		do i = 1, ncells
			write(iunit, '(i0)', advance='no') (cell_nodes_ptr(i+1) - cell_nodes_ptr(i))
			do j = cell_nodes_ptr(i), cell_nodes_ptr(i+1)-1
				write(iunit, '(1x,i0)', advance='no') cell_nodes(j) - 1
			end do
			write(iunit, *)
		end do
		
		!=== CELL TYPES: ===
		write(iunit, '(a, i0)') 'CELL_TYPES ', ncells
		do i = 1, ncells
			write(iunit, '(i0)') 7
		end do

		if (present(cell_type)) then
			write(iunit, '(a, i0)') 'CELL_DATA ', ncells
			write(iunit, '(a)') 'SCALARS cell_type int 1'
			write(iunit, '(a)') 'LOOKUP_TABLE default'
			do i = 1, ncells
				write(iunit, '(i0)') cell_type(i)
			end do
		end if

		close(iunit)
	end subroutine

	subroutine write_mesh_cell3d_vtk(filename, dim, nnodes, node_coords,&
									 ncells, cell_faces_ptr, cell_faces,&
									 nfaces, face_nodes_ptr, face_nodes)
		implicit none
		character(len=*), intent(in) :: filename
		integer, intent(in) :: dim, nnodes, ncells, nfaces
		real(8), intent(in) :: node_coords(:,:)
		integer, intent(in) :: cell_faces_ptr(:), cell_faces(:)
		integer, intent(in) :: face_nodes_ptr(:), face_nodes(:)

		integer :: i, j, iunit, cell_idx, face_idx
		integer :: nfaces_cell, face_id, nverts, pos
		integer :: total_cell_size, total_data
		integer, allocatable :: cell_buffer(:)
		integer :: buffer_pos

		open(newunit=iunit, file=filename, status='replace', action='write')
		write(iunit, '(a)') '# vtk DataFile Version 3.0'
		write(iunit, '(a)') 'Coarse mesh with polyhedra'
		write(iunit, '(a)') 'ASCII'
		write(iunit, '(a)') 'DATASET UNSTRUCTURED_GRID'

		!=== OUTPUT PINTS: ===
		write(iunit, '(a, i0, a)') 'POINTS ', nnodes, ' double'
		do i = 1, nnodes
			if (dim == 2) then
				write(iunit, '(3(1x,e23.15))') node_coords(1,i), node_coords(2,i), 0.0d0
			else
				write(iunit, '(3(1x,e23.15))') node_coords(1,i), node_coords(2,i), node_coords(3,i)
			end if
		end do

		!=== TOTAL CELL CONNECTIVITIES: ===
		total_cell_size = 0
		do cell_idx = 1, ncells
			nfaces_cell = cell_faces_ptr(cell_idx+1) - cell_faces_ptr(cell_idx)
			total_cell_size = total_cell_size + 1
			do i = cell_faces_ptr(cell_idx), cell_faces_ptr(cell_idx+1)-1
				face_idx = cell_faces(i)
				nverts = face_nodes_ptr(face_idx+1) - face_nodes_ptr(face_idx)
				total_cell_size = total_cell_size + 1
				total_cell_size = total_cell_size + nverts
			end do
		end do

		!=== OUTPUT CELLS: ===
		write(iunit, '(a, i0, a, i0)') 'CELLS ', ncells, ' ', total_cell_size
		allocate(cell_buffer(total_cell_size))
		buffer_pos = 1
		do cell_idx = 1, ncells
			nfaces_cell = cell_faces_ptr(cell_idx+1) - cell_faces_ptr(cell_idx)
			cell_buffer(buffer_pos) = nfaces_cell
			buffer_pos = buffer_pos + 1
			do i = cell_faces_ptr(cell_idx), cell_faces_ptr(cell_idx+1)-1
				face_idx = cell_faces(i)
				nverts = face_nodes_ptr(face_idx+1) - face_nodes_ptr(face_idx)
				cell_buffer(buffer_pos) = nverts
				buffer_pos = buffer_pos + 1
				do j = face_nodes_ptr(face_idx), face_nodes_ptr(face_idx+1)-1
					cell_buffer(buffer_pos) = face_nodes(j) - 1
					buffer_pos = buffer_pos + 1
				end do
			end do
		end do
		do i = 1, total_cell_size
			write(iunit, '(i0)', advance='no') cell_buffer(i)
			if (mod(i, 10) == 0 .or. i == total_cell_size) then
				write(iunit, *)
			else
				write(iunit, '(" ")', advance='no')
			end if
		end do
		deallocate(cell_buffer)

		!=== CELL TYPES: ===
		write(iunit, '(a, i0)') 'CELL_TYPES ', ncells
		do cell_idx = 1, ncells
			write(iunit, '(i0)') 42
		end do

		close(iunit)
	end subroutine




end module
