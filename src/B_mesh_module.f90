module mesh_module
use mesh_connectivity_module
use mesh_ordering_module
use string_utils_module
use metric2D_computation_module
use metric3D_computation_module
implicit none
private
!=======================================================================
!======================= MESH DATA STRUCTURE ===========================
!=======================================================================
type, public :: mesh_t
    integer :: dim = 3													!=== MESH DIMENSION ===													
	logical :: cache_misses_minimizer = .true.							!=== FLAG FOR REORDERING ALL DATA IN COMPACT MANNER ===

    !=== NODES: ===
    integer :: nnodes													!=== NUMBER OF NODES ===
    real(8), allocatable :: node_coords(:,:)   							!=== COORDINATES OF MESH NODES, DIM = (dim, nnodes) ===


    !=== FACES: ===													
    integer :: nfaces													!=== NUMBER OF FACES ===    
	!== FACES NODES (*face-node connectivity): ==
    integer, allocatable :: face_nodes(:)       						!=== ALL FACE NODES IN CSR FORMAT, DIM = (sum_i(number of nodes for face_i)) ===
    integer, allocatable :: face_nodes_ptr(:)   						!=== FACE NODES POINTER, DIM = (nfaces+1) ===
    !== FACES OWNERS: ==
    integer, allocatable :: face_left_cell(:)        					!=== FACE LEFT CELLS, DIM = (nfaces) ===
    integer, allocatable :: face_right_cell(:)       					!=== FACE RIGHT CELLS, DIM = (nfaces) ===
    !== FACE ZONE DATA: ==
    integer, allocatable :: face_zone(:)        						!=== FACE ZONE INFORMATION, DIM = (nfaces) === *used for boundary conditions
    integer, allocatable :: face_type(:)        						!=== FACE TYPE INFORMTAION, DIM = (nfaces) === *not used (boundary conditions in fluent)
	integer, allocatable :: face_bidx(:)        						!=== FACE LOCAL BOUNDARY INDEX, DIM = (nfaces) === *for internal faces bidx=-1
	character(len=128), allocatable :: zone_names(:), zone_subnames(:)  !=== FACE ZONES NAMES, DIM = (1:max(face_zone)) ===


    !=== CELLS: ===
    integer :: ncells													!=== NUMBER OF CELLS ===
    !== CELL FACES (*cell-face connectivity): ==
    integer, allocatable :: cell_faces(:)       						!=== ALL CELL FACES IN CSR FORMAT, DIM = (sum_i(number of faces for cell_i)) ===
    integer, allocatable :: cell_faces_ptr(:) 							!=== CELL FACES POINTER, DIM = (ncells+1) ===
    !== CELL NODES (*cell-node connectivity): ==
    integer, allocatable :: cell_nodes(:)       						!=== ALL CELL NODES IN CSR FORMAT, DIM = (sum_i(number of nodes for cell_i)) ===
    integer, allocatable :: cell_nodes_ptr(:) 							!=== CELL NODES POINTER, DIM = (ncells+1) ===
	!== AUXILIARY CELL DATA: ==
	integer, allocatable :: cell_type(:)								!=== ELEMENTS TYPE, DIM = (ncells) === *indexation according to fluent .cas file
	
	
    !=== BOUNDARY FACES: ===
    integer :: nbfaces													!=== NUMBER OF BOUNDARY FACES === 

    !=== MESH METRIC: ===	
    real(8), allocatable :: face_center(:,:)    						!=== FACE CENTERS COORDINATES, DIM = (dim, nfaces) ===
    real(8), allocatable :: face_normal(:,:)    						!=== FACE NORMALS, DIM = (dim, nfaces) ===
    real(8), allocatable :: face_area(:)    							!=== FACE AREAS, DIM = (nfaces) ===
    real(8), allocatable :: face_weight(:)      						!=== FACE WEIGHT FACTORS, DIM = (nfaces) ===
    real(8), allocatable :: cell_center(:,:)    						!=== CELL CENTERS COORDINATES, DIM = (dim, 1:ncells) ===
    real(8), allocatable :: cell_volume(:)      						!=== CELL VOLUMES VALUES, DIM = (ncells) ===

    integer :: verbosity = 0
contains
    procedure :: read_fluent => mesh_read_fluent
    procedure :: write_fluent => mesh_write_fluent
    procedure :: compute_metric => mesh_compute_metric
end type




contains
!=======================================================================
!================== MESH INPUT/OUTPUT SUBROUTINES ======================
!=======================================================================
	subroutine mesh_read_fluent(this, filename)
		implicit none
		class(mesh_t), intent(inout) :: this
		character(len=*), intent(in) :: filename
		
		integer :: iunit, ios, stat, i
		integer :: dim, nnodes, nfaces, ncells, nbfaces, max_face_zone_id
		integer, allocatable :: nnodes_per_face(:)
				
		if (this%verbosity > 0) then
			write(*,*) ' '
			write(*,*) ' Reading mesh file, inputfile: ', filename
		end if
		
		open(newunit=iunit, file=filename, status='old', action='read', iostat=ios)
		if (ios /= 0) then
			write(*,*) '<<< Mesh file opening error ! >>>'
			stop
		end if
		
		!=== READING MESH DIMENSIONS: ===
		call read_mesh_dimensions(iunit, stat, dim, nnodes, nfaces, ncells,&
								  max_face_zone_id, nnodes_per_face)
		if (stat .ne. 0) then
			write(*,*) '!!! Dimension reading error !!!'
			stop
		end if
		
		this%dim = dim
		this%nnodes = nnodes
		this%nfaces = nfaces
		this%ncells = ncells
		
		allocate(this%zone_names(max_face_zone_id))
		allocate(this%zone_subnames(max_face_zone_id))
		
		!=== BASIC ARRAYS ALLOCATION: ===
		allocate(this%node_coords(dim, nnodes),&
				 this%face_left_cell(nfaces),&
				 this%face_right_cell(nfaces),&
				 this%face_zone(nfaces),&
				 this%face_type(nfaces),&
				 this%cell_type(ncells))
				 
				 
		!=== FACE-NODE CONNECTIVITY: ===
		call face_node_connectivity(this%nfaces, nnodes_per_face,&
									this%face_nodes_ptr, this%face_nodes)
		
		!=== MESH DATA READING: ===
		call read_mesh_data(this, iunit)
		if (allocated(nnodes_per_face)) deallocate(nnodes_per_face)
		close(iunit)
		
		!=== REODERING CELL TO MINIMIZE CELL DIST IN CACHE USING RMC ===
		! ALGHORYTHM: 												 ===
		if (this%cache_misses_minimizer) then
			call reorder_cells_rcm(this%ncells, this%nfaces,&
								   this%face_left_cell, this%face_right_cell,&
								   this%face_nodes_ptr, this%face_nodes,&
								   this%cell_type)
		else
			call reorder_left_right_cells(this%nfaces, this%ncells,&
										  this%face_left_cell, this%face_right_cell,&
									      this%face_nodes_ptr, this%face_nodes)
		end if
		
		!=== FACES REORDERING (faces = [internal_block, external_block]): ===
		if (this%cache_misses_minimizer) then
			call reorder_faces_advanced(this%ncells, this%nfaces, this%nbfaces,&
										this%face_left_cell, this%face_right_cell,&
										this%face_zone, this%face_type,&
										this%face_nodes_ptr, this%face_nodes)
		else
			call reorder_faces(this%ncells, this%nfaces, this%nbfaces,&
							   this%face_left_cell, this%face_right_cell,&
							   this%face_zone, this%face_type,&
							   this%face_nodes_ptr, this%face_nodes)
		end if
		
		!=== CELL-FACE CONNECTIVITY: ===
		call cell_face_connectivity(this%ncells, this%nfaces,&
									this%face_left_cell, this%face_right_cell,&
									this%cell_faces_ptr, this%cell_faces)
		
		!=== CELL-NODE CONNECTIVITY: ===
		call cell_node_connectivity(this%nnodes, this%ncells,&
									this%face_nodes_ptr, this%face_nodes,&
									this%cell_faces_ptr, this%cell_faces,&
									this%cell_nodes_ptr, this%cell_nodes)
		
		!=== BOUNDARY FACES LOCAL NUMERATION: ===
		call boundary_faces_connectivity(this%nfaces, this%nbfaces, this%face_bidx)
		
		!=== ELEMENTS ORDERING: ===
		if (this%dim == 2) then
			call order_cell_nodes_2d(this%ncells, this%node_coords,&
									 this%cell_nodes_ptr, this%cell_nodes)
		else if (this%dim == 3) then
			call order_face_nodes_3d(this%nfaces, this%node_coords,&
									 this%face_nodes_ptr, this%face_nodes)
		end if
		
		!=== GHOST CELLS/BOUNDARY FACES-CELLS NUMERATION ===
		do i = this%nfaces-this%nbfaces+1, this%nfaces
			this%face_right_cell(i) = this%ncells+this%face_bidx(i)
		end do
		
		!=== MESH STATISTICS: ===
		if (this%verbosity > 0) then
			write(*, '(a, i0, a)') '  <<< Problem Dimension: ', this%dim, ' >>> '
			write(*, '(a, i0, a)') '  <<< Number of nodes: ', this%nnodes, ' >>> '
			write(*, '(a, i0, a)') '  <<< Number of cells: ', this%ncells, ' >>> '
			write(*, '(a, i0, a)') '  <<< Number of faces: ', this%nfaces, ' >>> '
			write(*, '(a, i0, a)') '  <<< Number of boundary faces: ', this%nbfaces, ' >>> '
			print*, ' '
		end if
		
		contains	
		!=============================================================== 
		!====================== READING SUBROUTINES ====================
		!=============================================================== 
		subroutine read_mesh_dimensions(iunit, stat, dim, nnodes, nfaces, ncells,&
										max_face_zone_id, temp_num_nodes_per_face)
			implicit none
			integer, intent(in) :: iunit
			integer, intent(inout) :: stat, dim, nnodes, nfaces, ncells, max_face_zone_id
			integer, allocatable, intent(inout) :: temp_num_nodes_per_face(:)
			
			integer :: ios, pos1, pos2, dim1, dim2
			integer :: face_nodes_num, face_zone_id
			integer :: zone_id, first_idx, last_idx, face_topology, i
			character(len=:), allocatable :: token, aux_line
			character(len=256) :: line
			
			stat = 0
			token = ''; aux_line = ''
			dim = -1; nnodes = -1; ncells = -1; nfaces = -1; max_face_zone_id = -1
			face_nodes_num = 0
			
			do
				read(iunit, '(a)', iostat=ios) line
				if (ios /= 0) then
					exit
				end if
				
				!=== DIMENSION: ===
				if (index(line, 'Dimension:') > 0) then
					read(iunit, '(a)') line
					pos1 = index(trim(line), '(')
					pos2 = index(trim(line), ')')
					
					aux_line = line(pos1+1:pos2-1)
					
					read(aux_line, *) dim1, dim2
					dim = max(dim1, dim2)
				end if			
				
				!=== NUMBER OF NODES: ===
				if (index(line, '(10') == 1) then
					pos1 = index(trim(line(2:)), '(') + 1
					pos2 = index(trim(line), ')')
					
					aux_line = line(pos1 + 1:pos2 - 1)
					call get_token(aux_line, 1, token)
					zone_id = hex_to_int(token)
					
					if (zone_id == 0) then
						call get_token(aux_line, 3, token)
						nnodes = hex_to_int(token)
					end if
				end if		
				
				!=== NUMBER OF CELLS: ===
				if (index(line, '(12') == 1) then
					pos1 = index(trim(line(2:)), '(') + 1
					pos2 = index(trim(line), ')')
					
					aux_line = line(pos1 + 1:pos2 - 1)
					call get_token(aux_line, 1, token)
					zone_id = hex_to_int(token)
					
					if (zone_id == 0) then
						call get_token(aux_line, 3, token)
						ncells = hex_to_int(token)
					end if
				end if	
				
				!=== NUMBER OF FACES: ===
				if (index(line, '(13') == 1) then
					pos1 = index(trim(line(2:)), '(') + 1
					pos2 = index(trim(line), ')')
					
					aux_line = line(pos1 + 1:pos2 - 1)
					call get_token(aux_line, 1, token)
					zone_id = hex_to_int(token)
					
					if (zone_id == 0) then
						call get_token(aux_line, 3, token)
						nfaces = hex_to_int(token)
						if (.not. allocated(temp_num_nodes_per_face)) &
							allocate(temp_num_nodes_per_face(nfaces))
						
					else 
						call get_token(aux_line, 2, token)
						first_idx = hex_to_int(token)
						
						call get_token(aux_line, 3, token)
						last_idx = hex_to_int(token)
						
						call get_token(aux_line, 5, token)
						face_topology = hex_to_int(token)
						
						if (face_topology == 0) then
							do i = first_idx, last_idx
								read(iunit, '(a)', iostat=ios) line
								
								call get_token(trim(line), 1, token)
								face_topology = hex_to_int(token)
								
								temp_num_nodes_per_face(i) = face_topology							
							end do
							
						else
							do i = first_idx, last_idx				
								temp_num_nodes_per_face(i) = face_topology	
							end do
							
						end if
					
					end if
				end if	
				
				!=== FACE ZONE NAMES: ===
				if (index(line, '(45') == 1) then
					pos1 = index(trim(line(2:)), '(') + 1
					pos2 = index(trim(line), ')')
					
					aux_line = line(pos1 + 1:pos2 - 1)
					read(aux_line, *) face_zone_id
					max_face_zone_id = max(max_face_zone_id, face_zone_id)
				end if
			end do
			
			rewind(iunit)
			
			if ( (dim <= 0) .or.&
				 (nnodes <= 0) .or.&
				 (ncells <= 0) .or.&
				 (nfaces <= 0) .or.&
				 max_face_zone_id < 0) then
				stat = -1
			end if
			
		end subroutine 
		
		subroutine read_mesh_data(this, iunit)
			implicit none
			integer, intent(in) :: iunit
			type(mesh_t), intent(inout) :: this
			
			integer :: zone_id, first_idx, last_idx, face_type, face_topology, cell_type
			integer :: i, j, ios, pos1, pos2, counter
			integer :: num_tokens
			integer :: num_nodes, node_idx, face_zone_id
			integer :: left_cell, right_cell
			integer, allocatable :: cell_type_arr(:)
			character(len=256) :: line, data_line
			character(len=128) :: zone_name, zone_subname
			character(len=:), allocatable :: token, aux_line
			
			
			token = ''; aux_line = ''
			
			do
				read(iunit, '(a)', iostat=ios) line
				if (ios /= 0) then
					exit
				end if
				
				!=== NODES DATA: ===
				if (index(line, '(10') == 1) then
					pos1 = index(trim(line(2:)), '(') + 1
					pos2 = index(trim(line), ')')
					
					aux_line = line(pos1 + 1:pos2 - 1)
					call get_token(aux_line, 1, token)
					zone_id = hex_to_int(token)
					
					if (zone_id == 0) cycle
					
					call get_token(aux_line, 2, token)
					first_idx = hex_to_int(token)
					call get_token(aux_line, 3, token)
					last_idx = hex_to_int(token)
					
					do i = first_idx, last_idx
						read(iunit, '(a)') data_line
						read(data_line, *) this%node_coords(:, i)
					end do
				end if	
				
				!=== FACES DATA: ===
				if (index(line, '(13') == 1) then
					pos1 = index(trim(line(2:)), '(') + 1
					pos2 = index(trim(line), ')')
					
					aux_line = line(pos1 + 1:pos2 - 1)
					
					call get_token(aux_line, 1, token)
					zone_id = hex_to_int(token)
					
					if (zone_id == 0) cycle
					
					call get_token(aux_line, 2, token)
					first_idx = hex_to_int(token)
					call get_token(aux_line, 3, token)
					last_idx = hex_to_int(token)
					
					call get_token(aux_line, 4, token)
					face_type = hex_to_int(token)
					call get_token(aux_line, 5, token)
					face_topology = hex_to_int(token)
										
					if (face_topology == 0) then
						do i = first_idx, last_idx
							read(iunit, '(a)') data_line
							num_tokens = get_num_tokens(data_line)
							
							call get_token(data_line, 1, token)
							num_nodes = hex_to_int(token)
							
							do j = 1, num_nodes
								call get_token(data_line, j + 1, token) 
								node_idx = hex_to_int(token)
								
								this%face_nodes(this%face_nodes_ptr(i) + j - 1) = node_idx								
							end do
							
							call get_token(data_line, num_nodes + 2, token)  
							left_cell = hex_to_int(token)
							call get_token(data_line, num_nodes + 3, token)  
							right_cell = hex_to_int(token)
							
							this%face_left_cell(i) = left_cell
							this%face_right_cell(i) = right_cell
							this%face_zone(i) = zone_id
							this%face_type(i) = face_type
						end do
					
					else 
						num_nodes = face_topology
						do i = first_idx, last_idx
							read(iunit, '(a)') data_line
							num_tokens = get_num_tokens(data_line)
														
							do j = 1, num_nodes
								call get_token(data_line, j, token) 
								node_idx = hex_to_int(token)
								
								this%face_nodes(this%face_nodes_ptr(i) + j - 1) = node_idx	
							end do
							
							call get_token(data_line, num_nodes + 1, token)  
							left_cell = hex_to_int(token)
							call get_token(data_line, num_nodes + 2, token)  
							right_cell = hex_to_int(token)
							
							this%face_left_cell(i) = left_cell
							this%face_right_cell(i) = right_cell
							this%face_zone(i) = zone_id
							this%face_type(i) = face_type					
						end do
					end if 
				end if
				
				!=== CELL DATA: ===
				if (index(line, '(12') == 1) then
					pos1 = index(trim(line(2:)), '(') + 1
					pos2 = index(trim(line), ')')
					
					aux_line = line(pos1 + 1:pos2 - 1)
					call get_token(aux_line, 1, token)
					zone_id = hex_to_int(token)
					
					if (zone_id == 0) cycle
					
					call get_token(aux_line, 2, token)
					first_idx = hex_to_int(token)
					call get_token(aux_line, 3, token)
					last_idx = hex_to_int(token)
					
					call get_token(aux_line, 5, token)
					cell_type = hex_to_int(token)
					
					if (cell_type .ne. 0) then
						do i = first_idx, last_idx
							this%cell_type(i) = cell_type
						end do
						
					else
						allocate(cell_type_arr(1:1 + (last_idx - first_idx)))
						read(iunit, *) (cell_type_arr(i), i = 1, 1 + (last_idx - first_idx))
						
						counter = 1
						do i = first_idx, last_idx
							this%cell_type(i) = cell_type_arr(counter)
							counter = counter + 1
						end do
						deallocate(cell_type_arr)
					end if
					
					
				end if
				
				!=== FACE ZONE NAMES: ===
				if (index(line, '(45') == 1) then
					pos1 = index(trim(line(2:)), '(') + 1
					pos2 = index(trim(line), ')')
					
					aux_line = line(pos1 + 1:pos2 - 1)
					read(aux_line, *) face_zone_id, zone_name, zone_subname
					this%zone_names(face_zone_id) = zone_name
					this%zone_subnames(face_zone_id) = zone_subname
				end if
				
			end do
			
			rewind(iunit)
		end subroutine
		
	end subroutine

	subroutine mesh_write_fluent(this, filename)
		implicit none
		class(mesh_t), intent(in) :: this
		character(len=*), intent(in) :: filename
		
		integer :: i, j, k, iunit, ios, ierr, stat
		integer :: num_zones, counter, current_zone
		integer :: zones_id(1:32), zones_type(1:32), zones_indexes(1:2, 1:32)
		
		open(newunit=iunit, file=filename, status='replace', action='write', iostat=ios)
		if (ios /= 0) then
			write(*,*) '<<< Mesh file opening error ! >>>'
			stop
		end if
		
		!=== WRITING MESH DIMENSION: ===
		write(iunit, '(a)') '(0 "GENERATED FROM CODE")'
		write(iunit, '(a)') '(0 "Dimension:")'
		write(iunit, '(a, i0, a, i0, a)') '(', 2, ' ', this%dim, ')'
		
		!=== WRITING MESH NODES: ===
		write(iunit, *) ' '
		write(iunit, '(a)') '(0 "Nodes:")'
		write(iunit, '(a, Z0, a, Z0, a)') '(10 (0 1 ', this%nnodes, ' 1 ', size(this%node_coords(:, 1)), '))'
		write(iunit, '(a, Z0, a, Z0, a)') '(10 (1 1 ', this%nnodes, ' 1 ', size(this%node_coords(:, 1)), ')('
		do i = 1, this%nnodes
			write(iunit, *) this%node_coords(:, i)
		end do
		write(iunit, '(a)') '))'
		
		!=== WRITING MESH FACES: ===
		write(iunit, *) ' '
		write(iunit, '(a)') '(0 "Faces:")'
		write(iunit, '(a, Z0, a)') '(13 (0 1 ', this%nfaces, ' 0))'
		
		zones_id(:) = 0
		num_zones = 0
		do i = 1, this%nfaces
			current_zone = this%face_zone(i)
			if (any(zones_id == current_zone)) cycle
			num_zones = num_zones + 1
			zones_id(num_zones) = current_zone
			zones_type(num_zones) = this%face_type(i)
		end do
		
		zones_indexes(1, :) = huge(1); zones_indexes(2, :) = -huge(1)	
		do j = 1, num_zones
			do i = 1, this%nfaces
				current_zone = this%face_zone(i)
				if (current_zone .ne. zones_id(j)) cycle
				
				zones_indexes(1, j) = min(zones_indexes(1, j), i)
				zones_indexes(2, j) = max(zones_indexes(2, j), i)
			end do
		end do
		
	
		do j = 1, num_zones
			write(iunit, '(a, Z0, a, Z0, a, Z0, a, Z0, a)') '(13 (', zones_id(j), ' ', zones_indexes(1, j), ' ', zones_indexes(2, j), ' ', zones_type(j), ' 0)('
			do i = zones_indexes(1, j), zones_indexes(2, j)
			
				write(iunit, '(Z0, a)', advance='no') this%face_nodes_ptr(i + 1) - this%face_nodes_ptr(i), ' '
    
				do k = this%face_nodes_ptr(i), this%face_nodes_ptr(i + 1) - 1
					write(iunit, '(Z0, a)', advance='no') this%face_nodes(k), ' '
				end do
				
				write(iunit, '(Z0, a, Z0)') this%face_left_cell(i), ' ', this%face_right_cell(i)
			end do
					
			write(iunit, '(a)') '))'
		end do
		
		!=== WRITING MESH CELLS: ===
		write(iunit, *) ' '
		write(iunit, '(a)') '(0 "Cells:")'
		write(iunit, '(a, Z0, a)') '(12 (0 1 ', this%ncells, ' 0))'
		write(iunit, '(a, Z0, a)') '(12 (2 1 ', this%ncells, ' 1 0)('
		write(iunit, '(Z0)') (this%cell_type(i), i = 1, this%ncells)
		write(iunit, '(a)') '))'
		
		
		write(iunit, *) ' '
		write(iunit, '(a)') '(0 "Zones:")'
		write(iunit, '(a)') '(45 (2 fluid fluid)())'
		do j = 1, num_zones
			write(iunit, '(a, i0, a, a, a, a, a, a)') '(45 (', zones_id(j), ' ', trim(adjustl(this%zone_names(zones_id(j)))), ' ', trim(adjustl(this%zone_subnames(zones_id(j)))), ')())'
		end do
		
		
		close(iunit)
	end subroutine



!=======================================================================
!================ MESH METRIC CALCULATION SUBROUTINES ==================
!=======================================================================
	subroutine mesh_compute_metric(this)
		implicit none
		class(mesh_t), intent(inout) :: this
		
		allocate(this%cell_center(this%dim, this%ncells + this%nbfaces),&
				 this%cell_volume(this%ncells),&
				 this%face_weight(this%nfaces),&
				 this%face_center(this%dim, this%nfaces),&
				 this%face_normal(this%dim, this%nfaces),&
				 this%face_area(this%nfaces))
		
		if (this%dim == 2) then 
			call compute_metric_2d(this%ncells, this%nfaces, this%nbfaces, this%dim,&
								   this%face_left_cell, this%face_right_cell,&
							       this%face_nodes_ptr, this%face_nodes,&
							       this%cell_nodes_ptr, this%cell_nodes,&
								   this%cell_faces_ptr, this%cell_faces,&
								   this%node_coords,&
								   this%face_center, this%face_normal,&
								   this%face_area, this%face_weight,&
								   this%cell_center, this%cell_volume)
		else if (this%dim == 3) then
			call compute_metric_3d(this%ncells, this%nfaces, this%nbfaces, this%dim,&
								   this%face_left_cell, this%face_right_cell,&
							       this%face_nodes_ptr, this%face_nodes,&
							       this%cell_nodes_ptr, this%cell_nodes,&
								   this%cell_faces_ptr, this%cell_faces,&
								   this%node_coords,&
								   this%face_center, this%face_normal,&
								   this%face_area, this%face_weight,&
								   this%cell_center, this%cell_volume)
		end if
	end subroutine




end module
