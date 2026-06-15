module structured_mesh_module
use mesh_ordering_module
use mesh_connectivity_module
implicit none
!=======================================================================
!=============== MESH BOUNDARY (EDGE) DATA STRUCTURE ===================
!=======================================================================
type boundary_t
	real(8), allocatable :: x(:), y(:)									!=== ORDERED COORINATES ===
	real(8), allocatable :: t(:) 										!=== CURVE LENGTH PARAMETER ===
	integer :: n														!=== NUMBER OF INITIAL POINTS ===
	integer :: face_zone = -1
	
	logical :: is_interface = .false.
	integer :: neighbor_block_id = -1
contains
	procedure :: init => init_boundary
	procedure :: get_point => get_point_at_t
	procedure :: read_file => read_boundary_file
end type

!=======================================================================
!===================== MESH BLOCK DATA STRUCTURE =======================
!=======================================================================
type block_t
    integer :: id
    type(boundary_t) :: edges(4) 										!=== EDGES = [Bottom, Top, Left, Right]
    integer :: Ni, Nj													!=== BLOCK SIZES ===
    real(8), allocatable :: coords(:,:,:)								!=== BLOCK COORDINATES ===
    real(8) :: delta_u, delta_v 										!=== BIAS FACTORS ===
    character(len=32) :: bias_type_u = 'exp', bias_type_v = 'exp'		!=== BIAS FACTOR TYPE ===
contains
	procedure :: initialize => initialize_block
end type

!=======================================================================
!================== STRUCTURED MESH DATA STRUCTURE =====================
!=======================================================================
type connectivity_t
    integer :: b1, e1
    integer :: b2, e2
    logical :: reverse = .false.
end type

type structured_mesh_t
    type(block_t), allocatable :: blocks(:)
    type(connectivity_t), allocatable :: connections(:)
    
    integer, allocatable :: node_offsets(:)
    integer, allocatable :: cell_offsets(:)
    
    character(len=128), allocatable :: zone_names(:), zone_subnames(:)
contains
    procedure :: read_map => read_mesh_map
    procedure :: assemble => assemble_to_unstructured
end type

contains
!=======================================================================
!=================== MESH BOUNDARY (LINE) METHODS ======================
!=======================================================================
    subroutine init_boundary(this, x_raw, y_raw)
		implicit none
        class(boundary_t), intent(inout) :: this
        real(8), intent(in) :: x_raw(:), y_raw(:)
        
        integer :: i
        real(8) :: dist, total_dist
        
        this%n = size(x_raw)
        if (allocated(this%x)) deallocate(this%x)
        if (allocated(this%y)) deallocate(this%y)
        if (allocated(this%t)) deallocate(this%t)
        allocate(this%x(this%n), this%y(this%n), this%t(this%n))
        this%x = x_raw
        this%y = y_raw
        
        !=== FINDING CURVE LEGNTH: ===
        this%t(1) = 0.0d0
        total_dist = 0.0d0
        do i = 2, this%n
            dist = sqrt((this%x(i)-this%x(i-1))**2 + (this%y(i)-this%y(i-1))**2)
            total_dist = total_dist + dist
            this%t(i) = total_dist
        end do
		
		!=== CURVE LENGTH PARAMETR: ===
        this%t = this%t/total_dist
    end subroutine

    function get_point_at_t(this, t_val) result(res)
		implicit none
		class(boundary_t), intent(in) :: this
		real(8), intent(in) :: t_val
		real(8) :: res(2)
		
		integer :: i, low, high, mid
		real(8) :: local_t
		
		!=== FINDING SEGMENT: ===
		low = 1
		high = this%n - 1
		
		do while (low <= high)
			mid = (low + high)/2
			if (this%t(mid+1) < t_val) then
				low = mid + 1
			else if (this%t(mid) > t_val) then
				high = mid - 1
			else
				i = mid
				exit
			end if
		end do
		
		!=== LINEAR INTERPOLATION: P(t) = LINEAR(P(t_i), P(t_i+1): ===
        local_t = (t_val - this%t(i))/(this%t(i+1) - this%t(i))
        res(1) = this%x(i) + local_t*(this%x(i+1) - this%x(i))
        res(2) = this%y(i) + local_t*(this%y(i+1) - this%y(i))
	end function
	
	subroutine read_boundary_file(this, filename)
		implicit none
		class(boundary_t), intent(inout) :: this
		character(len=*), intent(in) :: filename
		
		real(8), allocatable :: x(:), y(:)
		integer :: n, i, ios, unit_num
		
		open(newunit=unit_num, file=filename, status='old', iostat=ios)
		
		!=== NUMBER OF ROWS: ===
		n = 0
		do
			read(unit_num, *, iostat=ios)
			if (ios /= 0) exit
			n = n + 1
		end do
		
		allocate(x(n), y(n))
		rewind(unit_num)
		
		!=== READING DATA: ===	
		do i = 1, n
			read(unit_num, *) x(i), y(i)
		end do
				
		call this%init(x, y)
		close(unit_num)
		deallocate(x, y)
	end subroutine


!=======================================================================
!====================== MESH BLOCK METHODS =============================
!=======================================================================
	subroutine initialize_block(this, filename)
		implicit none
		class(block_t), intent(inout) :: this
		character(len=*), intent(in) :: filename
		
		logical :: flag
		integer :: iunit, ios
		integer :: pos1, pos2
		character(len=128) line
		character(len=:), allocatable :: aux_line
		real(8) :: C00(2), C10(2), C11(2), C01(2)
		
		if (allocated(this%coords)) deallocate(this%coords)
		
		flag = .false.
		
		open(newunit=iunit, file=filename, status='old')
		!=== READING RECTANGULAR BLOCK DATA: ===
		do 
			read(iunit, '(a)', iostat=ios) line
			if (ios /= 0) exit
			
			!=== BLOCK DIMENSION: ===
			if (index(line, 'Dimension:') > 0) then
				read(iunit, '(a)') line

				read(line, *) this%Ni, this%Nj
				allocate(this%coords(2, this%Ni, this%Nj))
			end if
			
			!=== BIAS FACTORS: ===
			if (index(line, 'Bias:') > 0) then
				!bias types:
				read(iunit, '(a)') line
				read(line, *) this%bias_type_u, this%bias_type_v
				
				!bias values:
				read(iunit, '(a)') line
				read(line, *) this%delta_u, this%delta_v
			end if
			
				!data reading(optional): 
			!=== 1. BLOCK COORDINATES (ORDERED BY ANGEL): ===
			if (index(line, 'Coordinates:') > 0) then
				call read_nodes(C00)
				call read_nodes(C10)
				call read_nodes(C11)
				call read_nodes(C01)
			end if
			
			!=== 2. BLOCK COORDINATES (CURVES FILENAMES): ===			![bottom top left right]
			if (index(line, 'Edges:') > 0) then
				read(iunit, '(a)') line; call this%edges(1)%read_file(trim(adjustl(line)))
				read(iunit, '(a)') line; call this%edges(2)%read_file(trim(adjustl(line)))
				read(iunit, '(a)') line; call this%edges(3)%read_file(trim(adjustl(line)))
				read(iunit, '(a)') line; call this%edges(4)%read_file(trim(adjustl(line)))
				flag = .true.
			end if
		end do
		close(iunit)
		
		!=== CREATING EDGES (if not created): ===
		if (.not. flag) then
			call this%edges(1)%init([C00(1), C10(1)], [C00(2), C10(2)])	!bottom
			call this%edges(2)%init([C01(1), C11(1)], [C01(2), C11(2)])	!top
			call this%edges(3)%init([C00(1), C01(1)], [C00(2), C01(2)])	!left
			call this%edges(4)%init([C10(1), C11(1)], [C10(2), C11(2)])	!right
		end if
		
		!=== BLOCK MESH GENERATION: ===
		call compute_tfi_block(this%edges(1), this%edges(2), this%edges(3), this%edges(4),&
							   this%Ni, this%Nj, this%coords,&
							   this%bias_type_u, this%bias_type_v,&
							   this%delta_u, this%delta_v)
		
		contains
		subroutine read_nodes(Cij)
			real(8), intent(inout) :: Cij(2)
			read(iunit, '(a)') line
			pos1 = index(trim(line), '(')
			pos2 = index(trim(line), ')')
			aux_line = line(pos1+1:pos2-1)
			read(aux_line, *) Cij(1), Cij(2)
		end subroutine
	end subroutine


!=======================================================================
!==================== STRUCTURED MESH METHODS ==========================
!=======================================================================
	subroutine read_mesh_map(this, filename)
		implicit none
		class(structured_mesh_t), intent(inout) :: this
		character(len=*), intent(in) :: filename
		
		integer :: iunit, ios
		integer :: counter
		integer :: nblocks, nconn, b1, b2, e1_idx, e2_idx, zone_id, i
		character(len=256) :: line, tmp, b_name, e_name, b_name2, e_name2
		integer :: p1, p2, p3, p4
		integer :: max_zone_id
		real(8) :: p1_b1(2), p1_b2(2), p2_b2(2), dist_start, dist_end
		
		open(newunit=iunit, file=filename, status='old')
		
		!=== NUMBER OF BLOCKS: ===
		do
			read(iunit, '(a)', iostat=ios) line
			if (ios /= 0) exit
			if (index(line, '#Blocks_Count:') > 0) then
				read(iunit, '(a)') line
				read(line, *) nblocks
				allocate(this%blocks(nblocks))
				exit
			end if
		end do
		rewind(iunit)
		
		!=== BLOCKS INITIALIZATION: ===
		counter = 1
		do
			read(iunit, '(a)', iostat=ios) line
			if (ios /= 0) exit
			
			if (index(line, '#Block:') > 0) then
				read(iunit, '(a)') line
				read(line, *) tmp, b_name
				call this%blocks(counter)%initialize(b_name)
				this%blocks(counter)%id = counter
				counter = counter + 1
			end if
		end do
		rewind(iunit)
		
		!=== NUMBER OF CONNECTIVITIES: ===
		nconn = 0 
		do
			read(iunit, '(a)', iostat=ios) line
			if (ios /= 0) exit
			if (index(line, '#Connectivity:') > 0) then
				do
					read(iunit, '(a)', iostat=ios) line
					if (index(line, '}') > 0 .or. ios /= 0) exit
					nconn = nconn + 1
				end do
			end if
		end do
		rewind(iunit)
		allocate(this%connections(nconn))
				
		!=== BLOCKS CONNECTIVITIES: ===
		nconn = 0
		do
			read(iunit, '(a)', iostat=ios) line
			if (ios /= 0) exit
			
			if (index(line, '#Connectivity:') > 0) then
				do
					read(iunit, '(a)', iostat=ios) line
					if (index(line, '}') > 0 .or. ios /= 0) exit
					
					p1 = index(line, 'Block:'); p2 = index(line, 'Edge:')
					p3 = index(line, '->'); p4 = index(line, 'Edge:', back=.true.)
					
					read(line(p1+6:p2-2), *) b1
					e_name = line(p2+5:index(line, ')')-1)
					
					read(line(index(line, 'Block:', back=.true.)+6:p4-2), *) b2
					e_name2 = line(p4+5:index(line, ')', back=.true.)-1)
					
					e1_idx = edge_name_to_idx(e_name)
					e2_idx = edge_name_to_idx(e_name2)
					
					this%blocks(b1)%edges(e1_idx)%is_interface = .true.
					this%blocks(b1)%edges(e1_idx)%neighbor_block_id = b2
					
					this%blocks(b2)%edges(e2_idx)%is_interface = .true.
					this%blocks(b2)%edges(e2_idx)%neighbor_block_id = b1
					
					nconn = nconn + 1
					if (b1 < b2) then
						this%connections(nconn)%b1 = b1
						this%connections(nconn)%e1 = e1_idx
						this%connections(nconn)%b2 = b2
						this%connections(nconn)%e2 = e2_idx
					else
						this%connections(nconn)%b1 = b2
						this%connections(nconn)%e1 = e2_idx
						this%connections(nconn)%b2 = b1
						this%connections(nconn)%e2 = e1_idx
					end if
					
					
					p1_b1 = this%blocks(b1)%edges(e1_idx)%get_point(0.0d0)
					p1_b2 = this%blocks(b2)%edges(e2_idx)%get_point(0.0d0)
					p2_b2 = this%blocks(b2)%edges(e2_idx)%get_point(1.0d0)

					dist_start = norm2(p1_b1 - p1_b2)
					dist_end = norm2(p1_b1 - p2_b2)

					if (dist_end < dist_start) then
						this%connections(nconn)%reverse = .true.
					else
						this%connections(nconn)%reverse = .false.
					end if
					
				end do
			end if
		end do
		rewind(iunit)
		
		max_zone_id = -1
		!=== BOUNDARY ZONES INDEXATION: ===
		do
			read(iunit, '(a)', iostat=ios) line
			if (ios /= 0) exit
			
			if (index(line, '#Boundaries:') > 0) then
				do
					read(iunit, '(a)', iostat=ios) line
					if (index(line, '}') > 0 .or. ios /= 0) exit
					
					p1 = index(line, 'Block:'); p2 = index(line, 'Edge:')
					p3 = index(line, '->'); p4 = index(line, 'Zone:')
					
					read(line(p1+6:p2-2), *) b1
					e_name = line(p2+5:index(line, ')')-1)
					e1_idx = edge_name_to_idx(e_name)
					read(line(p4+5:), *) zone_id		
					
					max_zone_id = max(max_zone_id, zone_id)		
					
					this%blocks(b1)%edges(e1_idx)%face_zone = zone_id
					
				end do
			end if
		end do
		rewind(iunit)
		
		!=== BOUNDARY ZONES NAMES: ===
		allocate(this%zone_names(max_zone_id + 2))
		allocate(this%zone_subnames(max_zone_id + 2))
		this%zone_names(1) = 'interior' 
		this%zone_subnames(1) = 'default-interior' 
		do
			read(iunit, '(a)', iostat=ios) line
			if (ios /= 0) exit
			
			if (index(line, '#ZoneNames:') > 0) then
				do
					read(iunit, '(a)', iostat=ios) line
					if (index(line, '}') > 0 .or. ios /= 0) exit
					
					p1 = index(line, 'Zone:'); p2 = index(line, 'Name:')
					p3 = index(line, '->')
					
					read(line(p1+5:p3-2), *) zone_id
					zone_id = zone_id + 2
					
					this%zone_names(zone_id) = 'wall'
					this%zone_subnames(zone_id) = line(p2+5:)
					
				end do
			end if
		end do
		rewind(iunit)
		
		
		close(iunit)
		
		
		!=== MESH GLOBAL INDICES: ===
		allocate(this%node_offsets(nblocks), this%cell_offsets(nblocks))
		this%node_offsets(1) = 0
		this%cell_offsets(1) = 0
		
		do i = 2, nblocks
			this%node_offsets(i) = this%node_offsets(i-1) + this%blocks(i-1)%Ni*this%blocks(i-1)%Nj
			this%cell_offsets(i) = this%cell_offsets(i-1) + (this%blocks(i-1)%Ni-1)*(this%blocks(i-1)%Nj-1)
		end do
		
		contains 
		function edge_name_to_idx(name) result(idx)
			implicit none
			character(len=*), intent(in) :: name
			integer :: idx
			select case(trim(adjustl(name)))
				case('bottom'); idx = 1
				case('top'); idx = 2
				case('left'); idx = 3
				case('right'); idx = 4
				case default; idx = 0
			end select
		end function

	end subroutine
	
	subroutine assemble_to_unstructured(this, ncells, nfaces, nnodes, nbfaces,&
										node_coords, face_nodes, face_nodes_ptr,&
										face_left_cell, face_right_cell,&
										face_zone, face_type, cell_type,&
										zone_names, zone_subnames)
		implicit none
		class(structured_mesh_t), intent(inout) :: this
		integer, intent(inout) :: ncells, nfaces, nnodes, nbfaces
		integer, allocatable, intent(inout) :: face_nodes(:), face_nodes_ptr(:),&
											   face_zone(:), face_type(:), cell_type(:),&
											   face_left_cell(:), face_right_cell(:)
		character(len=128), allocatable, intent(inout) :: zone_names(:), zone_subnames(:)
		real(8), allocatable, intent(inout) :: node_coords(:, :)
		
		real(8) :: p1(2), p2(2)
		integer :: b, i, j, k, n_total_nodes, n_total_cells
		integer :: id1, id2, root1, root2, glob_idx
		integer, allocatable :: node_remap(:), global_node_id(:)
		integer :: max_faces_guess, max_k
		integer :: edge_idx, face_idx
		integer :: Ni, Nj
		integer :: nodes(2)
		integer, allocatable :: interface_face_map(:, :), tmp_int(:)
		
		!=== ASSEMBLING NAMES: ===
		zone_names = this%zone_names
		zone_subnames = this%zone_subnames
		
		!=== ASSMEBLING NODES: ===
		if (allocated(node_coords)) deallocate(node_coords)
		!=== NUMBER OF RAW NODES: ===
		n_total_nodes = this%node_offsets(size(this%blocks)) +&
						this%blocks(size(this%blocks))%Ni*this%blocks(size(this%blocks))%Nj
		
		!=== NUMBER OF CELLS: ===
		n_total_cells = this%cell_offsets(size(this%blocks)) +&
						(this%blocks(size(this%blocks))%Ni-1)*(this%blocks(size(this%blocks))%Nj-1)
		ncells = n_total_cells

		!=== REMAP INITIALIZATION (NON UNIQUE NODES): ===
		allocate(node_remap(n_total_nodes))
		do i = 1, n_total_nodes
			node_remap(i) = i
		end do

		!=== LINKING INTERFACE NODES: ===
		do k = 1, size(this%connections)
			call link_interface_nodes(this%connections(k), node_remap)
		end do

		!=== GENERATION UNIQUE NODES INDICES: ===
		allocate(global_node_id(n_total_nodes))
		nnodes = 0
		global_node_id = -1
		do i = 1, n_total_nodes
			if (node_remap(i) == i) then
				nnodes = nnodes + 1
				global_node_id(i) = nnodes
			end if
		end do
		
		!=== REMAP UPDATE: ===
		do i = 1, n_total_nodes
			global_node_id(i) = global_node_id(find_root(i, node_remap))
		end do

		!=== FINAL MESH NODE COORDINATES: ===
		allocate(node_coords(2, nnodes))
		do b = 1, size(this%blocks)
			do j = 1, this%blocks(b)%Nj
				do i = 1, this%blocks(b)%Ni
					glob_idx = get_glb_idx(b, i, j)
					node_coords(:, global_node_id(glob_idx)) = this%blocks(b)%coords(:, i, j)
				end do
			end do
		end do
		
		
		
		!=== ASSEMBLING FACES: ===
		if (allocated(face_nodes)) deallocate(face_nodes)
		if (allocated(face_left_cell)) deallocate(face_left_cell)
		if (allocated(face_right_cell)) deallocate(face_right_cell)
		if (allocated(face_zone)) deallocate(face_zone)		
			
		max_faces_guess = 0
		do b = 1, size(this%blocks)
			max_faces_guess = max_faces_guess +&
							  this%blocks(b)%Ni*(this%blocks(b)%Nj-1) +&
							  this%blocks(b)%Nj*(this%blocks(b)%Ni-1)
		end do
		
		allocate(face_nodes(2*max_faces_guess))
		allocate(face_left_cell(max_faces_guess))
		allocate(face_right_cell(max_faces_guess))
		allocate(face_zone(max_faces_guess))
		
		!=== INTERFACE EDGES MAP: ===
		max_k = 0
		do b = 1, size(this%blocks)
			max_k = max(max_k, max(this%blocks(b)%Ni, this%blocks(b)%Nj))
		end do
		allocate(interface_face_map(max_k, size(this%connections)))
		interface_face_map = 0
		
		face_idx = 0
		nbfaces = 0

		!=== LOOP OVER BLOCKS: ===
		do b = 1, size(this%blocks)
			Ni = this%blocks(b)%Ni; Nj = this%blocks(b)%Nj

			!=== I-DIRECTION FACES: ===
			do j = 1, Nj - 1
				do i = 1, Ni
					edge_idx = 0
					if (i == 1)  edge_idx = 3 							!left
					if (i == Ni) edge_idx = 4							!right
					
					call process_face_logic(b, edge_idx, j,&
										    get_glb_idx(b, i, j), get_glb_idx(b, i, j+1),&
										    get_cell_idx(b, i-1, j), get_cell_idx(b, i, j))
				end do
			end do

			!=== J-DIRECTION FACES: ===
			do j = 1, Nj
				do i = 1, Ni - 1
					edge_idx = 0
					if (j == 1)  edge_idx = 1 							!bottom
					if (j == Nj) edge_idx = 2 							!top
					
					call process_face_logic(b, edge_idx, i,&
										    get_glb_idx(b, i, j), get_glb_idx(b, i+1, j),&
										    get_cell_idx(b, i, j-1), get_cell_idx(b, i, j))
				end do
			end do
		end do

		nfaces = face_idx
		!=== FACE NODES PTR FILLING: ===
		if (allocated(face_nodes_ptr)) deallocate(face_nodes_ptr)
		allocate(face_nodes_ptr(nfaces+1))
		face_nodes_ptr(1) = 1
		do i = 2, nfaces+1
			face_nodes_ptr(i) = face_nodes_ptr(i-1) + 2
		end do
		
		!=== FACE NODES FILLING: ===
		allocate(tmp_int(face_nodes_ptr(nfaces+1) - 1))
		tmp_int = face_nodes(1:face_nodes_ptr(nfaces+1) - 1)
		call move_alloc(tmp_int, face_nodes)
		
		!=== FACE LEFT/RIGHT CELLS FILLING: ===
		allocate(tmp_int(nfaces))
		tmp_int = face_left_cell(1:nfaces)
		call move_alloc(tmp_int, face_left_cell)
		
		allocate(tmp_int(nfaces))
		tmp_int = face_right_cell(1:nfaces)
		call move_alloc(tmp_int, face_right_cell)
		
		!=== FACE ZONE FILLING: ===
		allocate(tmp_int(nfaces))
		tmp_int = face_zone(1:nfaces)
		call move_alloc(tmp_int, face_zone)
		
		if (allocated(face_type)) deallocate(face_type)
		allocate(face_type(nfaces)) 
		face_type = 3
		
		
		call reorder_faces(ncells, nfaces, nbfaces, face_left_cell, face_right_cell,&
						   face_zone, face_type, face_nodes_ptr, face_nodes)
		call order_boundary_faces_by_zone(2, ncells, nfaces, nbfaces,&
                                          face_left_cell, face_right_cell,&
                                          face_zone, face_type,&
                                          face_nodes_ptr = face_nodes_ptr,&
                                          face_nodes = face_nodes)
		
		!=== FACE NODES ORDERING: ===
		block
			integer, allocatable :: cell_nodes(:), cell_nodes_ptr(:)
			integer, allocatable :: cell_faces(:), cell_faces_ptr(:)
			integer :: left_cell, right_cell, tmp
			integer :: node_start, node_end, n_nodes_in_cell
			integer :: current_cell_node, next_cell_node
			logical :: correct_order
			
			call cell_face_connectivity(ncells, nfaces,&
									    face_left_cell, face_right_cell,&
									    cell_faces_ptr, cell_faces)
			
			call cell_node_connectivity(nnodes, ncells,&
										face_nodes_ptr, face_nodes,&
										cell_faces_ptr, cell_faces,&
										cell_nodes_ptr, cell_nodes)
										
			call order_cell_nodes_2d(ncells, node_coords,&
									 cell_nodes_ptr, cell_nodes) 
										
			do i = 1, nfaces
				left_cell = face_left_cell(i)

				id1 = face_nodes(face_nodes_ptr(i))
				id2 = face_nodes(face_nodes_ptr(i)+1)

				node_start = cell_nodes_ptr(left_cell)
				node_end = cell_nodes_ptr(left_cell + 1) - 1
				n_nodes_in_cell = node_end - node_start + 1

				correct_order = .false.

				do j = 0, n_nodes_in_cell - 1
					current_cell_node = cell_nodes(node_start + j)
					
					if (current_cell_node == id1) then
						next_cell_node = cell_nodes(node_start + mod(j + 1, n_nodes_in_cell))
						
						if (next_cell_node == id2) then
							correct_order = .true.
						end if
						exit
					end if
				end do

				if (.not. correct_order) then
					tmp = face_nodes(face_nodes_ptr(i))
					face_nodes(face_nodes_ptr(i)) = face_nodes(face_nodes_ptr(i)+1)
					face_nodes(face_nodes_ptr(i)+1) = tmp
				end if
			end do
		
		end block	
		
		if (allocated(cell_type)) deallocate(cell_type)
		allocate(cell_type(ncells))
		cell_type = 3
				
		contains
		!=== RAW GLOBAL NODE INDEX: ===
		function get_glb_idx(b_idx, i, j) result(res)
			integer, intent(in) :: b_idx, i, j
			integer :: res
	
			res = this%node_offsets(b_idx) + (j - 1)*this%blocks(b_idx)%Ni + i
		end function

		!=== FINDING ROOT NODE INDEX: ===
		recursive function find_root(id, remap) result(root)
			integer, intent(in) :: id
			integer, intent(in) :: remap(:)
			integer :: root
			
			if (remap(id) == id) then
				root = id
			else
				root = find_root(remap(id), remap)
			end if
		end function

		!=== INTERFACE NODES LINKING SUBROUTINE: ===
		subroutine link_interface_nodes(conn, remap)
			type(connectivity_t), intent(in) :: conn
			integer, intent(inout) :: remap(:)
			
			integer :: b1, b2, e1, e2, k, k_eff, n_pts
			integer :: i1, j1, i2, j2
			integer :: id1, id2, root1, root2

			b1 = conn%b1; e1 = conn%e1
			b2 = conn%b2; e2 = conn%e2
			
			!=== FINDING NUMBER OF INTERFACES EDGE POINTS: ===
			if (e1 <= 2) then
				!bottom/top
				n_pts = this%blocks(b1)%Ni
			else
				!left/right
				n_pts = this%blocks(b1)%Nj
			end if


			!=== INTERFACE NODES LOOP: ===
			do k = 1, n_pts
				!=== OWNER: ===
				select case(e1)
					case(1); i1 = k; j1 = 1                     		!bottom
					case(2); i1 = k; j1 = this%blocks(b1)%Nj     		!top
					case(3); i1 = 1; j1 = k                      		!left
					case(4); i1 = this%blocks(b1)%Ni; j1 = k     		!right
				end select

				!=== NEIGHBOR: ===
				if (conn%reverse) then
					k_eff = n_pts - k + 1
				else
					k_eff = k
				end if
				select case(e2)
					case(1); i2 = k_eff; j2 = 1
					case(2); i2 = k_eff; j2 = this%blocks(b2)%Nj
					case(3); i2 = 1; j2 = k_eff
					case(4); i2 = this%blocks(b2)%Ni; j2 = k_eff
				end select

				id1 = get_glb_idx(b1, i1, j1)
				id2 = get_glb_idx(b2, i2, j2)
				
				!=== FINDING ROOT NODE (WITH MIN IDX): ===
				root1 = find_root(id1, remap)
				root2 = find_root(id2, remap)
				
				if (root1 /= root2) then
					remap(root2) = root1
				end if
			end do
		end subroutine
		
		subroutine process_face_logic(b_id, e_idx, k_idx, n1_raw, n2_raw, c_left_loc, c_right_loc)
			integer, intent(in) :: b_id, e_idx, k_idx, n1_raw, n2_raw, c_left_loc, c_right_loc
			integer :: n1, n2, cl, cr, conn_id, f_id, k_eff, n_intervals
			
			!=== GETTING GLOBAL NODE INDICES: ===
			n1 = global_node_id(n1_raw)
			n2 = global_node_id(n2_raw)
			
			!=== GETTING GLOBAL CELL INDICES: ===
			cl = 0; cr = 0
			if (c_left_loc  > 0) cl = c_left_loc  + this%cell_offsets(b_id)
			if (c_right_loc > 0) cr = c_right_loc + this%cell_offsets(b_id)

			!=== INTERNAL FACE: ===
			if (e_idx == 0) then
				face_idx = face_idx + 1
				face_left_cell(face_idx) = cl
				face_right_cell(face_idx) = cr
				nodes(1) = n1; nodes(2) = n2
				face_nodes(2*face_idx-1:2*face_idx) = nodes
				face_zone(face_idx) = 1
				return
			end if

			!=== EXTERNAL FACE: ===
			if (.not. this%blocks(b_id)%edges(e_idx)%is_interface) then
				face_idx = face_idx + 1
				nbfaces = nbfaces + 1
				face_left_cell(face_idx) = max(cl, cr)
				face_right_cell(face_idx) = 0
				nodes(1) = n1; nodes(2) = n2
				face_nodes(2*face_idx-1:2*face_idx) = nodes
				face_zone(face_idx) = this%blocks(b_id)%edges(e_idx)%face_zone + 2
				return
			end if

			!=== INTERFACE FACE: ===
			conn_id = find_conn_id(b_id, e_idx)
			k_eff = k_idx
			
			if (this%connections(conn_id)%b2 == b_id .and. this%connections(conn_id)%reverse) then
				if (e_idx <= 2) then
					n_intervals = this%blocks(b_id)%Ni - 1
				else
					n_intervals = this%blocks(b_id)%Nj - 1
				end if
				k_eff = n_intervals - k_idx + 1
			end if

			f_id = interface_face_map(k_eff, conn_id)
			
			if (f_id == 0) then
				!=== NEW FACE: ===
				face_idx = face_idx + 1
				interface_face_map(k_eff, conn_id) = face_idx
				face_left_cell(face_idx) = max(cl, cr)
				nodes(1) = n1; nodes(2) = n2
				face_nodes(2*face_idx-1:2*face_idx) = nodes
				face_zone(face_idx) = 1
			else
				face_right_cell(f_id) = max(cl, cr)
			end if
			
		end subroutine
		
		!=== FINDING CONNECTIVITY ID: ===
		function find_conn_id(b, e) result(res)
			integer, intent(in) :: b, e
			integer :: res, k
			res = 0
			do k = 1, size(this%connections)
				if ((this%connections(k)%b1 == b .and. this%connections(k)%e1 == e) .or. &
					(this%connections(k)%b2 == b .and. this%connections(k)%e2 == e)) then
					res = k
					return
				end if
			end do
		end function
		
		!=== CELL LOCAL INDEX: ===
		function get_cell_idx(b, i, j) result(res)
			integer, intent(in) :: b, i, j
			integer :: res
			if (i < 1 .or. i >= this%blocks(b)%Ni .or. j < 1 .or. j >= this%blocks(b)%Nj) then
				res = 0
			else
				res = i + (j-1)*(this%blocks(b)%Ni-1)
			end if
		end function

	end subroutine


!=======================================================================
!===================== TRANSFINITE INTERPOLATION =======================
!=======================================================================
	subroutine compute_tfi_block(bot, top, left, right, Ni, Nj, coords,&
								 bias_type_u, bias_type_v, delta_u, delta_v)
		implicit none
		type(boundary_t), intent(in) :: bot, top, left, right			!=== BLOCK CURVES ===
		real(8), intent(in) :: delta_u, delta_v							!=== BIAS FACTORS ===
		character(len=*), intent(in) :: bias_type_u, bias_type_v		!=== BIAS TYPE ===
		integer, intent(in) :: Ni, Nj 									!=== NUMBER OF NODES ===
		real(8), allocatable, intent(inout) :: coords(:, :, :) 			!=== NODES COORDINATES ===
		
		integer :: i, j
		real(8) :: u, u_linear, v, v_linear
		real(8) :: P_bot(2), P_top(2), P_left(2), P_right(2)
		real(8) :: C00(2), C10(2), C01(2), C11(2)
		real(8) :: B_c(2), H_c(2), V_c(2)
		
		if (allocated(coords)) deallocate(coords)
		allocate(coords(2, Ni, Nj))
		

		!=== CORNER NODES COORDINATES: ===
		C00 = bot%get_point(0.0d0) 										!=== LEFT+BOTTOM ===
		C10 = bot%get_point(1.0d0) 										!=== RIGHT+BOTTOM ===
		C01 = top%get_point(0.0d0) 										!=== LEFT+TOP ===
		C11 = top%get_point(1.0d0) 										!=== RIGHT+TOP ===

		!=== FINDING GRID NODES: ===
		do j = 1, Nj
			v_linear = real(j-1)/real(Nj-1)
			v = stretch_param(v_linear, delta_v, bias_type_v)
			
			P_left  = left%get_point(v)
			P_right = right%get_point(v)

			do i = 1, Ni
				u_linear = real(i-1)/real(Ni-1)
				u = stretch_param(u_linear, delta_u, bias_type_u)
				
				P_bot = bot%get_point(u)
				P_top = top%get_point(u)

				!=== TRANSFINITE INTERPOLATION: ===
				H_c = (1.0d0 - v)*P_bot + v*P_top
				V_c = (1.0d0 - u)*P_left + u*P_right
				
				B_c = (1.0d0 - u)*(1.0d0 - v)*C00 + u*(1.0d0 - v)*C10 +&
					  (1.0d0 - u)*v*C01 + u*v*C11

				coords(:, i, j) = H_c + V_c - B_c
			end do
		end do
		
		contains
		function stretch_param(t, delta, bias_type) result(t_new)
			implicit none
			real(8), intent(in) :: t, delta
			character(len=*), intent(in) :: bias_type
			
			real(8) :: t_new
			
			t_new = t
			if (dabs(delta) < 1.0d-6) return
			
			select case(trim(adjustl(bias_type)))
			case('tanh')
				t_new = 0.5d0*(1.0d0 + tanh(delta*(2.0d0*t - 1.0d0))/tanh(delta))
			case('exp')
				if (delta > 0.0d0) then
					t_new = (dexp(delta*t) - 1.0d0)/(dexp(delta) - 1.0d0)
				else
					t_new = 1.0d0 - (dexp(dabs(delta)*(1.0d0 - t)) - 1.0d0)/(dexp(dabs(delta)) - 1.0d0)
				end if
			end select
			
	
		end function
	
	end subroutine


!=======================================================================
!=================== RECOMPUTE USTRUCTURED MESH ========================
!=======================================================================
	subroutine recompute_unstrd_nodes(coords, nnodes, node_coords)
		implicit none
		real(8), intent(in) :: coords(:, :, :)
		real(8), allocatable, intent(inout) :: node_coords(:, :)
		integer, intent(inout) :: nnodes
		
		integer :: Ni, Nj, i, j
		integer :: global_idx
		
		Ni = size(coords, dim=2); Nj = size(coords, dim=3)
		nnodes = Ni*Nj
		if (allocated(node_coords)) deallocate(node_coords)
		allocate(node_coords(2, nnodes))
		
		do j = 1, Nj
			do i = 1, Ni
				global_idx = i + (j - 1)*Ni
				
				node_coords(:, global_idx) = coords(:, i, j)
			end do
		end do
		
	
	end subroutine

	subroutine recompute_unstrd_faces(coords, ncells, nfaces, face_nodes, face_nodes_ptr,&
									  face_left_cell, face_right_cell, face_zone)
		implicit none
		real(8), intent(in) :: coords(:, :, :)
		integer, intent(inout) :: ncells, nfaces
		integer, allocatable :: face_nodes(:), face_nodes_ptr(:), face_zone(:),&
								face_left_cell(:), face_right_cell(:)
		
		integer :: i, j, Ni, Nj, counter
		integer :: node1_idx, node2_idx, left_cell, right_cell, face_idx
		
		Ni = size(coords, dim=2); Nj = size(coords, dim=3)
		nfaces = (Ni-1)*Nj + (Nj-1)*Ni
		ncells = (Ni-1)*(Nj-1)
		
		if (allocated(face_nodes)) deallocate(face_nodes)
		if (allocated(face_nodes_ptr)) deallocate(face_nodes_ptr)
		if (allocated(face_zone)) deallocate(face_zone)
		if (allocated(face_left_cell)) deallocate(face_left_cell)
		if (allocated(face_right_cell)) deallocate(face_right_cell)
		
		allocate(face_nodes(2*nfaces), face_nodes_ptr(nfaces+1), face_zone(nfaces),&
			     face_left_cell(nfaces), face_right_cell(nfaces))
		
		face_nodes_ptr(1) = 1
		do i = 1, nfaces
			face_nodes_ptr(i + 1) = face_nodes_ptr(i) + 2
		end do
		
		face_idx = 0
		
		!=== I-DIRECTION FACES: ===
		do j = 1, Nj - 1
			do i = 1, Ni
				node1_idx = i + (j - 1)*Ni
				node2_idx = i + j*Ni
				
				left_cell = (i - 1) + (j - 1)*(Ni - 1)
				right_cell = i + (j - 1)*(Ni - 1)
				
				face_idx = face_idx + 1
				face_nodes(face_nodes_ptr(face_idx)) = node1_idx
				face_nodes(face_nodes_ptr(face_idx) + 1) = node2_idx
				
				face_left_cell(face_idx) = left_cell
				face_right_cell(face_idx) = right_cell
				
				face_zone(face_idx) = 0
				
				if (i == 1) then
					face_zone(face_idx) = -1
					face_left_cell(face_idx) = 0
				end if

				if (i == Ni) then
					face_zone(face_idx) = -1
					face_right_cell(face_idx) = 0
				end if
			end do
		end do
		
		!=== J-DIRECTION FACES: ===
		do j = 1, Nj
			do i = 1, Ni - 1
				node1_idx = i + (j - 1)*Ni
				node2_idx = i + 1 + (j - 1)*Ni
				
				left_cell = i + (j - 2)*(Ni - 1)
				right_cell = i + (j - 1)*(Ni - 1)
				
				face_idx = face_idx + 1
				face_nodes(face_nodes_ptr(face_idx)) = node1_idx
				face_nodes(face_nodes_ptr(face_idx) + 1) = node2_idx
				
				face_left_cell(face_idx) = left_cell
				face_right_cell(face_idx) = right_cell
				
				face_zone(face_idx) = 0
				if (j == 1) then
					face_zone(face_idx) = -1
					face_left_cell(face_idx) = 0
				end if

				if (j == Nj) then
					face_zone(face_idx) = -1
					face_right_cell(face_idx) = 0
				end if
			end do
		end do
	
	
	end subroutine


end module
