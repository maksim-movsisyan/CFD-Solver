module linear_solver_module
use linear_operator_module
use preconditioner_module
implicit none
!=======================================================================
!================= LINEAR SOLVER BASIC DATA STRUCTURE ==================
!=======================================================================
type, abstract :: linear_solver_t
	class(linear_operator_t), pointer :: A => null()
	class(preconditioner_t), pointer :: M => null()
	
	integer :: verbosity = -1 
	real(8) :: rel_tol = 1.0d-2
	real(8) :: abs_tol = 1.0d-6
	integer :: max_iter = 100
	integer :: iterations = 0
	real(8) :: init_residual = 0.0d0
	real(8) :: final_residual = 0.0d0
	real(8) :: residual_drop = 0.0d0
	logical :: converged = .false.
contains
	procedure(solve_interface), deferred :: solve
	procedure :: set_operator => solver_set_operator
	procedure :: set_preconditioner => solver_set_preconditioner
end type

abstract interface
	subroutine solve_interface(this, b, x)
		import :: linear_solver_t
		class(linear_solver_t), intent(inout) :: this
		real(8), intent(in), contiguous :: b(:)
		real(8), intent(inout), contiguous :: x(:)
	end subroutine
end interface



!=======================================================================
!========================= GMRES DATA STRUCTURE ========================
!=======================================================================
type, extends(linear_solver_t) :: gmres_solver_t
	integer :: restart = 30
	real(8) :: beta = 0.0d0
    integer :: outer_iter = 0
    
	real(8), allocatable :: V(:,:)    									! === Krylov Basis, DIM = (n x restart+1) ===
	real(8), allocatable :: H(:,:)    									! === Hessenbergs matrix, DIM = (restart+1 x restart) ===
	real(8), allocatable :: g(:)      									! === RSH in Krylovs Basis, DIM = (restart+1) ===
	real(8), allocatable :: c(:)      									! === Givenss rotantin - cos, DIM = (restart) ===
	real(8), allocatable :: s(:)                                        ! === Givenss rotantin - sin, DIM = (restart) ===
	real(8), allocatable :: y(:)      									! === Solution in Krylovs Basis, DIM = (restart) ===
	real(8), allocatable :: w(:)      									! === Work array, DIM = (n) ===
	real(8), allocatable :: w2(:)										! === Work array2, DIM = (n) ===
	real(8), allocatable :: r(:)      									! === Residual, DIM = (n) ===
	real(8), allocatable :: z(:)										! === Temp vector for peconditioning, DIM = (n) ===
contains
	procedure :: solve => gmres_solve
	procedure :: solve_no_restart => gmres_no_restart
	procedure :: allocate_workspace => gmres_allocate_workspace
	procedure :: deallocate_workspace => gmres_deallocate_workspace
	procedure :: print_stats => gmres_print_stats
end type



!=======================================================================
!========================= BiCGSTAB DATA STRUCTURE =====================
!=======================================================================
type, extends(linear_solver_t) :: bicgstab_solver_t


	real(8), allocatable :: p(:) 
	real(8), allocatable :: v(:) 
	real(8), allocatable :: t(:) 
	real(8), allocatable :: s(:) 
	real(8), allocatable :: r(:) 
	real(8), allocatable :: r_hat(:) 
	real(8), allocatable :: y_p(:) 
	real(8), allocatable :: y_s(:) 
	contains
	procedure :: solve => bicgstab_solve
	procedure :: allocate_workspace => bicgstab_allocate_workspace
	procedure :: deallocate_workspace => bicgstab_deallocate_workspace
	procedure :: print_stats => bicgstab_print_stats
	
end type



!=======================================================================
!======================== JACOBI DATA STRUCTURE ========================
!=======================================================================
type, extends(linear_solver_t) :: jacobi_solver_t
	real(8) :: omega = 1.d0
	real(8), allocatable :: r(:)
	real(8), allocatable :: diag(:)
	
	contains
	procedure :: solve => jacobi_solve	
	procedure :: allocate_workspace => jacobi_allocate_workspace
	procedure :: deallocate_workspace => jacobi_deallocate_workspace
end type



contains
!=======================================================================
!===================== LINEAR SOLVER BASIC METHODS =====================
!=======================================================================
	subroutine solver_set_operator(this, A)
		implicit none
		class(linear_solver_t), intent(inout) :: this
		class(linear_operator_t), target, intent(in) :: A
		
		this%A => A
	end subroutine

	subroutine solver_set_preconditioner(this, M)
		implicit none
		class(linear_solver_t), intent(inout) :: this
		class(preconditioner_t), target, intent(in) :: M
		
		this%M => M
	end subroutine




!=======================================================================
!============================== GMRES METHODS ==========================
!=======================================================================
	subroutine gmres_solve(this, b, x)
		implicit none
        class(gmres_solver_t), intent(inout) :: this
        real(8), intent(in), contiguous :: b(:)
        real(8), intent(inout), contiguous :: x(:)
        
        integer :: n, m, i, j, k, outer, info
        real(8) :: res, res0, res_prev
        real(8) :: hh, h1, h2, cs, sn, rho, delta
        real(8) :: time_start, time_end
        
        call cpu_time(time_start)
        
        !=== INITIALIZATION: ===
        n = size(b)
        m = this%restart
        
        if (.not. allocated(this%V)) then
            call this%allocate_workspace(n, m)
        end if
                
        this%converged = .false.
        this%iterations = 0
        this%outer_iter = 0
        
		!=== INITIAL RESIDUAL: ===
		call this%A%apply(x, this%r)
		do k = 1, n
			this%r(k) = b(k) - this%r(k)
		end do
        res0 = norm2(this%r)
        this%init_residual = res0
        
        res = 1.d0
		this%final_residual = res0
        this%residual_drop = res
        
        
        if (res0 < this%abs_tol) then
            this%converged = .true.
            if (this%verbosity > 0) then
                print *, 'GMRES: Converged without iterations'
            end if
            call cpu_time(time_end)
            if (this%verbosity > 2) then
                print *, 'GMRES: Time =', time_end - time_start, 'seconds'
            end if
            return
        end if
        
        !=== OUTER LOOP: ===
        outer_loop: do outer = 1, this%max_iter/m
            this%outer_iter = outer
            
            !=== KRTLOVS SUBSPACE INITIALIZATION ===
            do k = 1, n
				this%V(k, 1) = this%r(k)/res0
			end do
            this%g = 0.0d0
            this%g(1) = res0
            
            !=== HESSENBERGS MATRIX INITIALIZATION ===
            this%H = 0.0d0
            this%c = 0.0d0
            this%s = 0.0d0
            
            !=== INNER LOOP: ===
            do j = 1, m
                !=== PRECONDITIONING: w = M^{-1}*V(:,j) ===
                if (associated(this%M)) then
                    call this%M%apply(this%V(:, j), this%w)
                else
                    this%w = this%V(:, j)
                end if
                
                !=== MATRIX-VECOTR PRODUCT: w = A*w ===
                call this%A%apply(this%w, this%w2)
                this%w = this%w2					
                
                !=== MODIFIED GRAM-SCHMIDT ORTHOGONALIZATION: ===
                do i = 1, j
					this%H(i, j) = 0.d0
					do k = 1, n
						this%H(i, j) = this%H(i, j) + this%w(k)*this%V(k, i)
					end do

                    do k = 1, n
						this%w(k) = this%w(k) - this%H(i, j)*this%V(k, i)
                    end do
                end do
                
                this%H(j+1, j) = norm2(this%w)
                
                !=== BREAKDOWN CHECKING: ===
                if (abs(this%H(j+1, j)) < 1.0d-14) then
                    if (this%verbosity > 1) then
                        print *, 'GMRES: Arnoldi breakdown at iteration', j
                    end if
                    m = j
                    exit
                end if
                
                !=== NEXT BASIS VECTOR: ===
                do k = 1, n
					this%V(k, j+1) = this%w(k)/this%H(j+1, j)
				end do
                
                !=== GIVENSS ROTATION: ===
                do i = 1, j-1
                    h1 = this%H(i, j)
                    h2 = this%H(i+1, j)
                    this%H(i, j) = this%c(i)*h1 + this%s(i)*h2
                    this%H(i+1, j) = -this%s(i)*h1 + this%c(i)*h2
                end do
                
                hh = sqrt(this%H(j, j)**2 + this%H(j+1, j)**2)
                
                if (abs(hh) < 1.0d-14) then
                    this%c(j) = 1.0d0
                    this%s(j) = 0.0d0
                else
                    this%c(j) = this%H(j, j)/hh
                    this%s(j) = this%H(j+1, j)/hh
                end if
                
                this%H(j, j) = hh
                this%H(j+1, j) = 0.0d0
                
                h1 = this%g(j)
                h2 = 0.0d0  
                this%g(j) = this%c(j)*h1 + this%s(j)*h2
                this%g(j+1) = -this%s(j)*h1 + this%c(j)*h2
                
                !=== RESIDUAL EVALUTAION: ===
                res = abs(this%g(j+1))/res0
                this%final_residual = abs(this%g(j+1))
                this%iterations = (outer-1)*m + j
                
                !=== CONVEGENCE INFORMATION: ===
                if (this%verbosity > 1) then
                    if (mod(j, 10) == 0) then
                        print *, 'GMRES: Outer', outer, 'Inner', j, 'Residual drop =', res
                    end if
                end if
                
                !=== CONVEGENCE CHECKING: ===
                if (res < this%rel_tol) then
                    call backward_substitution(this%H(1:j, 1:j), &
                                              this%g(1:j), this%y(1:j), j)
                    
                    this%z = 0.0d0
					do i = 1, j
						this%z = this%z + this%y(i)*this%V(:, i)
					end do
                
                    if (associated(this%M)) then
						!=== w = M^{-1}*z ===
						call this%M%apply(this%z, this%w)
						do k = 1, n
							x(k) = x(k) + this%w(k)
						end do
					else
						do k = 1, n
							x(k) = x(k) + this%z(k)
						end do
					end if
					
                    
                    this%converged = .true.
                    
                    if (this%verbosity > 0) then
                        print *, 'GMRES: Converged in', this%iterations, 'iterations'
                        print *, 'GMRES: Final residual drop=', res
                    end if
                    
                    call cpu_time(time_end)
                    if (this%verbosity > 2) then
                        print *, 'GMRES: Time =', time_end - time_start, 'seconds'
                    end if
                    
                    return
                end if
            end do
            
            call backward_substitution(this%H(1:m, 1:m), this%g(1:m), &
                                      this%y(1:m), m)
            
            this%z = 0.0d0
			do i = 1, m
				do k = 1, n
					this%z(k) = this%z(k) + this%y(i)*this%V(k, i)
				end do
			end do
		
			if (associated(this%M)) then
				!=== w = M^{-1}*z ===
				call this%M%apply(this%z, this%w)
				do k = 1, n
					x(k) = x(k) + this%w(k)
				end do
			else
				do k = 1, n
					x(k) = x(k) + this%z(k)
				end do
			end if
					
            call this%A%apply(x, this%r)
            do k = 1, n
				this%r(k) = b(k) - this%r(k)
			end do
            
            res0 = norm2(this%r)
            res = res0/this%init_residual
            this%final_residual = res0
            
            if (res < this%rel_tol) then
                this%converged = .true.
                this%iterations = outer*m
                
                if (this%verbosity > 0) then
                    print *, 'GMRES: Converged after restart', outer
                    print *, 'GMRES: Final residual drop=', res
                end if
                
                call cpu_time(time_end)
                if (this%verbosity > 2) then
                    print *, 'GMRES: Time =', time_end - time_start, 'seconds'
                end if
                
                return
            end if
            
            if (this%iterations >= this%max_iter) then
                if (this%verbosity > 0) then
                    print *, 'GMRES: Maximum iterations reached'
                    print *, 'GMRES: Final residual drop=', res
                end if
                exit outer_loop
            end if
            
            res_prev = res
            
        end do outer_loop
        
        if (.not. this%converged) then
            if (this%verbosity > 0) then
                print *, 'GMRES: Failed to converge in', this%max_iter, 'iterations'
                print *, 'GMRES: Final residual drop=', res
            end if
        end if
        
        call cpu_time(time_end)
        if (this%verbosity > 2) then
            print *, 'GMRES: Total time =', time_end - time_start, 'seconds'
        end if
        
    end subroutine
    
    
    !=================================================================== 
    !=================== AUXILIARY GMRES SUBROUTINES =================== 
    !=================================================================== 
    subroutine backward_substitution(U, b, x, n)
        implicit none
        real(8), intent(in) :: U(n, n), b(n)
        real(8), intent(out) :: x(n)
        integer, intent(in) :: n
        
        integer :: i, j
        real(8) :: sum_val
        
        do i = n, 1, -1
            sum_val = b(i)
            do j = i+1, n
                sum_val = sum_val - U(i, j)*x(j)
            end do
            
            x(i) = sum_val/U(i, i)
        end do
        
    end subroutine
    
    subroutine generate_givens_rotation(a, b, c, s)
        implicit none
        real(8), intent(in) :: a, b
        real(8), intent(out) :: c, s
        
        real(8) :: r, t
        
        if (abs(b) < 1.0d-14) then
            c = 1.0d0
            s = 0.0d0
        else
            if (abs(b) > abs(a)) then
                t = a/b
                r = sqrt(1.0d0 + t*t)
                s = 1.0d0/r
                c = t * s
            else
                t = b/a
                r = sqrt(1.0d0 + t*t)
                c = 1.0d0/r
                s = t*c
            end if
        end if
        
    end subroutine
    
    subroutine apply_givens_rotation(c, s, x, y)
        implicit none
        real(8), intent(in) :: c, s
        real(8), intent(inout) :: x, y
        
        real(8) :: temp
        
        temp = c*x + s*y
        y = -s*x + c*y
        x = temp
    end subroutine
    
    
    !=================================================================== 
    !=================== AUXILIARY GMRES SUBROUTINES 2 =================
    !=================================================================== 
    subroutine gmres_allocate_workspace(this, n, m)
		implicit none
        class(gmres_solver_t), intent(inout) :: this
        integer, intent(in) :: n   
        integer, intent(in) :: m   
        
        call this%deallocate_workspace()
        
        allocate(this%V(n, m+1))       
        allocate(this%H(m+1, m))       
        allocate(this%g(m+1))          
        allocate(this%c(m))          
        allocate(this%s(m))            
        allocate(this%y(m))            
        allocate(this%w(n))
        allocate(this%w2(n)) 
        allocate(this%r(n))       
        allocate(this%z(n))       
        
        this%V = 0.0d0
        this%H = 0.0d0
        this%g = 0.0d0
        this%c = 0.0d0
        this%s = 0.0d0
        this%y = 0.0d0
        this%w = 0.0d0
        this%w2 = 0.0d0
        this%r = 0.0d0
        this%z = 0.0d0
    end subroutine
    
    subroutine gmres_deallocate_workspace(this)
		implicit none
        class(gmres_solver_t), intent(inout) :: this
        
        if (allocated(this%V)) deallocate(this%V)
        if (allocated(this%H)) deallocate(this%H)
        if (allocated(this%g)) deallocate(this%g)
        if (allocated(this%c)) deallocate(this%c)
        if (allocated(this%s)) deallocate(this%s)
        if (allocated(this%y)) deallocate(this%y)
        if (allocated(this%w)) deallocate(this%w)
        if (allocated(this%w2)) deallocate(this%w2)
        if (allocated(this%r)) deallocate(this%r)
        if (allocated(this%z)) deallocate(this%z)
        
    end subroutine
    
    subroutine gmres_print_stats(this)
		implicit none
        class(gmres_solver_t), intent(in) :: this
        
        print *, "========================================="
        print *, "GMRES Solver Statistics:"
        print *, "  Converged: ", this%converged
        print *, "  Iterations: ", this%iterations
        print *, "  Outer restarts: ", this%outer_iter
        print *, "  Initial residual: ", this%init_residual
        print *, "  Final residual: ", this%final_residual
        print *, "  Residual drop: ", this%residual_drop
        print *, "  Target relative tolerance: ", this%rel_tol
        print *, "  Target absolute tolerance: ", this%rel_tol
        print *, "  Restart size: ", this%restart
        print *, "  Max iterations: ", this%max_iter
        print *, "========================================="
        
    end subroutine
    
    
	!===================================================================
	!=================== GMRES WITH NO RESTARTERS ====================== 
	!===================================================================
    subroutine gmres_no_restart(this, b, x)
        implicit none
        class(gmres_solver_t), intent(inout) :: this
        real(8), intent(in), contiguous :: b(:)
        real(8), intent(inout), contiguous :: x(:)
        
        integer :: n, i, j, k, info
        real(8) :: res, res0, hh, h1, h2, cs, sn
        
        !=== INITIALIZATION: ===
        this%converged = .false.
        this%iterations = 0
        
        n = size(b)
        
        if (.not. allocated(this%V)) then
            allocate(this%V(n, this%max_iter+1))
            allocate(this%H(this%max_iter+1, this%max_iter))
            allocate(this%g(this%max_iter+1))
            allocate(this%c(this%max_iter))
            allocate(this%s(this%max_iter))
            allocate(this%y(this%max_iter))
            allocate(this%r(n), this%w(n), this%w2(n))
        end if
        
        
        !=== INITIAL RESIDUAL: ===
        call this%A%apply(x, this%r)
        this%r = b - this%r
        res0 = norm2(this%r)
        this%init_residual = res0
        res = 1.d0
        
        
        this%V(:, 1) = this%r/res0
        this%g = 0.0d0
        this%g(1) = res0
        
        this%H = 0.0d0
        
        do j = 1, this%max_iter
            !=== PRECONDITIONING: w = M^{-1}*V(:,j) ===
            if (associated(this%M)) then
                call this%M%apply(this%V(:, j), this%w)
            else
                this%w = this%V(:, j)
            end if
            
            !=== MATRIX-VECOTR PRODUCT: w = A*w ===
            call this%A%apply(this%w, this%w2)
            this%w = this%w2
            
            !=== MODIFIED GRAM-SCHMIDT ORTHOGONALIZATION: ===
            do i = 1, j
				do k = 1, n
					this%H(i, j) = this%H(i, j) + this%w(k)*this%V(k, i)
				end do

				do k = 1, n
					this%w(k) = this%w(k) - this%H(i, j)*this%V(k, i)
				end do
            end do
            
            this%H(j+1, j) = norm2(this%w)
            
            if (abs(this%H(j+1, j)) > 1.0d-14) then
                this%V(:, j+1) = this%w/this%H(j+1, j)
            else
                if (this%verbosity > 1) then
                    print *, 'GMRES: Breakdown at iteration', j
                end if
                exit
            end if
            
            !=== GIVENSS ROTATION: ===
            do i = 1, j-1
                call apply_givens_rotation(this%c(i), this%s(i), &
                                          this%H(i, j), this%H(i+1, j))
            end do
            
            
            call generate_givens_rotation(this%H(j, j), this%H(j+1, j), &
                                         this%c(j), this%s(j))
            
            call apply_givens_rotation(this%c(j), this%s(j), &
                                      this%H(j, j), this%H(j+1, j))
            this%H(j+1, j) = 0.0d0
            
            
            call apply_givens_rotation(this%c(j), this%s(j), &
                                      this%g(j), this%g(j+1))
            
            !=== CONVERGENCE CONTROL: ===
            res = abs(this%g(j+1))/res0
            this%iterations = j
            
            if (res < this%rel_tol) then
               
                call backward_substitution(this%H(1:j, 1:j), &
                                          this%g(1:j), this%y(1:j), j)
                
                this%z = 0.0d0
				do i = 1, j
					this%z = this%z + this%y(i)*this%V(:, i)
				end do
			
				if (associated(this%M)) then
					!=== w = M^{-1}*z ===
					call this%M%apply(this%z, this%w)
					x = x + this%w
				else
					x = x + this%z
				end if
                
                if (this%verbosity > 0) then
					print *, 'GMRES: Converged in', this%iterations, 'iterations'
					print *, 'GMRES: Final residual drop=', res
				end if
                    
                this%converged = .true.
                this%final_residual = abs(this%g(j+1))
                this%residual_drop = res
                return
            end if
        end do
        
        if (this%H(j, j) < 1e-14) j = j - 1
        
        call backward_substitution(this%H(1:j, 1:j), this%g(1:j), this%y(1:j), j)
        
        this%z = 0.0d0
		do i = 1, j
			this%z = this%z + this%y(i)*this%V(:, i)
		end do
	
		if (associated(this%M)) then
			!=== w = M^{-1}*z ===
			call this%M%apply(this%z, this%w)
			x = x + this%w
		else
			x = x + this%z
		end if
			        
		                

        this%final_residual = abs(this%g(j+1))
        this%residual_drop = res
        
        if (.not. this%converged) then
            if (this%verbosity > 0) then
                print *, 'GMRES: Failed to converge in', this%max_iter, 'iterations'
                print *, 'GMRES: Final residual drop=', res
            end if
        end if
        
    end subroutine




!=======================================================================
!=========================== BiCGSTAB METHODS ==========================
!=======================================================================
	subroutine bicgstab_solve(this, b, x)
		implicit none
		class(bicgstab_solver_t), intent(inout) :: this
		real(8), intent(in), contiguous :: b(:)
		real(8), intent(inout), contiguous :: x(:)
		
		integer :: n, j, k
		real(8) :: rho, rho_prev, alpha, omega, beta, res0
		real(8) :: norm_s, norm_r
		real(8) :: time_start, time_end
        
        call cpu_time(time_start)
		
		!=== INITIALIZATION: ===
		n = size(b)
        if (.not. allocated(this%p)) then
            call this%allocate_workspace(n)
        end if
        
        this%converged = .false.
        this%iterations = 0
        
		
		!=== INITIAL RESIDUAL: ===
        call this%A%apply(x, this%r)
        do k = 1, n
			this%r(k) = b(k) - this%r(k)
			this%r_hat(k) = this%r(k)
        end do
        res0 = norm2(this%r)
        this%init_residual = res0

		rho_prev = 1.d0
		alpha = 1.d0; omega = 1.d0
		this%v = 0.d0; this%p = 0.d0
		
		if (res0 < this%abs_tol) then
            this%converged = .true.
            if (this%verbosity > 0) then
                print *, 'BiCGSTAB: Converged without iterations'
            end if
            call cpu_time(time_end)
            if (this%verbosity > 2) then
                print *, 'BiCGSTAB: Time =', time_end - time_start, 'seconds'
            end if
            return
        end if
		
		do j = 1, this%max_iter
			rho = DOT_PRODUCT(this%r_hat, this%r)
			beta = (rho/rho_prev)*(alpha/omega)
			do k = 1, n
				this%p(k) = this%r(k) + beta*(this%p(k) - omega*this%v(k))
			end do
			
			!=== PRECONDITIONING: y_p = M^{-1}*p ===
			if (associated(this%M)) then
				call this%M%apply(this%p, this%y_p)
			else
				do k = 1, n
					this%y_p(k) = this%p(k)
				end do
			end if
			
			!=== MATRIX-VECOTR PRODUCT: v = A*y_p ===		
			call this%A%apply(this%y_p, this%v)		
   		
			alpha = rho/DOT_PRODUCT(this%r_hat, this%v)
			do k = 1, n
				this%s(k) = this%r(k) - alpha*this%v(k)
			end do
			
			
			!=== 1ST EXIT CONDITION: ===
			norm_s = norm2(this%s)
			if (norm_s/res0 < this%rel_tol .or. norm_s < this%abs_tol) then
				do k = 1, n
					x(k) = x(k) + alpha*this%y_p(k)
				end do
				this%converged = .true.
				this%final_residual = norm_s
				this%residual_drop = norm_s/res0
				this%iterations = j
				
				if (this%verbosity > 0) then
					print *, 'BiCGSTAB: Converged in', this%iterations, 'iterations'
					print *, 'BiCGSTAB: Final residual drop =', this%residual_drop
				end if
				
				call cpu_time(time_end)
				if (this%verbosity > 2) then
					print *, 'BiCGSTAB: Time =', time_end - time_start, 'seconds'
				end if
                    
				return
			end if
				
				
				
			!=== PRECONDITIONING: y_s = M^{-1}*s ===
			if (associated(this%M)) then
				call this%M%apply(this%s, this%y_s)
			else
				do k = 1, n
					this%y_s(k) = this%s(k)
				end do
			end if
			
			!=== MATRIX-VECOTR PRODUCT: t = A*y_s ===		
			call this%A%apply(this%y_s, this%t)			
           
			
			omega = DOT_PRODUCT(this%t, this%s)/DOT_PRODUCT(this%t, this%t)
			
			do k = 1, n
				x(k) = x(k) + alpha*this%y_p(k) + omega*this%y_s(k)
				this%r(k) = this%s(k) - omega*this%t(k)
			end do
			
			!=== 2ND EXIT CONDITION: ===
			norm_r = norm2(this%r)
			if (norm_r/res0 < this%rel_tol .or. norm_r < this%abs_tol) then
				this%converged = .true.
				this%final_residual = norm_r
				this%residual_drop = norm_r/res0
				this%iterations = j
				
				if (this%verbosity > 0) then
					print *, 'BiCGSTAB: Converged in', this%iterations, 'iterations'
					print *, 'BiCGSTAB: residual drop =', this%residual_drop
				end if
				
				call cpu_time(time_end)
				if (this%verbosity > 2) then
					print *, 'BiCGSTAB: Time =', time_end - time_start, 'seconds'
				end if
				
				return
			end if
			
			!=== CONVEGENCE INFORMATION: ===
			if (this%verbosity > 1) then
				if (mod(j, 10) == 0) then
					print *, 'BiCGSTAb: ', j, 'Residual drop =', norm_r/res0
				end if
			end if
			
			this%final_residual = norm_r
            this%iterations = j
			
			rho_prev = rho
		end do		
		
		
		if (.not. this%converged) then
            if (this%verbosity > 0) then
                print *, 'BiCGSTAb: Failed to converge in', this%max_iter, 'iterations'
                print *, 'BiCGSTAb: Final residual drop=', norm_r/res0
            end if
        end if
		
		call cpu_time(time_end)
        if (this%verbosity > 2) then
            print *, 'BiCGSTAB: Total time =', time_end - time_start, 'seconds'
        end if
        
	end subroutine
	
	
	!=================================================================== 
    !=================== AUXILIARY BiCGSTAB SUBROUTINES ================
    !=================================================================== 
    subroutine bicgstab_allocate_workspace(this, n)
		implicit none
        class(bicgstab_solver_t), intent(inout) :: this
        integer, intent(in) :: n   
        
        call this%deallocate_workspace()
        
        allocate(this%p(n),&
				 this%v(n),&
				 this%t(n),&
				 this%s(n),&
				 this%r(n),&
				 this%r_hat(n),&
				 this%y_p(n),&
				 this%y_s(n))
		
		this%p = 0.d0
		this%v = 0.d0
		this%t = 0.d0
		this%s = 0.d0
		this%r = 0.d0
		this%r_hat = 0.d0
		this%y_p = 0.d0
		this%y_s = 0.d0
    end subroutine
    
    subroutine bicgstab_deallocate_workspace(this)
		implicit none
        class(bicgstab_solver_t), intent(inout) :: this
        
        if (allocated(this%p)) deallocate(this%p)
        if (allocated(this%v)) deallocate(this%v)
        if (allocated(this%t)) deallocate(this%t)
        if (allocated(this%s)) deallocate(this%s)
        if (allocated(this%r)) deallocate(this%r)
        if (allocated(this%r_hat)) deallocate(this%r_hat)
        if (allocated(this%y_p)) deallocate(this%y_p)
        if (allocated(this%y_s)) deallocate(this%y_s)
    end subroutine
    
    subroutine bicgstab_print_stats(this)
		implicit none
        class(bicgstab_solver_t), intent(in) :: this
        
        print *, "========================================="
        print *, "BiCGSTAB Solver Statistics:"
        print *, "  Converged: ", this%converged
        print *, "  Iterations: ", this%iterations
        print *, "  Initial residual: ", this%init_residual
        print *, "  Final residual: ", this%final_residual
        print *, "  Residual drop: ", this%residual_drop
        print *, "  Target relative tolerance: ", this%rel_tol
        print *, "  Target absolute tolerance: ", this%rel_tol
        print *, "  Max iterations: ", this%max_iter
        print *, "========================================="
        
    end subroutine
    



!=======================================================================
!============================= JACOBI METHODS ==========================
!=======================================================================
	subroutine jacobi_solve(this, b, x)
		implicit none
		class(jacobi_solver_t), intent(inout) :: this
		real(8), intent(in), contiguous :: b(:)
		real(8), intent(inout), contiguous :: x(:)
		
		real(8) :: res0, res, ccheck
		integer :: i, j, n
	
		!=== INITIALIZATION: ===
		n = size(b)
        if (.not. allocated(this%r)) then
            call this%allocate_workspace(n)
        end if
        
        this%converged = .false.
        this%iterations = 0
        
        !=== INITIAL RESIDUAL: ===
        call this%A%apply(x, this%r)
		this%r = b - this%r
		
		res0 = norm2(this%r)
		this%init_residual = res0
		
        
        do j = 1, this%max_iter
			call this%A%get_diagonal(this%diag)
			
			
			do i = 1, n
				this%r(i) = this%r(i)/this%diag(i)
			end do
			
			x = x + this%omega*this%r
			
			call this%A%apply(x, this%r)
			this%r = b - this%r
			this%iterations = j
			
			ccheck = norm2(this%r)/res0
			if (ccheck < this%rel_tol) then
				this%residual_drop = ccheck
				this%converged = .true.
				if (this%verbosity > 0) then
					print *, 'JACOBI: Converged in', this%iterations, 'iterations'
					print *, 'JACOBI: residual drop =', this%residual_drop
				end if
				return
			end if
			
        end do
        
	end subroutine
	
	!=================================================================== 
    !=================== AUXILIARY JACOBI SUBROUTINES ==================
    !=================================================================== 
    subroutine jacobi_allocate_workspace(this, n)
		implicit none
        class(jacobi_solver_t), intent(inout) :: this
        integer, intent(in) :: n   
        
        call this%deallocate_workspace()
        
        allocate(this%r(n), this%diag(n))
		
		this%r = 0.d0
		this%diag = 1.d0
    end subroutine
    
    subroutine jacobi_deallocate_workspace(this)
		implicit none
		class(jacobi_solver_t), intent(inout) :: this
		
		if (allocated(this%r)) deallocate(this%r)
		if (allocated(this%diag)) deallocate(this%diag)
    end subroutine
    
end module
