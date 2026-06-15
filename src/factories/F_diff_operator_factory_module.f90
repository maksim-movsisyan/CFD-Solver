module diff_operator_factory_module
use diff_operator_module
use mesh_module, only: mesh_t
implicit none

!=======================================================================
!=============== DIFF OPERATOR FACTORY DATA STRUCTURE ==================
!=======================================================================
type :: gradient_factory_t
contains
	procedure :: create => factory_create_gradient
end type

type(gradient_factory_t) :: gradient_factory

contains
!=======================================================================
!====================== DIFF OPERATOR FACTORY METHODS ==================
!=======================================================================
	function create_gradient(gradient_type) result(gradient_ptr)
		implicit none
        character(len=*), intent(in) :: gradient_type
        class(gradient_t), pointer :: gradient_ptr
        
        nullify(gradient_ptr)
        
        select case (trim(adjustl(gradient_type)))
        case ('Green-Gauss', 'GG', 'GREEN-GAUSS', 'green-gauss',&
			  'Green_Gauss', 'G_G', 'GREEN_GAUSS', 'green_gauss')
            allocate(GG_gradient_t::gradient_ptr)
			select type (gradient_ptr)
			type is (GG_gradient_t)
				gradient_ptr%name = trim(adjustl(gradient_type))
			end select
            
        case ('Least-Squares', 'LS', 'LSQ', 'LEAST-SQUARES',&
			  'Least_Squares', 'L_S', 'L_SQ', 'LEAST_SQUARES')
            allocate(LSQ_gradient_t::gradient_ptr)
			select type (gradient_ptr)
			type is (LSQ_gradient_t)
				gradient_ptr%name = trim(adjustl(gradient_type))
			end select

            
        case default
            write(*,*) 'ERROR: Unknown gradient type: ', trim(adjustl(gradient_type))
            error stop 'Gradient operator factory error'
        end select
            
    end function
	
	function factory_create_gradient(this, mesh, gradient_type) result(gradient_ptr)
		implicit none
        class(gradient_factory_t), intent(in) :: this
        type(mesh_t), target, intent(in) :: mesh
        character(len=*), optional, intent(in) :: gradient_type
        
        class(gradient_t), pointer :: gradient_ptr
		character(len=128) :: grad_type
       
		grad_type = 'GG'
		if (present(gradient_type)) grad_type = gradient_type
       
        gradient_ptr => create_gradient(grad_type)
        
        call gradient_ptr%set_mesh(mesh)
    end function
    

end module
