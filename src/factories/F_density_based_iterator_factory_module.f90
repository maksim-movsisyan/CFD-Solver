module density_based_iterator_factory_module
use density_based_iterators_module
implicit none

integer, parameter, public :: SCHEME_EXPLICIT = 0
integer, parameter, public :: SCHEME_IMPLICIT = 1
integer, parameter, public :: SCHEME_JACOBIANFREE = 2
integer, parameter, public :: SCHEME_ALGHDIFF = 3
integer, parameter, public :: SCHEME_BLUSGS = 4
!=======================================================================
!============ DENSITY-BASED ITERATOR FACTORY DATA STRUCTURE ============
!=======================================================================
type :: db_iterator_factory_t
contains
	procedure, nopass :: create => dbi_create
end type

type(db_iterator_factory_t) :: db_iterator_factory

contains

    function dbi_create(SCHEME) result(dbi_ptr)
		implicit none
        integer, intent(in) :: SCHEME
        class(density_based_iterator_t), pointer :: dbi_ptr

        select case (SCHEME)
        case (SCHEME_EXPLICIT)
            allocate(db_explicit_iterator_t :: dbi_ptr)
        case(SCHEME_IMPLICIT)
			allocate(db_implicit_iterator_t :: dbi_ptr)
		case(SCHEME_JACOBIANFREE)
			allocate(db_jacfree_iterator_t :: dbi_ptr)
		case(SCHEME_ALGHDIFF)
			allocate(db_algdiff_iterator_t :: dbi_ptr)
		case(SCHEME_BLUSGS)
			allocate(db_blusgs_iterator_t :: dbi_ptr)
        case default
             write(*,*) 'ERROR: Unknown iterator ID: ', SCHEME
			 error stop 'Density-based iterator factory error'
        end select
    end function
end module
