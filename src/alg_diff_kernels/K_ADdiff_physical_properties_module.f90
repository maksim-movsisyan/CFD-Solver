module ADdiff_physical_properties_module
implicit none
contains
!  Differentiation of mu_air in forward (tangent) mode:
!   variations   of useful results: mu_air
!   with respect to varying inputs: t
!   RW status of diff variables: t:in mu_air:out
  REAL*8 FUNCTION MU_AIR_D(t, td, mu_air)
    IMPLICIT NONE
    REAL*8, PARAMETER :: mu_ref=0.00001716
    REAL*8, PARAMETER :: t_ref=273.15d0
    REAL*8, INTENT(IN) :: t
    REAL*8, INTENT(IN) :: td
    REAL*8 :: temp
    REAL*8 :: mu_air
    temp = (t/t_ref)**1.5d0/(t+111.d0)
    mu_air_d = mu_ref*(t_ref+111.d0)*(1.5d0*(t/t_ref)**0.5D0/t_ref-temp)&
&     *td/(t+111.d0)
    mu_air = mu_ref*(t_ref+111.d0)*temp
  END FUNCTION
  
!  Differentiation of lambda_air in forward (tangent) mode:
!   variations   of useful results: lambda_air
!   with respect to varying inputs: t
!   RW status of diff variables: t:in lambda_air:out
  REAL*8 FUNCTION LAMBDA_AIR_D(t, td, pr, cp, lambda_air)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: t, pr, cp
    REAL*8, INTENT(IN) :: td
    REAL*8 :: result1
    REAL*8 :: result1d
    REAL*8 :: lambda_air
    result1d = MU_AIR_D(t, td, result1)
    lambda_air_d = cp*result1d/pr
    lambda_air = result1*cp/pr
  END FUNCTION  

end module
