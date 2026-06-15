module physical_properties_module
implicit none
contains
pure real(8) function MU_air(T)
    implicit none
    real(8), parameter :: MU_ref = 0.00001716
    real(8), parameter :: T_ref = 273.15d0
    real(8), intent(in) :: T

    MU_air = MU_ref*(T/T_ref)**(1.5d0)*(T_ref + 111.d0)/(T + 111.d0)
end function

pure real(8) function lambda_air(T, Pr, cp)
    implicit none
    real(8), intent(in) :: T, Pr, cp

    lambda_air = MU_air(T)*cp/Pr
end function

end module
