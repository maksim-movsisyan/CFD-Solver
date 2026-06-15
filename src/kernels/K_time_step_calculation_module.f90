module time_step_calculation_module
use physical_properties_module
implicit none
contains
!=======================================================================
!============== PSEUDO TIME STEP CALCULATION SUBROUTINES ===============
!=======================================================================
	subroutine compute_ptime_step(dim, CFL, P, V, T, k, R_gas, dx, dtau_inv, dtau_visc)
		implicit none
		integer, intent(in) :: dim
		real(8), intent(in) :: CFL
		real(8), intent(in) :: P, V(dim), T
		real(8), intent(in) :: k, R_gas
		real(8), intent(in) :: dx
		real(8), intent(inout) :: dtau_inv, dtau_visc
		
		real(8) :: mu, Ro, A
		
		Ro = P/(R_gas*T)
		A = dsqrt(k*R_gas*T)
		mu = Mu_air(T)
		
		dtau_inv = CFL*dx/(A + Norm2(V))
		dtau_visc = CFL*Ro*dx**2/(2.d0*mu)
	end subroutine
	
end module
