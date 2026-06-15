module reiman_solver_module
implicit none
!=======================================================================
!=============== APPROXIMATE REIMAN SOLVERS INTERFACE ==================
!=======================================================================
interface
	subroutine reiman_solver_interface(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, P_R, U_R, V_R, W_R, T_R, n, fluxes)
		implicit none
		integer, intent(in) :: dim                       
		real(8), intent(in) :: k, R0, cp, cv
		real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L      	
		real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
		real(8), intent(in) :: n(dim)                       
		real(8), intent(out) :: fluxes(dim+2)
	end subroutine
end interface

contains
!=======================================================================
!===================== APPROXIMATE REIMAN SOLVERS ======================
!=======================================================================
pure subroutine HLL_reiman_solver(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, P_R, U_R, V_R, W_R, T_R, n, fluxes)
    implicit none
    integer, intent(in) :: dim                       
    real(8), intent(in) :: k, R0, cp, cv
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L      	
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)                       
    real(8), intent(out) :: fluxes(dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, E_L, E_R
    real(8) :: Vn_L, Vn_R, V2_L, V2_R
    real(8) :: S_L, S_R, S_plus, S_minus, sqrtRoL, sqrtRoR, invSumSqrt, invSumS
    real(8) :: H_tilda, V_tilda(3), Vn_tilda, A_tilda
    real(8) :: F_L(5), F_R(5), W_vecL(5), W_vecR(5)
    integer :: d


	!=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L)
    A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2*(dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    
    E_L = cv*T_L + 0.5d0*V2_L
    H_L = cp*T_L + 0.5d0*V2_L

    !=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R)
    A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    
    E_R = cv*T_R + 0.5d0*V2_R
    H_R = cp*T_R + 0.5d0*V2_R

    !=== ROE AVERAGING: ===
    sqrtRoL = dsqrt(Ro_L)
    sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)

    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    A_tilda = dsqrt((k - 1.0d0)*(H_tilda - 0.5d0*(V_tilda(1)**2 + V_tilda(2)**2 + V_tilda(3)**2*(dim-2))))

    !=== WAVES VELOCITY: ===
    S_L = min(Vn_L - A_L, Vn_tilda - A_tilda)
    S_R = max(Vn_R + A_R, Vn_tilda + A_tilda)

    !=== LEFT STATE FLUXES: ===
    F_L(1) = Ro_L*Vn_L
    F_L(2) = Ro_L*U_L*Vn_L + P_L*n(1)
    F_L(3) = Ro_L*V_L*Vn_L + P_L*n(2)
    F_L(4) = Ro_L*W_L*Vn_L + P_L*n(dim)
    F_L(dim+2) = Ro_L*H_L*Vn_L

    W_vecL(1) = Ro_L
    W_vecL(2) = Ro_L*U_L
    W_vecL(3) = Ro_L*V_L
    W_vecL(4) = Ro_L*W_L
    W_vecL(dim+2) = Ro_L*E_L
	
	!=== RIGHT STATE FLUXES: ===
    F_R(1) = Ro_R*Vn_R
    F_R(2) = Ro_R*U_R*Vn_R + P_R*n(1)
    F_R(3) = Ro_R*V_R*Vn_R + P_R*n(2)
    F_R(4) = Ro_R*W_R*Vn_R + P_R*n(dim)
    F_R(dim+2) = Ro_R*H_R*Vn_R

    W_vecR(1) = Ro_R
    W_vecR(2) = Ro_R*U_R
    W_vecR(3) = Ro_R*V_R
    W_vecR(4) = Ro_R*W_R
    W_vecR(dim+2) = Ro_R*E_R

    !=== FACE FLUXES: ===
    !if (S_L >= 0.0d0) then
    !    fluxes(1:dim+2) = F_L(1:dim+2)
    !else if (S_R <= 0.0d0) then
    !    fluxes(1:dim+2) = F_R(1:dim+2)
    !else
    !    fluxes(1:dim+2) = (S_R*F_L(1:dim+2) - S_L*F_R(1:dim+2) + S_L*S_R*(W_vecR(1:dim+2) - W_vecL(1:dim+2)))/(S_R - S_L)    
    !end if
	
	S_plus  = max(0.0d0, S_R)
	S_minus = min(0.0d0, S_L)
	invSumS = 1.d0/(S_plus - S_minus)
	fluxes(1:dim+2) = (S_plus*F_L(1:dim+2) - S_minus*F_R(1:dim+2) + S_plus*S_minus*(W_vecR(1:dim+2) - W_vecL(1:dim+2)))*invSumS
	
end subroutine

pure subroutine HLLC_reiman_solver(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, P_R, U_R, V_R, W_R, T_R, n, fluxes)				!IF ISIDE :( 
    implicit none
    integer, intent(in) :: dim                       
    real(8), intent(in) :: k, R0, cp, cv
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L      	
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)                       
    real(8), intent(out) :: fluxes(dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, E_L, E_R
    real(8) :: Vn_L, Vn_R, V2_L, V2_R
    real(8) :: S_L, S_R, S_star, sqrtRoL, sqrtRoR, invSumSqrt
    real(8) :: H_tilda, V_tilda(3), Vn_tilda, A_tilda
    real(8) :: F_L(5), F_R(5), W_vecL(5), W_vecR(5), D_vec(5)
    real(8) :: num, den

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L)
    A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2*(dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    
    E_L = cv*T_L + 0.5d0*V2_L
    H_L = cp*T_L + 0.5d0*V2_L

    !=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R)
    A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    
    E_R = cv*T_R + 0.5d0*V2_R
    H_R = cp*T_R + 0.5d0*V2_R

    !=== ROE AVERAGING: ===
    sqrtRoL = dsqrt(Ro_L)
    sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)
    
    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    A_tilda = dsqrt((k - 1.0d0)*(H_tilda - 0.5d0*(V_tilda(1)**2 + V_tilda(2)**2 + V_tilda(3)**2*(dim-2))))

    !=== WAVE VELOCITY: ===
    S_L = min(Vn_L - A_L, Vn_tilda - A_tilda)
    S_R = max(Vn_R + A_R, Vn_tilda + A_tilda)

    !=== CONTACT DISCONTINUITY WAVE VELOCITY: ===
    num = P_R - P_L + Ro_L*Vn_L*(S_L - Vn_L) - Ro_R*Vn_R*(S_R - Vn_R)
    den = Ro_L*(S_L - Vn_L) - Ro_R*(S_R - Vn_R)
    S_star = num/den

    !=== LEFT STATE FLUXES: ===
    F_L(1) = Ro_L*Vn_L
    F_L(2) = Ro_L*U_L*Vn_L + P_L*n(1)
    F_L(3) = Ro_L*V_L*Vn_L + P_L*n(2)
    F_L(4) = Ro_L*W_L*Vn_L + P_L*n(dim)
    F_L(dim+2) = Ro_L*H_L*Vn_L

    W_vecL(1) = Ro_L
    W_vecL(2) = Ro_L*U_L
    W_vecL(3) = Ro_L*V_L
    W_vecL(4) = Ro_L*W_L
    W_vecL(dim+2) = Ro_L*E_L
	
	!=== RIGHT STATE FLUXES: ===
    F_R(1) = Ro_R*Vn_R
    F_R(2) = Ro_R*U_R*Vn_R + P_R*n(1)
    F_R(3) = Ro_R*V_R*Vn_R + P_R*n(2)
    F_R(4) = Ro_R*W_R*Vn_R + P_R*n(dim)
    F_R(dim+2) = Ro_R*H_R*Vn_R

    W_vecR(1) = Ro_R
    W_vecR(2) = Ro_R*U_R
    W_vecR(3) = Ro_R*V_R
    W_vecR(4) = Ro_R*W_R
    W_vecR(dim+2) = Ro_R*E_R

    !=== FACE FLUXES: ===
    if (S_L >= 0.0d0) then
        fluxes(1:dim+2) = F_L(1:dim+2)
    else if (S_star >= 0.0d0) then
        D_vec(1) = 0.d0
        D_vec(2:dim+1) = n(1:dim)
        D_vec(dim+2) = S_star 
        
        fluxes(1:dim+2) = (S_star*(S_L*W_vecL(1:dim+2) - F_L(1:dim+2)) + &
                           S_L*(P_L + Ro_L*(S_L - Vn_L)*(S_star - Vn_L))*D_vec(1:dim+2))/(S_L - S_star)
    else if (S_R >= 0.0d0) then
        D_vec(1) = 0.d0
        D_vec(2:dim+1) = n(1:dim)
        D_vec(dim+2) = S_star
        fluxes(1:dim+2) = (S_star*(S_R*W_vecR(1:dim+2) - F_R(1:dim+2)) + &
                           S_R*(P_R + Ro_R*(S_R - Vn_R)*(S_star - Vn_R))*D_vec(1:dim+2))/(S_R - S_star)
    else
        fluxes(1:dim+2) = F_R(1:dim+2)
    end if
	
end subroutine

pure subroutine ROE_reiman_solver(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, P_R, U_R, V_R, W_R, T_R, n, fluxes)
    implicit none
    integer, intent(in) :: dim                       
    real(8), intent(in) :: k, R0, cp, cv
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L      	
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)                       
    real(8), intent(out) :: fluxes(dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, V2_L, V2_R, Vn_L, Vn_R
    real(8) :: Ro_tilda, V_tilda(3), H_tilda, A_tilda, V2_tilda, Vn_tilda
    real(8) :: dP, dVn, dRo, dV(3)
    real(8) :: lambda(dim+2), alpha(dim+2)
    real(8) :: F_L(5), F_R(5), Diss(5)
    real(8) :: sqrtRoL, sqrtRoR, invSumSqrt, gam1
    integer :: d

    gam1 = k - 1.0d0

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2*(dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L

    !=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2 * (dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R

    !=== ROE AVERAGING: ===
    sqrtRoL = dsqrt(Ro_L)
    sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)

    Ro_tilda = dsqrt(Ro_L*Ro_R)
    
    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    V2_tilda = V_tilda(1)**2 + V_tilda(2)**2 + (V_tilda(3)**2*(dim-2))
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    A_tilda = dsqrt(gam1*(H_tilda - 0.5d0*V2_tilda))

    !=== VARIABLES DIFFERENCES: ===
    dP = P_R - P_L
    dRo = Ro_R - Ro_L
    dV(1) = U_R - U_L
    dV(2) = V_R - V_L
    dV(3) = W_R - W_L
    dVn = dV(1)*n(1) + dV(2)*n(2) + (dV(3)*n(dim)*(dim-2))

    !=== ABS EIGEN VALUES: ===
    lambda(1) = abs(Vn_tilda - A_tilda)
    lambda(2:dim+1) = abs(Vn_tilda)
    lambda(dim+2) = abs(Vn_tilda + A_tilda)
    
    
    !=== WAVE AMPLITUDES: ===
    alpha(1) = (dP - Ro_tilda*A_tilda*dVn)/(2.0d0*A_tilda**2)
    alpha(dim+2) = (dP + Ro_tilda*A_tilda*dVn)/(2.0d0*A_tilda**2)
    alpha(2) = dRo - dP/(A_tilda**2)


    Diss = 0.0d0
    
    !=== WAVE 1: V_n - A: ===
    Diss(1) = lambda(1)*alpha(1)
    Diss(2:dim+1) = lambda(1)*alpha(1)*(V_tilda(1:dim) - A_tilda*n(1:dim))
    Diss(dim+2) = lambda(1)*alpha(1)*(H_tilda - A_tilda*Vn_tilda)

    !=== WAVE 2: V_n: ===
    Diss(1) = Diss(1) + lambda(2)*alpha(2)
    Diss(2:dim+1) = Diss(2:dim+1) + lambda(2)*alpha(2)*V_tilda(1:dim)
    Diss(dim+2) = Diss(dim+2) + lambda(2)*alpha(2)*0.5d0*V2_tilda

    !=== WAVE 3, 4: V_n: ===
    do d = 1, dim
        Diss(1+d) = Diss(1+d) + lambda(1+d)*Ro_tilda*(dV(d) - dVn*n(d))
        Diss(dim+2) = Diss(dim+2) + lambda(1+d)*Ro_tilda*(V_tilda(d)*dV(d))
    end do
    Diss(dim+2) = Diss(dim+2) - lambda(2)*Ro_tilda*(Vn_tilda*dVn)
    
    !=== WAVE 5: V_n + A: ===
    Diss(1) = Diss(1) + lambda(dim+2)*alpha(dim+2)
    Diss(2:dim+1) = Diss(2:dim+1) + lambda(dim+2)*alpha(dim+2)*(V_tilda(1:dim) + A_tilda*n(1:dim))
    Diss(dim+2) = Diss(dim+2) + lambda(dim+2)*alpha(dim+2)*(H_tilda + A_tilda*Vn_tilda)

    !=== LEFT/RIGHT FLUXES: ===
    F_L(1) = Ro_L*Vn_L
    F_L(2) = Ro_L*U_L*Vn_L + P_L*n(1)
    F_L(3) = Ro_L*V_L*Vn_L + P_L*n(2)
    F_L(4) = Ro_L*W_L*Vn_L + P_L*n(dim)
    F_L(dim+2) = Ro_L*H_L*Vn_L

    F_R(1) = Ro_R*Vn_R
    F_R(2) = Ro_R*U_R*Vn_R + P_R*n(1)
    F_R(3) = Ro_R*V_R*Vn_R + P_R*n(2)
    F_R(4) = Ro_R*W_R*Vn_R + P_R*n(dim)
    F_R(dim+2) = Ro_R*H_R*Vn_R

    !=== FACES FLUXES: ===
    fluxes(1:dim+2) = 0.5d0*(F_L(1:dim+2) + F_R(1:dim+2) - Diss(1:dim+2))

end subroutine

pure subroutine AUSM_reiman_solver(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, P_R, U_R, V_R, W_R, T_R, n, fluxes)
    implicit none
    integer, intent(in) :: dim                       
    real(8), intent(in) :: k, R0, cp, cv
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L      	
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)                       
    real(8), intent(out) :: fluxes(dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, V2_L, V2_R
    real(8) :: M_L, M_R, M_f, P_f, A_f
    real(8) :: Fi_L(5), Fi_R(5), F_p(5)
    real(8) :: Vn_L, Vn_R
    real(8) :: wL, wR

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L)
    A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2*(dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    
    H_L = cp*T_L + 0.5d0*V2_L

    !=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R)
    A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    
    H_R = cp*T_R + 0.5d0*V2_R

    !=== LEFT/RIGHT MACH NUMBERS: ===
    A_f = 0.5d0*(A_L + A_R)
    M_L = Vn_L/A_f
    M_R = Vn_R/A_f

    !=== CONVECTION-PRESSURE SPLITTING: ===
    M_f = M_plus(M_L) + M_minus(M_R)
    P_f = P_plus(M_L)*P_L + P_minus(M_R)*P_R

    !=== LEFT/RIGHT STATE CONVECTION VECTORS: ===
    Fi_L(1) = 1.0d0
    Fi_L(2) = U_L
    Fi_L(3) = V_L
    Fi_L(4) = W_L
    Fi_L(dim+2) = H_L

    Fi_R(1) = 1.0d0
    Fi_R(2) = U_R
    Fi_R(3) = V_R
    Fi_R(4) = W_R
    Fi_R(dim+2) = H_R

    !=== PRESSURE PART: ===
    F_p = 0.0d0
    F_p(2:dim+1) = P_f*n(1:dim)

    !=== FACE FLUXES: ===
    !if (M_f >= 0.0d0) then
    !    fluxes(1:dim+2) = A_f*M_f*Ro_L*Fi_L(1:dim+2) + F_p(1:dim+2)
    !else
    !    fluxes(1:dim+2) = A_f*M_f*Ro_R*Fi_R(1:dim+2) + F_p(1:dim+2)
    !end if
	
	wL = 0.5d0*(1.0d0 + sign(1.0d0, M_f))
	wR = 1.0d0 - wL
	fluxes(1:dim+2) = A_f*M_f*(wL*Ro_L*Fi_L(1:dim+2) + wR*Ro_R*Fi_R(1:dim+2)) + F_p(1:dim+2)
contains
    pure real(8) function M_plus(M)
        real(8), intent(in) :: M
        if (abs(M) >= 1.0d0) then
            M_plus = 0.5d0*(M + abs(M))
        else
            M_plus = 0.25d0*(M + 1.0d0)**2
        end if
    end function

    pure real(8) function M_minus(M)
        real(8), intent(in) :: M
        if (abs(M) >= 1.0d0) then
            M_minus = 0.5d0*(M - abs(M))
        else
            M_minus = -0.25d0*(M - 1.0d0)**2
        end if
    end function

    pure real(8) function P_plus(M)
        real(8), intent(in) :: M
        if (abs(M) >= 1.0d0) then
            P_plus = 0.5d0*(1.0d0 + sign(1.0d0, M))
        else
            P_plus = 0.25d0*(M + 1.0d0)**2 * (2.0d0 - M)
        end if
    end function

    pure real(8) function P_minus(M)
        real(8), intent(in) :: M
        if (abs(M) >= 1.0d0) then
            P_minus = 0.5d0*(1.0d0 - sign(1.0d0, M))
        else
            P_minus = 0.25d0*(M - 1.0d0)**2 * (2.0d0 + M)
        end if
    end function

end subroutine

pure subroutine AUSMplus_reiman_solver(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, P_R, U_R, V_R, W_R, T_R, n, fluxes)
    implicit none
    integer, intent(in) :: dim                       
    real(8), intent(in) :: k, R0, cp, cv
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L      	
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)                       
    real(8), intent(out) :: fluxes(dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, V2_L, V2_R
    real(8) :: M_L, M_R, M_f, P_f, A_f
    real(8) :: Fi_L(5), Fi_R(5), F_p(5)
    real(8) :: Vn_L, Vn_R
    real(8) :: wL, wR

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L)
    A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2*(dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    
    H_L = cp*T_L + 0.5d0*V2_L

    !=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R)
    A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    
    H_R = cp*T_R + 0.5d0*V2_R

    !=== LEFT/RIGHT MACH NUMBERS: ===
    A_f = dsqrt(A_L*A_R) 												!A_f = min(A_L**2/max(A_L, abs(Vn_L)), A_R**2/max(A_R, abs(Vn_R))) 
    M_L = Vn_L/A_f
    M_R = Vn_R/A_f

    !=== CONVECTION-PRESSURE SPLITTING: ===
    M_f = M_plus_plus(M_L) + M_minus_minus(M_R)
    P_f = P_plus_plus(M_L)*P_L + P_minus_minus(M_R)*P_R

    !=== LEFT/RIGHT STATE CONVECTION VECTORS: ===
    Fi_L(1) = 1.0d0
    Fi_L(2) = U_L
    Fi_L(3) = V_L
    Fi_L(4) = W_L
    Fi_L(dim+2) = H_L

    Fi_R(1) = 1.0d0
    Fi_R(2) = U_R
    Fi_R(3) = V_R
    Fi_R(4) = W_R
    Fi_R(dim+2) = H_R

    !=== PRESSURE PATR: ===
    F_p = 0.0d0
    F_p(2:dim+1) = P_f*n(1:dim)
	
	!=== FACE FLUXES: ===
    !if (M_f >= 0.0d0) then
    !    fluxes(1:dim+2) = A_f*M_f*Ro_L*Fi_L(1:dim+2) + F_p(1:dim+2)
    !else
    !    fluxes(1:dim+2) = A_f*M_f*Ro_R*Fi_R(1:dim+2) + F_p(1:dim+2)
    !end if
    
    wL = 0.5d0*(1.0d0 + sign(1.0d0, M_f))
	wR = 1.0d0 - wL
	fluxes(1:dim+2) = A_f*M_f*(wL*Ro_L*Fi_L(1:dim+2) + wR*Ro_R*Fi_R(1:dim+2)) + F_p(1:dim+2)

contains
    pure real(8) function M_plus_plus(M)
        real(8), intent(in) :: M
        if (M >= 1.0d0) then
            M_plus_plus = M
        else if (M <= -1.0d0) then
            M_plus_plus = 0.0d0
        else
            M_plus_plus = 0.25d0*(M + 1.0d0)**2 + 0.125d0*(M**2 - 1.0d0)**2
        end if
    end function

    pure real(8) function M_minus_minus(M)
        real(8), intent(in) :: M
        if (M >= 1.0d0) then
            M_minus_minus = 0.0d0
        else if (M <= -1.0d0) then
            M_minus_minus = M
        else
            M_minus_minus = -0.25d0*(M - 1.0d0)**2 - 0.125d0*(M**2 - 1.0d0)**2
        end if
    end function

    pure real(8) function P_plus_plus(M)
        real(8), intent(in) :: M
        if (M >= 1.0d0) then
            P_plus_plus = 1.0d0
        else if (M <= -1.0d0) then
            P_plus_plus = 0.0d0
        else
            P_plus_plus = 0.25d0*(M + 1.0d0)**2 * (2.0d0 - M) + 0.1875d0*M*(M**2 - 1.0d0)**2
        end if
    end function

    pure real(8) function P_minus_minus(M)
        real(8), intent(in) :: M
        if (M >= 1.0d0) then
            P_minus_minus = 0.0d0
        else if (M <= -1.0d0) then
            P_minus_minus = 1.0d0
        else
            P_minus_minus = 0.25d0*(M - 1.0d0)**2 * (2.0d0 + M) - 0.1875d0*M*(M**2 - 1.0d0)**2
        end if
    end function

end subroutine

end module

