module stability_operator_module
implicit none
interface
	subroutine stability_operator_interface(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, &
										P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(out) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)
	end subroutine
end interface

contains
!=======================================================================
!========================== JACOBIAN MATRICES ==========================
!=======================================================================
pure subroutine get_Aq_matrix(dim, R0, k, Vel, n, Ro, T, Aq)
	!===================================================================
	!==== JACOBIAN MATRIX Aq = dQ/dU*dF/dU*dU/dQ = dQ/dU*dF/dQ		====
	!==== Q - PRIMITIVE VARIABES									====
	!==== U - CONSERVATIVE VARIABLES								====
	!==== F - CONSERVATIVER FLUXES 									====
	!===================================================================
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: R0, k, Ro, T
    real(8), intent(in) :: Vel(3), n(dim)
    real(8), intent(inout) :: Aq(dim+2, dim+2)
    
    real(8) :: Vn, kRT
    integer :: d

    Vn = dot_product(Vel(1:dim), n)
    kRT = k*R0*T
    Aq = 0.0d0

    !=== CONTINUITIE EQUATION: ===
    Aq(1, 1) = Vn
    do d = 1, dim
        Aq(1, 1+d) = Ro*kRT*n(d) 
    end do

    !=== MOMENTUM EQUATION: ===
    do d = 1, dim
        Aq(1+d, 1) = n(d)/Ro
        Aq(1+d, 1+d) = Vn
    end do

    !=== ENERGY EQUATION: ===
    do d = 1, dim
        Aq(dim+2, 1+d) = (k - 1.0d0)*T*n(d)
    end do
    Aq(dim+2, dim+2) = Vn
end subroutine

pure subroutine get_Aq_visc_matrix(dim, cp, mu, lambda, RL, Vel, Ro, T, A, Aq)
	!===================================================================
	!==== JACOBIAN MATRIX APPROXATION FOR VISCOUS PART Aq=dQ/dU*dF/dQ ==
	!==== Q - PRIMITIVE VARIABES									====
	!==== U - CONSERVATIVE VARIABLES								====
	!==== F - CONSERVATIVER FLUXES 									====
	!===================================================================
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: cp, mu, lambda, RL
    real(8), intent(in) :: Vel(3), Ro, T, A
    real(8), intent(inout) :: Aq(dim+2, dim+2)
    
    real(8) :: teta, temp
    integer :: d
    
    teta = 1.d0/A**2
	temp = 1.d0/(RL*(T*cp*teta - 1.d0))
	
    Aq = 0.d0
	do d = 1, dim
		Aq(1, d+1) = -mu*Vel(d)*temp
		Aq(d+1, d+1) = mu/(RL*Ro)
		Aq(dim+2, d+1) = -mu*T*Vel(d)*teta*temp/Ro
	end do
	Aq(1, dim+2) = lambda*temp
	Aq(dim+2, dim+2) = lambda*T*teta*temp/Ro
	
end subroutine

pure subroutine get_Au_matrix(dim, k, Vel, n, Ro, E, P, Au)
	!===================================================================
	!==== JACOBIAN MATRIX Au = dF/dU								====
	!==== U - CONSERVATIVE VARIABLES								====
	!==== F - CONSERVATIVER FLUXES 									====
	!===================================================================
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, Ro, E, P
    real(8), intent(in) :: Vel(3), n(dim)
    real(8), intent(out) :: Au(dim+2, dim+2)

    real(8) :: Vn, V2, phi2, H, G1, G2
    integer :: i, j
	
	Vn = dot_product(Vel(1:dim), n)
    V2 = dot_product(Vel(1:dim), Vel(1:dim))
    H = E + P/Ro
    G1 = k - 1.0d0
    phi2 = 0.5d0*G1*V2

    Au = 0.0d0

    !=== CONTINUITIE EQUATION: ===
    do i = 1, dim
		Au(1, 1+i) = n(i)
    end do

    !=== MOMENTUM EQUATION: ===
    do i = 1, dim
        G2 = n(i)*G1
    
        Au(1+i, 1) = n(i)*phi2 - Vel(i)*Vn

        do j = 1, dim
            Au(1+i, 1+j) = Vel(i)*n(j) - G2*Vel(j)
        end do
        Au(1+i, 1+i) = Au(1+i, 1+i) + Vn
       
        Au(1+i, dim+2) = G2
    end do

    !=== ENERGY EQUATION: ===
    Au(dim+2, 1) = Vn*(phi2 - H)
    do j = 1, dim
        Au(dim+2, 1+j) = n(j)*H - G1*Vn*Vel(j)
    end do
    Au(dim+2, dim+2) = k*Vn

end subroutine

pure subroutine get_Au_visc_matrix(dim, cp, mu, lambda, RL, Vel, Ro, T, A, H, Au)
	!===================================================================
	!==== JACOBIAN MATRIX APPROXATION FOR VISCOUS PART Au=dF/dU		====
	!==== U - CONSERVATIVE VARIABLES								====
	!==== F - CONSERVATIVER FLUXES 									====
	!===================================================================
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: cp, mu, lambda, RL
    real(8), intent(in) :: Vel(3), Ro, T, A, H
    real(8), intent(inout) :: Au(dim+2, dim+2)
	
	real(8) :: teta, temp, C_mu, C_lam, V2
	integer :: d
    
    V2 = 0.d0
    do d = 1, dim
		V2 = V2 + Vel(d)**2
    end do
    
    C_mu  = mu/RL
    C_lam = lambda/RL
    teta = 1.d0/A**2
    temp = 1.d0/(Ro*(T*cp*teta - 1.d0))
	
	Au = 0.d0
	do d = 1, dim
		Au(d+1, 1) = -C_mu*Vel(d)/Ro
		Au(d+1, d+1) =  C_mu/Ro
        
        Au(dim+2, d+1) = -C_lam*(T*teta*Vel(d))*temp
	end do
	Au(dim+2, 1) = C_lam*T*(teta*(V2 - H) + 1.d0)*temp
    
	Au(dim+2, dim+2) = C_lam*(T*teta)*temp	

end subroutine

pure subroutine get_J_mtrx(dim, cp, k, A, Vel, H, Ro, T, S)
	!===================================================================
	!==== JACOBIAN MATRIX dQ/dU										====
	!==== Q - PRIMITIVE VARIABES 									====
	!==== U - CONSERVATIVE VARIABLES								====
	!===================================================================
	implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: cp, k, A, Ro, T, H
    real(8), intent(in) :: Vel(3)
    real(8), intent(inout) :: S(dim+2, dim+2)
    
    integer :: d
    real(8) :: den_base, den, V2, teta
    
	S = 0.d0
	
	den_base = k*T*cp - A**2.d0
	V2 = Vel(1)**2 + Vel(2)**2 + Vel(3)**2*(dim-2)
	
	!=== CONTINUITIE EQUATION: ===
	S(1, 1) = 0.5d0*V2*A**2/den_base
	do d = 1, dim
		S(1, 1+d) = -A**2*Vel(d)/den_base
	end do
	S(1, dim+2) = A**2/den_base
	
	!=== MOMENTUM EQUATION: ===
	do d = 1, dim
		S(1+d, 1+d) = 1.d0/Ro
		S(1+d, 1) = -Vel(d)/Ro
	end do
	
	!=== ENERGY EQUATION: ===
	S(dim+2, 1) = T*(A**2.d0 + k*(V2 - H))/(Ro*den_base)
	
	den = -k*T/(Ro*den_base)
	do d = 1, dim
		S(dim+2, 1+d) = den*Vel(d)
	end do
	S(dim+2, dim+2) = k*T/(Ro*den_base)
end subroutine


!=======================================================================
!=============== STABILITY OPERATORS PRIMITIVE VARIABLES ===============
!=======================================================================
pure subroutine HLL_stability_universal(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L,&
										P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, Vn_L, Vn_R, V2_L, V2_R
    real(8) :: S_L, S_R, sqrtRoL, sqrtRoR, invSumSqrt
    real(8) :: H_tilda, V_tilda(3), Vn_tilda, A_tilda, V2_tilda, Ro_tilda, T_tilda
    real(8) :: Ro_avg, T_avg, V_avg(3)
    real(8) :: Aq(dim+2, dim+2), coeff_A, coeff_E
    integer :: i

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L); A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2 * (dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L

	!=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R); A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R

    !=== ROE AVERAGING: ===
    Ro_tilda = dsqrt(Ro_L*Ro_R)
    sqrtRoL = dsqrt(Ro_L); sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)
    
    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    V2_tilda = V_tilda(1)**2 + V_tilda(2)**2 + (V_tilda(3)**2*(dim-2))
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    A_tilda = dsqrt((k-1.d0)*(H_tilda - 0.5d0*V2_tilda))
	T_tilda = A_tilda**2/(k*R0)
	
    S_L = min(Vn_L - A_L, Vn_tilda - A_tilda)
    S_R = max(Vn_R + A_R, Vn_tilda + A_tilda)

    !=== AVG MATRIX ON FACE: ===
    !V_avg = (1.d0 - g)*[U_L, V_L, W_L] + g*[U_R, V_R, W_R]
    !Ro_avg = (1.d0 - g)*Ro_L + g*Ro_R
    !T_avg = (1.d0 - g)*T_L + g*T_R
    call get_Aq_matrix(dim, R0, k, V_tilda, n, Ro_tilda, T_tilda, Aq)

    !=== STABILITY OPERATORS: ===
    Mtrx_L = 0.0d0
    Mtrx_R = 0.0d0

    if (S_L >= 0.0d0) then
        Mtrx_L = Aq
        
    else if (S_R <= 0.0d0) then
        Mtrx_R = Aq
        
    else
        coeff_A = 0.5d0*(1.0d0 + (S_R + S_L)/(S_R - S_L))
        coeff_E = S_R*S_L/(S_R - S_L)
        
        Mtrx_L = coeff_A*Aq
        Mtrx_R = (1.0d0 - coeff_A)*Aq
        
        do i = 1, dim+2
            Mtrx_L(i, i) = Mtrx_L(i, i) - coeff_E
            Mtrx_R(i, i) = Mtrx_R(i, i) + coeff_E
        end do
    end if
end subroutine

pure subroutine HLLC_stability_universal(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, &
										 P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, Vn_L, Vn_R, V2_L, V2_R
    real(8) :: S_L, S_R, S_star, sqrtRoL, sqrtRoR, invSumSqrt
    real(8) :: H_tilda, V_tilda(3), Vn_tilda, A_tilda, V2_tilda, Ro_tilda, T_tilda
    real(8) :: Ro_avg, T_avg, V_avg(3)
    real(8) :: Aq(dim+2, dim+2), coeff_A, coeff_E
    integer :: i, j

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L); A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2 * (dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L
	
	!=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R); A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R

    !=== ROE AVERAGING: ===
    Ro_tilda = dsqrt(Ro_L*Ro_R)
    sqrtRoL = dsqrt(Ro_L); sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)
    
    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    V2_tilda = V_tilda(1)**2 + V_tilda(2)**2 + (V_tilda(3)**2*(dim-2))
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    A_tilda = dsqrt((k-1.d0)*(H_tilda - 0.5d0*V2_tilda))
	T_tilda = A_tilda**2/(k*R0)
	
    S_L = min(Vn_L - A_L, Vn_tilda - A_tilda)
    S_R = max(Vn_R + A_R, Vn_tilda + A_tilda)

    S_star = ((P_R - P_L) + Ro_L*Vn_L*(S_L - Vn_L) - Ro_R*Vn_R*(S_R - Vn_R))/&
             (Ro_L*(S_L - Vn_L) - Ro_R*(S_R - Vn_R))

    !=== AVG MATRIX ON FACE: ===
    !V_avg = (1.d0 - g)*[U_L, V_L, W_L] + g*[U_R, V_R, W_R]
    !Ro_avg = (1.d0 - g)*Ro_L + g*Ro_R
    !T_avg = (1.d0 - g)*T_L + g*T_R
    call get_Aq_matrix(dim, R0, k, V_tilda, n, Ro_tilda, T_tilda, Aq)

    !=== STABILITY OPERATORS: ===
    Mtrx_L = 0.0d0
    Mtrx_R = 0.0d0

    if (S_L >= 0.0d0) then
        Mtrx_L = Aq
        
    else if (S_R <= 0.0d0) then
        Mtrx_R = Aq
        
    else if (S_star >= 0.0d0) then
        coeff_A = -S_star/(S_L - S_star)
        coeff_E = (S_star*S_L)/(S_L - S_star)
        Mtrx_L = coeff_A*Aq
        do i = 1, dim+2
            Mtrx_L(i, i) = Mtrx_L(i, i) + coeff_E
        end do
        
    else
        coeff_A = -S_star/(S_R - S_star)
        coeff_E = (S_star*S_R)/(S_R - S_star)
        Mtrx_R = coeff_A*Aq
        do i = 1, dim+2
            Mtrx_R(i, i) = Mtrx_R(i, i) + coeff_E
        end do
        
    end if

end subroutine

pure subroutine ROE_stability_universal(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, &
										P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)

    real(8) :: Ro_L, Ro_R, H_L, H_R, Vn_L, Vn_R, V2_L, V2_R
    real(8) :: Ro_tilda, V_tilda(3), H_tilda, A_tilda, V2_tilda, Vn_tilda, T_tilda
    real(8) :: Aq(dim+2, dim+2), T2_plus(dim+2, dim+2), T2_minus(dim+2, dim+2)
    real(8) :: r1(5), r5(5), l1(5), l5(5), lambda1, lambda5, inv_den
    real(8) :: sqrtRoL, sqrtRoR, invSumSqrt
    integer :: i, j

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2 * (dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L
	
	!=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R

    !=== ROE AVERAGING: ===
    Ro_tilda = dsqrt(Ro_L*Ro_R)
    sqrtRoL = dsqrt(Ro_L); sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)
    
    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    V2_tilda = V_tilda(1)**2 + V_tilda(2)**2 + (V_tilda(3)**2*(dim-2))
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    A_tilda = dsqrt((k-1.d0)*(H_tilda - 0.5d0*V2_tilda))
    T_tilda = A_tilda**2/(k*R0)

    !=== AVG MATRIX ON FACE: ===
	!call get_Aq_matrix(dim, R0, k, (1.d0-g)*[U_L, V_L, W_L] + g*[U_R, V_R, W_R], &
    !                   n, (1.d0-g)*Ro_L + g*Ro_R, (1.d0-g)*T_L + g*T_R, Aq)
    call get_Aq_matrix(dim, R0, k, V_tilda, n, Ro_tilda, T_tilda, Aq)

    !=== EIGEN VALUES: ===
    lambda1 = Vn_tilda + A_tilda
    lambda5 = Vn_tilda - A_tilda
    inv_den = 1.0d0/(Ro_tilda*A_tilda**2*(lambda5 - lambda1))

    !=== 1ST RIGHT EIGN VECTOR: ===
    r1 = 0.d0
    r1(1) = Ro_tilda*A_tilda**2
    r1(2:dim+1) = n(1:dim)*(Vn_tilda - lambda5)
    r1(dim+2) = A_tilda**2*(k - 1.d0)/(R0*k)
    r1 = -inv_den*r1

    !=== 5TH RIGHT EIGN VECTOR: ===
    r5 = 0.d0
    r5(1) = Ro_tilda*A_tilda**2
    r5(2:dim+1) = n(1:dim)*(Vn_tilda - lambda1)
    r5(dim+2) = A_tilda**2*(k - 1.d0)/(R0 * k)
    r5 = inv_den*r5

    !=== 1ST LEFT EIGN VECTOR: ===
    l1 = 0.d0
    l1(1) = lambda1 - Vn_tilda
    l1(2:dim+1) = Ro_tilda*A_tilda**2*n(1:dim)
	
	!=== 5TH LEFT EIGN VECTOR: ===
    l5 = 0.d0
    l5(1) = lambda5 - Vn_tilda
    l5(2:dim+1) = Ro_tilda*A_tilda**2*n(1:dim)

    !=== STABILITY OPERATORS: ===
    T2_plus = 0.d0
    T2_minus = 0.d0

    if (Vn_tilda >= 0.d0) then
        do j = 1, dim+2
            do i = 1, dim+2
                T2_plus(i, j) = r1(i)*(lambda1 - Vn_tilda)*l1(j) - &
                                r5(i)*min(Vn_tilda, Vn_tilda - lambda5)*l5(j)
                T2_minus(i, j) = r5(i)*min(0.d0, lambda5)*l5(j)
            end do
            T2_plus(j, j) = T2_plus(j, j) + Vn_tilda
        end do
    else
        do j = 1, dim+2
            do i = 1, dim+2
                T2_plus(i, j) = r1(i)*max(0.d0, lambda1)*l1(j)
                T2_minus(i, j) = r5(i)*(lambda5 - Vn_tilda)*l5(j) - &
                                 r1(i)*max(Vn_tilda, Vn_tilda - lambda1)*l1(j)
            end do
            T2_minus(j, j) = T2_minus(j, j) + Vn_tilda
        end do
    end if

    Mtrx_L = 0.5d0*Aq + 0.5d0*(T2_plus - T2_minus)
    Mtrx_R = 0.5d0*Aq - 0.5d0*(T2_plus - T2_minus)

end subroutine

pure subroutine ROE2_stability_universal(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, &
										 P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, Vn_L, Vn_R, V2_L, V2_R
    real(8) :: Ro_tilda, V_tilda(3), H_tilda, A_tilda, V2_tilda, Vn_tilda, T_tilda
    real(8) :: Aq(dim+2, dim+2), AbsA(dim+2, dim+2)
    real(8) :: r(5,5), l(5,5), lam(5), alpha(5)
    real(8) :: sqrtRoL, sqrtRoR, invSumSqrt, gam1, eps
    integer :: i, j, m

    gam1 = k - 1.0d0
    eps = 1.0d-14

    Ro_L = P_L/(R0*T_L); A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2 * (dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L

    Ro_R = P_R/(R0*T_R); A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R

    sqrtRoL = dsqrt(Ro_L); sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)
    
    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    V2_tilda = V_tilda(1)**2 + V_tilda(2)**2 + (V_tilda(3)**2*(dim-2))
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    A_tilda = dsqrt(max(eps, gam1*(H_tilda - 0.5d0*V2_tilda)))
    Ro_tilda = dsqrt(Ro_L*Ro_R)
    T_tilda = A_tilda**2/(k*R0)

    call get_Aq_matrix(dim, R0, k, V_tilda, n, Ro_tilda, T_tilda, Aq)

    lam(1) = Vn_tilda - A_tilda
    lam(2:dim+1) = Vn_tilda
    lam(dim+2) = Vn_tilda + A_tilda
    
    do i = 1, dim+2
        lam(i) = abs(lam(i))
        if (lam(i) < 0.1d0*A_tilda) lam(i) = (lam(i)**2 + (0.1d0*A_tilda)**2)/(0.2d0*A_tilda)
    end do

    r = 0.d0; l = 0.d0
    
    r(1,1) = 1.d0; r(2,1) = V_tilda(1)-A_tilda*n(1); r(3,1) = V_tilda(2)-A_tilda*n(2); r(4,1) = V_tilda(3)-A_tilda*n(dim)*(dim-2); r(dim+2,1) = H_tilda-A_tilda*Vn_tilda
    r(1,2) = 1.d0; r(2,2) = V_tilda(1); r(3,2) = V_tilda(2); r(4,2) = V_tilda(3); r(dim+2,2) = 0.5d0*V2_tilda
    r(2,3) = n(2); r(3,3) = -n(1); r(dim+2,3) = V_tilda(1)*n(2)-V_tilda(2)*n(1)
    if(dim==3) then
       r(2,4) = n(3); r(4,4) = -n(1); r(5,4) = V_tilda(1)*n(3)-V_tilda(3)*n(1)
    endif
    r(1,dim+2) = 1.d0; r(2,dim+2) = V_tilda(1)+A_tilda*n(1); r(3,dim+2) = V_tilda(2)+A_tilda*n(2); r(4,dim+2) = V_tilda(3)+A_tilda*n(dim)*(dim-2); r(dim+2,dim+2) = H_tilda+A_tilda*Vn_tilda

    l(1,1) = 0.5d0*(gam1*V2_tilda/(2.d0*A_tilda**2) + Vn_tilda/A_tilda)
    l(1,2) = -0.5d0*(gam1*V_tilda(1)/A_tilda**2 + n(1)/A_tilda)
    l(1,3) = -0.5d0*(gam1*V_tilda(2)/A_tilda**2 + n(2)/A_tilda)
    if(dim==3) l(1,4) = -0.5d0*(gam1*V_tilda(3)/A_tilda**2 + n(3)/A_tilda)
    l(1,dim+2) = 0.5d0*gam1/A_tilda**2

    l(2,1) = 1.d0 - gam1*V2_tilda/(2.d0*A_tilda**2)
    l(2,2) = gam1*V_tilda(1)/A_tilda**2
    l(2,3) = gam1*V_tilda(2)/A_tilda**2
    if(dim==3) l(2,4) = gam1*V_tilda(3)/A_tilda**2
    l(2,dim+2) = -gam1/A_tilda**2

    l(3,1) = -(V_tilda(1)*n(2)-V_tilda(2)*n(1))
    l(3,2) = n(2); l(3,3) = -n(1)

    if(dim==3) then
       l(4,1) = -(V_tilda(1)*n(3)-V_tilda(3)*n(1))
       l(4,2) = n(3); l(4,4) = -n(1)
    endif

    l(dim+2,1) = 0.5d0*(gam1*V2_tilda/(2.d0*A_tilda**2) - Vn_tilda/A_tilda)
    l(dim+2,2) = -0.5d0*(gam1*V_tilda(1)/A_tilda**2 - n(1)/A_tilda)
    l(dim+2,3) = -0.5d0*(gam1*V_tilda(2)/A_tilda**2 - n(2)/A_tilda)
    if(dim==3) l(dim+2,4) = -0.5d0*(gam1*V_tilda(3)/A_tilda**2 - n(3)/A_tilda)
    l(dim+2,dim+2) = 0.5d0*gam1/A_tilda**2

    AbsA = 0.d0
    do m = 1, dim+2
        do j = 1, dim+2
            do i = 1, dim+2
                AbsA(i,j) = AbsA(i,j) + lam(m) * r(i,m) * l(m,j)
            end do
        end do
    end do

    Mtrx_L = 0.5d0*(Aq + AbsA)
    Mtrx_R = 0.5d0*(Aq - AbsA)

end subroutine

pure subroutine AUSM_stability_universal(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, &
                                    P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g 
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, Vn_L, Vn_R, V2_L, V2_R
    real(8) :: M_L, M_R, M_f, A_f, Vel_i(3)
    real(8) :: tmp1, tmp2(3), tmp5, T_sqrt_const, J_mtrx(dim+2, dim+2)
    integer :: d, i, j, m
    real(8) :: M_tmp(dim+2, dim+2), val

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L); A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2 * (dim-2))
    Vn_L = (U_L*n(1) + V_L*n(2) + W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L

	!=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R); A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2 * (dim-2))
    Vn_R = (U_R*n(1) + V_R*n(2) + W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R
	
	!=== MACH NUMBERS: ===
    A_f = 0.5d0*(A_L + A_R)
    M_L = Vn_L/A_f
    M_R = Vn_R/A_f
    M_f = m_plus(M_L) + m_minus(M_R)

	!=== STABILITY OPERATORS: ===
    Mtrx_L = 0.0d0
    Mtrx_R = 0.0d0
    T_sqrt_const = dsqrt(k*R0)

    if (M_f >= 0.d0) then	
		Vel_i(1) = U_L; Vel_i(2) = V_L; Vel_i(3) = W_L
		
        !=== LEFT CELL OPERATOR: ===
        tmp1 = A_f*M_f/(R0*T_L)
        do d = 1, dim
			tmp2(d) = A_f*Ro_L*teta_L(A_f, n(d), M_L)
		end do
        tmp5 = Ro_L/4.d0*T_sqrt_const/dsqrt(T_L)*(M_f + A_f*gama_L(A_f, M_L)) - A_f*M_f*Ro_L/T_L

			!=== CONTINUITIE EQUATION: ===
        Mtrx_L(1, 1) = tmp1
        Mtrx_L(1, 2:dim+1) = tmp2(1:dim)
        Mtrx_L(1, dim+2) = tmp5

			!=== MOMENTUM EQUATION: ===
        do i = 1, dim
            Mtrx_L(1+i, 1) = n(i)*p_plus(M_L) + tmp1*Vel_i(i)
            do j = 1, dim
                Mtrx_L(1+i, 1+j) = n(i)*psi_L(A_f, M_L, n(j))*P_L + tmp2(j)*Vel_i(i)
            end do
            Mtrx_L(1+i, 1+i) = Mtrx_L(1+i, 1+i) + A_f*M_f*Ro_L
            Mtrx_L(1+i, dim+2) = n(i)*dd_L(A_f, M_L)*P_L + tmp5*Vel_i(i)
        end do

			!=== ENERGY EQUIATION: ===
        Mtrx_L(dim+2, 1) = tmp1*H_L
        do j = 1, dim
            Mtrx_L(dim+2, 1+j) = tmp2(j)*H_L + A_f*M_f*Ro_L*Vel_i(j)
        end do
        Mtrx_L(dim+2, dim+2) = tmp5*H_L + A_f*M_f*Ro_L*cp


		
		
		
		!=== RIGHT CELL OPERATOR: ===
        tmp1 = 0.d0
        do d = 1, dim
			tmp2(d) = A_f*Ro_L*teta_R(A_f, n(d), M_R)
		end do
        tmp5 = Ro_L/4.d0*T_sqrt_const/dsqrt(T_R)*(M_f + A_f*gama_R(A_f, M_R))

			!=== CONTINUITIE EQUATION: ===
        Mtrx_R(1, 1) = tmp1
        Mtrx_R(1, 2:dim+1) = tmp2(1:dim)
        Mtrx_R(1, dim+2) = tmp5

			!=== MOMENTUM EQUATION: ===
        do i = 1, dim
            Mtrx_R(1+i, 1) = n(i)*p_minus(M_R)
            do j = 1, dim
                Mtrx_R(1+i, 1+j) = n(i)*psi_r(A_f, M_R, n(j))*P_R + tmp2(j)*Vel_i(i)
            end do
            Mtrx_R(1+i, dim+2) = n(i)*dd_R(A_f, M_R)*P_R + tmp5*Vel_i(i)
        end do

			!=== ENERGY EQUIATION: ===
        Mtrx_R(dim+2, 1) = tmp1
        do j = 1, dim
            Mtrx_R(dim+2, 1+j) = tmp2(j)*H_L
        end do
        Mtrx_R(dim+2, dim+2) = tmp5*H_L
        
    else
        Vel_i(1) = U_R; Vel_i(2) = V_R; Vel_i(3) = W_R
		
        !=== LEFT CELL OPERATOR: ===
        tmp1 = 0.d0
        do d = 1, dim
			tmp2(d) = A_f*Ro_R*teta_L(A_f, n(d), M_L)
		end do
        tmp5 = Ro_R/4.d0*T_sqrt_const/dsqrt(T_L)*(M_f + A_f*gama_L(A_f, M_L))

			!=== CONTINUITIE EQUATION: ===
        Mtrx_L(1, 1) = tmp1
        Mtrx_L(1, 2:dim+1) = tmp2(1:dim)
        Mtrx_L(1, dim+2) = tmp5

			!=== MOMENTUM EQUATION: ===
        do i = 1, dim
            Mtrx_L(1+i, 1) = n(i)*p_plus(M_L)
            do j = 1, dim
                Mtrx_L(1+i, 1+j) = n(i)*psi_L(A_f, M_L, n(j))*P_L + tmp2(j)*Vel_i(i)
            end do
            Mtrx_L(1+i, 1+i) = Mtrx_L(1+i, 1+i)
            Mtrx_L(1+i, dim+2) = n(i)*dd_L(A_f, M_L)*P_L + tmp5*Vel_i(i)
        end do

			!=== ENERGY EQUIATION: ===
        Mtrx_L(dim+2, 1) = tmp1
        do j = 1, dim
            Mtrx_L(dim+2, 1+j) = tmp2(j)*H_R
        end do
        Mtrx_L(dim+2, dim+2) = tmp5*H_R


		
		
		
		!=== RIGHT CELL OPERATOR: ===
        tmp1 = A_f*M_f/(R0*T_R)
        do d = 1, dim
			tmp2(d) = A_f*Ro_R*teta_R(A_f, n(d), M_R)
		end do
        tmp5 = Ro_R/4.d0*T_sqrt_const/dsqrt(T_R)*(M_f + A_f*gama_R(A_f, M_R)) - A_f*M_f*Ro_R/T_R

			!=== CONTINUITIE EQUATION: ===
        Mtrx_R(1, 1) = tmp1
        Mtrx_R(1, 2:dim+1) = tmp2(1:dim)
        Mtrx_R(1, dim+2) = tmp5

			!=== MOMENTUM EQUATION: ===
        do i = 1, dim
            Mtrx_R(1+i, 1) = n(i)*p_minus(M_R) + tmp1*Vel_i(i)
            do j = 1, dim
                Mtrx_R(1+i, 1+j) = n(i)*psi_R(A_f, M_R, n(j))*P_R + tmp2(j)*Vel_i(i)
            end do
            Mtrx_R(1+i, 1+i) = Mtrx_R(1+i, 1+i) + A_f*M_f*Ro_R
            Mtrx_R(1+i, dim+2) = n(i)*dd_R(A_f, M_R)*P_R + tmp5*Vel_i(i)
        end do

			!=== ENERGY EQUIATION: ===
        Mtrx_R(dim+2, 1) = tmp1*H_R
        do j = 1, dim
            Mtrx_R(dim+2, 1+j) = tmp2(j)*H_R + A_f*M_f*Ro_R*Vel_i(j)
        end do
        Mtrx_R(dim+2, dim+2) = tmp5*H_R + A_f*M_f*Ro_R*cp
    end if
	
	!=== CONS TO PRIM: ===
	Vel_i(1) = U_L; Vel_i(2) = V_L; Vel_i(3) = W_L
	call get_J_mtrx(dim, cp, k, A_L, Vel_i, H_L, Ro_L, T_L, J_mtrx)
	do j = 1, dim + 2
		M_tmp(:, j) = 0.0d0
		do m = 1, dim + 2
			val = Mtrx_L(m, j)
			do i = 1, dim + 2
				M_tmp(i, j) = M_tmp(i, j) + J_mtrx(i, m)*val
			end do
		end do
	end do
	Mtrx_L = M_tmp
	
	Vel_i(1) = U_R; Vel_i(2) = V_R; Vel_i(3) = W_R
    call get_J_mtrx(dim, cp, k, A_R, Vel_i, H_R, Ro_R, T_R, J_mtrx)
    do j = 1, dim + 2
		M_tmp(:, j) = 0.0d0
		do m = 1, dim + 2
			val = Mtrx_R(m, j)
			do i = 1, dim + 2
				M_tmp(i, j) = M_tmp(i, j) + J_mtrx(i, m)*val
			end do
		end do
	end do
	Mtrx_R = M_tmp
    
	contains
	pure real(8) function teta_L(A_f, n_i, M)
		real(8), intent(in) :: A_f, n_i, M
		if (abs(M) > 1.0d0) then
			teta_L = max(0.0d0, sign(1.0d0, M))*(n_i/A_f)
		else
			teta_L = 0.5d0*(M + 1.0d0)*(n_i/A_f)
		end if
	end function

	pure real(8) function teta_R(A_f, n_i, M)
		real(8), intent(in) :: A_f, n_i, M
		if (abs(M) > 1.0d0) then
			teta_R = max(0.0d0, sign(1.0d0, -M))*(n_i/A_f)
		else
			teta_R = -0.5d0*(M - 1.0d0)*(n_i/A_f)
		end if
	end function

	pure real(8) function gama_L(A_f, M)
		real(8), intent(in) :: A_f, M
		if (abs(M) > 1.0d0) then
			gama_L = -max(0.0d0, sign(1.0d0, M))*(M/A_f)
		else
			gama_L = -0.5d0*(M + 1.0d0)*M/A_f
		end if
	end function

	pure real(8) function gama_R(A_f, M)
		real(8), intent(in) :: A_f, M
		if (abs(M) > 1.0d0) then
			gama_R = -max(0.0d0, sign(1.0d0, -M))*(M/A_f)
		else
			gama_R = 0.5d0*(M - 1.0d0)*M/A_f
		end if
	end function

	pure real(8) function psi_L(A_f, M, n_i)
		real(8), intent(in) :: A_f, M, n_i
		psi_L = merge(0.0d0, 0.75d0 * (1.0d0 - M**2)*(n_i/A_f), abs(M) > 1.0d0)
	end function

	pure real(8) function psi_R(A_f, M, n_i)
		real(8), intent(in) :: A_f, M, n_i
		psi_R = merge(0.0d0, 0.75d0 * (M**2 - 1.0d0)*(n_i/A_f), abs(M) > 1.0d0)
	end function

	pure real(8) function dd_L(A_f, M)
		real(8), intent(in) :: A_f, M
		dd_L = merge(0.0d0, 0.75d0*(M**3 - M)/A_f, abs(M) > 1.0d0)
	end function

	pure real(8) function dd_R(A_f, M)
		real(8), intent(in) :: A_f, M
		dd_R = merge(0.0d0, 0.75d0*(M**3 - M)/A_f, abs(M) > 1.0d0)
	end function

	pure real(8) function m_plus(M)
		real(8), intent(in) :: M
		m_plus = merge(0.5d0*(M + abs(M)), 0.25d0*(M + 1.0d0)**2, abs(M) >= 1.0d0)
	end function

	pure real(8) function m_minus(M)
		real(8), intent(in) :: M
		m_minus = merge(0.5d0*(M - abs(M)), -0.25d0*(M - 1.0d0)**2, abs(M) >= 1.0d0)
	end function

	pure real(8) function p_plus(M)
		real(8), intent(in) :: M
		p_plus = merge(0.5d0 * (1.0d0 + sign(1.0d0, M)), 0.25d0*(M + 1.0d0)**2*(2.0d0 - M), abs(M) >= 1.0d0)
	end function

	pure real(8) function p_minus(M)
		real(8), intent(in) :: M
		p_minus = merge(0.5d0*(1.0d0 - sign(1.0d0, M)), 0.25d0*(M - 1.0d0)**2*(2.0d0 + M), abs(M) >= 1.0d0)
	end function
end subroutine

pure subroutine AUSM_plus_stability_universal(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, &
										P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g 
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, Vn_L, Vn_R, V2_L, V2_R
    real(8) :: M_L, M_R, M_f, A_f, Vel_i(3)
    real(8) :: tmp1, tmp2(3), tmp5, T_sqrt_const, J_mtrx(dim+2, dim+2)
    integer :: d, i, j, m
    real(8) :: M_tmp(dim+2, dim+2), val

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L); A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2 * (dim-2))
    Vn_L = (U_L*n(1) + V_L*n(2) + W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L

	!=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R); A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2 * (dim-2))
    Vn_R = (U_R*n(1) + V_R*n(2) + W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R
	
	!=== MACH NUMBERS: ===
    A_f = dsqrt(A_L*A_R)
	M_L = Vn_L/A_f
	M_R = Vn_R/A_f
	M_f = m_plus(M_L) + m_minus(M_R)

	!=== STABILITY OPERATORS: ===
    Mtrx_L = 0.0d0
    Mtrx_R = 0.0d0
    T_sqrt_const = dsqrt(k*R0)

    if (M_f >= 0.d0) then	
		Vel_i(1) = U_L; Vel_i(2) = V_L; Vel_i(3) = W_L
		
        !=== LEFT CELL OPERATOR: ===
        tmp1 = A_f*M_f/(R0*T_L)
        do d = 1, dim
			tmp2(d) = A_f*Ro_L*teta_L(A_f, n(d), M_L)
		end do
        tmp5 = Ro_L/4.d0*T_sqrt_const/dsqrt(T_L*A_L/A_R)*(M_f + A_f*gama_L(A_f, M_L)) - A_f*M_f*Ro_L/T_L

			!=== CONTINUITIE EQUATION: ===
        Mtrx_L(1, 1) = tmp1
        Mtrx_L(1, 2:dim+1) = tmp2(1:dim)
        Mtrx_L(1, dim+2) = tmp5

			!=== MOMENTUM EQUATION: ===
        do i = 1, dim
            Mtrx_L(1+i, 1) = n(i)*p_plus(M_L) + tmp1*Vel_i(i)
            do j = 1, dim
                Mtrx_L(1+i, 1+j) = n(i)*psi_L(A_f, M_L, n(j))*P_L + tmp2(j)*Vel_i(i)
            end do
            Mtrx_L(1+i, 1+i) = Mtrx_L(1+i, 1+i) + A_f*M_f*Ro_L
            Mtrx_L(1+i, dim+2) = n(i)*dd_L(A_f, M_L)*P_L + tmp5*Vel_i(i)
        end do

			!=== ENERGY EQUIATION: ===
        Mtrx_L(dim+2, 1) = tmp1*H_L
        do j = 1, dim
            Mtrx_L(dim+2, 1+j) = tmp2(j)*H_L + A_f*M_f*Ro_L*Vel_i(j)
        end do
        Mtrx_L(dim+2, dim+2) = tmp5*H_L + A_f*M_f*Ro_L*cp


		
		
		
		!=== RIGHT CELL OPERATOR: ===
        tmp1 = 0.d0
        do d = 1, dim
			tmp2(d) = A_f*Ro_L*teta_R(A_f, n(d), M_R)
		end do
        tmp5 = Ro_L/4.d0*T_sqrt_const/dsqrt(T_R*A_R/A_L)*(M_f + A_f*gama_R(A_f, M_R))

			!=== CONTINUITIE EQUATION: ===
        Mtrx_R(1, 1) = tmp1
        Mtrx_R(1, 2:dim+1) = tmp2(1:dim)
        Mtrx_R(1, dim+2) = tmp5

			!=== MOMENTUM EQUATION: ===
        do i = 1, dim
            Mtrx_R(1+i, 1) = n(i)*p_minus(M_R)
            do j = 1, dim
                Mtrx_R(1+i, 1+j) = n(i)*psi_R(A_f, M_R, n(j))*P_R + tmp2(j)*Vel_i(i)
            end do
            Mtrx_R(1+i, dim+2) = n(i)*dd_R(A_f, M_R)*P_R + tmp5*Vel_i(i)
        end do

			!=== ENERGY EQUIATION: ===
        Mtrx_R(dim+2, 1) = tmp1
        do j = 1, dim
            Mtrx_R(dim+2, 1+j) = tmp2(j)*H_L
        end do
        Mtrx_R(dim+2, dim+2) = tmp5*H_L
        
    else
        Vel_i(1) = U_R; Vel_i(2) = V_R; Vel_i(3) = W_R
		
        !=== LEFT CELL OPERATOR: ===
        tmp1 = 0.d0
        do d = 1, dim
			tmp2(d) = A_f*Ro_R*teta_L(A_f, n(d), M_L)
		end do
        tmp5 = Ro_R/4.d0*T_sqrt_const/dsqrt(T_L*A_L/A_R)*(M_f + A_f*gama_L(A_f, M_L))

			!=== CONTINUITIE EQUATION: ===
        Mtrx_L(1, 1) = tmp1
        Mtrx_L(1, 2:dim+1) = tmp2(1:dim)
        Mtrx_L(1, dim+2) = tmp5

			!=== MOMENTUM EQUATION: ===
        do i = 1, dim
            Mtrx_L(1+i, 1) = n(i)*p_plus(M_L)
            do j = 1, dim
                Mtrx_L(1+i, 1+j) = n(i)*psi_L(A_f, M_L, n(j))*P_L + tmp2(j)*Vel_i(i)
            end do
            Mtrx_L(1+i, 1+i) = Mtrx_L(1+i, 1+i)
            Mtrx_L(1+i, dim+2) = n(i)*dd_L(A_f, M_L)*P_L + tmp5*Vel_i(i)
        end do

			!=== ENERGY EQUIATION: ===
        Mtrx_L(dim+2, 1) = tmp1
        do j = 1, dim
            Mtrx_L(dim+2, 1+j) = tmp2(j)*H_R
        end do
        Mtrx_L(dim+2, dim+2) = tmp5*H_R


		
		
		
		!=== RIGHT CELL OPERATOR: ===
        tmp1 = A_f*M_f/(R0*T_R)
        do d = 1, dim
			tmp2(d) = A_f*Ro_R*teta_R(A_f, n(d), M_R)
		end do
        tmp5 = Ro_R/4.d0*T_sqrt_const/dsqrt(T_R*A_R/A_L)*(M_f + A_f*gama_R(A_f, M_R)) - A_f*M_f*Ro_R/T_R

			!=== CONTINUITIE EQUATION: ===
        Mtrx_R(1, 1) = tmp1
        Mtrx_R(1, 2:dim+1) = tmp2(1:dim)
        Mtrx_R(1, dim+2) = tmp5

			!=== MOMENTUM EQUATION: ===
        do i = 1, dim
            Mtrx_R(1+i, 1) = n(i)*p_minus(M_R) + tmp1*Vel_i(i)
            do j = 1, dim
                Mtrx_R(1+i, 1+j) = n(i)*psi_R(A_f, M_R, n(j))*P_R + tmp2(j)*Vel_i(i)
            end do
            Mtrx_R(1+i, 1+i) = Mtrx_R(1+i, 1+i) + A_f*M_f*Ro_R
            Mtrx_R(1+i, dim+2) = n(i)*dd_R(A_f, M_R)*P_R + tmp5*Vel_i(i)
        end do

			!=== ENERGY EQUIATION: ===
        Mtrx_R(dim+2, 1) = tmp1*H_R
        do j = 1, dim
            Mtrx_R(dim+2, 1+j) = tmp2(j)*H_R + A_f*M_f*Ro_R*Vel_i(j)
        end do
        Mtrx_R(dim+2, dim+2) = tmp5*H_R + A_f*M_f*Ro_R*cp
    end if


    !=== CONS TO PRIM: ===
    Vel_i(1) = U_L; Vel_i(2) = V_L; Vel_i(3) = W_L
	call get_J_mtrx(dim, cp, k, A_L, Vel_i, H_L, Ro_L, T_L, J_mtrx)
	do j = 1, dim + 2
		M_tmp(:, j) = 0.0d0
		do m = 1, dim + 2
			val = Mtrx_L(m, j)
			do i = 1, dim + 2
				M_tmp(i, j) = M_tmp(i, j) + J_mtrx(i, m)*val
			end do
		end do
	end do
	Mtrx_L = M_tmp
	
	Vel_i(1) = U_R; Vel_i(2) = V_R; Vel_i(3) = W_R
    call get_J_mtrx(dim, cp, k, A_R, Vel_i, H_R, Ro_R, T_R, J_mtrx)
    do j = 1, dim + 2
		M_tmp(:, j) = 0.0d0
		do m = 1, dim + 2
			val = Mtrx_R(m, j)
			do i = 1, dim + 2
				M_tmp(i, j) = M_tmp(i, j) + J_mtrx(i, m)*val
			end do
		end do
	end do
	Mtrx_R = M_tmp
    
	contains
	pure real(8) function poly_base(M)
		real(8), intent(in) :: M
		poly_base = (M**2 - 1.0d0)**2
	end function

	pure real(8) function teta_L(A_f, n_i, M)
		real(8), intent(in) :: A_f, n_i, M
		real(8) :: res
		if (M > 1.0d0) then
			res = 1.0d0
		else if (M < -1.0d0) then
			res = 0.0d0
		else
			res = (0.5d0*(M + 1.0d0) + 0.5d0*M*(M**2 - 1.0d0))
		end if
		teta_L = res*(n_i/A_f)
	end function

	pure real(8) function teta_R(A_f, n_i, M)
		real(8), intent(in) :: A_f, n_i, M
		real(8) :: res
		if (M > 1.0d0) then
			res = 0.0d0
		else if (M < -1.0d0) then
			res = 1.0d0
		else
			res = (-0.5d0*(M - 1.0d0) - 0.5d0*M*(M**2 - 1.0d0))
		end if
		teta_R = res*(n_i/A_f)
	end function

	pure real(8) function gama_L(A_f, M)
		real(8), intent(in) :: A_f, M
		real(8) :: res
		if (M > 1.0d0) then
			res = -M
		else if (M < -1.0d0) then
			res = 0.0d0
		else
			res = -0.5d0*M*(M + 1.0d0) - 0.5d0*(M**2)*(M**2 - 1.0d0)
		end if
		gama_L = res/A_f
	end function

	pure real(8) function gama_R(A_f, M)
		real(8), intent(in) :: A_f, M
		real(8) :: res
		if (M > 1.0d0) then
			res = 0.0d0
		else if (M < -1.0d0) then
			res = -M
		else
			res = 0.5d0*M*(M - 1.0d0) + 0.5d0*(M**2)*(M**2 - 1.0d0)
		end if
		gama_R = res/A_f
	end function

	pure real(8) function psi_L(A_f, M, n_i)
		real(8), intent(in) :: A_f, M, n_i
		real(8) :: dP
		dP = merge(0.0d0, 0.75d0*(1.0d0-M**2) + 0.1875d0*(5.0d0*M**4 - 6.0d0*M**2 + 1.0d0), abs(M) > 1.0d0)
		psi_L = dP*(n_i/A_f)
	end function

	pure real(8) function psi_R(A_f, M, n_i)
		real(8), intent(in) :: A_f, M, n_i
		real(8) :: dP
		dP = merge(0.0d0, 0.75d0*(M**2-1.0d0) - 0.1875d0*(5.0d0*M**4 - 6.0d0*M**2 + 1.0d0), abs(M) > 1.0d0)
		psi_R = dP*(n_i/A_f)
	end function

	pure real(8) function dd_L(A_f, M)
		real(8), intent(in) :: A_f, M
		real(8) :: dP
		dP = merge(0.0d0, 0.75d0*(M**3 - M) + 0.1875d0*(-5.0d0*M**5 + 6.0d0*M**3 - M), abs(M) > 1.0d0)
		dd_L = dP/A_f
	end function

	pure real(8) function dd_R(A_f, M)
		real(8), intent(in) :: A_f, M
		real(8) :: dP
		dP = merge(0.0d0, 0.75d0*(M**3 - M) - 0.1875d0*(-5.0d0*M**5 + 6.0d0*M**3 - M), abs(M) > 1.0d0)
		dd_R = dP/A_f
	end function

	pure real(8) function m_plus(M)
		real(8), intent(in) :: M
		if (M > 1.0d0) then
			m_plus = M
		else if (M < -1.0d0) then
			m_plus = 0.0d0
		else
			m_plus = 0.25d0*(M + 1.0d0)**2 + 0.125d0*(M**2 - 1.0d0)**2
		end if
	end function

	pure real(8) function m_minus(M)
		real(8), intent(in) :: M
		if (M > 1.0d0) then
			m_minus = 0.0d0
		else if (M < -1.0d0) then
			m_minus = M
		else
			m_minus = -0.25d0*(M - 1.0d0)**2 - 0.125d0*(M**2 - 1.0d0)**2
		end if
	end function

	pure real(8) function p_plus(M)
		real(8), intent(in) :: M
		if (M >= 1.0d0) then
			p_plus = 1.0d0
		else if (M <= -1.0d0) then
			p_plus = 0.0d0
		else
			p_plus = 0.25d0*(M + 1.0d0)**2 * (2.0d0 - M) + 0.1875d0*M*(M**2 - 1.0d0)**2
		end if
	end function

	pure real(8) function p_minus(M)
		real(8), intent(in) :: M
		if (M >= 1.0d0) then
			p_minus = 0.0d0
		else if (M <= -1.0d0) then
			p_minus = 1.0d0
		else
			p_minus = 0.25d0*(M - 1.0d0)**2 * (2.0d0 + M) - 0.1875d0*M*(M**2 - 1.0d0)**2
		end if
	end function

end subroutine


!=======================================================================
!=============== STABILITY OPERATORS CONSERVATIVE VARIABLES ============
!=======================================================================
pure subroutine HLL_stability_conservative(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L,&
                                            P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)
	
	real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, Vn_L, Vn_R, V2_L, V2_R
    real(8) :: S_L, S_R, sqrtRoL, sqrtRoR, invSumSqrt
    real(8) :: H_tilda, V_tilda(3), Vn_tilda, A_tilda, V2_tilda, Ro_tilda, T_tilda, E_tilda, P_tilda
    real(8) :: Ro_avg, T_avg, V_avg(3)
    real(8) :: Au(dim+2, dim+2), coeff_A, coeff_E
    integer :: i
    
    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L); A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2 * (dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L

	!=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R); A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R

    !=== ROE AVERAGING: ===
    Ro_tilda = dsqrt(Ro_L*Ro_R)
    sqrtRoL = dsqrt(Ro_L); sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)
    
    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    V2_tilda = V_tilda(1)**2 + V_tilda(2)**2 + (V_tilda(3)**2*(dim-2))
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    A_tilda = dsqrt((k-1.d0)*(H_tilda - 0.5d0*V2_tilda))
	T_tilda = A_tilda**2/(k*R0)
	E_tilda = cv*T_tilda + 0.5d0*V2_tilda
	P_tilda = (k - 1.d0)/k*(H_tilda - 0.5d0*V2_tilda)*Ro_tilda
	
    S_L = min(Vn_L - A_L, Vn_tilda - A_tilda)
    S_R = max(Vn_R + A_R, Vn_tilda + A_tilda)
    
    !=== AVG MATRIX ON FACE: ===
    call get_Au_matrix(dim, k, V_tilda, n, Ro_tilda, E_tilda, P_tilda, Au)

    !=== STABILITY OPERATORS: ===
    if (S_L >= 0.0d0) then
        Mtrx_L = Au
        Mtrx_R = 0.0d0
    else if (S_R <= 0.0d0) then
        Mtrx_L = 0.0d0
        Mtrx_R = Au
    else
        coeff_A = 0.5d0*(1.0d0 + (S_R + S_L)/(S_R - S_L))
        coeff_E = S_R*S_L/(S_R - S_L)
        
        Mtrx_L = coeff_A*Au
        Mtrx_R = (1.0d0 - coeff_A)*Au
        
        do i = 1, dim+2
            Mtrx_L(i, i) = Mtrx_L(i, i) - coeff_E
            Mtrx_R(i, i) = Mtrx_R(i, i) + coeff_E
        end do
    end if
end subroutine

pure subroutine HLLC_stability_conservative(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, &
										 P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)

    real(8) :: Ro_L, Ro_R, A_L, A_R, H_L, H_R, Vn_L, Vn_R, V2_L, V2_R
    real(8) :: S_L, S_R, S_star, sqrtRoL, sqrtRoR, invSumSqrt
    real(8) :: H_tilda, V_tilda(3), Vn_tilda, A_tilda, V2_tilda, Ro_tilda, T_tilda, E_tilda, P_tilda
    real(8) :: Ro_avg, T_avg, V_avg(3)
    real(8) :: Au(dim+2, dim+2), coeff_A, coeff_E
    integer :: i, j

    !=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L); A_L = dsqrt(k*R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2 * (dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L
	
	!=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R); A_R = dsqrt(k*R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R

    !=== ROE AVERAGING: ===
    Ro_tilda = dsqrt(Ro_L*Ro_R)
    sqrtRoL = dsqrt(Ro_L); sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)
    
    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    V2_tilda = V_tilda(1)**2 + V_tilda(2)**2 + (V_tilda(3)**2*(dim-2))
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    A_tilda = dsqrt((k-1.d0)*(H_tilda - 0.5d0*V2_tilda))
	T_tilda = A_tilda**2/(k*R0)
	E_tilda = cv*T_tilda + 0.5d0*V2_tilda
	P_tilda = (k - 1.d0)/k*(H_tilda - 0.5d0*V2_tilda)*Ro_tilda
	
    S_L = min(Vn_L - A_L, Vn_tilda - A_tilda)
    S_R = max(Vn_R + A_R, Vn_tilda + A_tilda)

    S_star = ((P_R - P_L) + Ro_L*Vn_L*(S_L - Vn_L) - Ro_R*Vn_R*(S_R - Vn_R))/&
             (Ro_L*(S_L - Vn_L) - Ro_R*(S_R - Vn_R))

    !=== AVG MATRIX ON FACE: ===
    call get_Au_matrix(dim, k, V_tilda, n, Ro_tilda, E_tilda, P_tilda, Au)

    !=== STABILITY OPERATORS: ===
    Mtrx_L = 0.0d0
    Mtrx_R = 0.0d0

    if (S_L >= 0.0d0) then
        Mtrx_L = Au
        
    else if (S_R <= 0.0d0) then
        Mtrx_R = Au
        
    else if (S_star >= 0.0d0) then
        coeff_A = -S_star/(S_L - S_star)
        coeff_E = (S_star*S_L)/(S_L - S_star)
        Mtrx_L = coeff_A*Au
        do i = 1, dim+2
            Mtrx_L(i, i) = Mtrx_L(i, i) + coeff_E
        end do
        
    else
        coeff_A = -S_star/(S_R - S_star)
        coeff_E = (S_star*S_R)/(S_R - S_star)
        Mtrx_R = coeff_A*Au
        do i = 1, dim+2
            Mtrx_R(i, i) = Mtrx_R(i, i) + coeff_E
        end do
        
    end if

end subroutine

pure subroutine ROE_stability_conservative(dim, k, R0, cp, cv, P_L, U_L, V_L, W_L, T_L, &
                                         P_R, U_R, V_R, W_R, T_R, g, n, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, cp, cv, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)

    real(8) :: Ro_L, Ro_R, H_L, H_R, Vn_L, Vn_R, V2_L, V2_R
    real(8) :: Ro_tilda, V_tilda(3), H_tilda, A_tilda, V2_tilda, Vn_tilda, E_tilda, P_tilda
    real(8) :: Au(dim+2, dim+2), T2_plus(dim+2, dim+2), T2_minus(dim+2, dim+2)
    real(8) :: r1(5), r5(5), l1(5), l5(5), lambda1, lambda5
    real(8) :: sqrtRoL, sqrtRoR, invSumSqrt, phi2, inv2A2, k1, eps, delta_ij
    integer :: i, j

	!=== LEFT STATE: ===
    Ro_L = P_L/(R0*T_L)
    V2_L = U_L**2 + V_L**2 + (W_L**2 * (dim-2))
    Vn_L = U_L*n(1) + V_L*n(2) + (W_L*n(dim)*(dim-2))
    H_L = cp*T_L + 0.5d0*V2_L
	
	!=== RIGHT STATE: ===
    Ro_R = P_R/(R0*T_R)
    V2_R = U_R**2 + V_R**2 + (W_R**2*(dim-2))
    Vn_R = U_R*n(1) + V_R*n(2) + (W_R*n(dim)*(dim-2))
    H_R = cp*T_R + 0.5d0*V2_R
    
	!=== ROE AVERAGING: ===
    Ro_tilda = dsqrt(Ro_L*Ro_R)
    sqrtRoL = dsqrt(Ro_L); sqrtRoR = dsqrt(Ro_R)
    invSumSqrt = 1.0d0/(sqrtRoL + sqrtRoR)
    
    V_tilda(1) = (sqrtRoL*U_L + sqrtRoR*U_R)*invSumSqrt
    V_tilda(2) = (sqrtRoL*V_L + sqrtRoR*V_R)*invSumSqrt
    V_tilda(3) = (sqrtRoL*W_L + sqrtRoR*W_R)*invSumSqrt
    V2_tilda = V_tilda(1)**2 + V_tilda(2)**2 + (V_tilda(3)**2*(dim-2))
    Vn_tilda = V_tilda(1)*n(1) + V_tilda(2)*n(2) + (V_tilda(3)*n(dim)*(dim-2))
    
    H_tilda = (sqrtRoL*H_L + sqrtRoR*H_R)*invSumSqrt
    A_tilda = dsqrt((k-1.d0)*(H_tilda - 0.5d0*V2_tilda))
	E_tilda = H_tilda/k + (k-1.d0)/k*0.5d0*V2_tilda 
    P_tilda = (k-1.d0)/k*(H_tilda - 0.5d0*V2_tilda)*Ro_tilda
    
    !=== AVG MATRIX ON FACE: ===
	!call get_Au_matrix(dim, R0, k, (1.d0-g)*[U_L, V_L, W_L] + g*[U_R, V_R, W_R], &
    !                   n, (1.d0-g)*Ro_L + g*Ro_R, (1.d0-g)*T_L + g*T_R, Aq)
    call get_Au_matrix(dim, k, V_tilda, n, Ro_tilda, E_tilda, P_tilda, Au)
    
    
    !=== EIGEN VALUES: ===
    eps = 0.05d0*A_tilda
    lambda1 = Vn_tilda - A_tilda
    lambda5 = Vn_tilda + A_tilda
    
    k1 = (k - 1.0d0)
    phi2 = 0.5d0*k1*V2_tilda											!phi2 = k1*V2_tilda
    inv2A2 = 0.5d0/(A_tilda**2)
		
	!=== 1ST RIGHT EIGEN VECTOR: ===
    r1(1) = 1.0d0
    r1(2:dim+1) = V_tilda(1:dim) - A_tilda*n(1:dim)
    r1(dim+2) = H_tilda - Vn_tilda*A_tilda
    
    !=== 5TH RIGHT EIGEN VECTOR: ===
    r5(1) = 1.0d0
    r5(2:dim+1) = V_tilda(1:dim) + A_tilda*n(1:dim)
    r5(dim+2) = H_tilda + Vn_tilda*A_tilda
    
    !=== 1ST LEFT EIGEN VECTOR: ===
    l1(1) = (phi2 + A_tilda*Vn_tilda)*inv2A2							
    l1(2:dim+1) = -(k1*V_tilda(1:dim) + A_tilda*n(1:dim))*inv2A2		!l1(2:dim+1) = -(2.d0*k1*V_tilda(1:dim) + A_tilda*n(1:dim))*inv2A2
    l1(dim+2) = k1*inv2A2												!l1(dim+2) = 2.d0*k1*inv2A2
    
    !=== 5TH LEFT EIGEN VECTOR: ===
    l5(1) = (phi2 - A_tilda*Vn_tilda)*inv2A2
    l5(2:dim+1) = -(k1*V_tilda(1:dim) - A_tilda*n(1:dim))*inv2A2		!l5(2:dim+1) = -(2.d0*k1*V_tilda(1:dim) - A_tilda*n(1:dim))*inv2A2
    l5(dim+2) = k1*inv2A2												!l5(dim+2) = 2.d0*k1*inv2A2

    !=== STABILITY OPERATORS: ===
    T2_plus  = 0.d0
	T2_minus = 0.d0
	
	!if (Vn_tilda >= 0.d0) then
	!	do j = 1, dim+2
	!		do i = 1, dim+2
	!			T2_plus(i, j) = r1(i)*(lambda1 - Vn_tilda)*l1(j) - &
	!							r5(i)*min(Vn_tilda, Vn_tilda - lambda5)*l5(j)
	!			T2_minus(i, j) = r5(i)*min(0.d0, lambda5)*l5(j)
	!		end do
	!		T2_plus(j, j) = T2_plus(j, j) + Vn_tilda
	!	end do
	!else
	!	do j = 1, dim+2
	!		do i = 1, dim+2
	!			T2_plus(i, j) = r1(i)*max(0.d0, lambda1)*l1(j)
	!			T2_minus(i, j) = r5(i)*(lambda5 - Vn_tilda)*l5(j) - &
	!							 r1(i)*max(Vn_tilda, Vn_tilda - lambda1)*l1(j)
	!		end do
	!		T2_minus(j, j) = T2_minus(j, j) + Vn_tilda
	!	end do
	!end if
	    
    do j = 1, dim+2
        do i = 1, dim+2
            T2_plus(i, j) = r1(i)*max(0.d0, lambda1)*l1(j) + r5(i)*max(0.d0, lambda5)*l5(j)
            T2_minus(i, j) = r1(i)*min(0.d0, lambda1)*l1(j) + r5(i)*min(0.d0, lambda5)*l5(j)
        end do
    end do
    
    do j = 1, dim+2
        do i = 1, dim+2
            delta_ij = 0.d0
            if (i == j) delta_ij = 1.d0
            
            T2_plus(i, j)  = T2_plus(i, j) + max(0.d0, Vn_tilda)*(delta_ij - (r1(i)*l1(j) + r5(i)*l5(j)))
            T2_minus(i, j) = T2_minus(i, j) + min(0.d0, Vn_tilda)*(delta_ij - (r1(i)*l1(j) + r5(i)*l5(j)))
        end do
    end do
    
    Mtrx_L = 0.5d0*(Au + (T2_plus - T2_minus))
    Mtrx_R = 0.5d0*(Au - (T2_plus - T2_minus))

end subroutine





!=======================================================================
!=========== VISCOUS STABILITY OPERATORS PRIMITIVE VARIABLES ===========
!=======================================================================
pure subroutine VISCOUS_stability_universal(dim, k, R0, mu, lambda, cp, P_L, U_L, V_L, W_L, T_L, &
											P_R, U_R, V_R, W_R, T_R, g, n, ksi, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, mu, lambda, cp, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim), ksi(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)
	
	real(8) :: P_face, Vel_face(3), T_face, A_face, Ro_face
	real(8) :: RL
	integer :: d
	
	!=== FACE VALUES: ===
	P_face = (1.d0 - g)*P_L + g*P_R
	Vel_face(1) = (1.d0 - g)*U_L + g*U_R
	Vel_face(2) = (1.d0 - g)*V_L + g*V_R
	Vel_face(3) = (1.d0 - g)*W_L + g*W_R
	T_face = (1.d0 - g)*T_L + g*T_R
	A_face = dsqrt(k*R0*T_face)
	Ro_face = P_face/(T_face*R0)
	
	RL = norm2(ksi)
		
	!=== STABILITY OPERATOR: ===
	call get_Aq_visc_matrix(dim, cp, mu, lambda, RL, Vel_face, Ro_face, T_face, A_face, Mtrx_L)	
	
	Mtrx_R = -Mtrx_L
	
end subroutine


!=======================================================================
!========= VISCOUS STABILITY OPERATORS CONSERVATIVE VARIABLES ==========
!=======================================================================
pure subroutine VISCOUS_stability_conservative(dim, k, R0, mu, lambda, cp, P_L, U_L, V_L, W_L, T_L, &
											   P_R, U_R, V_R, W_R, T_R, g, n, ksi, Mtrx_L, Mtrx_R)
    implicit none
    integer, intent(in) :: dim
    real(8), intent(in) :: k, R0, mu, lambda, cp, g
    real(8), intent(in) :: P_L, U_L, V_L, W_L, T_L
    real(8), intent(in) :: P_R, U_R, V_R, W_R, T_R
    real(8), intent(in) :: n(dim), ksi(dim)
    real(8), intent(inout) :: Mtrx_L(dim+2, dim+2), Mtrx_R(dim+2, dim+2)
	
	real(8) :: P_face, Vel_face(3), T_face, A_face, Ro_face, H_face, V2_face
	real(8) :: RL
	integer :: d
	
	!=== FACE VALUES: ===
	P_face = (1.d0 - g)*P_L + g*P_R
	Vel_face(1) = (1.d0 - g)*U_L + g*U_R
	Vel_face(2) = (1.d0 - g)*V_L + g*V_R
	Vel_face(3) = (1.d0 - g)*W_L + g*W_R
	T_face = (1.d0 - g)*T_L + g*T_R
	A_face = dsqrt(k*R0*T_face)
	Ro_face = P_face/(T_face*R0)
    V2_face = Vel_face(1)**2 + Vel_face(2)**2 + Vel_face(3)**2*(dim-2)
    H_face = cp*T_face + 0.5d0*V2_face
	
	RL = norm2(ksi)
		
	!=== STABILITY OPERATOR: ===
	call get_Au_visc_matrix(dim, cp, mu, lambda, RL, Vel_face, Ro_face, T_face, A_face, H_face, Mtrx_L)
		
	Mtrx_R = -Mtrx_L	
end subroutine

end module

