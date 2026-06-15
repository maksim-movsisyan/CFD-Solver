module ADdiff_reiman_solver_module
implicit none
contains
!=======================================================================
!================ ALG-DIFF APPROXIMATE REIMAN SOLVERS ==================
!=======================================================================
!  Differentiation of hll_reiman_solver in forward (tangent) mode:
!   variations   of useful results: fluxes
!   with respect to varying inputs: w_l w_r t_l t_r v_l v_r p_l
!                u_l p_r u_r
!   RW status of diff variables: w_l:in w_r:in t_l:in t_r:in v_l:in
!                v_r:in p_l:in fluxes:out u_l:in p_r:in u_r:in
  PURE SUBROUTINE HLL_REIMAN_SOLVER_D(dim, k, r0, cp, cv, p_l, p_ld, u_l, &
&   u_ld, v_l, v_ld, w_l, w_ld, t_l, t_ld, p_r, p_rd, u_r, u_rd, v_r, &
&   v_rd, w_r, w_rd, t_r, t_rd, n, fluxesd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: k, r0, cp, cv
    REAL*8, INTENT(IN) :: p_l, u_l, v_l, w_l, t_l
    REAL*8, INTENT(IN) :: p_ld, u_ld, v_ld, w_ld, t_ld
    REAL*8, INTENT(IN) :: p_r, u_r, v_r, w_r, t_r
    REAL*8, INTENT(IN) :: p_rd, u_rd, v_rd, w_rd, t_rd
    REAL*8, INTENT(IN) :: n(dim)
    REAL*8, INTENT(OUT) :: fluxesd(dim+2)
    REAL*8 :: ro_l, ro_r, a_l, a_r, h_l, h_r, e_l, e_r
    REAL*8 :: ro_ld, ro_rd, a_ld, a_rd, h_ld, h_rd, e_ld, e_rd
    REAL*8 :: vn_l, vn_r, v2_l, v2_r
    REAL*8 :: vn_ld, vn_rd, v2_ld, v2_rd
    REAL*8 :: s_l, s_r, s_plus, s_minus, sqrtrol, sqrtror, invsumsqrt, &
&   invsums
    REAL*8 :: s_ld, s_rd, s_plusd, s_minusd, sqrtrold, sqrtrord, &
&   invsumsqrtd, invsumsd
    REAL*8 :: h_tilda, v_tilda(3), vn_tilda, a_tilda
    REAL*8 :: h_tildad, v_tildad(3), vn_tildad, a_tildad
    REAL*8 :: f_l(5), f_r(5), w_vecl(5), w_vecr(5)
    REAL*8 :: f_ld(5), f_rd(5), w_vecld(5), w_vecrd(5)
    INTEGER :: d
    INTRINSIC DSQRT
    INTRINSIC MIN
    INTRINSIC MAX
    REAL*8 :: arg1
    REAL*8 :: arg1d
    REAL*8 :: temp
    DOUBLE PRECISION :: temp0
    REAL*8, DIMENSION(dim+2) :: temp1
    
!=== LEFT STATE: ===
    temp = p_l/(r0*t_l)
    ro_ld = (p_ld-temp*r0*t_ld)/(r0*t_l)
    ro_l = temp
    arg1d = k*r0*t_ld
    arg1 = k*r0*t_l
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_ld = 0.0_8
    ELSE
      a_ld = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_l = temp0
    v2_ld = 2*u_l*u_ld + 2*v_l*v_ld + (dim-2)*2*w_l*w_ld
    v2_l = u_l**2 + v_l**2 + w_l**2*(dim-2)
    vn_ld = n(1)*u_ld + n(2)*v_ld + n(dim)*(dim-2)*w_ld
    vn_l = u_l*n(1) + v_l*n(2) + w_l*n(dim)*(dim-2)
    e_ld = cv*t_ld + 0.5d0*v2_ld
    e_l = cv*t_l + 0.5d0*v2_l
    h_ld = cp*t_ld + 0.5d0*v2_ld
    h_l = cp*t_l + 0.5d0*v2_l
    
!=== RIGHT STATE: ===
    temp = p_r/(r0*t_r)
    ro_rd = (p_rd-temp*r0*t_rd)/(r0*t_r)
    ro_r = temp
    arg1d = k*r0*t_rd
    arg1 = k*r0*t_r
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_rd = 0.0_8
    ELSE
      a_rd = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_r = temp0
    v2_rd = 2*u_r*u_rd + 2*v_r*v_rd + (dim-2)*2*w_r*w_rd
    v2_r = u_r**2 + v_r**2 + w_r**2*(dim-2)
    vn_rd = n(1)*u_rd + n(2)*v_rd + n(dim)*(dim-2)*w_rd
    vn_r = u_r*n(1) + v_r*n(2) + w_r*n(dim)*(dim-2)
    e_rd = cv*t_rd + 0.5d0*v2_rd
    e_r = cv*t_r + 0.5d0*v2_r
    h_rd = cp*t_rd + 0.5d0*v2_rd
    h_r = cp*t_r + 0.5d0*v2_r
    
!=== ROE AVERAGING: ===
    temp0 = DSQRT(ro_l)
    IF (ro_l .EQ. 0.0) THEN
      sqrtrold = 0.0_8
    ELSE
      sqrtrold = ro_ld/(2.D0*DSQRT(ro_l))
    END IF
    sqrtrol = temp0
    temp0 = DSQRT(ro_r)
    IF (ro_r .EQ. 0.0) THEN
      sqrtrord = 0.0_8
    ELSE
      sqrtrord = ro_rd/(2.D0*DSQRT(ro_r))
    END IF
    sqrtror = temp0
    temp = 1.0/(sqrtrol+sqrtror)
    invsumsqrtd = -(temp*(sqrtrold+sqrtrord)/(sqrtrol+sqrtror))
    invsumsqrt = temp
    v_tildad = 0.0_8
    temp = sqrtrol*u_l + sqrtror*u_r
    v_tildad(1) = invsumsqrt*(u_l*sqrtrold+sqrtrol*u_ld+u_r*sqrtrord+&
&     sqrtror*u_rd) + temp*invsumsqrtd
    v_tilda(1) = temp*invsumsqrt
    temp = sqrtrol*v_l + sqrtror*v_r
    v_tildad(2) = invsumsqrt*(v_l*sqrtrold+sqrtrol*v_ld+v_r*sqrtrord+&
&     sqrtror*v_rd) + temp*invsumsqrtd
    v_tilda(2) = temp*invsumsqrt
    temp = sqrtrol*w_l + sqrtror*w_r
    v_tildad(3) = invsumsqrt*(w_l*sqrtrold+sqrtrol*w_ld+w_r*sqrtrord+&
&     sqrtror*w_rd) + temp*invsumsqrtd
    v_tilda(3) = temp*invsumsqrt
    temp = sqrtrol*h_l + sqrtror*h_r
    h_tildad = invsumsqrt*(h_l*sqrtrold+sqrtrol*h_ld+h_r*sqrtrord+&
&     sqrtror*h_rd) + temp*invsumsqrtd
    h_tilda = temp*invsumsqrt
    vn_tildad = n(1)*v_tildad(1) + n(2)*v_tildad(2) + n(dim)*(dim-2)*&
&     v_tildad(3)
    vn_tilda = v_tilda(1)*n(1) + v_tilda(2)*n(2) + v_tilda(3)*n(dim)*(&
&     dim-2)
    arg1d = (k-1.0d0)*(h_tildad-0.5d0*(2*v_tilda(1)*v_tildad(1)+2*&
&     v_tilda(2)*v_tildad(2)+(dim-2)*2*v_tilda(3)*v_tildad(3)))
    arg1 = (k-1.0d0)*(h_tilda-0.5d0*(v_tilda(1)**2+v_tilda(2)**2+v_tilda&
&     (3)**2*(dim-2)))
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_tildad = 0.0_8
    ELSE
      a_tildad = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_tilda = temp0
    IF (vn_l - a_l .GT. vn_tilda - a_tilda) THEN
      s_ld = vn_tildad - a_tildad
      s_l = vn_tilda - a_tilda
    ELSE
      s_ld = vn_ld - a_ld
      s_l = vn_l - a_l
    END IF
    IF (vn_r + a_r .LT. vn_tilda + a_tilda) THEN
      s_rd = vn_tildad + a_tildad
      s_r = vn_tilda + a_tilda
    ELSE
      s_rd = vn_rd + a_rd
      s_r = vn_r + a_r
    END IF
    
!=== LEFT STATE FLUXES: ===
    f_ld = 0.0_8
    f_ld(1) = vn_l*ro_ld + ro_l*vn_ld
    f_l(1) = ro_l*vn_l
    f_ld(2) = vn_l*(u_l*ro_ld+ro_l*u_ld) + ro_l*u_l*vn_ld + n(1)*p_ld
    f_l(2) = ro_l*u_l*vn_l + p_l*n(1)
    f_ld(3) = vn_l*(v_l*ro_ld+ro_l*v_ld) + ro_l*v_l*vn_ld + n(2)*p_ld
    f_l(3) = ro_l*v_l*vn_l + p_l*n(2)
    f_ld(4) = vn_l*(w_l*ro_ld+ro_l*w_ld) + ro_l*w_l*vn_ld + n(dim)*p_ld
    f_l(4) = ro_l*w_l*vn_l + p_l*n(dim)
    f_ld(dim+2) = vn_l*(h_l*ro_ld+ro_l*h_ld) + ro_l*h_l*vn_ld
    f_l(dim+2) = ro_l*h_l*vn_l
    w_vecld = 0.0_8
    w_vecld(1) = ro_ld
    w_vecl(1) = ro_l
    w_vecld(2) = u_l*ro_ld + ro_l*u_ld
    w_vecl(2) = ro_l*u_l
    w_vecld(3) = v_l*ro_ld + ro_l*v_ld
    w_vecl(3) = ro_l*v_l
    w_vecld(4) = w_l*ro_ld + ro_l*w_ld
    w_vecl(4) = ro_l*w_l
    w_vecld(dim+2) = e_l*ro_ld + ro_l*e_ld
    w_vecl(dim+2) = ro_l*e_l
    
!=== RIGHT STATE FLUXES: ===
    f_rd = 0.0_8
    f_rd(1) = vn_r*ro_rd + ro_r*vn_rd
    f_r(1) = ro_r*vn_r
    f_rd(2) = vn_r*(u_r*ro_rd+ro_r*u_rd) + ro_r*u_r*vn_rd + n(1)*p_rd
    f_r(2) = ro_r*u_r*vn_r + p_r*n(1)
    f_rd(3) = vn_r*(v_r*ro_rd+ro_r*v_rd) + ro_r*v_r*vn_rd + n(2)*p_rd
    f_r(3) = ro_r*v_r*vn_r + p_r*n(2)
    f_rd(4) = vn_r*(w_r*ro_rd+ro_r*w_rd) + ro_r*w_r*vn_rd + n(dim)*p_rd
    f_r(4) = ro_r*w_r*vn_r + p_r*n(dim)
    f_rd(dim+2) = vn_r*(h_r*ro_rd+ro_r*h_rd) + ro_r*h_r*vn_rd
    f_r(dim+2) = ro_r*h_r*vn_r
    w_vecrd = 0.0_8
    w_vecrd(1) = ro_rd
    w_vecr(1) = ro_r
    w_vecrd(2) = u_r*ro_rd + ro_r*u_rd
    w_vecr(2) = ro_r*u_r
    w_vecrd(3) = v_r*ro_rd + ro_r*v_rd
    w_vecr(3) = ro_r*v_r
    w_vecrd(4) = w_r*ro_rd + ro_r*w_rd
    w_vecr(4) = ro_r*w_r
    w_vecrd(dim+2) = e_r*ro_rd + ro_r*e_rd
    w_vecr(dim+2) = ro_r*e_r
    
    IF (0.0d0 .LT. s_r) THEN
      s_plusd = s_rd
      s_plus = s_r
    ELSE
      s_plus = 0.0d0
      s_plusd = 0.0_8
    END IF
    IF (0.0d0 .GT. s_l) THEN
      s_minusd = s_ld
      s_minus = s_l
    ELSE
      s_minus = 0.0d0
      s_minusd = 0.0_8
    END IF
    temp = 1.0/(s_plus-s_minus)
    invsumsd = -(temp*(s_plusd-s_minusd)/(s_plus-s_minus))
    invsums = temp
    temp1 = s_plus*f_l(1:dim+2) - s_minus*f_r(1:dim+2) + s_plus*s_minus*&
&     (w_vecr(1:dim+2)-w_vecl(1:dim+2))
    fluxesd(1:dim+2) = invsums*(f_l(1:dim+2)*s_plusd+s_plus*f_ld(1:dim+2&
&     )-f_r(1:dim+2)*s_minusd-s_minus*f_rd(1:dim+2)+(w_vecr(1:dim+2)-&
&     w_vecl(1:dim+2))*(s_minus*s_plusd+s_plus*s_minusd)+s_plus*s_minus*&
&     (w_vecrd(1:dim+2)-w_vecld(1:dim+2))) + temp1*invsumsd
  END SUBROUTINE

!  Differentiation of hllc_reiman_solver in forward (tangent) mode:
!   variations   of useful results: fluxes
!   with respect to varying inputs: w_l w_r t_l t_r v_l v_r p_l
!                u_l p_r u_r
!   RW status of diff variables: w_l:in w_r:in t_l:in t_r:in v_l:in
!                v_r:in p_l:in fluxes:out u_l:in p_r:in u_r:in
  PURE SUBROUTINE HLLC_REIMAN_SOLVER_D(dim, k, r0, cp, cv, p_l, p_ld, u_l, &
&   u_ld, v_l, v_ld, w_l, w_ld, t_l, t_ld, p_r, p_rd, u_r, u_rd, v_r, &
&   v_rd, w_r, w_rd, t_r, t_rd, n, fluxesd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: k, r0, cp, cv
    REAL*8, INTENT(IN) :: p_l, u_l, v_l, w_l, t_l
    REAL*8, INTENT(IN) :: p_ld, u_ld, v_ld, w_ld, t_ld
    REAL*8, INTENT(IN) :: p_r, u_r, v_r, w_r, t_r
    REAL*8, INTENT(IN) :: p_rd, u_rd, v_rd, w_rd, t_rd
    REAL*8, INTENT(IN) :: n(dim)
    REAL*8, INTENT(OUT) :: fluxesd(dim+2)
    REAL*8 :: ro_l, ro_r, a_l, a_r, h_l, h_r, e_l, e_r
    REAL*8 :: ro_ld, ro_rd, a_ld, a_rd, h_ld, h_rd, e_ld, e_rd
    REAL*8 :: vn_l, vn_r, v2_l, v2_r
    REAL*8 :: vn_ld, vn_rd, v2_ld, v2_rd
    REAL*8 :: s_l, s_r, s_star, sqrtrol, sqrtror, invsumsqrt
    REAL*8 :: s_ld, s_rd, s_stard, sqrtrold, sqrtrord, invsumsqrtd
    REAL*8 :: h_tilda, v_tilda(3), vn_tilda, a_tilda
    REAL*8 :: h_tildad, v_tildad(3), vn_tildad, a_tildad
    REAL*8 :: f_l(5), f_r(5), w_vecl(5), w_vecr(5), d_vec(5)
    REAL*8 :: f_ld(5), f_rd(5), w_vecld(5), w_vecrd(5), d_vecd(5)
    REAL*8 :: num, den
    REAL*8 :: numd, dend
    INTRINSIC DSQRT
    INTRINSIC MIN
    INTRINSIC MAX
    REAL*8 :: arg1
    REAL*8 :: arg1d
    REAL*8 :: temp
    DOUBLE PRECISION :: temp0
    REAL*8, DIMENSION(dim+2) :: temp1
    REAL*8, DIMENSION(dim+2) :: temp2
!=== LEFT STATE: ===
    temp = p_l/(r0*t_l)
    ro_ld = (p_ld-temp*r0*t_ld)/(r0*t_l)
    ro_l = temp
    arg1d = k*r0*t_ld
    arg1 = k*r0*t_l
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_ld = 0.0_8
    ELSE
      a_ld = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_l = temp0
    v2_ld = 2*u_l*u_ld + 2*v_l*v_ld + (dim-2)*2*w_l*w_ld
    v2_l = u_l**2 + v_l**2 + w_l**2*(dim-2)
    vn_ld = n(1)*u_ld + n(2)*v_ld + n(dim)*(dim-2)*w_ld
    vn_l = u_l*n(1) + v_l*n(2) + w_l*n(dim)*(dim-2)
    e_ld = cv*t_ld + 0.5d0*v2_ld
    e_l = cv*t_l + 0.5d0*v2_l
    h_ld = cp*t_ld + 0.5d0*v2_ld
    h_l = cp*t_l + 0.5d0*v2_l
!=== RIGHT STATE: ===
    temp = p_r/(r0*t_r)
    ro_rd = (p_rd-temp*r0*t_rd)/(r0*t_r)
    ro_r = temp
    arg1d = k*r0*t_rd
    arg1 = k*r0*t_r
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_rd = 0.0_8
    ELSE
      a_rd = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_r = temp0
    v2_rd = 2*u_r*u_rd + 2*v_r*v_rd + (dim-2)*2*w_r*w_rd
    v2_r = u_r**2 + v_r**2 + w_r**2*(dim-2)
    vn_rd = n(1)*u_rd + n(2)*v_rd + n(dim)*(dim-2)*w_rd
    vn_r = u_r*n(1) + v_r*n(2) + w_r*n(dim)*(dim-2)
    e_rd = cv*t_rd + 0.5d0*v2_rd
    e_r = cv*t_r + 0.5d0*v2_r
    h_rd = cp*t_rd + 0.5d0*v2_rd
    h_r = cp*t_r + 0.5d0*v2_r
!=== ROE AVERAGING: ===
    temp0 = DSQRT(ro_l)
    IF (ro_l .EQ. 0.0) THEN
      sqrtrold = 0.0_8
    ELSE
      sqrtrold = ro_ld/(2.D0*DSQRT(ro_l))
    END IF
    sqrtrol = temp0
    temp0 = DSQRT(ro_r)
    IF (ro_r .EQ. 0.0) THEN
      sqrtrord = 0.0_8
    ELSE
      sqrtrord = ro_rd/(2.D0*DSQRT(ro_r))
    END IF
    sqrtror = temp0
    temp = 1.0/(sqrtrol+sqrtror)
    invsumsqrtd = -(temp*(sqrtrold+sqrtrord)/(sqrtrol+sqrtror))
    invsumsqrt = temp
    v_tildad = 0.0_8
    temp = sqrtrol*u_l + sqrtror*u_r
    v_tildad(1) = invsumsqrt*(u_l*sqrtrold+sqrtrol*u_ld+u_r*sqrtrord+&
&     sqrtror*u_rd) + temp*invsumsqrtd
    v_tilda(1) = temp*invsumsqrt
    temp = sqrtrol*v_l + sqrtror*v_r
    v_tildad(2) = invsumsqrt*(v_l*sqrtrold+sqrtrol*v_ld+v_r*sqrtrord+&
&     sqrtror*v_rd) + temp*invsumsqrtd
    v_tilda(2) = temp*invsumsqrt
    temp = sqrtrol*w_l + sqrtror*w_r
    v_tildad(3) = invsumsqrt*(w_l*sqrtrold+sqrtrol*w_ld+w_r*sqrtrord+&
&     sqrtror*w_rd) + temp*invsumsqrtd
    v_tilda(3) = temp*invsumsqrt
    temp = sqrtrol*h_l + sqrtror*h_r
    h_tildad = invsumsqrt*(h_l*sqrtrold+sqrtrol*h_ld+h_r*sqrtrord+&
&     sqrtror*h_rd) + temp*invsumsqrtd
    h_tilda = temp*invsumsqrt
    vn_tildad = n(1)*v_tildad(1) + n(2)*v_tildad(2) + n(dim)*(dim-2)*&
&     v_tildad(3)
    vn_tilda = v_tilda(1)*n(1) + v_tilda(2)*n(2) + v_tilda(3)*n(dim)*(&
&     dim-2)
    arg1d = (k-1.0d0)*(h_tildad-0.5d0*(2*v_tilda(1)*v_tildad(1)+2*&
&     v_tilda(2)*v_tildad(2)+(dim-2)*2*v_tilda(3)*v_tildad(3)))
    arg1 = (k-1.0d0)*(h_tilda-0.5d0*(v_tilda(1)**2+v_tilda(2)**2+v_tilda&
&     (3)**2*(dim-2)))
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_tildad = 0.0_8
    ELSE
      a_tildad = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_tilda = temp0
    IF (vn_l - a_l .GT. vn_tilda - a_tilda) THEN
      s_ld = vn_tildad - a_tildad
      s_l = vn_tilda - a_tilda
    ELSE
      s_ld = vn_ld - a_ld
      s_l = vn_l - a_l
    END IF
    IF (vn_r + a_r .LT. vn_tilda + a_tilda) THEN
      s_rd = vn_tildad + a_tildad
      s_r = vn_tilda + a_tilda
    ELSE
      s_rd = vn_rd + a_rd
      s_r = vn_r + a_r
    END IF
!=== CONTACT DISCONTINUITY WAVE VELOCITY: ===
    numd = p_rd - p_ld + (s_l-vn_l)*(vn_l*ro_ld+ro_l*vn_ld) + ro_l*vn_l*&
&     (s_ld-vn_ld) - (s_r-vn_r)*(vn_r*ro_rd+ro_r*vn_rd) - ro_r*vn_r*(&
&     s_rd-vn_rd)
    num = p_r - p_l + ro_l*vn_l*(s_l-vn_l) - ro_r*vn_r*(s_r-vn_r)
    dend = (s_l-vn_l)*ro_ld + ro_l*(s_ld-vn_ld) - (s_r-vn_r)*ro_rd - &
&     ro_r*(s_rd-vn_rd)
    den = ro_l*(s_l-vn_l) - ro_r*(s_r-vn_r)
    s_stard = (numd-num*dend/den)/den
    s_star = num/den
!=== LEFT STATE FLUXES: ===
    f_ld = 0.0_8
    f_ld(1) = vn_l*ro_ld + ro_l*vn_ld
    f_l(1) = ro_l*vn_l
    f_ld(2) = vn_l*(u_l*ro_ld+ro_l*u_ld) + ro_l*u_l*vn_ld + n(1)*p_ld
    f_l(2) = ro_l*u_l*vn_l + p_l*n(1)
    f_ld(3) = vn_l*(v_l*ro_ld+ro_l*v_ld) + ro_l*v_l*vn_ld + n(2)*p_ld
    f_l(3) = ro_l*v_l*vn_l + p_l*n(2)
    f_ld(4) = vn_l*(w_l*ro_ld+ro_l*w_ld) + ro_l*w_l*vn_ld + n(dim)*p_ld
    f_l(4) = ro_l*w_l*vn_l + p_l*n(dim)
    f_ld(dim+2) = vn_l*(h_l*ro_ld+ro_l*h_ld) + ro_l*h_l*vn_ld
    f_l(dim+2) = ro_l*h_l*vn_l
    w_vecld = 0.0_8
    w_vecld(1) = ro_ld
    w_vecl(1) = ro_l
    w_vecld(2) = u_l*ro_ld + ro_l*u_ld
    w_vecl(2) = ro_l*u_l
    w_vecld(3) = v_l*ro_ld + ro_l*v_ld
    w_vecl(3) = ro_l*v_l
    w_vecld(4) = w_l*ro_ld + ro_l*w_ld
    w_vecl(4) = ro_l*w_l
    w_vecld(dim+2) = e_l*ro_ld + ro_l*e_ld
    w_vecl(dim+2) = ro_l*e_l
!=== RIGHT STATE FLUXES: ===
    f_rd = 0.0_8
    f_rd(1) = vn_r*ro_rd + ro_r*vn_rd
    f_r(1) = ro_r*vn_r
    f_rd(2) = vn_r*(u_r*ro_rd+ro_r*u_rd) + ro_r*u_r*vn_rd + n(1)*p_rd
    f_r(2) = ro_r*u_r*vn_r + p_r*n(1)
    f_rd(3) = vn_r*(v_r*ro_rd+ro_r*v_rd) + ro_r*v_r*vn_rd + n(2)*p_rd
    f_r(3) = ro_r*v_r*vn_r + p_r*n(2)
    f_rd(4) = vn_r*(w_r*ro_rd+ro_r*w_rd) + ro_r*w_r*vn_rd + n(dim)*p_rd
    f_r(4) = ro_r*w_r*vn_r + p_r*n(dim)
    f_rd(dim+2) = vn_r*(h_r*ro_rd+ro_r*h_rd) + ro_r*h_r*vn_rd
    f_r(dim+2) = ro_r*h_r*vn_r
    w_vecrd = 0.0_8
    w_vecrd(1) = ro_rd
    w_vecr(1) = ro_r
    w_vecrd(2) = u_r*ro_rd + ro_r*u_rd
    w_vecr(2) = ro_r*u_r
    w_vecrd(3) = v_r*ro_rd + ro_r*v_rd
    w_vecr(3) = ro_r*v_r
    w_vecrd(4) = w_r*ro_rd + ro_r*w_rd
    w_vecr(4) = ro_r*w_r
    w_vecrd(dim+2) = e_r*ro_rd + ro_r*e_rd
    w_vecr(dim+2) = ro_r*e_r
!=== FACE FLUXES: ===
    IF (s_l .GE. 0.0d0) THEN
      fluxesd(1:dim+2) = f_ld(1:dim+2)
    ELSE IF (s_star .GE. 0.0d0) THEN
      d_vec(1) = 0.d0
      d_vec(2:dim+1) = n(1:dim)
      d_vecd = 0.0_8
      d_vecd(dim+2) = s_stard
      d_vec(dim+2) = s_star
      temp1 = s_l*w_vecl(1:dim+2) - f_l(1:dim+2)
      temp = p_l + ro_l*(s_l-vn_l)*(s_star-vn_l)
      temp2 = (s_star*temp1+temp*s_l*d_vec(1:dim+2))/(s_l-s_star)
      fluxesd(1:dim+2) = (temp1*s_stard+s_star*(w_vecl(1:dim+2)*s_ld+s_l&
&       *w_vecld(1:dim+2)-f_ld(1:dim+2))+s_l*d_vec(1:dim+2)*(p_ld+(&
&       s_star-vn_l)*((s_l-vn_l)*ro_ld+ro_l*(s_ld-vn_ld))+ro_l*(s_l-vn_l&
&       )*(s_stard-vn_ld))+temp*(d_vec(1:dim+2)*s_ld+s_l*d_vecd(1:dim+2)&
&       )-temp2*(s_ld-s_stard))/(s_l-s_star)
    ELSE IF (s_r .GE. 0.0d0) THEN
      d_vec(1) = 0.d0
      d_vec(2:dim+1) = n(1:dim)
      d_vecd = 0.0_8
      d_vecd(dim+2) = s_stard
      d_vec(dim+2) = s_star
      temp2 = s_r*w_vecr(1:dim+2) - f_r(1:dim+2)
      temp = p_r + ro_r*(s_r-vn_r)*(s_star-vn_r)
      temp1 = (s_star*temp2+temp*s_r*d_vec(1:dim+2))/(s_r-s_star)
      fluxesd(1:dim+2) = (temp2*s_stard+s_star*(w_vecr(1:dim+2)*s_rd+s_r&
&       *w_vecrd(1:dim+2)-f_rd(1:dim+2))+s_r*d_vec(1:dim+2)*(p_rd+(&
&       s_star-vn_r)*((s_r-vn_r)*ro_rd+ro_r*(s_rd-vn_rd))+ro_r*(s_r-vn_r&
&       )*(s_stard-vn_rd))+temp*(d_vec(1:dim+2)*s_rd+s_r*d_vecd(1:dim+2)&
&       )-temp1*(s_rd-s_stard))/(s_r-s_star)
    ELSE
      fluxesd(1:dim+2) = f_rd(1:dim+2)
    END IF
  END SUBROUTINE

!  Differentiation of roe_reiman_solver in forward (tangent) mode:
!   variations   of useful results: fluxes
!   with respect to varying inputs: w_l w_r t_l t_r v_l v_r p_l
!                u_l p_r u_r
!   RW status of diff variables: w_l:in w_r:in t_l:in t_r:in v_l:in
!                v_r:in p_l:in fluxes:out u_l:in p_r:in u_r:in
  PURE SUBROUTINE ROE_REIMAN_SOLVER_D(dim, k, r0, cp, cv, p_l, p_ld, u_l, &
&   u_ld, v_l, v_ld, w_l, w_ld, t_l, t_ld, p_r, p_rd, u_r, u_rd, v_r, &
&   v_rd, w_r, w_rd, t_r, t_rd, n, fluxesd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: k, r0, cp, cv
    REAL*8, INTENT(IN) :: p_l, u_l, v_l, w_l, t_l
    REAL*8, INTENT(IN) :: p_ld, u_ld, v_ld, w_ld, t_ld
    REAL*8, INTENT(IN) :: p_r, u_r, v_r, w_r, t_r
    REAL*8, INTENT(IN) :: p_rd, u_rd, v_rd, w_rd, t_rd
    REAL*8, INTENT(IN) :: n(dim)
    REAL*8, INTENT(OUT) :: fluxesd(dim+2)
    REAL*8 :: ro_l, ro_r, a_l, a_r, h_l, h_r, v2_l, v2_r, vn_l, vn_r
    REAL*8 :: ro_ld, ro_rd, h_ld, h_rd, v2_ld, v2_rd, vn_ld, vn_rd
    REAL*8 :: ro_tilda, v_tilda(3), h_tilda, a_tilda, v2_tilda, vn_tilda
    REAL*8 :: ro_tildad, v_tildad(3), h_tildad, a_tildad, v2_tildad, &
&   vn_tildad
    REAL*8 :: dp, dvn, dro, dv(3)
    REAL*8 :: dpd, dvnd, drod, dvd(3)
    REAL*8 :: lambda(dim+2), alpha(dim+2)
    REAL*8 :: lambdad(dim+2), alphad(dim+2)
    REAL*8 :: f_l(5), f_r(5), diss(5)
    REAL*8 :: f_ld(5), f_rd(5), dissd(5)
    REAL*8 :: sqrtrol, sqrtror, invsumsqrt, gam1
    REAL*8 :: sqrtrold, sqrtrord, invsumsqrtd
    INTEGER :: d
    INTRINSIC DSQRT
    INTRINSIC ABS
    REAL*8 :: arg1
    REAL*8 :: arg1d
    REAL*8 :: temp
    DOUBLE PRECISION :: temp0
    REAL*8, DIMENSION(dim) :: temp1
    REAL*8 :: temp2
    gam1 = k - 1.0d0
!=== LEFT STATE: ===
    temp = p_l/(r0*t_l)
    ro_ld = (p_ld-temp*r0*t_ld)/(r0*t_l)
    ro_l = temp
    v2_ld = 2*u_l*u_ld + 2*v_l*v_ld + (dim-2)*2*w_l*w_ld
    v2_l = u_l**2 + v_l**2 + w_l**2*(dim-2)
    vn_ld = n(1)*u_ld + n(2)*v_ld + n(dim)*(dim-2)*w_ld
    vn_l = u_l*n(1) + v_l*n(2) + w_l*n(dim)*(dim-2)
    h_ld = cp*t_ld + 0.5d0*v2_ld
    h_l = cp*t_l + 0.5d0*v2_l
!=== RIGHT STATE: ===
    temp = p_r/(r0*t_r)
    ro_rd = (p_rd-temp*r0*t_rd)/(r0*t_r)
    ro_r = temp
    v2_rd = 2*u_r*u_rd + 2*v_r*v_rd + (dim-2)*2*w_r*w_rd
    v2_r = u_r**2 + v_r**2 + w_r**2*(dim-2)
    vn_rd = n(1)*u_rd + n(2)*v_rd + n(dim)*(dim-2)*w_rd
    vn_r = u_r*n(1) + v_r*n(2) + w_r*n(dim)*(dim-2)
    h_rd = cp*t_rd + 0.5d0*v2_rd
    h_r = cp*t_r + 0.5d0*v2_r
!=== ROE AVERAGING: ===
    temp0 = DSQRT(ro_l)
    IF (ro_l .EQ. 0.0) THEN
      sqrtrold = 0.0_8
    ELSE
      sqrtrold = ro_ld/(2.D0*DSQRT(ro_l))
    END IF
    sqrtrol = temp0
    temp0 = DSQRT(ro_r)
    IF (ro_r .EQ. 0.0) THEN
      sqrtrord = 0.0_8
    ELSE
      sqrtrord = ro_rd/(2.D0*DSQRT(ro_r))
    END IF
    sqrtror = temp0
    temp = 1.0/(sqrtrol+sqrtror)
    invsumsqrtd = -(temp*(sqrtrold+sqrtrord)/(sqrtrol+sqrtror))
    invsumsqrt = temp
    temp0 = DSQRT(ro_l*ro_r)
    IF (ro_l*ro_r .EQ. 0.0) THEN
      ro_tildad = 0.0_8
    ELSE
      ro_tildad = (ro_r*ro_ld+ro_l*ro_rd)/(2.D0*DSQRT(ro_l*ro_r))
    END IF
    ro_tilda = temp0
    v_tildad = 0.0_8
    temp = sqrtrol*u_l + sqrtror*u_r
    v_tildad(1) = invsumsqrt*(u_l*sqrtrold+sqrtrol*u_ld+u_r*sqrtrord+&
&     sqrtror*u_rd) + temp*invsumsqrtd
    v_tilda(1) = temp*invsumsqrt
    temp = sqrtrol*v_l + sqrtror*v_r
    v_tildad(2) = invsumsqrt*(v_l*sqrtrold+sqrtrol*v_ld+v_r*sqrtrord+&
&     sqrtror*v_rd) + temp*invsumsqrtd
    v_tilda(2) = temp*invsumsqrt
    temp = sqrtrol*w_l + sqrtror*w_r
    v_tildad(3) = invsumsqrt*(w_l*sqrtrold+sqrtrol*w_ld+w_r*sqrtrord+&
&     sqrtror*w_rd) + temp*invsumsqrtd
    v_tilda(3) = temp*invsumsqrt
    v2_tildad = 2*v_tilda(1)*v_tildad(1) + 2*v_tilda(2)*v_tildad(2) + (&
&     dim-2)*2*v_tilda(3)*v_tildad(3)
    v2_tilda = v_tilda(1)**2 + v_tilda(2)**2 + v_tilda(3)**2*(dim-2)
    vn_tildad = n(1)*v_tildad(1) + n(2)*v_tildad(2) + n(dim)*(dim-2)*&
&     v_tildad(3)
    vn_tilda = v_tilda(1)*n(1) + v_tilda(2)*n(2) + v_tilda(3)*n(dim)*(&
&     dim-2)
    temp = sqrtrol*h_l + sqrtror*h_r
    h_tildad = invsumsqrt*(h_l*sqrtrold+sqrtrol*h_ld+h_r*sqrtrord+&
&     sqrtror*h_rd) + temp*invsumsqrtd
    h_tilda = temp*invsumsqrt
    arg1d = gam1*(h_tildad-0.5d0*v2_tildad)
    arg1 = gam1*(h_tilda-0.5d0*v2_tilda)
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_tildad = 0.0_8
    ELSE
      a_tildad = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_tilda = temp0
!=== VARIABLES DIFFERENCES: ===
    dpd = p_rd - p_ld
    dp = p_r - p_l
    drod = ro_rd - ro_ld
    dro = ro_r - ro_l
    dvd = 0.0_8
    dvd(1) = u_rd - u_ld
    dv(1) = u_r - u_l
    dvd(2) = v_rd - v_ld
    dv(2) = v_r - v_l
    dvd(3) = w_rd - w_ld
    dv(3) = w_r - w_l
    dvnd = n(1)*dvd(1) + n(2)*dvd(2) + n(dim)*(dim-2)*dvd(3)
    dvn = dv(1)*n(1) + dv(2)*n(2) + dv(3)*n(dim)*(dim-2)
    IF (vn_tilda - a_tilda .GE. 0.) THEN
      lambdad = 0.0_8
      lambdad(1) = vn_tildad - a_tildad
      lambda(1) = vn_tilda - a_tilda
    ELSE
      lambdad = 0.0_8
      lambdad(1) = a_tildad - vn_tildad
      lambda(1) = -(vn_tilda-a_tilda)
    END IF
    IF (vn_tilda .GE. 0.) THEN
      lambdad(2:dim+1) = vn_tildad
      lambda(2:dim+1) = vn_tilda
    ELSE
      lambdad(2:dim+1) = -vn_tildad
      lambda(2:dim+1) = -vn_tilda
    END IF
    IF (vn_tilda + a_tilda .GE. 0.) THEN
      lambdad(dim+2) = vn_tildad + a_tildad
      lambda(dim+2) = vn_tilda + a_tilda
    ELSE
      lambdad(dim+2) = -vn_tildad - a_tildad
      lambda(dim+2) = -(vn_tilda+a_tilda)
    END IF
!=== WAVE AMPLITUDES: ===
    alphad = 0.0_8
    temp0 = 2.0d0*(a_tilda*a_tilda)
    temp = (dp-ro_tilda*a_tilda*dvn)/temp0
    alphad(1) = (dpd-dvn*(a_tilda*ro_tildad+ro_tilda*a_tildad)-ro_tilda*&
&     a_tilda*dvnd-temp*2.0d0*2*a_tilda*a_tildad)/temp0
    alpha(1) = temp
    temp0 = 2.0d0*(a_tilda*a_tilda)
    temp = (dp+ro_tilda*a_tilda*dvn)/temp0
    alphad(dim+2) = (dpd+dvn*(a_tilda*ro_tildad+ro_tilda*a_tildad)+&
&     ro_tilda*a_tilda*dvnd-temp*2.0d0*2*a_tilda*a_tildad)/temp0
    alpha(dim+2) = temp
    temp = dp/(a_tilda*a_tilda)
    alphad(2) = drod - (dpd-temp*2*a_tilda*a_tildad)/a_tilda**2
    alpha(2) = dro - temp
    diss = 0.0d0
!=== WAVE 1: V_n - A: ===
    dissd = 0.0_8
    dissd(1) = alpha(1)*lambdad(1) + lambda(1)*alphad(1)
    diss(1) = lambda(1)*alpha(1)
    temp1 = v_tilda(1:dim) - n(1:dim)*a_tilda
    dissd(2:dim+1) = temp1*(alpha(1)*lambdad(1)+lambda(1)*alphad(1)) + &
&     lambda(1)*alpha(1)*(v_tildad(1:dim)-n(1:dim)*a_tildad)
    diss(2:dim+1) = lambda(1)*alpha(1)*temp1
    dissd(dim+2) = (h_tilda-a_tilda*vn_tilda)*(alpha(1)*lambdad(1)+&
&     lambda(1)*alphad(1)) + lambda(1)*alpha(1)*(h_tildad-vn_tilda*&
&     a_tildad-a_tilda*vn_tildad)
    diss(dim+2) = lambda(1)*alpha(1)*(h_tilda-a_tilda*vn_tilda)
!=== WAVE 2: V_n: ===
    dissd(1) = dissd(1) + alpha(2)*lambdad(2) + lambda(2)*alphad(2)
    diss(1) = diss(1) + lambda(2)*alpha(2)
    dissd(2:dim+1) = dissd(2:dim+1) + v_tilda(1:dim)*(alpha(2)*lambdad(2&
&     )+lambda(2)*alphad(2)) + lambda(2)*alpha(2)*v_tildad(1:dim)
    diss(2:dim+1) = diss(2:dim+1) + lambda(2)*alpha(2)*v_tilda(1:dim)
    dissd(dim+2) = dissd(dim+2) + 0.5d0*(alpha(2)*(v2_tilda*lambdad(2)+&
&     lambda(2)*v2_tildad)+lambda(2)*v2_tilda*alphad(2))
    diss(dim+2) = diss(dim+2) + lambda(2)*alpha(2)*0.5d0*v2_tilda
!=== WAVE 3, 4: V_n: ===
    DO d=1,dim
      temp = dv(d) - n(d)*dvn
      temp2 = lambda(1+d)*ro_tilda
      dissd(1+d) = dissd(1+d) + temp*(ro_tilda*lambdad(1+d)+lambda(1+d)*&
&       ro_tildad) + temp2*(dvd(d)-n(d)*dvnd)
      diss(1+d) = diss(1+d) + temp2*temp
      temp2 = lambda(1+d)*dv(d)
      dissd(dim+2) = dissd(dim+2) + ro_tilda*v_tilda(d)*(dv(d)*lambdad(1&
&       +d)+lambda(1+d)*dvd(d)) + temp2*(v_tilda(d)*ro_tildad+ro_tilda*&
&       v_tildad(d))
      diss(dim+2) = diss(dim+2) + temp2*(ro_tilda*v_tilda(d))
    END DO
    dissd(dim+2) = dissd(dim+2) - ro_tilda*vn_tilda*(dvn*lambdad(2)+&
&     lambda(2)*dvnd) - lambda(2)*dvn*(vn_tilda*ro_tildad+ro_tilda*&
&     vn_tildad)
    diss(dim+2) = diss(dim+2) - lambda(2)*ro_tilda*(vn_tilda*dvn)
!=== WAVE 5: V_n + A: ===
    dissd(1) = dissd(1) + alpha(dim+2)*lambdad(dim+2) + lambda(dim+2)*&
&     alphad(dim+2)
    diss(1) = diss(1) + lambda(dim+2)*alpha(dim+2)
    temp1 = v_tilda(1:dim) + n(1:dim)*a_tilda
    temp2 = lambda(dim+2)*alpha(dim+2)
    dissd(2:dim+1) = dissd(2:dim+1) + temp1*(alpha(dim+2)*lambdad(dim+2)&
&     +lambda(dim+2)*alphad(dim+2)) + temp2*(v_tildad(1:dim)+n(1:dim)*&
&     a_tildad)
    diss(2:dim+1) = diss(2:dim+1) + temp2*temp1
    temp2 = lambda(dim+2)*alpha(dim+2)
    dissd(dim+2) = dissd(dim+2) + (h_tilda+a_tilda*vn_tilda)*(alpha(dim+&
&     2)*lambdad(dim+2)+lambda(dim+2)*alphad(dim+2)) + temp2*(h_tildad+&
&     vn_tilda*a_tildad+a_tilda*vn_tildad)
    diss(dim+2) = diss(dim+2) + temp2*(h_tilda+a_tilda*vn_tilda)
!=== LEFT/RIGHT FLUXES: ===
    f_ld = 0.0_8
    f_ld(1) = vn_l*ro_ld + ro_l*vn_ld
    f_l(1) = ro_l*vn_l
    f_ld(2) = vn_l*(u_l*ro_ld+ro_l*u_ld) + ro_l*u_l*vn_ld + n(1)*p_ld
    f_l(2) = ro_l*u_l*vn_l + p_l*n(1)
    f_ld(3) = vn_l*(v_l*ro_ld+ro_l*v_ld) + ro_l*v_l*vn_ld + n(2)*p_ld
    f_l(3) = ro_l*v_l*vn_l + p_l*n(2)
    f_ld(4) = vn_l*(w_l*ro_ld+ro_l*w_ld) + ro_l*w_l*vn_ld + n(dim)*p_ld
    f_l(4) = ro_l*w_l*vn_l + p_l*n(dim)
    f_ld(dim+2) = vn_l*(h_l*ro_ld+ro_l*h_ld) + ro_l*h_l*vn_ld
    f_l(dim+2) = ro_l*h_l*vn_l
    f_rd = 0.0_8
    f_rd(1) = vn_r*ro_rd + ro_r*vn_rd
    f_r(1) = ro_r*vn_r
    f_rd(2) = vn_r*(u_r*ro_rd+ro_r*u_rd) + ro_r*u_r*vn_rd + n(1)*p_rd
    f_r(2) = ro_r*u_r*vn_r + p_r*n(1)
    f_rd(3) = vn_r*(v_r*ro_rd+ro_r*v_rd) + ro_r*v_r*vn_rd + n(2)*p_rd
    f_r(3) = ro_r*v_r*vn_r + p_r*n(2)
    f_rd(4) = vn_r*(w_r*ro_rd+ro_r*w_rd) + ro_r*w_r*vn_rd + n(dim)*p_rd
    f_r(4) = ro_r*w_r*vn_r + p_r*n(dim)
    f_rd(dim+2) = vn_r*(h_r*ro_rd+ro_r*h_rd) + ro_r*h_r*vn_rd
    f_r(dim+2) = ro_r*h_r*vn_r
!=== FACES FLUXES: ===
    fluxesd(1:dim+2) = 0.5d0*(f_ld(1:dim+2)+f_rd(1:dim+2)-dissd(1:dim+2)&
&     )
  END SUBROUTINE

!  Differentiation of ausm_reiman_solver in forward (tangent) mode:
!   variations   of useful results: fluxes
!   with respect to varying inputs: w_l w_r t_l t_r v_l v_r p_l
!                u_l p_r u_r
!   RW status of diff variables: w_l:in w_r:in t_l:in t_r:in v_l:in
!                v_r:in p_l:in fluxes:out u_l:in p_r:in u_r:in
  SUBROUTINE AUSM_REIMAN_SOLVER_D(dim, k, r0, cp, cv, p_l, p_ld, u_l, &
&   u_ld, v_l, v_ld, w_l, w_ld, t_l, t_ld, p_r, p_rd, u_r, u_rd, v_r, &
&   v_rd, w_r, w_rd, t_r, t_rd, n, fluxesd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: k, r0, cp, cv
    REAL*8, INTENT(IN) :: p_l, u_l, v_l, w_l, t_l
    REAL*8, INTENT(IN) :: p_ld, u_ld, v_ld, w_ld, t_ld
    REAL*8, INTENT(IN) :: p_r, u_r, v_r, w_r, t_r
    REAL*8, INTENT(IN) :: p_rd, u_rd, v_rd, w_rd, t_rd
    REAL*8, INTENT(IN) :: n(dim)
    REAL*8, INTENT(OUT) :: fluxesd(dim+2)
    REAL*8 :: ro_l, ro_r, a_l, a_r, h_l, h_r, v2_l, v2_r
    REAL*8 :: ro_ld, ro_rd, a_ld, a_rd, h_ld, h_rd, v2_ld, v2_rd
    REAL*8 :: m_l, m_r, m_f, p_f, a_f
    REAL*8 :: m_ld, m_rd, m_fd, p_fd, a_fd
    REAL*8 :: fi_l(5), fi_r(5), f_p(5)
    REAL*8 :: fi_ld(5), fi_rd(5), f_pd(5)
    REAL*8 :: vn_l, vn_r
    REAL*8 :: vn_ld, vn_rd
    REAL*8 :: wl, wr
    INTRINSIC DSQRT
    INTRINSIC SIGN
    REAL*8 :: arg1
    REAL*8 :: arg1d
    REAL*8 :: result1
    REAL*8 :: result1d
    REAL*8 :: result2
    REAL*8 :: result2d
    REAL*8 :: temp
    DOUBLE PRECISION :: temp0
    REAL*8, DIMENSION(dim+2) :: temp1
!=== LEFT STATE: ===
    temp = p_l/(r0*t_l)
    ro_ld = (p_ld-temp*r0*t_ld)/(r0*t_l)
    ro_l = temp
    arg1d = k*r0*t_ld
    arg1 = k*r0*t_l
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_ld = 0.0_8
    ELSE
      a_ld = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_l = temp0
    v2_ld = 2*u_l*u_ld + 2*v_l*v_ld + (dim-2)*2*w_l*w_ld
    v2_l = u_l**2 + v_l**2 + w_l**2*(dim-2)
    vn_ld = n(1)*u_ld + n(2)*v_ld + n(dim)*(dim-2)*w_ld
    vn_l = u_l*n(1) + v_l*n(2) + w_l*n(dim)*(dim-2)
    h_ld = cp*t_ld + 0.5d0*v2_ld
    h_l = cp*t_l + 0.5d0*v2_l
!=== RIGHT STATE: ===
    temp = p_r/(r0*t_r)
    ro_rd = (p_rd-temp*r0*t_rd)/(r0*t_r)
    ro_r = temp
    arg1d = k*r0*t_rd
    arg1 = k*r0*t_r
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_rd = 0.0_8
    ELSE
      a_rd = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_r = temp0
    v2_rd = 2*u_r*u_rd + 2*v_r*v_rd + (dim-2)*2*w_r*w_rd
    v2_r = u_r**2 + v_r**2 + w_r**2*(dim-2)
    vn_rd = n(1)*u_rd + n(2)*v_rd + n(dim)*(dim-2)*w_rd
    vn_r = u_r*n(1) + v_r*n(2) + w_r*n(dim)*(dim-2)
    h_rd = cp*t_rd + 0.5d0*v2_rd
    h_r = cp*t_r + 0.5d0*v2_r
!=== LEFT/RIGHT MACH NUMBERS: ===
    a_fd = 0.5d0*(a_ld+a_rd)
    a_f = 0.5d0*(a_l+a_r)
    m_ld = (vn_ld-vn_l*a_fd/a_f)/a_f
    m_l = vn_l/a_f
    m_rd = (vn_rd-vn_r*a_fd/a_f)/a_f
    m_r = vn_r/a_f
!=== CONVECTION-PRESSURE SPLITTING: ===
    result1d = M_PLUS_D(m_l, m_ld, result1)
    result2d = M_MINUS_D(m_r, m_rd, result2)
    m_fd = result1d + result2d
    m_f = result1 + result2
    result1d = P_PLUS_D(m_l, m_ld, result1)
    result2d = P_MINUS_D(m_r, m_rd, result2)
    p_fd = p_l*result1d + result1*p_ld + p_r*result2d + result2*p_rd
    p_f = result1*p_l + result2*p_r
!=== LEFT/RIGHT STATE CONVECTION VECTORS: ===
    fi_l(1) = 1.0d0
    fi_ld = 0.0_8
    fi_ld(2) = u_ld
    fi_l(2) = u_l
    fi_ld(3) = v_ld
    fi_l(3) = v_l
    fi_ld(4) = w_ld
    fi_l(4) = w_l
    fi_ld(dim+2) = h_ld
    fi_l(dim+2) = h_l
    fi_r(1) = 1.0d0
    fi_rd = 0.0_8
    fi_rd(2) = u_rd
    fi_r(2) = u_r
    fi_rd(3) = v_rd
    fi_r(3) = v_r
    fi_rd(4) = w_rd
    fi_r(4) = w_r
    fi_rd(dim+2) = h_rd
    fi_r(dim+2) = h_r
!=== PRESSURE PART: ===
    f_p = 0.0d0
    f_pd = 0.0_8
    f_pd(2:dim+1) = n(1:dim)*p_fd
    f_p(2:dim+1) = p_f*n(1:dim)
    wl = 0.5d0*(1.0d0+SIGN(1.0d0, m_f))
    wr = 1.0d0 - wl
    temp1 = wl*ro_l*fi_l(1:dim+2) + wr*ro_r*fi_r(1:dim+2)
    fluxesd(1:dim+2) = temp1*(m_f*a_fd+a_f*m_fd) + a_f*m_f*(wl*(fi_l(1:&
&     dim+2)*ro_ld+ro_l*fi_ld(1:dim+2))+wr*(fi_r(1:dim+2)*ro_rd+ro_r*&
&     fi_rd(1:dim+2))) + f_pd(1:dim+2)

  CONTAINS
!  Differentiation of m_plus in forward (tangent) mode:
!   variations   of useful results: m_plus
!   with respect to varying inputs: m
    REAL*8 FUNCTION M_PLUS_D(m, md, m_plus)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      REAL*8, INTENT(IN) :: md
      INTRINSIC ABS
      REAL*8 :: abs0
      REAL*8 :: abs1
      REAL*8 :: abs1d
      REAL*8 :: m_plus
      IF (m .GE. 0.) THEN
        abs0 = m
      ELSE
        abs0 = -m
      END IF
      IF (abs0 .GE. 1.0d0) THEN
        IF (m .GE. 0.) THEN
          abs1d = md
          abs1 = m
        ELSE
          abs1d = -md
          abs1 = -m
        END IF
        m_plus_d = 0.5d0*(md+abs1d)
        m_plus = 0.5d0*(m+abs1)
      ELSE
        m_plus_d = 0.25d0*2*(m+1.0d0)*md
        m_plus = 0.25d0*(m+1.0d0)**2
      END IF
    END FUNCTION M_PLUS_D

    PURE REAL*8 FUNCTION M_PLUS(m)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      INTRINSIC ABS
      REAL*8 :: abs0
      REAL*8 :: abs1
      IF (m .GE. 0.) THEN
        abs0 = m
      ELSE
        abs0 = -m
      END IF
      IF (abs0 .GE. 1.0d0) THEN
        IF (m .GE. 0.) THEN
          abs1 = m
        ELSE
          abs1 = -m
        END IF
        m_plus = 0.5d0*(m+abs1)
      ELSE
        m_plus = 0.25d0*(m+1.0d0)**2
      END IF
    END FUNCTION M_PLUS

!  Differentiation of m_minus in forward (tangent) mode:
!   variations   of useful results: m_minus
!   with respect to varying inputs: m
    REAL*8 FUNCTION M_MINUS_D(m, md, m_minus)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      REAL*8, INTENT(IN) :: md
      INTRINSIC ABS
      REAL*8 :: abs0
      REAL*8 :: abs1
      REAL*8 :: abs1d
      REAL*8 :: m_minus
      IF (m .GE. 0.) THEN
        abs0 = m
      ELSE
        abs0 = -m
      END IF
      IF (abs0 .GE. 1.0d0) THEN
        IF (m .GE. 0.) THEN
          abs1d = md
          abs1 = m
        ELSE
          abs1d = -md
          abs1 = -m
        END IF
        m_minus_d = 0.5d0*(md-abs1d)
        m_minus = 0.5d0*(m-abs1)
      ELSE
        m_minus_d = -(0.25d0*2*(m-1.0d0)*md)
        m_minus = -(0.25d0*(m-1.0d0)**2)
      END IF
    END FUNCTION M_MINUS_D

    PURE REAL*8 FUNCTION M_MINUS(m)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      INTRINSIC ABS
      REAL*8 :: abs0
      REAL*8 :: abs1
      IF (m .GE. 0.) THEN
        abs0 = m
      ELSE
        abs0 = -m
      END IF
      IF (abs0 .GE. 1.0d0) THEN
        IF (m .GE. 0.) THEN
          abs1 = m
        ELSE
          abs1 = -m
        END IF
        m_minus = 0.5d0*(m-abs1)
      ELSE
        m_minus = -(0.25d0*(m-1.0d0)**2)
      END IF
    END FUNCTION M_MINUS

!  Differentiation of p_plus in forward (tangent) mode:
!   variations   of useful results: p_plus
!   with respect to varying inputs: m
    REAL*8 FUNCTION P_PLUS_D(m, md, p_plus)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      REAL*8, INTENT(IN) :: md
      INTRINSIC ABS
      INTRINSIC SIGN
      REAL*8 :: abs0
      REAL*8 :: p_plus
      IF (m .GE. 0.) THEN
        abs0 = m
      ELSE
        abs0 = -m
      END IF
      IF (abs0 .GE. 1.0d0) THEN
        p_plus = 0.5d0*(1.0d0+SIGN(1.0d0, m))
        p_plus_d = 0.0_8
      ELSE
        p_plus_d = 0.25d0*((2.0d0-m)*2*(m+1.0d0)-(m+1.0d0)**2)*md
        p_plus = 0.25d0*(m+1.0d0)**2*(2.0d0-m)
      END IF
    END FUNCTION P_PLUS_D

    PURE REAL*8 FUNCTION P_PLUS(m)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      INTRINSIC ABS
      INTRINSIC SIGN
      REAL*8 :: abs0
      IF (m .GE. 0.) THEN
        abs0 = m
      ELSE
        abs0 = -m
      END IF
      IF (abs0 .GE. 1.0d0) THEN
        p_plus = 0.5d0*(1.0d0+SIGN(1.0d0, m))
      ELSE
        p_plus = 0.25d0*(m+1.0d0)**2*(2.0d0-m)
      END IF
    END FUNCTION P_PLUS

!  Differentiation of p_minus in forward (tangent) mode:
!   variations   of useful results: p_minus
!   with respect to varying inputs: m
    REAL*8 FUNCTION P_MINUS_D(m, md, p_minus)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      REAL*8, INTENT(IN) :: md
      INTRINSIC ABS
      INTRINSIC SIGN
      REAL*8 :: abs0
      REAL*8 :: p_minus
      IF (m .GE. 0.) THEN
        abs0 = m
      ELSE
        abs0 = -m
      END IF
      IF (abs0 .GE. 1.0d0) THEN
        p_minus = 0.5d0*(1.0d0-SIGN(1.0d0, m))
        p_minus_d = 0.0_8
      ELSE
        p_minus_d = 0.25d0*((m+2.0d0)*2*(m-1.0d0)+(m-1.0d0)**2)*md
        p_minus = 0.25d0*(m-1.0d0)**2*(2.0d0+m)
      END IF
    END FUNCTION P_MINUS_D

    PURE REAL*8 FUNCTION P_MINUS(m)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      INTRINSIC ABS
      INTRINSIC SIGN
      REAL*8 :: abs0
      IF (m .GE. 0.) THEN
        abs0 = m
      ELSE
        abs0 = -m
      END IF
      IF (abs0 .GE. 1.0d0) THEN
        p_minus = 0.5d0*(1.0d0-SIGN(1.0d0, m))
      ELSE
        p_minus = 0.25d0*(m-1.0d0)**2*(2.0d0+m)
      END IF
    END FUNCTION P_MINUS

  END SUBROUTINE

!  Differentiation of ausmplus_reiman_solver in forward (tangent) mode:
!   variations   of useful results: fluxes
!   with respect to varying inputs: w_l w_r t_l t_r v_l v_r p_l
!                u_l p_r u_r
!   RW status of diff variables: w_l:in w_r:in t_l:in t_r:in v_l:in
!                v_r:in p_l:in fluxes:out u_l:in p_r:in u_r:in
  SUBROUTINE AUSMPLUS_REIMAN_SOLVER_D(dim, k, r0, cp, cv, p_l, p_ld, u_l&
&   , u_ld, v_l, v_ld, w_l, w_ld, t_l, t_ld, p_r, p_rd, u_r, u_rd, v_r, &
&   v_rd, w_r, w_rd, t_r, t_rd, n, fluxesd)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: dim
    REAL*8, INTENT(IN) :: k, r0, cp, cv
    REAL*8, INTENT(IN) :: p_l, u_l, v_l, w_l, t_l
    REAL*8, INTENT(IN) :: p_ld, u_ld, v_ld, w_ld, t_ld
    REAL*8, INTENT(IN) :: p_r, u_r, v_r, w_r, t_r
    REAL*8, INTENT(IN) :: p_rd, u_rd, v_rd, w_rd, t_rd
    REAL*8, INTENT(IN) :: n(dim)
    REAL*8, INTENT(OUT) :: fluxesd(dim+2)
    REAL*8 :: ro_l, ro_r, a_l, a_r, h_l, h_r, v2_l, v2_r
    REAL*8 :: ro_ld, ro_rd, a_ld, a_rd, h_ld, h_rd, v2_ld, v2_rd
    REAL*8 :: m_l, m_r, m_f, p_f, a_f
    REAL*8 :: m_ld, m_rd, m_fd, p_fd, a_fd
    REAL*8 :: fi_l(5), fi_r(5), f_p(5)
    REAL*8 :: fi_ld(5), fi_rd(5), f_pd(5)
    REAL*8 :: vn_l, vn_r
    REAL*8 :: vn_ld, vn_rd
    REAL*8 :: wl, wr
    INTRINSIC DSQRT
    INTRINSIC SIGN
    REAL*8 :: arg1
    REAL*8 :: arg1d
    REAL*8 :: result1
    REAL*8 :: result1d
    REAL*8 :: result2
    REAL*8 :: result2d
    REAL*8 :: temp
    DOUBLE PRECISION :: temp0
    REAL*8, DIMENSION(dim+2) :: temp1
!=== LEFT STATE: ===
    temp = p_l/(r0*t_l)
    ro_ld = (p_ld-temp*r0*t_ld)/(r0*t_l)
    ro_l = temp
    arg1d = k*r0*t_ld
    arg1 = k*r0*t_l
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_ld = 0.0_8
    ELSE
      a_ld = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_l = temp0
    v2_ld = 2*u_l*u_ld + 2*v_l*v_ld + (dim-2)*2*w_l*w_ld
    v2_l = u_l**2 + v_l**2 + w_l**2*(dim-2)
    vn_ld = n(1)*u_ld + n(2)*v_ld + n(dim)*(dim-2)*w_ld
    vn_l = u_l*n(1) + v_l*n(2) + w_l*n(dim)*(dim-2)
    h_ld = cp*t_ld + 0.5d0*v2_ld
    h_l = cp*t_l + 0.5d0*v2_l
!=== RIGHT STATE: ===
    temp = p_r/(r0*t_r)
    ro_rd = (p_rd-temp*r0*t_rd)/(r0*t_r)
    ro_r = temp
    arg1d = k*r0*t_rd
    arg1 = k*r0*t_r
    temp0 = DSQRT(arg1)
    IF (arg1 .EQ. 0.0) THEN
      a_rd = 0.0_8
    ELSE
      a_rd = arg1d/(2.D0*DSQRT(arg1))
    END IF
    a_r = temp0
    v2_rd = 2*u_r*u_rd + 2*v_r*v_rd + (dim-2)*2*w_r*w_rd
    v2_r = u_r**2 + v_r**2 + w_r**2*(dim-2)
    vn_rd = n(1)*u_rd + n(2)*v_rd + n(dim)*(dim-2)*w_rd
    vn_r = u_r*n(1) + v_r*n(2) + w_r*n(dim)*(dim-2)
    h_rd = cp*t_rd + 0.5d0*v2_rd
    h_r = cp*t_r + 0.5d0*v2_r
!=== LEFT/RIGHT MACH NUMBERS: ===
!A_f = min(A_L**2/max(A_L, abs(Vn_L)), A_R**2/max(A_R, abs(Vn_R))) 
    temp0 = DSQRT(a_l*a_r)
    IF (a_l*a_r .EQ. 0.0) THEN
      a_fd = 0.0_8
    ELSE
      a_fd = (a_r*a_ld+a_l*a_rd)/(2.D0*DSQRT(a_l*a_r))
    END IF
    a_f = temp0
    m_ld = (vn_ld-vn_l*a_fd/a_f)/a_f
    m_l = vn_l/a_f
    m_rd = (vn_rd-vn_r*a_fd/a_f)/a_f
    m_r = vn_r/a_f
!=== CONVECTION-PRESSURE SPLITTING: ===
    result1d = M_PLUS_PLUS_D(m_l, m_ld, result1)
    result2d = M_MINUS_MINUS_D(m_r, m_rd, result2)
    m_fd = result1d + result2d
    m_f = result1 + result2
    result1d = P_PLUS_PLUS_D(m_l, m_ld, result1)
    result2d = P_MINUS_MINUS_D(m_r, m_rd, result2)
    p_fd = p_l*result1d + result1*p_ld + p_r*result2d + result2*p_rd
    p_f = result1*p_l + result2*p_r
!=== LEFT/RIGHT STATE CONVECTION VECTORS: ===
    fi_l(1) = 1.0d0
    fi_ld = 0.0_8
    fi_ld(2) = u_ld
    fi_l(2) = u_l
    fi_ld(3) = v_ld
    fi_l(3) = v_l
    fi_ld(4) = w_ld
    fi_l(4) = w_l
    fi_ld(dim+2) = h_ld
    fi_l(dim+2) = h_l
    fi_r(1) = 1.0d0
    fi_rd = 0.0_8
    fi_rd(2) = u_rd
    fi_r(2) = u_r
    fi_rd(3) = v_rd
    fi_r(3) = v_r
    fi_rd(4) = w_rd
    fi_r(4) = w_r
    fi_rd(dim+2) = h_rd
    fi_r(dim+2) = h_r
!=== PRESSURE PATR: ===
    f_p = 0.0d0
    f_pd = 0.0_8
    f_pd(2:dim+1) = n(1:dim)*p_fd
    f_p(2:dim+1) = p_f*n(1:dim)
!=== FACE FLUXES: ===
    wl = 0.5d0*(1.0d0+SIGN(1.0d0, m_f))
    wr = 1.0d0 - wl
    temp1 = wl*ro_l*fi_l(1:dim+2) + wr*ro_r*fi_r(1:dim+2)
    fluxesd(1:dim+2) = temp1*(m_f*a_fd+a_f*m_fd) + a_f*m_f*(wl*(fi_l(1:&
&     dim+2)*ro_ld+ro_l*fi_ld(1:dim+2))+wr*(fi_r(1:dim+2)*ro_rd+ro_r*&
&     fi_rd(1:dim+2))) + f_pd(1:dim+2)

  CONTAINS
!  Differentiation of m_plus_plus in forward (tangent) mode:
!   variations   of useful results: m_plus_plus
!   with respect to varying inputs: m
    REAL*8 FUNCTION M_PLUS_PLUS_D(m, md, m_plus_plus)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      REAL*8, INTENT(IN) :: md
      REAL*8 :: m_plus_plus
      IF (m .GE. 1.0d0) THEN
        m_plus_plus_d = md
        m_plus_plus = m
      ELSE IF (m .LE. -1.0d0) THEN
        m_plus_plus = 0.0d0
        m_plus_plus_d = 0.0_8
      ELSE
        m_plus_plus_d = (0.25d0*2*(m+1.0d0)+0.125d0*2**2*(m**2-1.0d0)*m)&
&         *md
        m_plus_plus = 0.25d0*(m+1.0d0)**2 + 0.125d0*(m**2-1.0d0)**2
      END IF
    END FUNCTION M_PLUS_PLUS_D

    PURE REAL*8 FUNCTION M_PLUS_PLUS(m)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      IF (m .GE. 1.0d0) THEN
        m_plus_plus = m
      ELSE IF (m .LE. -1.0d0) THEN
        m_plus_plus = 0.0d0
      ELSE
        m_plus_plus = 0.25d0*(m+1.0d0)**2 + 0.125d0*(m**2-1.0d0)**2
      END IF
    END FUNCTION M_PLUS_PLUS

!  Differentiation of m_minus_minus in forward (tangent) mode:
!   variations   of useful results: m_minus_minus
!   with respect to varying inputs: m
    REAL*8 FUNCTION M_MINUS_MINUS_D(m, md, m_minus_minus)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      REAL*8, INTENT(IN) :: md
      REAL*8 :: m_minus_minus
      IF (m .GE. 1.0d0) THEN
        m_minus_minus = 0.0d0
        m_minus_minus_d = 0.0_8
      ELSE IF (m .LE. -1.0d0) THEN
        m_minus_minus_d = md
        m_minus_minus = m
      ELSE
        m_minus_minus_d = -((0.25d0*2*(m-1.0d0)+0.125d0*2**2*(m**2-1.0d0&
&         )*m)*md)
        m_minus_minus = -(0.25d0*(m-1.0d0)**2) - 0.125d0*(m**2-1.0d0)**2
      END IF
    END FUNCTION M_MINUS_MINUS_D

    PURE REAL*8 FUNCTION M_MINUS_MINUS(m)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      IF (m .GE. 1.0d0) THEN
        m_minus_minus = 0.0d0
      ELSE IF (m .LE. -1.0d0) THEN
        m_minus_minus = m
      ELSE
        m_minus_minus = -(0.25d0*(m-1.0d0)**2) - 0.125d0*(m**2-1.0d0)**2
      END IF
    END FUNCTION M_MINUS_MINUS

!  Differentiation of p_plus_plus in forward (tangent) mode:
!   variations   of useful results: p_plus_plus
!   with respect to varying inputs: m
    REAL*8 FUNCTION P_PLUS_PLUS_D(m, md, p_plus_plus)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      REAL*8, INTENT(IN) :: md
      REAL*8 :: temp
      REAL*8 :: p_plus_plus
      IF (m .GE. 1.0d0) THEN
        p_plus_plus = 1.0d0
        p_plus_plus_d = 0.0_8
      ELSE IF (m .LE. -1.0d0) THEN
        p_plus_plus = 0.0d0
        p_plus_plus_d = 0.0_8
      ELSE
        temp = m*m - 1.0d0
        p_plus_plus_d = (0.25d0*((2.0d0-m)*2*(m+1.0d0)-(m+1.0d0)**2)+&
&         0.1875d0*(temp**2+m**2*2**2*temp))*md
        p_plus_plus = 0.25d0*((m+1.0d0)*(m+1.0d0)*(2.0d0-m)) + 0.1875d0*&
&         (m*(temp*temp))
      END IF
    END FUNCTION P_PLUS_PLUS_D

    PURE REAL*8 FUNCTION P_PLUS_PLUS(m)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      IF (m .GE. 1.0d0) THEN
        p_plus_plus = 1.0d0
      ELSE IF (m .LE. -1.0d0) THEN
        p_plus_plus = 0.0d0
      ELSE
        p_plus_plus = 0.25d0*(m+1.0d0)**2*(2.0d0-m) + 0.1875d0*m*(m**2-&
&         1.0d0)**2
      END IF
    END FUNCTION P_PLUS_PLUS

!  Differentiation of p_minus_minus in forward (tangent) mode:
!   variations   of useful results: p_minus_minus
!   with respect to varying inputs: m
    REAL*8 FUNCTION P_MINUS_MINUS_D(m, md, p_minus_minus)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      REAL*8, INTENT(IN) :: md
      REAL*8 :: temp
      REAL*8 :: p_minus_minus
      IF (m .GE. 1.0d0) THEN
        p_minus_minus = 0.0d0
        p_minus_minus_d = 0.0_8
      ELSE IF (m .LE. -1.0d0) THEN
        p_minus_minus = 1.0d0
        p_minus_minus_d = 0.0_8
      ELSE
        temp = m*m - 1.0d0
        p_minus_minus_d = (0.25d0*((m+2.0d0)*2*(m-1.0d0)+(m-1.0d0)**2)-&
&         0.1875d0*(temp**2+m**2*2**2*temp))*md
        p_minus_minus = 0.25d0*((m-1.0d0)*(m-1.0d0)*(m+2.0d0)) - &
&         0.1875d0*(m*(temp*temp))
      END IF
    END FUNCTION P_MINUS_MINUS_D

    PURE REAL*8 FUNCTION P_MINUS_MINUS(m)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: m
      IF (m .GE. 1.0d0) THEN
        p_minus_minus = 0.0d0
      ELSE IF (m .LE. -1.0d0) THEN
        p_minus_minus = 1.0d0
      ELSE
        p_minus_minus = 0.25d0*(m-1.0d0)**2*(2.0d0+m) - 0.1875d0*m*(m**2&
&         -1.0d0)**2
      END IF
    END FUNCTION P_MINUS_MINUS

  END SUBROUTINE

end module

