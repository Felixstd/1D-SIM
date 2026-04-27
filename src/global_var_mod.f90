MODULE global_var
use size
! h    : average sea ice thickness  [m]
! uice : ice velocity               [m/s]

!   h(0)--u(1)--h(1)--u(2)--...u(n)--h(n)--u(n+1)--h(n+1)
!   h(0) and h(n+1) are used for open bcs

  IMPLICIT NONE
  DOUBLE PRECISION :: h(0:nx+1)
  DOUBLE PRECISION :: A(0:nx+1)
  DOUBLE PRECISION :: hn1(0:nx+1)
  DOUBLE PRECISION :: An1(0:nx+1)
  DOUBLE PRECISION :: hn2(0:nx+1)
  DOUBLE PRECISION :: An2(0:nx+1)
  DOUBLE PRECISION :: Pp_half(0:nx+1)
  DOUBLE PRECISION :: Tp_half(0:nx+1)
  DOUBLE PRECISION :: P_half(0:nx+1)
  DOUBLE PRECISION :: bathy(0:nx+1)
  DOUBLE PRECISION :: scaling(1:nx+1)

  ! Muphi
  DOUBLE PRECISION :: Peq(0:nx+1)
  DOUBLE PRECISION :: Pp(0:nx+1)
  DOUBLE PRECISION :: Pmax(0:nx+1)
  DOUBLE PRECISION :: Phi_I(0:nx+1)
  DOUBLE PRECISION :: shear_I(0:nx+1)
  DOUBLE PRECISION :: mu_I(0:nx+1)
  DOUBLE PRECISION :: Inertial(0:nx+1)
  DOUBLE PRECISION :: W_sigma(0:nx+1)


  ! Energy 
  DOUBLE PRECISION :: P_pot(0:nx+1)
  DOUBLE PRECISION :: P_fric_R(0:nx+1)
  DOUBLE PRECISION :: P_fric_S(0:nx+1)
  DOUBLE PRECISION :: P_fric_R_visc(0:nx+1)
  DOUBLE PRECISION :: P_fric_S_visc(0:nx+1)
  DOUBLE PRECISION :: P_fric_R_vZ_pE(0:nx+1)
  DOUBLE PRECISION :: P_fric_R_pZ_vE(0:nx+1)
  DOUBLE PRECISION :: P_fric_R_plas(0:nx+1)
  DOUBLE PRECISION :: P_fric_S_plas(0:nx+1)
  DOUBLE PRECISION :: P_lat(0:nx+1)
  DOUBLE PRECISION :: P_h(0:nx+1)
  DOUBLE PRECISION :: P_w(0:nx+1)


END MODULE global_var
