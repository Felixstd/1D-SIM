MODULE rheology
!
! C : ice strength parameter
! Pstar : ice compression strength parameter
! ell_2 : 1/ellipticity**2
! kt: for tensile strength, T = kt x P
  IMPLICIT NONE
  DOUBLE PRECISION :: C, Pstar, alpha, alpha2, kt
  DOUBLE PRECISION :: e, e_2, small2, denomin_P
  DOUBLE PRECISION :: l_scale, l_scale2

END MODULE rheology

MODULE muphi

! mu_0       : 
! mu_infty   : 
! c_phi      : 
! I_0        : 
! d_average  : 
! phi_0      :

  IMPLICIT NONE

  DOUBLE PRECISION :: mu_0, mu_infty 
  DOUBLE PRECISION :: c_phi
  DOUBLE PRECISION :: I_0
  DOUBLE PRECISION :: d_average
  DOUBLE PRECISION :: Phi_0
  DOUBLE PRECISION :: mu_b, D
  DOUBLE PRECISION :: n, eta_max, zeta_max


END MODULE muphi

MODULE properties

  IMPLICIT NONE
  DOUBLE PRECISION :: rho, rhowater, ge
  LOGICAL :: advection_mom

END MODULE properties

MODULE forcing

  IMPLICIT NONE
  DOUBLE PRECISION :: Cda, Cdw, small1, uwind
  double precision :: rhoair, Cdair, Cdwater

END MODULE forcing

MODULE resolution

  IMPLICIT NONE
  DOUBLE PRECISION :: Deltax, Deltax2
  DOUBLE PRECISION :: Deltat, DtoverDx, Deltate
  DOUBLE PRECISION :: T_tot

END MODULE resolution

MODULE numerical

  IMPLICIT NONE
  INTEGER :: N_sub, maxiteSOR, maxiteGMRES, iteSOR_pre
  INTEGER :: Nmax_OL, solver
  DOUBLE PRECISION :: T, smallA
  DOUBLE PRECISION :: omega, tol_SOR, dropini, gamma_nl

END MODULE numerical

