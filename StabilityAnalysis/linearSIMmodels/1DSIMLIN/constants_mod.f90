MODULE rheology
!
! C : ice strength parameter
! Pstar : ice compression strength parameter
! ell_2 : 1/ellipticity**2
! kt: for tensile strength, T = kt x P
  IMPLICIT NONE
  DOUBLE PRECISION :: C, Pstar, alpha, alpha2, kt
  DOUBLE PRECISION :: e_2, small2

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
  DOUBLE PRECISION :: mu_b


END MODULE muphi

MODULE properties

  IMPLICIT NONE
  DOUBLE PRECISION :: rho, rhowater, ge

END MODULE properties

MODULE resolution

  IMPLICIT NONE
  DOUBLE PRECISION :: Deltax, Deltax2
  DOUBLE PRECISION :: Deltat, DtoverDx, Deltate
  INTEGER :: T_tot

END MODULE resolution

MODULE numerical

  IMPLICIT NONE
  DOUBLE PRECISION :: T, smallA

END MODULE numerical