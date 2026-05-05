MODULE option
!
  IMPLICIT NONE
  logical :: linear_drag, linear_viscous, constant_wind, rampupwind
  logical :: rep_closure, oceanSIM, Asselin, DiagStress, mechenergy
  integer :: IMEX, BDF2, rheo
  character(LEN=20) :: adv_scheme, regularization, initcond, initcond_vel
  logical :: Pstart_change, P0_constant
  logical :: nonlocal

END MODULE option
