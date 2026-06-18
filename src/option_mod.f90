MODULE option
!
  IMPLICIT NONE
  logical :: linear_drag, linear_viscous, constant_wind, rampupwind
  logical :: diverging_winds
  logical :: rep_closure, oceanSIM, Asselin, DiagStress, mechenergy
  integer :: IMEX, BDF2, rheo
  character(LEN=20) :: adv_scheme, regularization, initcond, initcond_vel
  logical :: nonlocal

END MODULE option
