MODULE option
!
  IMPLICIT NONE
  logical :: rep_closure
  integer :: IMEX, BDF2, rheo
  character(LEN=20) :: adv_scheme, regularization

END MODULE option