subroutine ice_strength ( hin, Ain, uin)
  use size
  use rheology
  use resolution
  use global_var
  use muphi
  use option
  use properties

  implicit none

  integer :: i
  
  double precision, intent(in) :: hin(0:nx+1), Ain(0:nx+1), uin(0:nx+1)  


  Pp_half(0)    = 0d0 ! ! sea ice pressure / 2d0
  Pp_half(nx+1) = 0d0
  Tp_half(0)    = 0d0 ! ! sea ice pressure / 2d0
  Tp_half(nx+1) = 0d0


  if (pp_linear) then 
    do i = 1, nx
      Pp_half(i) = 0.5d0 * Pstar * hin(i) * Ain(i) 
      Tp_half(i) = kt*Pp_half(i)
    enddo
  else

    do i = 1, nx
      Pp_half(i) = 0.5d0 * Pstar * hin(i) * dexp(-C * ( 1d0 - Ain(i) ) )
      Tp_half(i) = kt*Pp_half(i)
    enddo
  endif

  ! print*, "Ice strength constant = ", C

!------- set p = 0 at open boundaries for proper care of open bc --------------
!                    see p.1241-1242 for details              
!--- set dh/dx, dA/dx = 0 at the outside cell when there is an open bc --------

  return
end subroutine ice_strength
      




