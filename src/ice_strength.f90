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
  double precision :: deno, P0, rho_prime, dudx
  double precision :: Pstar_prime(1:nx+1)
  
  rho_prime = ((rhowater - rho)/rhowater)*rho


  if (Pstart_change) then 
    do i = 1, nx ! for tracer points

			dudx = ( uin(i+1) - uin(i) ) / Deltax

      ! deno = alpha*sqrt( (dudx)**2d0 ) ! small2 is there to avoid div by zero
      deno = alpha2*(dudx)**2d0

      if (P0_constant) then 
        P0 = rho_prime*ge*20d0
      else
        P0 = rho_prime*ge*hin(i)
      endif

      Pstar_prime(i) = P0 + (Pstar - P0)*tanh(deno/sqrt(deno+denomin_P**2d0))
      ! Pstar_prime(i) = P0 + (Pstar - P0)*tanh(deno/denomin_P)
    enddo
  
  else  
    do i = 1, nx
      Pstar_prime(i) = Pstar
    enddo
  endif

  if (rheo .eq. 1) then

    Pp_half(0)    = 0d0 ! ! sea ice pressure / 2d0
    Pp_half(nx+1) = 0d0
    Tp_half(0)    = 0d0 ! ! sea ice pressure / 2d0
    Tp_half(nx+1) = 0d0

    do i = 1, nx
      Pp_half(i) = 0.5d0 * Pstar_prime(i) * hin(i) * dexp(-C * ( 1d0 - Ain(i) ) )
      Tp_half(i) = kt*Pp_half(i)
    enddo
  
  elseif ( rheo .eq. 2 ) then
    
    Pp_half(0)    = 0d0 ! ! sea ice pressure / 2d0
    Pp_half(nx+1) = 0d0
    Tp_half(0)    = 0d0 ! ! sea ice pressure / 2d0
    Tp_half(nx+1) = 0d0

    do i = 1, nx
      Pp_half(i) = Pstar * hin(i) * dexp(-C * ( 1d0 - Ain(i) ))
      Tp_half(i) = kt*Pp_half(i)
      ! Peq(i) = rho * hin(i) * (( d_average * shear_I(i) ) / ( Ain(i) - Phi_0 ))**2 

    enddo

  endif
    

!------- set p = 0 at open boundaries for proper care of open bc --------------
!                    see p.1241-1242 for details              
!--- set dh/dx, dA/dx = 0 at the outside cell when there is an open bc --------

  return
end subroutine ice_strength
      




