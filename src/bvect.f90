subroutine bvect(tauair, un1, Cw, b)
  use size
  use resolution
  use properties
  use numerical
  use global_var
  use shallow_water
  use option
  use rheology

  implicit none
      
  integer :: i
  double precision h_at_u, a_at_u
  double precision :: lap_P_half(1:nx+1)
  double precision, intent(in) :: tauair(1:nx+1), un1(1:nx+1), Cw(1:nx+1)
  double precision, intent(out):: b(1:nx+1)


  if (nonlocal) then 
    call laplacian(2d0*P_half, lap_P_half)
  endif

  b(1)    = 0d0 ! close bc
  b(nx+1) = 0d0 ! close bc 
  ! b(nx) = 0d0
  do i = 2, nx

     h_at_u = ( h(i) + h(i-1) ) / 2d0
     a_at_u = ( A(i) + A(i-1) ) / 2d0
     a_at_u=max(a_at_u, smallA)

     b(i) = a_at_u*tauair(i) + a_at_u*Cw(i)*uwn2(i) - &
            ( P_half(i) - P_half(i-1) ) / Deltax + ( rho * h_at_u * un1(i) ) / Deltat - &
            rho * h_at_u * ge * ( etawn1(i) - etawn1(i-1) ) / Deltax

    if (nonlocal) then 
        b(i) = b(i) + l_scale2*( lap_P_half(i) - lap_P_half(i-1) ) / Deltax
    endif

  enddo
   
  return
end subroutine bvect
    
