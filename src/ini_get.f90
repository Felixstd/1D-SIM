!****************************************************************************
!     initial thickness and velocity field (u is at t=0, h is at t=Deltat/2)
!****************************************************************************

! no restart possible for the moment jfl WATCHOUT
! no open bcs possible for the moment jfl WATCHOUT

subroutine ini_get (utp, restart, expres, ts_res)

  use size
  use global_var
  use shallow_water
  use MOMeqSW_output
  use option
  use muphi
  use resolution

  implicit none
     
  logical, intent(in) :: restart
  integer, intent(in) :: expres, ts_res
  integer :: i
  double precision, intent(inout)  :: utp(1:nx+1)
  double precision :: rdnb, small, ramp, eps, window, x, L
  double precision :: x1, x2, c1, c2, b1, b2, par1, par2, w, g
  double precision :: sharpness, i_prime

  character(LEN=30) filename  ! restart file name 

  small = 0.0001d0
  eps = 20000 !controls sharpness of the transition

  allocate(etaw(0:nx+1), etawn1(0:nx+1), etawn2(0:nx+1))
  allocate(uw(1:nx+1), uwn1(1:nx+1), uwn2(1:nx+1))
  
  if (oceanSIM) then 
     allocate(duwdt(1:nx+1), gedetawdx(1:nx+1), buw(1:nx+1))
     allocate(tauiw(1:nx+1), tauaw(1:nx+1))
  endif
  
  utp(1)    = 0d0 ! close bc
  utp(nx)   = 0d0
  utp(nx+1) = 0d0 ! close bc
  h(0)    = 0d0
  h(nx+1) = 0d0
  A(0)    = 0d0
  A(nx+1) = 0d0
  uw      = 0d0
  uwn1    = 0d0
  uwn2    = 0d0
  etaw    = 0d0
  etawn1  = 0d0
  etawn2  = 0d0
  
  scaling=1d0 ! initialize scaling field (only used for JFNK)

  if (restart) then 
     print *, 'Restart code should be verified'
     stop
     write (filename,'("output/h_",i3.3,".",i2.2)') ts_res, expres
     open (10, file = filename, status = 'old')
     
     write (filename,'("output/A_",i3.3,".",i2.2)') ts_res, expres
     open (11, file = filename, status = 'old')

     write (filename,'("output/u_",i3.3,".",i2.2)') ts_res, expres
     open (12, file = filename, status = 'old')

     read (10,*) ( h(i), i = 0, nx+1 )
     read (11,*) ( A(i), i = 0, nx+1 )
     read (12,*) ( utp(i), i = 1, nx+1 )

     close(10)
     close(11)
     close(12)

  else ! specify initial fields

  sharpness = 0.0000001**(-0.1)

  do i = 1, nx

      if ((initcond_vel .eq. 'GraysmoothConv') .and. (initcond .eq. 'stepsmooth')) then 
         A(i) = min(max(tanh((real(i)-real(nx)/3)/sharpness) &
                     - tanh((real(i)-2*real(nx)/3)/sharpness), 0d0)/5, 1d0)
            if (A(i) .gt. 0d0) then  
               h(i) = A(i)
            endif

      elseif (initcond .eq. 'step') then
         A(i) = min(max((real(10000d0)-(real(i)-real(nx)/2d0)**2)**(1/20d0), 0d0), 1d0)
         if (A(i) .gt. 0d0) then  
            h(i) = 1d0
         endif

      elseif (initcond .eq. 'stepsmooth') then
            A(i) = min(max(tanh((real(i)-real(nx)/3)/sharpness) &
                     - tanh((real(i)-2*real(nx)/3)/sharpness), 0d0)/2, 1d0)
            if (A(i) .gt. 0d0) then  
               h(i) = A(i)
            endif

      
      elseif (initcond .eq. 'gaussian') then
         A(i) = min(0.1d0*exp(-D*(real(i) - real(nx)/2d0)**2d0), 1d0)
         h(i) = min(0.1d0*exp(-D*(real(i) - real(nx)/2d0)**2d0), 1d0)

      elseif (initcond .eq. 'supergaussian') then
         A(i) = min(exp(-D*(real(i)-real(nx)/2d0)**(2*n)), 1d0)
         h(i) = min(exp(-D*(real(i)-real(nx)/2d0)**(2*n)), 1d0)

      elseif (initcond .eq. 'gaussianh') then
         A(i) = 1d0
         h(i) = min(0.1*exp(-D*(real(i) - real(nx)/2d0)**2d0), 1d0)

      elseif (initcond .eq. 'constantAsteph') then
         A(i) = 1d0
         h(i) = min(max((real(10000d0)-(real(i)-real(nx)/2d0)**2)**(1/20d0), 0d0), 1d0)

      elseif (initcond .eq. 'constanthstepA') then
         h(i) = 1d0
         A(i) = min(max((real(10000d0)-(real(i)-real(nx)/2d0)**2)**(1/20d0), 0d0), 1d0)

      elseif (initcond .eq. 'constants') then 
         A(i) = 1d0
         h(i) = 1d0


      endif

  enddo

    do i = 2, nx

      if (initcond_vel .eq. 'Gray') then 
      
         if (h(i) > 0d0) then 
            utp(i) = (i - nx/2d0)/1e5
         
         else 
            utp(i) = 0
         endif

      elseif (initcond_vel .eq. 'Graysmooth') then

         L = nx* Deltax
         x = (i-1) * Deltax
         ramp = (x - L/2)
         window = 0.5 * (tanh((ramp + 80e3)/eps) - tanh((ramp - 80e3)/eps))
         utp(i) = (ramp/1d8)*window

      elseif (initcond_vel .eq. 'GraysmoothConv') then
         
         L = nx* Deltax
         x = (i-1) * Deltax
         ramp = (x - L/2)
         window = 0.5 * (tanh((ramp + 80e3)/eps) - tanh((ramp - 80e3)/eps))
         utp(i) = -(ramp/1d8)*window
      
      elseif (initcond_vel .eq. 'parabolaConv') then 
         w = 10d0
         g = 100d0
         c1 = nx/2d0 - g/2d0
         c2 = nx/2d0 + g/2d0

         x1 = (i - c1)/w
         x2 = (i - c2)/w
         if (abs(x1) < 1) then
            b1 = exp(-1d0/(1-x1**2d0))
         else
            b1 = 0
         endif

         if (abs(x2) < 1) then
            b2 = exp(-1d0/(1-x2**2d0))
         else
            b2 = 0
         endif


         par1 = -(1d0-((i - c1)/w)**2d0)*b1/w
         par2 = (1d0-((i - c2)/w)**2d0)*b2/w

         utp(i) = par1 + par2+0.001


      elseif (initcond_vel .eq. 'zero') then 
         utp(i) = 0d0
      endif 
      
      uw(i)  = 0d0
  enddo

  
  if (oceanSIM) then 
     uwn1=uw
     uwn2=uw
     etawn1=etaw
     etawn2=etaw
     duwdt=0d0
     gedetawdx=0d0
     tauiw=0d0
     tauaw=0d0
     buw=0d0
  endif
  
!  do i = 11, nx-10
!     h(i) = 1d0
!     A(i) = 0.7d0
!     A(i) = i/(nx*1d0) - 0.5d0/(1d0*nx) ! 0 at West wall and 1 at East wall
!     h(i) = max(1d-06, h(i))
!     bathy(i)=100d0
!  enddo

!  do i = 1, 20
!     bathy(i)=10d0
!  enddo
 
!  do i = 1, 19
!     h(i) = 1d0
!     A(i) = 1d0
!  enddo
 
  endif

  return
end subroutine ini_get






