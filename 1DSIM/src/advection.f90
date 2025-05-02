subroutine FCT_continuity(utp, hn1in, An1in, Aout, hout)

    use size
    use rheology
    use muphi
    use option
    use resolution
    use global_var

    implicit none 
    
    integer                      :: i
    double precision, intent(in) :: utp(0:nx+1)
    double precision, intent(in) :: hn1in(0:nx+1), An1in(0:nx+1)
    double precision, intent(out) :: hout(0:nx+1), Aout(0:nx+1)
    double precision             :: hout_est(0:nx+1), Aout_est(0:nx+1)
    double precision             :: fluxphalf_low, fluxmhalf_low
    double precision             :: fluxphalf_high, fluxmhalf_high
    double precision             :: Afluxmhalf, Afluxphalf
    double precision             :: ACphalf, ACmhalf
    double precision             :: DxoverDt


    DxoverDt = Deltax/Deltat


    hout(0) = 0d0    ! closed b.c.s
    hout(nx+1) = 0d0
    Aout(0) = 0d0
    Aout(nx+1) = 0d0

    fluxphalf_low =  0d0
    fluxmhalf_low =  0d0
    fluxphalf_high = 0d0
    fluxmhalf_high = 0d0

    Aout_est(1) = 0d0
    Aout_est(nx) = 0d0

    do i = 2, nx-1

        ! hout(i) = hn1in(i)
        !-------- Continuity for Sea Ice Concentration --------!
        call calc_flux_continuity(utp(i), utp(i+1), utp(i-1), An1in(i-1), &
            An1in(i), An1in(i+1), fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high)

        Afluxphalf = fluxphalf_high - fluxphalf_low
        Afluxmhalf = fluxmhalf_high - fluxmhalf_low

        Aout_est(i) = An1in(i) - DtoverDx*(fluxphalf_low - fluxmhalf_low)

        Aout_est(i) = max(Aout_est(i), 0d0)
        Aout_est(i) = min(Aout_est(i), 1d0)


        !-------- Continuity for Sea Ice Thickness --------!
        call calc_flux_continuity(utp(i), utp(i+1), utp(i-1), hn1in(i-1), &
            hn1in(i), hn1in(i+1), fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high)
        
        Afluxphalf = fluxphalf_high - fluxphalf_low
        Afluxmhalf = fluxmhalf_high - fluxmhalf_low

        hout_est(i) = hn1in(i) - DtoverDx*(fluxphalf_low - fluxmhalf_low)

        hout_est(i) = max(hout_est(i), 0d0)


    enddo

    fluxphalf_low =  0d0
    fluxmhalf_low =  0d0
    fluxphalf_high = 0d0
    fluxmhalf_high = 0d0

    do i = 2, nx-1

        !-------- Continuity for Sea Ice Concentration --------!

        call calc_flux_continuity(utp(i), utp(i+1), utp(i-1), An1in(i-1), &
        An1in(i), An1in(i+1), fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high)

        Afluxphalf = fluxphalf_high - fluxphalf_low
        Afluxmhalf = fluxmhalf_high - fluxmhalf_low

        ACphalf = sign(1d0, Afluxphalf) * max(0d0, min(abs(Afluxphalf), &
                    sign(1d0, Afluxphalf)*(Aout_est(i+2) - Aout_est(i+1))*DxoverDt, &
                    sign(1d0, Afluxphalf)*(Aout_est(i) - Aout_est(i-1))*DxoverDt))

        ACmhalf = sign(1d0, Afluxmhalf) * max(0d0, min(abs(Afluxmhalf), &
                    sign(1d0, Afluxmhalf)*(Aout_est(i-2) - Aout_est(i-1))*DxoverDt, &
                    sign(1d0, Afluxmhalf)*(Aout_est(i+1) - Aout_est(i))*DxoverDt))

        Aout(i) = Aout_est(i) - DtoverDx*(ACphalf - ACmhalf)
        
        Aout(i) = max(Aout(i), 0d0)
        Aout(i) = min(Aout(i), 1d0) 


        !-------- Continuity for Sea Ice thickness --------!

        call calc_flux_continuity(utp(i), utp(i+1), utp(i-1), hn1in(i-1), &
             hn1in(i), hn1in(i+1), fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high)

        Afluxphalf = fluxphalf_high - fluxphalf_low
        Afluxmhalf = fluxmhalf_high - fluxmhalf_low

        ACphalf = sign(1d0, Afluxphalf) * max(0d0, min(abs(Afluxphalf), &
                    sign(1d0, Afluxphalf)*(hout_est(i+2) - hout_est(i+1))*DxoverDt, &
                    sign(1d0, Afluxphalf)*(hout_est(i) - hout_est(i-1))*DxoverDt))

        ACmhalf = sign(1d0, Afluxmhalf) * max(0d0, min(abs(Afluxmhalf), &
                    sign(1d0, Afluxmhalf)*(hout_est(i-2) - hout_est(i-1))*DxoverDt, &
                    sign(1d0, Afluxmhalf)*(hout_est(i+1) - hout_est(i))*DxoverDt))

        hout(i) = hout_est(i) - DtoverDx*(ACphalf - ACmhalf)
        
        hout(i) = max(Aout(i), 0d0)

    enddo

    Aout(1) = 0d0
    Aout(nx) = 0d0
    hout(1) = 0d0
    hout(nx) = 0d0


    return
end subroutine FCT_continuity


subroutine calc_flux_continuity(ui, uip1, uim1, Tim1, Ti, Tip1, fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high)


    use resolution

    implicit none

    double precision, intent(in)  :: ui, uip1, uim1, Tim1, Ti, Tip1 ! T=tracer
    double precision, intent(out) :: fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high
    double precision :: flux, aphalf, amhalf

    
    !------ Lower order scheme -------!
    !------    Upwind scheme   -------!

    if (ui .gt. Deltax/Deltat) print *, 'WARNING: u > dx/dt', ui, Deltax/Deltat

    if (ui .ge. 0d0) then ! left side of cell                                                                                     
        fluxmhalf_low = ui*Tim1
    else
        fluxmhalf_low = ui*Ti
    endif

    if (uip1 .ge. 0d0) then ! right side of cell                                                                                  
        fluxphalf_low = uip1*Ti
    else
        fluxphalf_low = uip1*Tip1
    endif

    ! print* ,fluxphalf_low, fluxmhalf_low

    !------ Higher order scheme -------!
    !------    Lax-Wendroff     -------!

    if (uip1 .ne. ui) then

        aphalf = (uip1*Tip1 - ui*Ti) / (uip1 - ui)
    
    else
        aphalf = Ti
    endif
    
    fluxphalf_high = 1 / 2 * ((uip1*Tip1 + ui*Ti) - aphalf**2 * DtoverDx *(uip1 - ui))

    if (uim1 .ne. ui) then

        amhalf = (ui*Ti - uim1*Tim1) / (ui - uim1)
    else
        amhalf = Ti
    endif

    fluxmhalf_high = 1/2*((ui*Ti + uim1*Tim1) - amhalf**2 * DtoverDx * (ui - uim1))

    return

end subroutine calc_flux_continuity


subroutine FCT_momentum(utp, hn1in, An1in, mn1in, mout)

    use size
    use rheology
    use muphi
    use option
    use properties
    use resolution
    use global_var

    implicit none 
    
    integer                      :: i
    double precision, intent(in) :: utp(0:nx+1)
    double precision, intent(in) :: hn1in(0:nx+1), An1in(0:nx+1), mn1in(0:nx+1)
    double precision, intent(out):: mout(0:nx+1)
    double precision             :: mout_est(0:nx+1), T(0:nx+1)
    double precision             :: fluxphalf_low, fluxmhalf_low
    double precision             :: fluxphalf_high, fluxmhalf_high
    double precision             :: Afluxmhalf, Afluxphalf
    double precision             :: ACphalf, ACmhalf
    double precision             :: DxoverDt


    DxoverDt = Deltax/Deltat


    ! utp(0) = 0d0    ! closed b.c.s
    ! utp(nx+1) = 0d0
    ! utp(0) = 0d0
    ! utp(nx+1) = 0d0

    fluxphalf_low =  0d0
    fluxmhalf_low =  0d0
    fluxphalf_high = 0d0
    fluxmhalf_high = 0d0

    do i = 0, nx+1
        T(i) = rho*hn1in(i)*utp(i)**2 + sig11(i)
        ! print*, T(i)
    enddo
    T(0) = 0d0
    T(1) = 0d0
    T(nx) = 0d0
    T(nx+1) = 0d0


    do i = 2, nx-1

        ! hout(i) = hn1in(i)
        !-------- Continuity for Sea Ice Concentration --------!
        call calc_flux_momentum(mn1in(i), mn1in(i+1), mn1in(i-1), T(i-1), &
            T(i), T(i+1), fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high)

        Afluxphalf = fluxphalf_high - fluxphalf_low
        Afluxmhalf = fluxmhalf_high - fluxmhalf_low

        mout_est(i) = mn1in(i) - DtoverDx*(fluxphalf_low - fluxmhalf_low)
        ! print*, 'est', mout_est(i)

    enddo

    fluxphalf_low =  0d0
    fluxmhalf_low =  0d0
    fluxphalf_high = 0d0
    fluxmhalf_high = 0d0

    do i = 2, nx-1

        !-------- Continuity for Sea Ice Concentration --------!

        call calc_flux_momentum(mn1in(i), mn1in(i+1), mn1in(i-1), T(i-1), &
        T(i), T(i+1), fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high)

        Afluxphalf = fluxphalf_high - fluxphalf_low
        Afluxmhalf = fluxmhalf_high - fluxmhalf_low

        ACphalf = sign(1d0, Afluxphalf) * max(0d0, min(abs(Afluxphalf), &
                    sign(1d0, Afluxphalf)*(mout_est(i+2) - mout_est(i+1))*DxoverDt, &
                    sign(1d0, Afluxphalf)*(mout_est(i) - mout_est(i-1))*DxoverDt))

        ACmhalf = sign(1d0, Afluxmhalf) * max(0d0, min(abs(Afluxmhalf), &
                    sign(1d0, Afluxmhalf)*(mout_est(i-2) - mout_est(i-1))*DxoverDt, &
                    sign(1d0, Afluxmhalf)*(mout_est(i+1) - mout_est(i))*DxoverDt))

        mout(i) = mout_est(i) - DtoverDx*(ACphalf - ACmhalf)
        ! print*, 'out', mout(i)
    enddo
    
    mout(0) = 0d0
    mout(nx+1) = 0d0
    mout(1) = 0d0
    mout(nx) = 0d0

    return
end subroutine FCT_momentum


subroutine calc_flux_momentum(ui, uip1, uim1, Tim1, Ti, Tip1, fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high)

    
    use resolution
    use properties
    use global_var


    implicit none

    double precision, intent(in)  :: ui, uip1, uim1, Tim1, Ti, Tip1 ! T=tracer
    double precision, intent(out) :: fluxphalf_low, fluxmhalf_low, fluxphalf_high, fluxmhalf_high
    double precision :: flux, aphalf, amhalf, uip1_local, uim1_local

    !------ Lower order scheme -------!
    !------    Upwind scheme   -------!

    uip1_local = (Tip1 - Ti) / (uip1 - ui + 1e-20)
    ! print*, uip1_local
    if (uip1_local .ge. 0d0) then
        fluxphalf_low = Ti
    else
        fluxphalf_low = Tip1
    endif

    uim1_local = (Ti - Tim1) / (ui - uim1+ 1e-20)
    ! print*, uim1_local
    if (uim1_local .ge. 0d0) then
        fluxmhalf_low = Tim1
    else
        fluxmhalf_low = Ti
    endif

    !------ Higher order scheme -------!
    !------    Lax-Wendroff scheme   -------!

    
    ! if (uip1 .ne. ui) then

    aphalf = (Tip1 - Ti) / (uip1 - ui)
    
    ! else
        ! aphalf = Ti
    ! endif
    
    fluxphalf_high = 1 / 2 * ((Tip1 + Ti) - aphalf**2 * DtoverDx *(Tip1 - Ti))

    ! if (uim1 .ne. ui) then

    amhalf = (Ti - Tim1) / (ui - uim1)
    ! else
    !     amhalf = Ti
    ! endif

    fluxmhalf_high = 1/2*((Ti + Tim1) - amhalf**2 * DtoverDx * (Ti - Tim1))


                                                            
end subroutine calc_flux_momentum
