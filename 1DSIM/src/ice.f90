program ice

    use size 
    use rheology
    use muphi
    use properties
    use resolution
    use global_var
    use numerical
    use option

    implicit none 

    logical :: restart
    integer :: i, ii, ts, tsini, nstep, tsfin, k, s, Nmax_OL, solver, idiag
    integer :: out_step(11), expnb, expres, ts_res
    integer :: time_out, n_output, output_time
    double precision :: e
    double precision :: u(0:nx+1), un1(0:nx+1), un2(0:nx+1)
    double precision :: zeta(0:nx+1), eta(0:nx+1), mn1(0:nx+1), m(0:nx+1)
    double precision :: meanvalue, time1, time2, timecrap
    double precision :: nbhr

    character filename*64

    out_step = 0

    !------------------------------------------------------------------------
!     Input by user
!------------------------------------------------------------------------
    expnb      = 1
    rheo           = 1
    restart        = .false.
    regularization = 'tanh' ! tanh, Kreyscher, capping (Hibler)
    adv_scheme     = 'upwind' ! upwind, upwindRK2, semilag

    !   T_tot      = 2*60*60
    T_tot      = 10*24*60*60
    Deltat     = 10 ! time step [s]
    !   nstep      = 1440     ! lenght of the run in nb of time steps
    nstep = T_tot/Deltat !lenght of the run in nb of time steps
    ! nstep = 1
    
    write (filename,'("output/time_run.",i2.2)') expnb
    open (10, file = filename, status = 'unknown')
    time_out = 0
    output_time = 0
    write(10, *) (output_time)
    out_step(1) = 0
    do n_output = 2, 11
        time_out = time_out + T_tot / 10
        output_time = output_time + nstep/10
        write(10, *) (output_time)
        out_step(n_output)=output_time   
    enddo
    close(10)

    
    expres     = 2
    ts_res     = 50 ! time level of restart (!!! watchout for Deltat !!!)
    !   out_step(1)=1
    !   out_step(2)=1440

    if (restart) then
        tsini = ts_res + 1
    else
        tsini = 1
    endif
  
    tsfin = tsini - 1 + nstep

    if ( nx .eq. 100 ) then 
        Deltax   =  20d03  ! grid size [m], the domain is always 2000 km 
    elseif  ( nx .eq. 200 ) then
        Deltax   =  10d03            
    elseif  (( nx .eq. 400 ) .or. (nx .eq. 600)) then
        Deltax   =  10d03       
    elseif  ( nx .eq. 500 ) then
        Deltax   =  4d03        
    else
        print *,  'Wrong grid size dimension', nx
        STOP
    endif

    Deltax2 = Deltax ** 2
    DtoverDx = Deltat / Deltax

    C          = 20d0         ! ice strength parameter (watchout no A for now)
    Pstar      = 27.5d03      ! ice compression strength parameter
    e          = 2d0          ! ratio long to short axis of ellipse
    e_2        = 1/(e**2d0)   !
    alpha      = sqrt(1d0 + e_2)
    alpha2     = 1d0 + e_2
    kt         = 0d0          ! T = kt * P (1.0 in Konig and Holland, 2010)
    rho        = 900d0        ! ice density
    ge         = 9.8d0        ! Earth's gravitional acceleration
    
    ! Mu-Phi Parameters 
    d_average  = 1d03
    mu_0       = 0.1
    mu_infty   = 0.8
    I_0        = 1e-3         
    mu_b       = 0.1
    Phi_0      = 1
    c_phi      = 1


    call ini_get (u, restart, expres, ts_res)
    call output_results(ts, expnb, u, zeta, eta)
    un1=u
    hn1=h
    An1=A

    do ts = tsini, tsfin

        nbhr = nbhr + Deltat / 3600d0
        print *, 'time level, cumulative time (h) =', ts, nbhr
        
        
        call cpu_time(timecrap)
        call cpu_time(time1)

        un1=u
        hn1=h
        An1=A

        do i = 0, nx+1
            mn1 = rho*hn1(i)*un1(i)
        enddo

        call shear(un1)
        call ice_strength (hn1, An1) ! standard approach no IMEX 
        call inertial_number()
        call angle_friction_mu()
        
        call viscouscoefficient(zeta, eta)
        call stress(un1, zeta, eta)

        call FCT_continuity(un1, hn1, An1, A, h)
        call FCT_momentum(un1, hn1, An1,mn1, m)

        do i = 0, nx+1
            ! print*, m(i)
            if (h(i) < 1e-5) then 
                u(i) = 0d0
            else
                u(i) = m(i)/(rho*h(i) + 1e-12)
            endif

            ! if (h(i) + 1e-12 > 1e-12)
            ! print*, h(i), u(i)
        enddo

!------------------------------------------------------------------------
!     output results
!------------------------------------------------------------------------

        if (any(ts == out_step(1:10))) then
            print *, 'outputting results'
            call output_results(ts, expnb, un1, zeta, eta)
            ! call output_file(e, gamma_nl, solver, expnb)
        endif

    enddo


end program ice