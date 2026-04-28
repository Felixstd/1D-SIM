!*************************************************************************
!     program ice:
!       1D model that calculate the ice thickness (h), concentration (A) 
!       and ice velocity (u).
!       
!       momentum equation:  rho*h(u^n-u^n-1)/Deltat = f(u^n,h^n-1,A^n-1)
!                           This equation is solved implicitly for u^n
!                           (u). u^n-1 is the previous time step solution.
!
!       continuity equation:h^n = f(h^n-1, u^n)  
!                           The new value of h (and A) is obtained by 
!                           advecting h^n-1 with u^n. 
!       
!       author: JF Lemieux
!       version: 1.0 (20 april 2012)
!
!************************************************************************
!------ no advection for A -------!
program ice

    use size
    use rheology
    use muphi
    use forcing
    use properties
    use resolution
    use global_var
    use shallow_water
    use MOMeqSW_output
    use numerical
    use option
    use diag_stress
  
    implicit none

    logical :: p_flag, restart
    integer :: i, ii, ts, tsini, nstep, tsfin, k, s, idiag, readnamelist
    integer :: out_step(21), expnb, expres, ts_res, fgmres_its, fgmres_per_ts
    integer, save :: Nfail, meanN ! nb of failures, mean Newton ite per ts
    double precision :: u(1:nx+1), un1(1:nx+1), un2(1:nx+1)
    double precision :: tauair(1:nx+1)    ! tauair
    double precision :: b(1:nx+1)         ! b vector
    double precision :: zeta(0:nx+1), eta(0:nx+1), sigma(0:nx+1), Cw(1:nx+1), Cb(1:nx+1)
    double precision :: F_uk1(1:nx+1), R_uk1(1:nx+1) ! could use F for R
    double precision :: meanvalue, time1, time2, timecrap
    double precision :: L2norm, nl_target, nbhr

    !-- FSTD --!
    logical :: output_diag
    integer :: time_out, n_output, output_time, count
    double precision :: P_pot_avg, P_fric_R_avg, P_fric_S_avg, P_lat_avg, P_h_avg, P_w_avg, P_fric_R_vZ_pE_avg, P_fric_R_pZ_vE_avg
    double precision :: P_fric_R_visc_avg, P_fric_S_visc_avg, P_fric_R_plas_avg, P_fric_S_plas_avg
    double precision, allocatable :: h_time(:), A_time(:), utp_time(:)
    double precision, allocatable :: P_pot_Tavg(:)
    double precision, allocatable :: P_fric_R_Tavg(:), P_fric_S_Tavg(:)
    double precision, allocatable :: P_fric_S_visc_Tavg(:), P_fric_R_visc_Tavg(:)
    double precision, allocatable :: P_fric_S_plas_Tavg(:), P_fric_R_plas_Tavg(:)
    double precision, allocatable :: P_lat_Tavg(:)
    double precision, allocatable :: P_h_Tavg(:)
    double precision, allocatable :: P_w_Tavg(:)
    double precision, allocatable :: P_fric_R_vZ_pE_Tavg(:), P_fric_R_pZ_vE_Tavg(:)
    
    character filename*64

    out_step = 0
    sigma    = 0d0 ! initial stresses are zero
    Nfail    = 0
    meanN    = 0

!------------------------------------------------------------------------
!     Read run information
!------------------------------------------------------------------------
    print *, 'Read namelist?'
    read *, readnamelist
    PRINT *, readnamelist

    print *, 'Use restart file for (h,A,u,v)?'
    read  *, restart

    print *, 'experiment #?'
    read  *, expnb
    PRINT *, expnb

!------------------------------------------------------------------------
!     Default settings and parameters
!------------------------------------------------------------------------

    call get_default          

!------------------------------------------------------------------------
!     Default settings and parameters
!------------------------------------------------------------------------

    nstep = T_tot/Deltat !lenght of the run in nb of time steps

    !-- writing the output time steps and outputting them --!
    write (filename,'("output/time_run.",i2.2)') expnb
    open (10, file = filename, status = 'unknown')
    
    time_out = 0
    output_time = 0
    write(10, *) (output_time)
    
    out_step(1) = 1
    
    do n_output = 2, 21
        time_out = time_out + T_tot / 20
        output_time = output_time + nstep/20
        print*, output_time
        write(10, *) (output_time)
        out_step(n_output)=output_time   
    
    enddo
    
    close(10)

    !-- Restart --!
    expres     = 2
    ts_res     = 50 ! time level of restart (!!! watchout for Deltat !!!)

!------------------------------------------------------------------------ 
! verify choice of solver and options
!------------------------------------------------------------------------ 

    if (solver .eq. 3 .or. solver .eq. 4) then
        if (IMEX .ne. 0 .and. BDF2 .ne. 0) then
            print *, 'set IMEX=0 and BDF2=0'
            stop
        endif
    endif
    
    if (BDF2 .eq. 1 .and. adv_scheme .ne. 'upwindRK2') then ! semilag should work (leapfrog...) check this
        print *, 'set adv_scheme = upwindRK2'
        stop
    endif
  
!------------------------------------------------------------------------ 
!     Set first time level depending on restart specifications                
!------------------------------------------------------------------------

    if (restart) then
        tsini = ts_res + 1
    else
        tsini = 1
    endif
  
    tsfin = tsini - 1 + nstep

    !--- Allocate variables base on time ---!
    allocate(h_time(tsfin))
    allocate(A_time(tsfin))
    allocate(utp_time(tsfin))

    allocate(P_lat_Tavg(tsfin))
    allocate(P_fric_R_Tavg(tsfin))
    allocate(P_fric_S_Tavg(tsfin))
    allocate(P_h_Tavg(tsfin))
    allocate(P_w_Tavg(tsfin))
    allocate(P_pot_Tavg(tsfin))
    allocate(P_fric_R_visc_Tavg(tsfin))
    allocate(P_fric_S_visc_Tavg(tsfin))
    allocate(P_fric_R_plas_Tavg(tsfin))
    allocate(P_fric_S_plas_Tavg(tsfin))
    allocate(P_fric_R_vZ_pE_Tavg(tsfin))
    allocate(P_fric_R_pZ_vE_Tavg(tsfin))

    
!------------------------------------------------------------------------
!     Define a flag for the precond (T) or solver (F)
!------------------------------------------------------------------------

    p_flag = .true.
    if (solver .eq. 1) p_flag = .false.
   
    output_diag = .false.

!------------------------------------------------------------------------
!     initial conditions
!------------------------------------------------------------------------

    call ini_get (u, restart, expres, ts_res)
    call output_results(ts, expnb, solver, u, zeta, eta)
    un1=u
    hn1=h
    An1=A
    tauair = 0d0 ! initialization (watchout for restart)
    nbhr = 0d0
    fgmres_per_ts = 0

!   print*, u
    count = 1
    do ts = tsini, tsfin
     
        nbhr = nbhr + Deltat / 3600d0
        print *, 'time level, cumulative time (h) =', ts, nbhr

        utp_time(count) = u(int(nx/2))
        h_time(count)   = h(int(nx/2))
        A_time(count)   = A(int(nx/2))

        call cpu_time(timecrap)
        call cpu_time(time1)

        if ( BDF2 .eq. 1 ) un2 = un1 ! BDF2 needs u at 3 time levels
                                    ! Attention not initialized the 1st time level.

        if ( adv_scheme .eq. 'semilag') then ! semilag is 3 time level scheme
            hn2 = hn1
            An2 = An1
        endif

!------------------------------------------------------------------------
!     update previous time level solutions
!------------------------------------------------------------------------

        un1=u
        hn1=h
        An1=A
        if (oceanSIM) then 
            uwn2   = uwn1 
            uwn1   = uw
            etawn2 = etawn1
            etawn1 = etaw
        endif
   
        if (IMEX .eq. 0) then
            
            if (rheo .eq. 2) then
                call ice_strength (hn1, An1) ! standard approach no IMEX 
            else 
                call ice_strength (hn1, An1, un1) 
            endif
        
        endif

!------- get wind forcing (independent of u) -----------------------------

        call wind_forcing (tauair, ts)
     
!------- Solves NL mom eqn at specific time step with solver1, 2 or 3
!        F(u) = A(u)u - b(u) = 0, u is the solution vector
!------- Begining of outer loop (OL) or Newton iterations ----------------
  
        if (solver .eq. 1 .or. solver .eq. 2 ) then ! implicit

            if (solver .eq. 2 ) call calc_scaling (An1) ! scaling=1 for other solvers
     
            do k = 1, Nmax_OL 
        
                if (IMEX .gt. 0) then ! IMEX method 1 or 2
                    call advection (un1, u, hn1, An1, hn2, An2, h, A) ! advect tracers
                    ! call ice_strength (h, A) ! Pp_half is Pp/2 where Pp is the ice strength (Tp_half: tensile strength)
                    ! call shear(un1)
                    call ice_strength (h, A, u) ! standard approach no IMEX 
                    ! call inertial_number()
                    ! call angle_friction_mu()
                endif


                call viscouscoefficient (u, zeta, eta) ! u is u^k-1
                call Cw_coefficient (u, Cw, Cb)            ! u is u^k-1
                call calc_R (u, zeta, eta, Cw, Cb, tauair, R_uk1)
                call Fu (u, un1, un2, h, R_uk1, F_uk1) 

                L2norm = sqrt(DOT_PRODUCT(F_uk1,F_uk1))
                ! print*, L2norm
         
                if (k .eq. 1) then  
                    nl_target = gamma_nl*L2norm
                    ! nl_target = 1d-05
                    !	  call output_ini_L2norm(ts,L2norm,expnb)
                endif

                if (L2norm .lt. nl_target .or. L2norm .lt. 1d-08) exit

                if (solver .eq. 1) then
                    print *, 'L2-norm after k ite=', ts, k-1, L2norm
                    
                    call bvect(tauair, un1, Cw, b)
                    ! call output_sor(ts,k,solver, expnb,b)
                    call SOR (b, u, un1, h, A, zeta, eta, Cw, Cb, p_flag, ts)
      !             call SOR_A (b, u, zeta, eta, Cw, k, ts)
                elseif (solver .eq. 2) then
                    call prepFGMRES_NK(u, h, A, F_uk1, zeta, eta, Cw, Cb, un1, un2, tauair, &
                                    L2norm, k, ts, fgmres_its)
      !           call SOR_J(u, F_uk1, zeta, eta, Cw, upts, tauair, k, ts)
                endif
                
                fgmres_per_ts = fgmres_per_ts + fgmres_its
                
                if (k .eq. Nmax_OL) Nfail = Nfail + 1

            enddo
    
            meanN = meanN + k-1


            if (output_diag) then 
                call output_residual(ts,k,expnb,F_uk1)
                call output_nb_ite (ts, k ,fgmres_per_ts, expnb)
                call output_ini_L2norm(ts,L2norm,expnb)
            endif

   !  call output_residual()

        elseif (solver .eq. 3 .or. solver .eq. 4) then ! explicit (EVP)
     
            call viscouscoefficient (u, zeta, eta) ! u is u^k-1
            call Cw_coefficient (u, Cw, Cb)            ! u is u^k-1
            call EVP2solver(tauair, u, ts)
     
        endif

!------------------------------------------------------------------------
!     Diagnostic of energy
!------------------------------------------------------------------------   
        if (mechenergy) then 
            call mechanical_energy(u, zeta, eta, Cw)
            call average_mech_energy(P_pot_avg, P_fric_R_avg, P_fric_S_avg, P_fric_R_visc_avg, P_fric_S_visc_avg, &
                                P_fric_R_plas_avg, P_fric_S_plas_avg,P_lat_avg, P_h_avg, P_w_avg, &
                                P_fric_R_vZ_pE_avg, P_fric_R_pZ_vE_avg)
            
            P_pot_Tavg(ts) = P_pot_avg
            P_fric_R_Tavg(ts) = P_fric_R_avg
            P_fric_S_Tavg(ts) = P_fric_S_avg
            P_fric_R_visc_Tavg(ts) = P_fric_R_visc_avg
            P_fric_S_visc_Tavg(ts) = P_fric_S_visc_avg
            P_fric_R_plas_Tavg(ts) = P_fric_R_plas_avg
            P_fric_S_plas_Tavg(ts) = P_fric_S_plas_avg
            P_lat_Tavg(ts) = P_lat_avg
            P_h_Tavg(ts) = P_h_avg
            P_w_Tavg(ts) = P_w_avg
            P_fric_R_pZ_vE_Tavg(ts) = P_fric_R_pZ_vE_avg
            P_fric_R_vZ_pE_Tavg(ts) = P_fric_R_vZ_pE_avg

        endif

!------------------------------------------------------------------------
!     Diagnostic of ice-ocean stress
!------------------------------------------------------------------------     
        if (DiagStress) then
            call Cw_coefficient (u, Cw, Cb)
            call calc_diag_stress (idiag, u, Cw, tauair)
        endif
!------------------------------------------------------------------------       
     
        call cpu_time(time2)
        print *, 'cpu time = ', time2-time1

        if (IMEX .eq. 0) call advection (un1, u, hn1, An1, hn2, An2, h, A)


        if (rheo .eq. 2) call energydissipation(u, zeta, eta, P_half)

!     call meantracer(h,meanvalue)

!------------------------------------------------------------------------
!     Shallow water model
!------------------------------------------------------------------------
  
        if (oceanSIM) then
       
            call advect_etaw (etaw)
            call Cw_coefficient (u, Cw, Cb)
        
            if (IMEX .eq. 0) then
                call momentum_uw (idiag, tauair, Cdair, Cw, An1, u)
            elseif (IMEX .gt. 0) then
                call momentum_uw (idiag, tauair, Cdair, Cw, A, u)
            endif
            
            if (Asselin) then
	            call Asselin_filter (etaw, uw)
            endif
        endif

!------------------------------------------------------------------------
!     output results
!------------------------------------------------------------------------

        if (any(ts == out_step(1:21))) then
            print *, 'outputting results'
            call output_results(ts, expnb, solver, u, zeta, eta)
            call output_file(expnb)
        endif

!------------------------------------------------------------------------
!     calculate diagnostics            
!------------------------------------------------------------------------

!     call check_neg_vel(u)
        call minmaxtracer(h,1,ts)
        call minmaxtracer(A,2,ts)
        call minmaxtracer(u,3,ts)
        
        if (oceanSIM) then
            call minmaxtracer(etaw,4,ts)
            call minmaxtracer(uw,5,ts)
            if (DiagStress) call output_diag_stress (ts, expnb, idiag)
        endif
        
        count = count+1

    enddo
  
    if (solver .eq. 1) then
        print *, 'Nb failures, mean ite of Picard: ', Nfail, (meanN*1d0)/(nstep*1d0)
    elseif (solver .eq. 2) then
        print *, 'Nb failures, mean ite of JFNK: ', Nfail, (meanN*1d0)/(nstep*1d0)
        print *, 'mean nb of fgmres it per time level: ', fgmres_per_ts/(nstep*1d0)
    endif

    call output_mech_energy(P_pot_Tavg, P_fric_R_Tavg, P_fric_S_Tavg,P_fric_R_visc_Tavg, &
                            P_fric_S_visc_Tavg, P_fric_R_plas_Tavg, P_fric_S_plas_Tavg, P_lat_Tavg, &
                            P_h_Tavg, P_w_Tavg, P_fric_R_vZ_pE_Tavg, P_fric_R_pZ_vE_Tavg, &
                            ts, expnb, solver, nstep)

    !call output_times(ts, expnb, solver, nstep, A_time, h_time, utp_time)
    deallocate(A_time, h_time, utp_time) 
    deallocate(etaw, etawn1, etawn2, uw, uwn1, uwn2) 
    if (oceanSIM) then
        deallocate(duwdt, gedetawdx, tauiw, tauaw, buw) 
    endif
  
end program ice