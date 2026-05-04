!--------------------------------------------------------------------
! Those are a collection of subroutines used to define the 
! model parameters. 
! 
! This file contains two subroutines:
! 1. get_default: default values for parameters
! 2. read_namelist: read the parameters defined in the namelist if provided. 
!
! Written by FSTD mostly taken from the SIM model. 
!--------------------------------------------------------------------


subroutine get_default

    use size
    use rheology
    use muphi
    use forcing
    use properties
    use resolution
    use global_var
    use shallow_water
    use numerical
    use option

    implicit none 

    rheo           = 1
    linear_drag    = .false.
    linear_viscous = .false. ! linear viscous instead of viscous-plastic
    constant_wind  = .true.  ! T: 10m/s, F: spat and temp varying winds
    rampupwind     = .false.
    uwind          = 10d0    ! uwind velocity for constant_wind
    rep_closure    = .false. ! replacement closure (see Kreysher et al. 2000)
    regularization = 'tanh'  ! tanh, Kreyscher, capping (Hibler)
    adv_scheme     = 'upwind'! upwind, upwindRK2, semilag
    oceanSIM       = .false. ! for shallow water model
    implicitDrag   = .false. ! for uwater mom eq.
    Asselin        = .false. ! Asselin filter for uw and etaw
    DiagStress     = .false. ! diagnostic for ice-ocean stress at i=idiag
    Agamma         = 1d-02   ! Asselin filter parameter
    T_tot          = 30000 ! Total length of the simulation [s]
    Deltat         = 0.1! time step [s]



    !---------------------------------------------------------
    ! Well-posedness and energy options (FSTD)
    !---------------------------------------------------------
    initcond       = 'stepsmooth' ! initial conditions for h and A
    initcond_vel   = 'Graysmooth' ! initial conditions for u 
    mechenergy     = .false.      ! compute the power dissipations
    advection_mom  = .false.      ! advection term in solver
    Pstart_change  = .false.      ! new P parametrization
    P0_constant    = .false.


    !---------------------------------------------------------
    ! Solver options 
    !---------------------------------------------------------
    solver     = 1    ! 1: Picard+SOR, 2: JFNK, 3: EVP, 4: EVP*
    IMEX       = 0    ! 0: no IMEX, 1: Jdu=-F(IMEX), 2: J(IMEX)du=-F(IMEX) 
    BDF2       = 0    ! 0: standard, 1: Backward difference formula (2nd order)

    Nmax_OL    = 1000

    N_sub = 900

    omega      = 1d0    ! relax parameter for SOR
    tol_SOR    = 1d-10  ! tol for SOR solver
    maxiteSOR  = 10000  ! max nb of ite for SOR
    iteSOR_pre = 10     ! nb of iterations for the SOR precond
    maxiteGMRES= 50     ! max nb of ite for GMRES
    gamma_nl = 1d-04
    dropini  = 1.5d0    ! defines initial drop in L2norm before gamma = 0.01
    small1   = 1d-10    ! to have a continuously diff water drag term
    small2   = 1d-22    ! to have a continuously diff rheology term
    smallA   = 1d-03    ! for num stab of Atw and Ata (in zones with ~no ice)


!------------------------------------------------------------------------
!     Define Deltax and check CFL based on input by user
!------------------------------------------------------------------------

    if ( nx .eq. 100 ) then 
        Deltax   =  20d03  ! grid size [m], the domain is always 2000 km 
    elseif (nx .eq. 50) then 
        Deltax = 10d03
    elseif  ( nx .eq. 200 ) then
        Deltax   =  1            
    elseif  ( nx .eq. 500 ) then
        Deltax   =  1000
    elseif  ( nx .eq. 625 ) then
        Deltax   =  800
    elseif  ( nx .eq. 800 ) then
        Deltax   =  625
    elseif  ( nx .eq. 1000 ) then
        Deltax   =  500
    elseif  ( nx .eq. 2000 ) then
        Deltax   =  250
    elseif  ( nx .eq. 4000 ) then
        Deltax   =  125
    elseif  ( nx .eq. 5000 ) then
        Deltax   =  100
    elseif  ( nx .eq. 8000 ) then
        Deltax   =  62.5d0
    
    else
        print *,  'Wrong grid size dimension', nx
        STOP
    endif

    Deltax2 = Deltax ** 2
    DtoverDx = Deltat / Deltax

!------------------------------------------------------------------------
!     Define constants
!------------------------------------------------------------------------

    C          = 20d0         ! ice strength parameter (watchout no A for now)
    Pstar      = 27.5d03! ice compression strength parameter
    e          = 2d0    ! ratio long to short axis of ellipse
    denomin_P  = 2d-09

    e_2        = 1/(e**2d0)   !
    alpha      = sqrt(1d0 + e_2)
    alpha2     = 1d0 + e_2
    kt         = 0d0          ! T = kt * P (1.0 in Konig and Holland, 2010)

    Cdwater    = 5.5d-30! water-ice drag coeffient
    Cdair      = 0d0    ! air-ice drag coeffient 
    Cdairw     = 0d0    ! air-water drag coeffient  
    rhoair     = 1.3d0        ! air density
    rho        = 900d0        ! ice density
    rhowater   = 1026d0       ! water density
    ge         = 9.8d0        ! Earth's gravitional acceleration
    Hw         = 2.5d0          ! mean water depth (for shallow water model)
    bw         = 0.0005d0/Hw  ! friction term for the uw momentum eq (always implicit).

    Cda        = rhoair   * Cdair
    Cdw        = rhowater * Cdwater
  
!------------------------------------------------------------------------
!     Define Constants specific to GC rheology (rheology = 2)
!------------------------------------------------------------------------
    d_average  = 1000
    mu_0       = 0.2
    mu_infty   = 0.8
    mu_b       = 1.05d0
    I_0        = 1e-3
    Phi_0      = 1
    c_phi      = 1
    D = 0.0001 !DEFAULT was -0.00001
    n = 2 !for super gaussian 
    zeta_max   = 1d9
    eta_max    = 1d9

end subroutine get_default


subroutine read_namelist 

    use size
    use rheology
    use muphi
    use forcing
    use properties
    use resolution
    use global_var
    use shallow_water
    use numerical
    use option

    implicit none

    
    integer :: nml_error, filenb
    character filename*32

    !---- namelist variables -------
            
    namelist /option_nml/ &
        rheo, linear_drag,                                        &
        linear_viscous, constant_wind, rampupwind,                &
        uwind, rep_closure, regularization, adv_scheme, oceanSIM, &
        implicitDrag, Asselin, DiagStress,  solver, IMEX, BDF2,   &
        initcond, initcond_vel, mechenergy, advection_mom, Pstart_change, &
        P0_constant
  

    namelist /numerical_param_nml/ &
        Deltat, T_tot, Agamma, gamma_nl, Nmax_OL, tol_SOR, maxiteSOR

    namelist /phys_param_nml/ &
        Pstar, C, e, rhoair, rho, rhowater, &
        Cdair, Cdwater, zeta_max, denomin_P

    filename ='namelistSIM'
    filenb = 10

    print *, 'Reading namelist values'
    
    open (filenb, file=filename, status='old',iostat=nml_error)
    if (nml_error /= 0) then
        nml_error = -1
    else
        nml_error =  1
    endif
        
    do while (nml_error > 0)
        print*,'Reading option_nml'
        read(filenb, nml=option_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading other_nml'
        read(filenb, nml=numerical_param_nml,iostat=nml_error)
        if (nml_error /= 0) exit
        print*,'Reading phys_param_nml'
        read(filenb, nml=phys_param_nml,iostat=nml_error)
        print *, nml_error
    enddo

    close(filenb)

    e_2        = 1/(e**2d0)   !
    alpha      = sqrt(1d0 + e_2)
    alpha2     = 1d0 + e_2

    Cda        = rhoair   * Cdair
    Cdw        = rhowater * Cdwater

end subroutine read_namelist