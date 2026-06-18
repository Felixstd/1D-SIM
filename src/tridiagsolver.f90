!--------------------------------------------------------
! Subroutine used to solved the ice strength helmholtz equation. 
! The equation is a nonlocal regularization of the VP rheology. 
!
! It is given by 
!               P_xx = 1/(l^2)(P-P_H)
! where P_H: Hibler and l = l_E*A. 
!
! This equation is solved using a Tridiagonal solver. 
!
! FSTD
!--------------------------------------------------------


subroutine solve_helmholtz_P(P_loc, Ain, P_nl)

    
    use size
    use resolution
    use rheology
    use option
    
    implicit none

    double precision, intent(in)  :: P_loc(0:nx+1)
    double precision, intent(out) :: P_nl(0:nx+1)
    double precision, intent(in) :: Ain(0:nx+1)

    double precision :: lower(1:nx), main_d(1:nx), upper(1:nx), rhs(1:nx)
    double precision :: w, Ai, l_scale2i
    integer :: i

    ! Coefficients of l_scale2*d2P/dx2 - P = -P_loc
    ! lower(i) =  l_scale2/dx^2
    ! main_d(i) = -(2*l_scale2/dx^2 + 1)
    ! upper(i) =  l_scale2/dx^2
    ! rhs(i)   = -P_loc(i)

    do i = 1, nx
        Ai = max(Ain(i),1d-10)

        if (l_c_exp) then 
            l_scale2i = l_scale2*exp(-C*(1-Ai))
        else
            l_scale2i = (l_scale2*Ai**2d0)
        endif
        lower(i)  =  l_scale2i / Deltax2
        main_d(i) = -(2d0 * l_scale2i / Deltax2 + 1d0)
        upper(i)  =  l_scale2i / Deltax2
        rhs(i)    = -P_loc(i)
    enddo

    ! Boundary conditions: P_nl = P_loc at boundaries (Dirichlet)
    ! Enforce by modifying first and last rows
    main_d(1)  = 1d0
    upper(1)   = 0d0
    rhs(1)     = P_loc(1)

    main_d(nx) = 1d0
    lower(nx)  = 0d0
    rhs(nx)    = P_loc(nx)

    ! Thomas algorithm (forward sweep)
    do i = 2, nx
        w         = lower(i) / main_d(i-1)
        main_d(i) = main_d(i) - w * upper(i-1)
        rhs(i)    = rhs(i)    - w * rhs(i-1)
    enddo

    ! Back substitution
    P_nl(nx) = rhs(nx) / main_d(nx)
    do i = nx-1, 1, -1
        P_nl(i) = (rhs(i) - upper(i) * P_nl(i+1)) / main_d(i)
    enddo

    ! Ghost points
    P_nl(0)    = 0d0
    P_nl(nx+1) = 0d0

end subroutine solve_helmholtz_P


subroutine solve_timeHelmholtz_P(Pn1, P_loc, Ain, P_nl)
    
        use size
    use resolution
    use rheology
    
    implicit none

    double precision, intent(in)  :: Pn1(0:nx+1)
    double precision, intent(in)  :: P_loc(0:nx+1)
    double precision, intent(out) :: P_nl(0:nx+1)
    double precision, intent(in) :: Ain(0:nx+1)

    double precision :: lower(1:nx), main_d(1:nx), upper(1:nx), rhs(1:nx)
    double precision :: w, Ai
    integer :: i

    ! Coefficients of l_scale2*d2P/dx2 - P = -P_loc
    ! lower(i) =  l_scale2/dx^2
    ! main_d(i) = -(2*l_scale2/dx^2 + 1)
    ! upper(i) =  l_scale2/dx^2
    ! rhs(i)   = -P_loc(i)

    do i = 1, nx
        Ai = max(Ain(i),1d-10)
        lower(i)  =  -(l_scale2*Ai**2d0) / Deltax2
        main_d(i) = t_scale/(Deltat)+(2d0 * (l_scale2*Ai**2d0) / Deltax2 + 1d0)
        upper(i)  =  -(l_scale2*Ai**2d0) / Deltax2
        rhs(i)    = t_scale/Deltat*Pn1(i)+P_loc(i)
    enddo

    ! Boundary conditions: P_nl = P_loc at boundaries (Dirichlet)
    ! Enforce by modifying first and last rows
    main_d(1)  = 1d0
    upper(1)   = 0d0
    rhs(1)     = P_loc(1)

    main_d(nx) = 1d0
    lower(nx)  = 0d0
    rhs(nx)    = P_loc(nx)

    ! Thomas algorithm (forward sweep)
    do i = 2, nx
        w         = lower(i) / main_d(i-1)
        main_d(i) = main_d(i) - w * upper(i-1)
        rhs(i)    = rhs(i)    - w * rhs(i-1)
    enddo

    ! Back substitution
    P_nl(nx) = rhs(nx) / main_d(nx)
    do i = nx-1, 1, -1
        P_nl(i) = (rhs(i) - upper(i) * P_nl(i+1)) / main_d(i)
    enddo

    ! Ghost points
    P_nl(0)    = 0d0
    P_nl(nx+1) = 0d0


end subroutine 
