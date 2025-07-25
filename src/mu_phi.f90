subroutine angle_friction_mu(Inertial_num, muI)

    use muphi
    use global_var
    use resolution

    implicit none

    integer i, j
    double precision eps, mu

    double precision, intent(in) :: Inertial_num
    double precision, intent(out) :: muI

    eps = 1d-12

    do i = 1, nx

        muI = max(mu_0 + ( mu_infty - mu_0 ) / ( I_0/Inertial_num + 1), mu_0)
        ! mu_I(i) = max(mu_0 + ( mu_infty - mu_0 ) / ( I_0/Inertial + 1), mu_0)
        ! mu_I(i) = mu_0 + ( mu_infty - mu_0 ) / ( I_0/Inertial(i) + 1)

        ! print*, mu_I

    enddo

! 

end subroutine angle_friction_mu


subroutine inertial_number(shear, htp, P, Inertial_num)

    use muphi
    use global_var
    use resolution
    use properties

    implicit none


    integer i, j
    double precision eps, mu
    double precision, intent(in) :: shear, htp, P
    double precision, intent(out) :: Inertial_num

    eps = 1d-12

    do i = 1, nx

        Inertial_num = min(d_average * shear * SQRT(rho * htp/P), 1d0)
        ! Inertial(i) = d_average * shear_I(i) * SQRT(rho * h(i)/Pp_half(i))
    enddo

    return

end subroutine inertial_number


subroutine shear(utp)
    
    use global_var
    use resolution
    use rheology

    implicit none

    integer i, j

    double precision dudx, dvdy, dudy, dvdx, land, lowA
    double precision, intent(in):: utp(0:nx+2)


    do i = 1, nx

        dudx = ( utp(i+1) - utp(i) ) / Deltax  
!----- stresses and strain rates at the grid center -------------------------   

        ! shear_I(i) = abs(dudx)
        shear_I(i) = sqrt( (dudx)**2d0+small2)

    enddo

    return
end subroutine shear