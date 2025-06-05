subroutine angle_friction_mu 

    use muphi
    use global_var
    use resolution

    implicit none

    integer i, j
    double precision eps, mu

    eps = 1d-12

    do i = 0, nx+1


        ! mu = mu_0 + ( mu_infty - mu_0 ) / ( I_0/(1-A(i) + 1))


        ! mu_I(i) = mu_infty * tanh(mu/mu_infty)

        mu_I(i) = max(mu_0 + ( mu_infty - mu_0 ) / ( I_0/Inertial(i) + 1), mu_0)
        ! mu_I(i) = mu_0 + mu_infty*(1-A(i))

    enddo

! 

end subroutine angle_friction_mu


subroutine inertial_number

    use muphi
    use global_var
    use resolution
    use properties

    implicit none


    integer i, j
    double precision eps, mu

    eps = 1d-12

    do i = 0, nx+1

        Inertial(i) = min(d_average * shear_I(i) * SQRT(rho * h(i)/Pp(i)), 1d0)
        ! Inertial(i) = 1

    enddo

    return

end subroutine inertial_number


subroutine shear(utp)
    
    use global_var
    use resolution

    implicit none


    ! include 'parameter.h'
    ! include 'CB_Dyndim.h'
    ! include 'CB_DynVariables.h'
    ! include 'CB_DynForcing.h'
    ! include 'CB_const.h'
    ! include 'CB_mask.h'
    ! include 'CB_options.h'

    integer i, j

    double precision dudx, dvdy, dudy, dvdx, land, lowA
    double precision, intent(in):: utp(0:nx+2)

    ! double precision shear(0:nx+1,0:ny+1)

    do i = 0, nx+1
        ! print*, utp(i)
        ! dudx = ( utp(i+1) - utp(i) ) / Deltax

        if (i .eq. nx+1) then
            dudx = ( utp(i) - utp(i-1) ) / Deltax
        else
            dudx = ( utp(i+1) - utp(i) ) / Deltax  
        endif  
!----- stresses and strain rates at the grid center -------------------------   

        shear_I(i) = abs(dudx)

    enddo

    return
end subroutine shear