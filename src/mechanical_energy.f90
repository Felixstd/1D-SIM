! This is a subroutine to compute the mechanical energy
! of sea ice. 
!
! We have : 
!   p_k = p_h + p_a + p_w + p_i + p_g
!
! 

subroutine mechanical_energy(utp, zeta, eta, Cw)
    
    use size
    use rheology
    use global_var
    use muphi
    use option
    use properties
    use resolution
    
    implicit none 

    integer :: i, count_v, count_p
    double precision, intent(in) :: utp(0:nx+1), zeta(0:nx+1), eta(0:nx+1), Cw(1:nx+1)
    double precision :: sig11(0:nx+1), div(0:nx+1)
    double precision :: shear, flux, a_at_u_i, a_at_u_ip1, calc_flux


    ! - Compute the stress tensor -!
    do i = 1, nx
        div(i) = ( utp(i+1) - utp(i) ) / Deltax
        sig11(i) = (eta(i)+zeta(i))*div(i) - P_half(i)
    enddo
    count_v = 0
    count_p = 0
    !-- Compute the different energies --!
    do i = 1, nx

        !- Only where theres ice
        if (A(i) > 0.0) then 
            P_pot(i) = P_half(i)*div(i)

            if (abs(zeta(i) - eta_max) < 1d-12) then
                ! print*, 'here'
                P_fric_R_visc(i) = -zeta(i)*div(i)**2d0
                P_fric_S_visc(i) = -eta(i)*sqrt(div(i)**2d0)**2d0
                count_v  = count_v + 1
            else
                P_fric_R_plas(i) = -zeta(i)*div(i)**2d0
                P_fric_S_plas(i) = -eta(i)*sqrt(div(i)**2d0)**2d0
                count_p  = count_p + 1
            endif
            
            P_fric_R(i) = -zeta(i)*div(i)**2d0
            P_fric_S(i) = -eta(i)*sqrt(div(i)**2d0)**2d0

            flux=calc_flux(utp(i),utp(i+1),sig11(i-1),sig11(i), sig11(i+1))
            P_lat(i) = flux/Deltax
            
            flux=calc_flux(utp(i),utp(i+1),h(i-1),h(i), h(i+1))
            P_h(i) = 1d0/2d0*rho*utp(i)**2d0*(-flux/Deltax)

            a_at_u_i = ( A(i) + A(i-1) ) / 2d0
            a_at_u_i = max(a_at_u_i, 1d-10)

            a_at_u_ip1 = ( A(i+1) + A(i) ) / 2d0
            a_at_u_ip1=max(a_at_u_ip1, 1d-10)

            P_w(i) = -(utp(i)*a_at_u_i*Cw(i) + utp(i+1)*a_at_u_ip1*Cw(i+1))/2
        ! else
        !     P_pot(i) = 0d0
        !     P_fric_R(i) = 0d0
        !     P_fric_S(i) = 0d0
        !     P_lat(i) = 0d0
        !     P_h(i) = 0d0
        !     P_w(i) = 0d0
        endif


    enddo

    print*, 'Number of viscous and plastic grid cells ', count_p, count_v 

end subroutine mechanical_energy

subroutine average_mech_energy(P_pot_avg, P_fric_R_avg, P_fric_S_avg, P_fric_R_visc_avg, P_fric_S_visc_avg, &
                                P_fric_R_plas_avg, P_fric_S_plas_avg, P_lat_avg, P_h_avg, P_w_avg)

    use size
    use rheology
    use global_var
    use muphi
    use option
    use properties
    use resolution
    
    implicit none 

    integer :: i, count, count_pfrv, count_pfrp, count_pfsv, count_pfsp
    double precision, intent(out) :: P_pot_avg, P_fric_R_avg, P_fric_S_avg, P_lat_avg, P_h_avg, P_w_avg
    double precision, intent(out) :: P_fric_R_visc_avg, P_fric_S_visc_avg, P_fric_R_plas_avg, P_fric_S_plas_avg
    ! Create logical mask for A > 0.15
    
    P_pot_avg        = 0.0
    P_fric_R_avg     = 0.0
    P_fric_S_avg     = 0.0
    P_lat_avg        = 0.0
    P_h_avg          = 0.0
    P_w_avg          = 0.0
    P_fric_R_visc_avg = 0.0
    P_fric_S_visc_avg = 0.0
    P_fric_R_plas_avg = 0.0
    P_fric_S_plas_avg = 0.0
    count_pfrv = 0
    count_pfrp = 0
    count_pfsv = 0
    count_pfsp = 0

    do i = 1, nx
        ! if (A(i) > 0.001) then
            count = count + 1
            P_pot_avg        = P_pot_avg        + P_pot(i)
            P_fric_R_avg     = P_fric_R_avg     + P_fric_R(i)
            P_fric_S_avg     = P_fric_S_avg     + P_fric_S(i)
            P_lat_avg        = P_lat_avg        + P_lat(i)
            P_h_avg          = P_h_avg          + P_h(i)
            P_w_avg          = P_w_avg          + P_w(i)
            P_fric_R_visc_avg = P_fric_R_visc_avg + P_fric_R_visc(i)
            P_fric_S_visc_avg = P_fric_S_visc_avg + P_fric_S_visc(i)
            P_fric_R_plas_avg = P_fric_R_plas_avg + P_fric_R_plas(i)    
            P_fric_S_plas_avg = P_fric_S_plas_avg + P_fric_S_plas(i)
            ! if (abs(P_fric_R_visc(i)) > 0d0) then
            !     P_fric_R_visc_avg = P_fric_R_visc_avg + P_fric_R_visc(i)
            !     count_pfrv = count_pfrv + 1
            ! endif
            
            ! if (abs(P_fric_S_visc(i)) > 0d0) then
            !     P_fric_S_visc_avg = P_fric_S_visc_avg + P_fric_S_visc(i)
            !     count_pfsv = count_pfsv + 1
            ! endif

            ! if (abs(P_fric_R_plas(i)) > 0d0) then
            !     P_fric_R_plas_avg = P_fric_R_plas_avg + P_fric_R_plas(i)
            !     count_pfrp = count_pfrp + 1
            ! endif

            ! if (abs(P_fric_S_plas(i)) > 0d0) then
            !     P_fric_S_plas_avg = P_fric_S_plas_avg + P_fric_S_plas(i)
            !     count_pfsp = count_pfsp + 1
            ! endif
        ! endif
    enddo


    ! if (count > 0) then
    !     P_pot_avg        = P_pot_avg        / count
    !     P_fric_R_avg     = P_fric_R_avg     / count
    !     P_fric_S_avg     = P_fric_S_avg     / count
    !     P_lat_avg        = P_lat_avg        / count
    !     P_h_avg          = P_h_avg          / count
    !     P_w_avg          = P_w_avg          / count
    !     print* ,'HERE', P_fric_R_visc_avg
    !     P_fric_R_visc_avg = P_fric_R_visc_avg / count_pfrv
    !     P_fric_S_visc_avg = P_fric_S_visc_avg / count_pfsv
    !     P_fric_R_plas_avg = P_fric_R_plas_avg / count_pfrp
    !     P_fric_S_plas_avg = P_fric_S_plas_avg / count_pfsp
    ! end if
    ! print*, P_fric_R_visc_avg

end subroutine average_mech_energy
