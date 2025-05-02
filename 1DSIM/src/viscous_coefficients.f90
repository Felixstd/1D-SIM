subroutine viscouscoefficient(zeta, eta)

    use size
    use rheology
    use muphi
    use global_var
    use option
    use resolution

    implicit none 

    integer                         :: i
    double precision, intent(out)   :: zeta(0:nx+1), eta(0:nx+1)
    double precision                :: dudx, deno, denonum, denomin
    double precision                :: eta_max, shearmax

    denomin=2d-09
    eta_max = 1d12

    do i = 0, nx+1

        if (rheo .eq. 1) then
            if (regularization .eq. 'tanh') then
            
                shearmax = max(shear_I(i), 1d-20)
                
                zeta(i) = 1d08
                ! zeta(i) = 2d08*Pp(i) * tanh((mu_b) / (shearmax * 2d08))
                ! eta(i)  = eta_max * tanh((mu_I(i)) * Pp(i) / (shearmax * eta_max))
                eta(i) = 1d08
            endif
            print*, eta(i), zeta(i)

            P(i) = Pstar
        
        elseif (rheo .eq. 2) then
      
            ! dudx = ( utp(i+1) - utp(i) ) / Deltax

            if ( regularization .eq. 'tanh' ) then

                deno = alpha*sqrt( (shear_I(i))**2d0 + 1d-12 ) ! small2 is there to avoid div by zero

                zeta(i) = ((Pp(i)/2)/denomin)*tanh(denomin*(1d0/deno))
            endif

            eta(i)  = zeta(i) * e_2

            if ( regularization .eq. 'tanh' ) then
                P(i) = (Pp(i)/2d0) * ( deno / denomin ) * tanh(denomin*(1d0/deno))
            endif
        ! P(i) = 0d0
        endif


    enddo

    return

end subroutine viscouscoefficient