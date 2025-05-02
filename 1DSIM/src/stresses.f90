subroutine stress(utp, zeta, eta)

    use rheology
    use muphi
    use global_var
    use option
    use resolution

    implicit none 

    integer :: i
    double precision, intent(in)  :: utp(0:nx+1), zeta(0:nx+1), eta(0:nx+1)
    double precision              :: dudx

    do i = 0, nx+1

        if (rheo .eq. 1) then

            if (i .eq. nx+1) then
                dudx = (utp(i) - utp(i-1)) / Deltax
            else 
                dudx = (utp(i+1) - utp(i)) / Deltax

            endif

            sig11(i) = (zeta(i) + eta(i) / 2) * dudx - P(i)
            ! print*,'sig11', sig11(i)
        
        elseif (rheo .eq. 2) then

            if (i .eq. nx+1) then
                dudx = (utp(i) - utp(i-1)) / Deltax
            else 
                dudx = (utp(i+1) - utp(i)) / Deltax
            endif

            sig11(i) = (eta(i) + zeta(i))*dudx - P(i)
        endif


    enddo

    return

end subroutine