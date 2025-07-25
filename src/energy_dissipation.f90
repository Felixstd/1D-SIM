subroutine energydissipation(utp, zeta, eta, P)

	use size
	use resolution
	use rheology
	use option
	use global_var
	use muphi
    use properties

    implicit none

	integer :: i

	double precision, intent(in) :: utp(1:nx+1), zeta(0:nx+1), eta(0:nx+1), P(0:nx+1)
	double precision :: dudx, deno, denonum, denomin

    do i = 1, nx
        dudx = ( utp(i+1) - utp(i) ) / Deltax

        W_sigma(i) = ((zeta(i) + eta(i))*dudx - P(i))*dudx
    
    enddo

end subroutine