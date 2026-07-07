subroutine u_linearization(utp, uk2, ul)

    use size
    use resolution
    use rheology

    implicit none 

    double precision, dimension(1:nx+1), intent(in) :: uk2, utp
    double precision, dimension(1:nx+1), intent(out) :: ul


    ul = (utp + uk2) / 2d0

end subroutine



