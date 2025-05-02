subroutine ice_strength(hin, Ain)

    use size
    use rheology
    use global_var

    implicit none 

    integer                      :: i
    double precision, intent(in) :: hin(0:nx+1), Ain(0:nx+1)


    do i = 0, nx+1

        Pp(i) = Pstar * h(i) * dexp(-C * ( 1d0 - A(i) ))


    enddo

    Pp(0) = 0d0
    Pp(nx+1) = 0d0
     
    return

end subroutine ice_strength